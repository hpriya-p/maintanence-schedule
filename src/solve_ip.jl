
using JSON
using JuMP, MathOptInterface, Cbc
const MOI = MathOptInterface

const ztol = 1e-6
const inf = 1e20

INTERVENTIONS = "Interventions"
NUM_SCENARIO = "Scenarios_number"
RESOURCES = "Resources"
SEASONS = "Seasons"
EXCLUSIONS = "Exclusions"
QUANTILE = "Quantile"
ALPHA = "Alpha"


function constructMILPfromFile(filename)
    parsed_file = JSON.parsefile(filename)

# Constants
    n = length(keys(get(parsed_file, INTERVENTIONS, [])))
    C = length(keys(get(parsed_file, RESOURCES, [])))
    T = get(parsed_file, "T", 0)
    St = get(parsed_file, NUM_SCENARIO, 0)
    S = max(St...)
    tau = get(parsed_file, QUANTILE, 0)
    alpha = get(parsed_file, ALPHA, 0)

    risk_arr = zeros(n, T, T, S) # (intervention, current time, start time, scenario)
# Populate risk_arr
    interventions_dict = Dict()
    ctr = 1
    for (intervention, int_values) in pairs(get(parsed_file, INTERVENTIONS, []))
        inter = ctr
            if(haskey(interventions_dict, intervention)) 
               inter = get(interventions_dict, intervention, ctr) 
            else 
            interventions_dict[intervention] = ctr
            end 
        for (curtime, ctime_values) in pairs(get(int_values, "risk", []))
            for (start_time, stime_values) in pairs(ctime_values)
                for s in 1:getindex(St, parse(Int64, curtime))
                    setindex!(risk_arr, stime_values[s], inter, parse(Int64, curtime), parse(Int64, start_time), s)
                end
            end
        end
        ctr += 1
    end

    M = max(risk_arr...) + 1

    r_arr = zeros(n, T, T, C) # (intervention, current time, start time, resource)
    resources_dict = Dict()
# Populate r_arr
    for (intervention, int_values) in pairs(get(parsed_file, INTERVENTIONS, []))
        ctr = 1
        for (resource, resource_values) in pairs(get(int_values, "workload", []))
            res = ctr
            if(haskey(resources_dict, resource)) 
               res = get(resources_dict, resource, ctr) 
            else 
            resources_dict[resource] = ctr
            end 

            for (curr_time, curr_time_values) in pairs(resource_values)
                for (start_time, start_time_load) in pairs(curr_time_values)
                    setindex!(r_arr, start_time_load, get(interventions_dict, intervention, 1), parse(Int64, curr_time), parse(Int64, start_time), res)
                end
            end
            ctr += 1
        end
    end

    u = zeros(C, T)
    for (res, res_values) in pairs(get(parsed_file, RESOURCES, []))
        setindex!(u, get(res_values, "max", []), get(resources_dict, res, 1), 1:T)
    end

    l = zeros(C, T)
    for (res, res_values) in pairs(get(parsed_file, RESOURCES, []))
        setindex!(l, get(res_values, "min", []), get(resources_dict, res, 1), 1:T)
    end

    deltas = zeros(n, T)
    for (intervention, int_values) in pairs(get(parsed_file, INTERVENTIONS, []))
        setindex!(deltas, get(int_values, "Delta", []), get(interventions_dict, intervention, 1), 1:T)
    end

# Construct JuMP Model
    j_model = Model()
    t[1:n] = @variable(j_model, t[1:n], Int)
    d[1:n] = @variable(j_model, d[1:n], Int)
    y[1:n, 1:T] = @variable(j_model, y[1:n, 1:T], Bin)
    z1[1:n, 1:T] = @variable(j_model, z1[1:n, 1:T], Bin)
    z2[1:n, 1:T] = @variable(j_model, z2[1:n, 1:T], Bin)
    a[1:n, 1:T] = @variable(j_model, a[1:n, 1:T], Bin)
    c[1:n, 1:T, 1:T] = @variable(j_model, c[1:n, 1:T, 1:T], Bin) # of the form c[intervention, current time, start time]
    risk[1:S, 1:T] = @variable(j_model, risk[1:S, 1:T])
    w[1:S, 1:T] = @variable(j_model, w[1:S, 1:T], Bin)
    q[1:T] = @variable(j_model, q[1:T])
    e[1:T] = @variable(j_model, e[1:T])
    mean_risk_t = [sum(risk[s, t1] for s in 1:getindex(St, t1)) / getindex(St, t1) for t1 in 1:T]
    mean_risk = sum(getindex(mean_risk_t, t1) for t1 in 1:T) / T
    quant_risk = [q[t1] - getindex(mean_risk_t, t1) for t1 in 1:T]
    vars = Dict("t"=>t, "y"=>y, "q"=>q, "w"=>w, "a"=>a, "c"=>c, "d"=>d, "z1"=>z1, "z2"=>z2, "risk"=>risk, "e"=>e, "quant_risk"=>quant_risk, "mean_risk"=>mean_risk_t)
    println("Adding constraints 1-25 (without 22)...")
    constr1 = [@constraint(j_model, t[i] >= time * y[i,time]) for i in 1:n for time in 1:T]
    constr2 = [@constraint(j_model, t[i] <= time + T * (1 - y[i, time])) for i in 1:n for time in 1:T]
    constr3a = [@constraint(j_model, z1[i, t1] <= sum(y[i, j] for j in 1:t1) + ztol)  for i in 1:n for t1 in 1:T]
    constr3b = [@constraint(j_model, z1[i, t1] >= sum(y[i, j] for j in 1:t1) - ztol)  for i in 1:n for t1 in 1:T]
    constr4 = [@constraint(j_model, d[i] >= t1 * z2[i, t1]) for i in 1:n for t1 in 1:T]
    constr5 = [@constraint(j_model, d[i] <= t1 + (T + 1) * z2[i, t1]) for i in 1:n for t1 in 1:T]
    constr6 = [@constraint(j_model, a[i, t1] >= z1[i, t1] + z2[i, t1] - 1) for i in 1:n for t1 in 1:T]
    constr7 = [@constraint(j_model, a[i, t1] <= z1[i, t1]) for i in 1:n for t1 in 1:T]
    constr8 = [@constraint(j_model, a[i, t1] <= z2[i, t1]) for i in 1:n for t1 in 1:T]
    constr9 = [@constraint(j_model, c[i, time1, time2] >= y[i, time2] + a[i, time1] - 1) for i in 1:n for time1 in 1:T for time2 in 1:T]
    constr10 = [@constraint(j_model, c[i, time1, time2] <= y[i, time2]) for i in 1:n for time1 in 1:T for time2 in 1:T]
    constr11 = [@constraint(j_model, c[i, time1, time2] <= a[i, time1]) for i in 1:n for time1 in 1:T for time2 in 1:T]
    constr12 = [@constraint(j_model, d[i] <= time + deltas[i, time] + T * (1 - y[i, time])) for i in 1:n for time in 1:T]
    constr13 = [@constraint(j_model, d[i] >= time + deltas[i, time] - T * (1 - y[i, time])) for i in 1:n for time in 1:T]
    constr14 = [@constraint(j_model, risk[s, t1] == sum(getindex(risk_arr, i, t1, t2, s) * c[i, t1, t2] for i in 1:n for t2 in 1:T)) for t1 in 1:T for s in 1:S]
    constr15 = [@constraint(j_model, q[t1] >= risk[s, t1] - M * (1 - w[s, t1])) for s in 1:S for t1 in 1:T]
    constr16 = [@constraint(j_model, q[t1] >= 0) for t1 in 1:T]
    constr17 = [@constraint(j_model, e[t1] >= 0) for t1 in 1:T]
    constr18 = [@constraint(j_model, e[t1] >= getindex(quant_risk, t1)) for t1 in 1:T]
    constr17 = [@constraint(j_model, d[i] <= T + 1) for i in 1:n]
    #constr18 = [@constraint(j_model, d[i] >= 0) for i in 1:n]
    constr19 = [@constraint(j_model, sum(y[i, time] for time in 1:T) == 1) for i in 1:n]
    constr20 = [@constraint(j_model, sum(getindex(r_arr, i, t1, t2, cj) * c[i, t1, t2] for i in 1:n for t2 in 1:T) >= l[cj, t1])  for t1 in 1:T for cj in 1:C]
    constr21 = [@constraint(j_model, sum(getindex(r_arr, i, t1, t2, cj) * c[i, t1, t2] for i in 1:n for t2 in 1:T) <= u[cj, t1])  for t1 in 1:T for cj in 1:C]
    constr23 = [@constraint(j_model, sum(w[s, t1] for s in 1:getindex(St, t1)) >= tau * getindex(St, t1)) for t1 in 1:T]
    constr24 = [@constraint(j_model, t[i] <= T) for i in 1:n]
    constr25 = [@constraint(j_model, t[i] >= 0) for i in 1:n]

    println("Adding Constraint 22...")
# Constraint 22
    for (exc, exc_values) in pairs(get(parsed_file, EXCLUSIONS, []))
        i1 = get(interventions_dict, getindex(exc_values, 1), 1)
        i2 = get(interventions_dict, getindex(exc_values, 2), 1)
        season = getindex(exc_values, 3)
        season_vals = [parse(Int64, x) for x in get(get(parsed_file, SEASONS, []), season, [])]
        [@constraint(j_model, a[i1, time] + a[i2, time] <= 1) for time in season_vals]
    end
    println("Adding objective...")
    @objective(j_model, Min, alpha * mean_risk + (1 - alpha) * sum(e[t1] for t1 in 1:T) / T)
    return j_model, vars, Dict("Interventions" => interventions_dict, "Resources" => resources_dict)
end

milp_model, var_dict, mappings = constructMILPfromFile("example1.json")
println("Constructed MILP Model")
set_optimizer(milp_model, Cbc.Optimizer)
optimize!(milp_model)
println(value.(get(var_dict, "t", [])))
