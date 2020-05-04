
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
        if (haskey(interventions_dict, intervention))
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
            if (haskey(resources_dict, resource))
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
    println(deltas)

# Construct JuMP Model
    j_model = Model()
    @variable(j_model, x[1:n, 1:T], Bin)
    @variable(j_model, y[1:n, 1:T], Bin)
    @variable(j_model, z[1:n, 1:T, 1:T], Bin) # of the form z[intervention, current time, start time]
    @variable(j_model, risk[1:S, 1:T])
    @variable(j_model, w[1:S, 1:T], Bin)
    @variable(j_model, q[1:T] >= 0)
    @variable(j_model, e[1:T] >= 0)
    mean_risk_t = [sum(risk[s, t1] for s in 1:getindex(St, t1)) / getindex(St, t1) for t1 in 1:T]
    mean_risk = sum(getindex(mean_risk_t, t1) for t1 in 1:T) / T
    quant_risk = [q[t1] - getindex(mean_risk_t, t1) for t1 in 1:T]
    println(string("T:", T, " n:", n, " S:", S, " C:", C))
    print("Constraint 1...")
    [@constraint(j_model, sum(y[i, t] for t in 1:T) == 1) for i in 1:n]
    println("Complete")
    print("Constraint 2...")
    [@constraint(j_model, sum((deltas[i, j] + j)* y[i, j] for j in 1:T) <= T+1) for i in 1:n]
    println("Complete")
    print("Constraint 3...")
    [@constraint(j_model, sum((deltas[i, j] + j)* y[i, j] for j in 1:T) >=(t+1)*x[i,t]) for t in 1:T for i in 1:n]
    println("Complete")
    print("Constraint 4...")
    [@constraint(j_model, sum((deltas[i, j] + j)* y[i, j] for j in 1:T) <=t  + (T+1)*x[i,t]) for t in 1:T for i in 1:n]
    println("Complete")
    print("Constraint 5...")
    [@constraint(j_model, z[i, t1, t2] >= y[i, t2] + sum(y[i, j] for j in 1:t1) + x[i, t1] - 2) for i in 1:n for t1 in 1:T for t2 in 1:T]
    println("Complete")
    print("Constraint 6...")
    [@constraint(j_model, 3* z[i, t1, t2] <= y[i, t2] + sum(y[i, j] for j in 1:t1) + x[i, t1]) for i in 1:n for t1 in 1:T for t2 in 1:T]
    println("Complete")
    print("Constraint 7...")
    [@constraint(j_model, sum(getindex(r_arr, i, t1, t2, c) * z[i, t1, t2] for i in 1:n for t2 in 1:T) >= l[c, t1])  for t1 in 1:T for c in 1:C]
    println("Complete")
    print("Constraint 8...")
    [@constraint(j_model, sum(getindex(r_arr, i, t1, t2, c) * z[i, t1, t2] for i in 1:n for t2 in 1:T) <= u[c, t1])  for t1 in 1:T for c in 1:C]
    println("Complete")
    print("Constraint 9...")
    for (exc, exc_values) in pairs(get(parsed_file, EXCLUSIONS, []))
        i1 = get(interventions_dict, getindex(exc_values, 1), 1)
        i2 = get(interventions_dict, getindex(exc_values, 2), 1)
        season = getindex(exc_values, 3)
        season_vals = get(get(parsed_file, SEASONS, []), season, [])
        if(length(season_vals) > 0 && season_vals[1] isa String)
            season_vals = [parse(Int64, x) for x in season_vals]
        end
        [@constraint(j_model, sum(z[i1, t, j] for j in 1:t) + sum(z[i2, t, k] for k in 1:t) <= 1) for t in season_vals]
    end
    println("Complete")
    print("Constraint 10...")
    [@constraint(j_model, sum(w[s, t1] for s in 1:getindex(St, t1)) >= tau * getindex(St, t1)) for t1 in 1:T]
    println("Complete")
    print("Constraint 11...")
    [@constraint(j_model, risk[s, t1] == sum(getindex(risk_arr, i, t1, t2, s) * z[i, t1, t2] for i in 1:n for t2 in 1:T)) for t1 in 1:T for s in 1:S]
    println("Complete")
    print("Constraint 12...")
    [@constraint(j_model, q[t1] >= risk[s, t1] - M * (1 - w[s, t1])) for s in 1:S for t1 in 1:T]
    println("Complete")
    print("Constraint 15...")
    [@constraint(j_model, e[t1] >= getindex(quant_risk, t1)) for t1 in 1:T]
    println("Complete")
    print("Objective...")
    @objective(j_model, Min, alpha * mean_risk + (1 - alpha) * sum(e[t1] for t1 in 1:T) / T)
    println("Complete")
    vars = Dict("delta" => deltas, "y" => y, "q" => q, "w" => w, "x" => x, "z" => z, "risk" => risk, "e" => e, "quant_risk" => quant_risk, "mean_risk" => mean_risk_t)
    return j_model, vars, Dict("Interventions" => interventions_dict, "Resources" => resources_dict)
end

function writeToFile(filename, y, T, map) 
    open(filename, "w") do f
        for (key, val) in pairs(get(map, "Interventions", []))
            write(f, string(key, " ",argmax(y[val, 1:T]), "\n"))
            println(string(key, " ",argmax(y[val, 1:T])))
        end
    end

end

function main(args)
    println(args)
milp_model, var_dict, mappings = constructMILPfromFile(args[1])
println("Constructed MILP Model")
set_optimizer(milp_model, Cbc.Optimizer)
optimize!(milp_model)
y = value.(get(var_dict, "y", []))
println(y)
writeToFile(args[2], y, 3, mappings)
end
main(ARGS)