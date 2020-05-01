######################################################################
# The problem is assumed to be in SEF. The tableau that is returned is
# for the original instance with slack variables (added by the solver). 
# That is, given an instance in the form 
# 	A x = b
# 	l <= x <= u
# The optimal tableau for the following augmented instance is returned:
# 	[A | I] [x // s] <= b
# 	l <= x <= u
# 	s >= 0
# Command line arguments: FILE_NAME
######################################################################


using JuMP, MathOptInterface
const MOI = MathOptInterface

const ztol = 1e-9
const inf = 1e20


function getP2Model(interventions, deltas, T)
	# Pass in a parameters Object! 

	"""
	Interventions: List<Int>

	"""
	#Constants
	n = length(interventions)

	#Construct JuMP Model 
	j_model = Model()
	t[1:n] = @variable(j_model, t[1:n], Int)
	d[1:n] = @variable(j_model, d[1:n], Int)
	y[1:n, 1:T] = @variable(j_model, y[1:n, 1:T], Bin)
	a[1:n, 1:T] = @variable(j_model, a[1:n, 1:T], Bin)
	c[1:n, 1:T, 1:T] = @variable(j_model, c[1:n, 1:T, 1:T], Bin)
	k[1:S, 1:T, 1:n] =  @variable(j_model, k[1:S, 1:T, 1:n]) #risk
	risk[1:S, 1:T] = @variable(j_model, risk[1:S, 1:T])
	r[1:C, 1:T, 1:n] = @variable(j_model, r[1:C, 1:T, 1:n])

	q[1:T] = @variable(j_model, q[1:T])	

	constr1 = [@constraint(j_model, t[i] >= time*y[i,time]) for i in 1:n for time in 1:T]
	constr2 = [@constraint(j_model, t[i] <= time + T*(1-y[i, time])) for i in 1:n for time in 1:T]	
	constr3 = [@constraint(j_model, t[i] <= time + (1- a[i, time])*T) for i in 1:n for time in 1:T]
	constr4 = [@constraint(j_model, d[i] >= time*a[i, time]) for i in 1:n for time in 1:T]
	constr5 = [@constraint(j_model, c[i, time1, time2] >= y[i, time1] + a[i, time2] - 1) for i in 1:n for time1 in 1:T for time2 in 1:T]
	constr6 = [@constraint(j_model, c[i, time1, time2] <= y[i, time1]) for i in 1:n for time1 in 1:T for time2 in 1:T]
	constr7 = [@constraint(j_model, c[i, time1, time2] <= a[i, time2]) for i in 1:n for time1 in 1:T for time2 in 1:T]
	#constr8 = [@constraint(j_model, d[i] <= time + deltas[i, time] + T * (1 - y[i, time])) for i in 1:n for time in 1:T]
	#constr9 = [@constraint(j_model, d[i] >= time + deltas[i, time] - T*(1 - y[i, time])) for i in 1:n for time in 1:T]
	constr10 = [@constraint(j_model, k[s, t1, i] <= risk[s, t1, i, t2] + M * (1 - c[i, t2, t1])) for i in 1:n for t1 in 1:T for t2 in 1:T]
	constr11 = [@constraint(j_model, k[s, t1, i] >= risk[s, t1, i, t2] - M * (1 - c[i, t2, t1])) for i in 1:n for t1 in 1:T for t2 in 1:T]
	constr12 = [@constraint(j_model, k[s, t1, i] <= M * c[i, t2, t1]) for i in 1:n for t1 in 1:T for t2 in 1:T]
	constr13 = [@constraint(j_model, k[s, t1, i] >= 0) for i in 1:n for t1 in 1:T]
	constr14 = [@constraint(j_model, risk[s, t1] == sum(k[s, t1, i] for i in 1:n))]
	constr15 = [@constraint(j_model, r[c1, t1, i] <= res[c1, t1, i, t2] + u[c1, t1]*(1-c[i, t2, t1])) for c1 in 1:C for t1 in 1:T for t2 in 1:T for i in 1:n]
	constr16 = [@constraint(j_model, r[c1, t1, i] >= res[c1, t1, i, t2] - u[c1, t1]*(1-c[i, t2, t1])) for c1 in 1:C for t1 in 1:T for t2 in 1:T for i in 1:n]
	constr17 = [@constraint(j_model, r[c1, t1, i] <= u[c1, t1] * c[i, t2, t1]) for c1 in 1:C for t1 in 1:T for t2 in 1:T for i in 1:n]
	constr18 = [@constraint(j_model, r[c1, t1, i] >= 0) for c1 in 1:C for t1 in 1:T for i in 1:n]
	constr19 = [@constraint(j_model, q[t1] >= )]
	constr21 = [@constraint(j_model, d[i] <= T+ 1) for i in 1:n] 
	constr23 = [@constraint(j_model, sum(y[i, time] for time in 1:T) == 1) for i in 1:n]
	
	print(j_model)
	return j_model
end
getP2Model([1, 2], 1, 3)