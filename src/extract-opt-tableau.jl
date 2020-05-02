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


using CPLEX, JuMP, MathProgBase, MathOptInterface, LinearAlgebra
const MOI = MathOptInterface
const MPB = MathProgBase

const ztol = 1e-9
const inf = 1e20


function getTableau(file_name)
	#Load information from file
	model = MPB.LinearQuadraticModel(CplexSolver())
	MPB.loadproblem!(model, file_name)
	A = MPB.getconstrmatrix(model)
	c = MPB.getobj(model)
	m, n = size(A)
	l = MPB.getvarLB(model)
	u = MPB.getvarUB(model)
	vtype = MPB.getvartype(model)
	b = MPB.getconstrLB(model)

	#Construct JuMP Model 
	j_model = Model(with_optimizer(CPLEX.Optimizer))
	@variable(j_model, x[1:n])
	vars = x
	single_var_constr = [@constraint(j_model, 
				x[i] in MOI.Interval(l[i], u[i])) for i in 1:n]
	
	#Considering LP-relaxation, so will not add integrality constraints
	## (The following code adds integrality constraints and hence is commented)	
		#for i in 1:n
		#	if(vtype[i] == :Bin || vtype[i] == :Int) 
		#		set_integer(x[i])
		#	end
		#end

	all_constr = A * x
	lst_of_constr = [@constraint(j_model, all_constr[i] == b[i]) for i in 1:m]
	if(MPB.getsense(model) == :Min)
		@objective(j_model, Min, sum(c[i]*x[i] for i=1:n ))
	else 
		@objective(j_model, Max, sum(c[i]*x[i] for i = 1:n))
	end 

	#Solve LP-Relaxation
	optimize!(j_model)
	xbar = [MOI.get(j_model, MOI.VariablePrimal(), vars[i]) for i in 1:n]
	sbar = b - A*xbar #Slack variables

	#get Basis
	status = CPLEX.get_basis(j_model.moi_backend.optimizer.model.inner)
	status_x = status[1]
	status_s = status[2]
	B_X = [i for i in 1:n if status_x[i] == :Basic]
	B_S = [n + i for i in 1:m if status_s[i] == :Basic]


	##The inverse of a sparse matrix is often a dense matrix; it is best practice to convert
	##to a dense matrix prior to computing the inverse

	A_aug = hcat(Matrix(A), Matrix{Int64}(I, m, m))
	B_aug = vcat(B_X, B_S)
	return inv(A_aug[:, B_aug])*A_aug
end


f = ARGS[1]
@show getTableau(f)
