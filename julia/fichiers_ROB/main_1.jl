using JuMP
using GLPK

m = 10
n = 10
p = 3 # rare species
q = 3 # common species
alpha = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

# Tableau de probabilités p[k, i, j]
proba = zeros(6, 10, 10)

# Remplissage des valeurs
data = [
    1 4 3 0.4; 1 5 3 0.3; 1 6 2 0.4; 1 8 6 0.3; 1 8 7 0.2; 1 9 5 0.2; 1 9 6 0.4;
    2 2 7 0.2; 2 3 7 0.4; 2 4 7 0.2; 2 5 9 0.4; 2 5 10 0.3; 2 7 2 0.5; 2 9 6 0.2; 2 9 7 0.2;
    3 2 4 0.2; 3 2 5 0.3; 3 2 7 0.4; 3 3 7 0.4; 3 4 7 0.5; 3 5 9 0.2; 3 5 10 0.4; 3 9 2 0.3;
    4 1 2 0.3; 4 1 4 0.3; 4 3 2 0.4; 4 4 3 0.4; 4 5 1 0.3; 4 5 3 0.2; 4 5 5 0.2; 4 6 2 0.2;
    4 6 3 0.4; 4 7 5 0.4; 4 7 9 0.3; 4 8 7 0.2; 4 8 9 0.5; 4 9 7 0.4;
    5 1 10 0.4; 5 2 1 0.3; 5 2 10 0.3; 5 3 5 0.5; 5 6 3 0.4; 5 6 6 0.2; 5 6 7 0.2;
    5 7 5 0.4; 5 7 9 0.5; 5 8 9 0.4; 5 9 2 0.5; 5 9 3 0.2; 5 9 4 0.4; 5 9 5 0.4;
    6 1 3 0.4; 6 1 4 0.4; 6 1 6 0.5; 6 1 7 0.3; 6 1 8 0.3; 6 2 9 0.2; 6 3 5 0.4;
    6 3 9 0.4; 6 3 10 0.2; 6 4 9 0.3; 6 5 5 0.4; 6 8 1 0.5; 6 9 1 0.2; 6 10 3 0.2
]

for i in 1:size(data, 1)
        proba[Int(data[i,1]), Int(data[i,2]), Int(data[i,3])] = data[i,4]
end

c = [6 6 6 4 4 4 4 8 8 8;
     6 6 6 4 4 4 4 8 8 8;
     6 6 6 4 4 4 4 8 8 8;
     5 5 5 3 3 3 3 7 7 7;
     5 5 5 3 3 3 3 7 7 7;
     5 5 5 3 3 3 3 7 7 7;
     5 5 5 3 3 3 3 7 7 7;
     4 4 4 6 6 6 6 5 5 5;
     4 4 4 6 6 6 6 5 5 5;
     4 4 4 6 6 6 6 5 5 5]



function solve(m::Int, n::Int, p::Int, q::Int, alpha::Vector{Float64}, proba::Array{Float64, 3}, c::Matrix{Int64})

    # Create a JuMP model
    model = Model(GLPK.Optimizer)
    # Limite le temps d'exécution à 60 secondes
    set_time_limit_sec(model, 60.0)

    # Define variables
    @variable(model, x[1:m, 1:n], Bin) #zone sélectionnée
    @variable(model, y[1:m, 1:n], Bin) #zone centrale

    # Define objective function
    @objective(model, Min, sum(c[k,j] * x[k, j] for k in 1:m for j in 1:n))

    for i in 1:m
        for j in 1:n
            @constraint(model, y[1,j]==0 )
            @constraint(model, y[m,j]==0 )
            @constraint(model, y[i,1]==0 )
            @constraint(model, y[i,n]==0 )
        end
    end

    # Define constraints
    for i in 2:m-1
        for j in 2:n-1
            @constraint(model, y[i,j]<=sum(x[r,s] for r in i-1:i+1 for s in j-1:j+1)/9)
        end
    end



    # Constraint for rare species
    for k in 1:p
        @constraint(model, sum(log(1-proba[k,i,j])*y[i,j] for i in 1:m for j in 1:n) <= log(1-alpha[k]))
    end
    # Constraint for common species
    for k in p+1:p+q
        @constraint(model, sum(log(1-proba[k,i,j])*x[i,j] for i in 1:m for j in 1:n) <= log(1-alpha[k]))
    end

    # Solve the optimization problem
    optimize!(model)

    println(termination_status(model) )
    # Check if the problem was solved to optimality
    if termination_status(model) == MOI.OPTIMAL
        # Get the optimal solution
        x_opt = JuMP.value.(x)
        y_opt = JuMP.value.(y)
        println(x_opt)
        println(y_opt)
        # println("Optimal solution: x = $x_opt, y = $y_opt")
        println("Objective value: ", objective_value(model))
        println("Objective bound: ", objective_bound(model))
    else
        println("Optimization problem could not be solved.")
        println("MOI termination status: ", termination_status(model))

    end
end

solve(m,n,p,q,alpha,proba,c)

