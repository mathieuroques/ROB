using JuMP
using CPLEX

N=8; #nombre d'individus dans la population
Nm=4; #nombre males
Nf=4; #nombre femmes
C=1; #nombre de paires de chromosomes
G=5;  #nombre de locus par chromosomes
A=2; #nombre d'allèles par gène 
T=50; #h
init=0.001; #theta 1

M = [
    [1 1 1 1 2;
     1 1 1 2 1;
     1 1 2 1 1;
     1 1 2 2 1;
     1 1 3 1 2;
     1 1 3 2 1;
     1 1 4 1 2;
     1 1 4 2 1;
     1 1 5 1 1;
     1 1 5 2 2],

    [2 1 1 1 2;
     2 1 1 2 2;
     2 1 2 1 2;
     2 1 2 2 1;
     2 1 3 1 1;
     2 1 3 2 2;
     2 1 4 1 1;
     2 1 4 2 2;
     2 1 5 1 1;
     2 1 5 2 1],

    [3 1 1 1 2;
     3 1 1 2 2;
     3 1 2 1 2;
     3 1 2 2 1;
     3 1 3 1 2;
     3 1 3 2 1;
     3 1 4 1 1;
     3 1 4 2 2;
     3 1 5 1 1;
     3 1 5 2 2],

    [4 1 1 1 2;
     4 1 1 2 2;
     4 1 2 1 1;
     4 1 2 2 1;
     4 1 3 1 1;
     4 1 3 2 2;
     4 1 4 1 2;
     4 1 4 2 2;
     4 1 5 1 2;
     4 1 5 2 2],

    [5 1 1 1 2;
     5 1 1 2 1;
     5 1 2 1 1;
     5 1 2 2 1;
     5 1 3 1 1;
     5 1 3 2 1;
     5 1 4 1 2;
     5 1 4 2 2;
     5 1 5 1 1;
     5 1 5 2 2],

    [6 1 1 1 1;
     6 1 1 2 2;
     6 1 2 1 1;
     6 1 2 2 1;
     6 1 3 1 2;
     6 1 3 2 2;
     6 1 4 1 2;
     6 1 4 2 1;
     6 1 5 1 2;
     6 1 5 2 1],

    [7 1 1 1 1;
     7 1 1 2 1;
     7 1 2 1 1;
     7 1 2 2 1;
     7 1 3 1 2;
     7 1 3 2 2;
     7 1 4 1 1;
     7 1 4 2 1;
     7 1 5 1 1;
     7 1 5 2 1],

    [8 1 1 1 2;
     8 1 1 2 2;
     8 1 2 1 1;
     8 1 2 2 1;
     8 1 3 1 2;
     8 1 3 2 2;
     8 1 4 1 2;
     8 1 4 2 1;
     8 1 5 1 2;
     8 1 5 2 1]
]

function solve(N::Int, Nm::Int, Nf::Int, C::Int, G::Int, A::Int, T::Int, init::Float64, M::Vector{Matrix{Int64}})

    # Create a JuMP model
    time_start = time()
    model = Model(CPLEX.Optimizer)

    # Define variables
    @variable(model, x[1:N]>=0, Int) # nombre de progénitures
    @variable(model, P[1:A,1:G]>=0) # proba que l'allèle disparaisse
    @variable(model, t[1:A,1:G]) # variable de linéarisation
    

    # Define objective function
    @objective(model, Min, (1/(G*A))*sum(P[j,i] for i in 1:G for j in 1:A))

    # Define constraints
    for i in 1:G
        for j in 1:A
            # Probabilité d'extinction
            @constraint(model, P[j, i] >= t[j, i] - sum(x[k] for k in 1:N if ((M[k][i*2-1,5] == j)&&(M[k][i*2,5] == j)) ))
            @constraint(model, P[j,i]<=1)
            # Relaxation
            for r in 1:T
                theta_r = init^((T-r)/(T-1))
                @constraint(model, log(theta_r) + (1/theta_r)*(t[j,i]-theta_r) >= sum(x[k] * log(0.5) for k in 1:N if ((M[k][i*2-1,5] == j)+(M[k][i*2,5] == j))==1))
            end
        end
    end

    # Progéniture max / individu
    for k in 1:N
        @constraint(model, x[k]==2)
    end

    @constraint(model, sum(x[k] for k in 1:N) == N)
    @constraint(model, sum(x[h] for h in 1:Nm) == sum(x[f] for f in Nm+1:N))


    # Solve the model
    optimize!(model)

    # Get the optimal solution
    solution = Dict()
    for k in 1:N
        solution[k] = value(x[k])
    end

    time_end = time()
    execution_time = time_end - time_start
    println("Temps d'execution : ", execution_time)

    return solution


    # println(termination_status(model) )
    # # Check if the problem was solved to optimality
    # if termination_status(model) == MOI.OPTIMAL
    #     # Get the optimal solution
    #     x_opt = JuMP.value.(x)
    #     y_opt = JuMP.value.(y)
    #     println(x_opt)
    #     println(y_opt)
    #     # println("Optimal solution: x = $x_opt, y = $y_opt")
    #     println("Objective value: ", objective_value(model))
    #     println("Objective bound: ", objective_bound(model))
    # else
    #     println("Optimization problem could not be solved.")
    #     println("MOI termination status: ", termination_status(model))

    # end


end

solve(N,Nm,Nf,C,G,A,T,init,M)

