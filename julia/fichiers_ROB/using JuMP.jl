using JuMP
using GLPK

N=8;
Nm=4;
Nf=4;
C=1;
G=5;
A=2;
T=50;
init=0.001;

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

function solve(N::Int, Nm::Int, Nf::Int, C::Int, G::Int, A::Int, T::Int, init::Float64, M::Vector{Array{Int64}})

    # Create a JuMP model
    model = Model(GLPK.Optimizer)
    # Limite le temps d'exécution à 60 secondes
    set_time_limit_sec(model, 60.0)

    # Define variables
    @variable(model, x[1:N], Int) # nombre de progénitures
    @variable(model, P[1:A,1:G]>=0) # proba que l'allèle disparaisse
    @variable(model, t[1:A,1:G]) # variable de linéarisation
    

    # Define objective function
    @objective(model, Min, (1/G*A)*sum(P[j,i] for i in 1:G for j in 1:A))

    # Define constraints
    for i in 1:G
        for j in 1:A
            @constraint(model, P[j, i] >= t[j, i] - sum(x[k] for k in 1:N if ((M[k][i*2-1,5] == j)&&(M[k][i*2,5] == j)) ))
            @constraint(model, P[j,i]<=1)
            for r in 1:T
                @constraint(model, log(init*((T-r)/(T-1))) + (1/(init*((T-r)/(T-1))))*(t[j,i]-init*((T-r)/(T-1))) >= sum(x[k] * log(0.5) for k in 1:N if ((M[k][i*2-1,5] == j)+(M[k][i*2,5] == j))==1))
            end
            
        end
    end

#     Jeanne, douce étoile dans le ciel, Ton nom résonne comme une merveille. Ton sourire illumine nos jours, Et ton amour remplit nos cœurs.

# Jeanne, tu es une fleur épanouie, Ta présence est une mélodie. Ton regard brille d'une lueur, Qui éclaire notre chemin avec douceur.

# Jeanne, tu es une force tranquille, Une âme noble et subtile. Ta bienveillance est un trésor, Qui nous guide vers un avenir meilleur.

# Jeanne, tu es une inspiration, Un rayon de soleil dans notre horizon. Que ta vie soit remplie de bonheur, Et que chaque jour soit une douceur.

    @constraint(model, sum(x[k] for k in 1:N) == N)
    @constraint(model, sum(x[h] for h in 1:Nm) == sum(x[f] for f in Nm+1:N))


    # Solve the model
    optimize!(model)

    # Get the optimal solution
    solution = Dict()
    for k in 1:N
        solution[k] = value(x[k])
    end


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

