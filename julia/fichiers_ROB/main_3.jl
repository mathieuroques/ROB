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

function solve_relax(N::Int, Nm::Int, Nf::Int, C::Int, G::Int, A::Int, T::Int, init::Float64, M::Vector{Matrix{Int64}})

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
        @constraint(model, x[k]<=3)
    end

    # Population size
    @constraint(model, c_taille_pop, sum(x[k] for k in 1:N) == 2*N)
    @constraint(model, c_femme_homme, sum(x[h] for h in 1:Nm) == sum(x[f] for f in Nm+1:N))

    # Solve the model
    optimize!(model)

    println("\n","Paramètres du modèle : ")
    println("h : ", T)
    println("N : ", N)
    println("G : ", G)
    println("A : ", A,"\n")

    # println(termination_status(model) )
    # Check if the problem was solved to optimality
    if termination_status(model) == MOI.OPTIMAL
        # Get the optimal solution
        x_opt = JuMP.value.(x)
        P = JuMP.value.(P)
        # println(x_opt)
        for i in size(P)
            if P[i]!=0
                println(i, " ", P[i])
            end
        end
        # println("Optimal solution: x = $x_opt, y = $y_opt")
        println("Borne inférieure : ", objective_bound(model))
    else
        println("Optimization problem could not be solved.")
        println("MOI termination status: ", termination_status(model))
    end

    time_end = time()
    execution_time = time_end - time_start
    println("Temps d'execution : ", execution_time)
    println("Nombre de noeuds explorés : ",JuMP.node_count(model))

    if termination_status(model) == MOI.OPTIMAL
        return x_opt
    end

end



# Solution réalisable injectée dans le problème initial
function solve_initial_pb(N::Int, Nm::Int, Nf::Int, C::Int, G::Int, A::Int, T::Int, init::Float64, M::Vector{Matrix{Int64}}, x_admissible::Vector{Float64})

    # Create a JuMP model
    time_start = time()
    model = Model(CPLEX.Optimizer)

    # Define variables
    @variable(model, P[1:A,1:G]>=0) # proba que l'allèle disparaisse    

    # Define objective function
    @objective(model, Min, (1/(G*A))*sum(P[j,i] for i in 1:G for j in 1:A))

    # Define constraints
    for i in 1:G
        for j in 1:A
            # Probabilité d'extinction
            @constraint(model, P[j, i] >= prod(0.5^x_admissible[k] for k in 1:N if ((M[k][i*2-1,5] == j)+(M[k][i*2,5] == j))==1)- sum(x_admissible[k] for k in 1:N if ((M[k][i*2-1,5] == j)&&(M[k][i*2,5] == j)) ))
            @constraint(model, P[j,i]<=1)
        end
    end

    # Solve the model
    optimize!(model)

    # println(termination_status(model) )
    # Check if the problem was solved to optimality
    if termination_status(model) == MOI.OPTIMAL
        # Get the optimal solution
        P = JuMP.value.(P)
        somme = 0

        for i in 1:10
            if P[i]!=0
                somme +=1
            end
        end  
        println("Probabilité d'extinction non nulle pour ", somme, " allèles.")
      

        # println("Optimal solution: x = $x_opt, y = $y_opt")
        println("Espérance du nombre d'allèles perdus : ", objective_value(model))
        # println("Objective bound: ", objective_bound(model))
    else
        println("Optimization problem could not be solved.")
        println("MOI termination status: ", termination_status(model))
    end

    time_end = time()
    execution_time = time_end - time_start
    println("Temps d'execution : ", execution_time)


end



function run(N::Int, Nm::Int, Nf::Int, C::Int, G::Int, A::Int, T::Int, init::Float64, M::Vector{Matrix{Int64}})
    println("                                                                     ")
    println("----------------------------MODELE RELACHE---------------------------")
    println("                                                                     ")

    x_admissible = solve_relax(N,Nm,Nf,C,G,A,T,init,M) #solution admissible pour la relaxation du problème à injecter dans le pb initial

    println("                                                                     ")
    println("----------------------------MODELE INITIAL---------------------------")
    println("                                                                     ")

    solve_initial_pb(N,Nm,Nf,C,G,A,T,init,M,x_admissible)
end

for h in [50,100,500]
    run(N,Nm,Nf,C,G,A,h,init,M)
end