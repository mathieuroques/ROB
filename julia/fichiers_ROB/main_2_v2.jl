using JuMP
using CPLEX

# Données
c = [7 3 10 10 2 8 6 4 5 5;
     7 7 10 5 2 8 6 3 9 9;
     7 3 4 6 3 2 4 9 7 8;
     6 2 7 6 4 7 5 10 7 8;
     2 4 3 4 9 6 4 9 8 4;
     7 5 2 9 8 9 5 6 10 10;
     5 2 3 7 9 9 4 9 6 3;
     5 2 9 4 2 8 6 9 3 4;
     9 6 5 4 5 6 8 9 6 6;
     8 8 7 7 3 5 8 3 9 9]

n = 10
Amin = 30
Amax = 35
B = 920
# Amin = 20, Amax = 21, B = 520
# Amin = 70, Amax = 75, B = 3500
lambda0 = 20
M = 1+sqrt(2)*n

function distance_entre_centres(i1, j1, i2, j2)
    # Calcule la distance euclidienne entre deux points (i1, j1) et (i2, j2)
    return sqrt((i1 - i2)^2 + (j1 - j2)^2)
end

function solve(n::Int,Amin::Int,Amax::Int,B::Int,lambda::Int,M::Float64,c::Matrix)
    model = Model(CPLEX.Optimizer)
    # Limite le temps d'exécution à 60 secondes
    # set_time_limit_sec(model, 60.0)

    ### Variables
    @variable(model, x[1:n, 1:n], Bin) #zone de coordonnées i,j sélectionnée
    @variable(model, y[1:n, 1:n, 1:n, 1:n], Bin) #zones de coordonnées i,j et k,l toutes les deux sélectionnées
    @variable(model, z[1:n, 1:n]) #linéarisation du min
    @variable(model, d_tilde[1:n, 1:n, 1:n, 1:n]) #distance entre les zones de coordonnées i,j et k;l

    ### Objectif
    @objective(model, Max, sum(z[i,j] for i in 1:n for j in 1:n) -lambda*sum(x[i,j] for i in 1:n for j in 1:n))

    ### Contraintes
    # Aires
    @constraint(model, sum(x[i,j] for i in 1:n for j in 1:n)<=Amax)
    @constraint(model, sum(x[i,j] for i in 1:n for j in 1:n)>=Amin)

    # Budget max
    @constraint(model, sum(10*x[i,j]*c[i,j] for i in 1:n for j in 1:n)<=B)

    # Linéarisations
    for i in 1:n
        for j in 1:n
            for k in 1:n
                for l in 1:n
                    # Linéarisation produit x[i,j]x[k,l]
                    @constraint(model,y[i,j,k,l]<=x[i,j])
                    @constraint(model,y[i,j,k,l]<=x[k,l])
                    @constraint(model,y[i,j,k,l]>=x[i,j]+x[k,l]-1)
                    @constraint(model,y[i,j,k,l]>=0)
                    d_ijkl = distance_entre_centres(i,j,k,l)
                    @constraint(model,d_tilde[i,j,k,l]==d_ijkl*y[i,j,k,l]+M*(1-y[i,j,k,l]))  
                    
                    #Linéarisation du min
                    if (i,j)!=(k,l)
                        @constraint(model,z[i,j]<=d_tilde[i,j,k,l])
                    end
                end
            end
        end
    end


    ### Résoudre le problème
    optimize!(model)

    # Affichage solution
    if termination_status(model) == MOI.OPTIMAL
        # Afficher la valeur de l'objectif
        objectiveValue = objective_value(model)
        println("Valeur de l'objectif : ", objectiveValue)

        println("Aires protégées : ")
        for i in 1:n
            for j in 1:n
                print(" ",JuMP.value(x[i,j]))
            end
            println("")
            println("_____________________________________________")

        end
        return objectiveValue, JuMP.value(lambda), model
    else
        println("Aucun solution trouvée dans le temps imparti.")
    end   
end


solve(n, Amin, Amax, B, lambda0, M, c)
println("Nombre de noeuds explorés : ",JuMP.node_count(model))




function dinkelbach_algorithm(lambda0)
    lambda = lambda0
    
    while true
        # Step 2: Calculate v(λ) and Find xλ such that v(λ) = f(xλ) - λg(xλ)
        v_lambda,new_lambda, model = solve(n, Amin, Amax, B, lambda0, M, c) 
        
        if v_lambda > 0 && v_lambda < 1e-6
            # Step 4: Update λ
            lambda = new_lambda
        else
            # xλ is an optimal solution of (P)
            return model
        end
    end
end

