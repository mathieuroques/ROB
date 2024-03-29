using JuMP
using CPLEX

C = 0 # nombre minimum d'aires qui doivent être non coupées




# Définition des paramètres
# Instance 1
# M = 10  
# N = 10 
# w1 = 1 # coefficient de pondération w1
# w2 = 5 # coefficient de pondération w2
# l =3 # longueur du côté de chaque parcelle
# g = 1.26157 

# t = [84 68 97 98 64 89 82 71 74 76;
# 87 83 98 75 60 90 78 67 92 94;
# 84 68 70 81 67 61 73 92 86 90;
# 79 62 86 79 73 84 76 98 84 90;
# 62 72 66 72 92 80 71 91 87 70;
# 85 77 63 93 90 94 76 81 99 98;
# 76 63 66 84 94 93 72 92 79 65;
# 76 63 92 69 60 88 79 93 66 73;
# 92 82 77 72 77 81 89 95 80 80;
# 88 89 83 86 69 78 91 64 94 92]

# Instance 2
# N = 5
# M = 5
# w1 = 2
# w2 = 1
# l = 3
# g = 1.26157

# t = [10 10 10 1 10;
# 10 10 1 1 10;
# 10 10 1 10 10;
# 1 10 10 10 10;
# 1 10 10 10 10]


# Définition de la matrice Aij
A = [Tuple{Int, Int}[] for i in 1:M, j in 1:N]

for i in 1:M
    for j in 1:N
        if i > 1 && j > 1 && i < M && j < N
            A[i, j] = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
        elseif i == 1 && j == 1
            A[i, j] = [(i+1, j), (i, j+1)]
        elseif i == 1 && j == N
            A[i, j] = [(i+1, j), (i, j-1)]
        elseif i == M && j == 1
            A[i, j] = [(i-1, j), (i, j+1)]
        elseif i == M && j == N
            A[i, j] = [(i-1, j), (i, j-1)]
        elseif i == 1 && j != 1 && j != N
            A[i, j] = [(i+1, j), (i, j-1), (i, j+1)]
        elseif i == M && j != 1 && j != N
            A[i, j] = [(i-1, j), (i, j-1), (i, j+1)]
        elseif j == 1 && i != 1 && i != M
            A[i, j] = [(i-1, j), (i+1, j), (i, j+1)]
        elseif j == N && i != 1 && i != M
            A[i, j] = [(i-1, j), (i+1, j), (i, j-1)]
        end
    end
end

###########################################################################
println("                                                                       ")
println("------------------------------MODELE (P1)------------------------------")
println("                                                                       ")

time_start_1 = time()
# Définition du modèle pour (P1)
model_1 = Model(CPLEX.Optimizer)
set_optimizer_attribute(model_1, "CPX_PARAM_MIPINTERVAL", 10)


# Déclaration des variables
@variable(model_1, x[i in 1:M, j in 1:N], Bin)  # Variable booléenne xij
@variable(model_1, d[i in 1:M, j in 1:N], Int)  # Variable booléenne xij
for i in 1:M
    for j in 1:N
        set_lower_bound(d[i, j], 0)
    end
end


# Fonction objective
@objective(model_1, Max, w1* sum(t[i, j] * (1 - x[i, j]) for i in 1:M, j in 1:N) + w2 * g * l *sum( (4 * x[i, j] - d[i, j]) for i in 1:M, j in 1:N))

# Contraintes

# #initialisation
# for i in 1:M
#     for j in 1:N
#         if i==1 && j==1
#             @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i+1,j) (i,j+1)]) - 4 * (1 - x[i, j]))
#         end
#         if i==1 && j==N
#             @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i+1,j) (i,j-1)]) - 4 * (1 - x[i, j]))
#         end
#         if i==M && j==1
#             @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i,j+1)]) - 4 * (1 - x[i, j]))
#         end
#         if i==M && j==N
#             @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i,j-1)]) - 4 * (1 - x[i, j]))
#         end
#         if i==1 && j!=1 && j!=N
#             @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i+1,j) (i,j-1) (i,j+1)]) - 4 * (1 - x[i, j]) )
#         end
#         if i==M && j!=1 && j!=N
#             @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i,j-1) (i,j+1)]) - 4 * (1 - x[i, j]) )
#         end
#         if j==1 && i!=1 && i!=M
#             @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i+1,j)  (i,j+1)]) - 4 * (1 - x[i, j]) )
#         end
#         if j==N && i!=1 && i!=M
#             @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i+1,j) (i,j-1) ]) - 4 * (1 - x[i, j]) )
#         end
#         if i!=1 && i!=M && j!=1 && j!=N
#             @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i+1,j) (i,j-1) (i,j+1)]) - 4 * (1 - x[i, j]) )
#         end
#     end
# end

# Contrainte sur le nombre de frontières non lisières
for i in 1:M
    for j in 1:N
        @constraint(model_1, d[i, j] >= sum(x[k, l] for (k,l) in A[i,j]) -  size(A[i,j], 1)* (1 - x[i, j]) )
    end
end

# Contrainte sur le nombre de zone non coupées minimal
# @constraint(model_1, sum(x[i, j] for i in 1:M, j in 1:N) >= C)


# Résolution du modèle
optimize!(model_1)

# Récupération des résultats
solution_x_1 = value.(x)

# Affichage solution
if termination_status(model_1) == MOI.OPTIMAL
    # Afficher la valeur de l'objectif
    objectiveValue_1 = objective_value(model_1)
    println("Valeur de l'objectif : ", objectiveValue_1)

    # println("Aires protégées : ")
    # for i in 1:M
    #     for j in 1:N
    #         print(" ",JuMP.value(x[i,j]))
    #     end
    #     println("")
    #     println("_____________________________________________")
    # end
else
    println("Aucun solution trouvée dans le temps imparti.")
end   
println("Nombre de noeuds explorés : ",JuMP.node_count(model_1))
time_end_1 = time()
println("Temps total d'exécution : ", time_end_1 - time_start_1)

####################################################################
println("                                                                       ")
println("------------------------------MODELE (P2)------------------------------")
println("                                                                       ")

time_start_2 = time()

# Définition du modèle
model_2 = Model(CPLEX.Optimizer)

# Déclaration des variables
@variable(model_2, x[1:M, 1:N], Bin)  # Variable booléenne xij
@variable(model_2, y[1:M, 1:N, 1:M, 1:N], Bin)  # Variable booléenne yijkl

# Fonction objective
@objective(model_2, Max, w1 * sum(t[i, j] * (1 - x[i, j]) for i in 1:M, j in 1:N) +
                      w2 * g * l * sum(x[i, j] * 4 - 
                                      sum(y[i, j, k, l] for (k,l) in A[i,j]) for i in 1:M, j in 1:N))

# Contraintes
for i in 1:M
    for j in 1:N
        for (k,l) in A[i, j]  # Utilisation de la matrice Aij définie précédemment
            @constraint(model_2, y[i, j, k, l] >= x[i, j] + x[k, l] - 1)  # Contrainte yijkl >= xi + xk - 1
            @constraint(model_2, y[i, j, k, l] >= 0)  # Contrainte yijkl >= 0
        end
    end
end

# Contrainte sur le nombre de zone non coupées minimal
# @constraint(model_2, sum(x[i, j] for i in 1:M, j in 1:N) >= C)


# Résolution du modèle
optimize!(model_2)

# Récupération des résultats
solution_x_2 = value.(x)
solution_y_2 = value.(y)

# Affichage solution
if termination_status(model_2) == MOI.OPTIMAL
    # Afficher la valeur de l'objectif
    objectiveValue_2 = objective_value(model_2)
    println("Valeur de l'objectif : ", objectiveValue_2)

    println("Aires protégées : ")
    for i in 1:M
        for j in 1:N
            print(" ",JuMP.value(x[i,j]))
        end
        println("")
        println("_____________________________________________")
    end
else
    println("Aucun solution trouvée dans le temps imparti.")
end   

println("Nombre de noeuds explorés : ",JuMP.node_count(model_2))
time_end_2 = time()
println("Temps total d'exécution : ", time_end_2 - time_start_2)
