using JuMP
using GLPK

# Définition des paramètres
M = 10  # Remplacez m par la valeur appropriée
N = 10  # Remplacez n par la valeur appropriée
w1 = 1 # Remplacez par le coefficient de pondération w1
w2 = 5 # Remplacez par le coefficient de pondération w2
l =3 # Remplacez par la longueur du côté de chaque parcelle
g = 1.26157 # 
P=1;

t = [1 84 68 97 98 64 89 82 71 74 76;
2 87 83 98 75 60 90 78 67 92 94;
3 84 68 70 81 67 61 73 92 86 90;
4 79 62 86 79 73 84 76 98 84 90;
5 62 72 66 72 92 80 71 91 87 70;
6 85 77 63 93 90 94 76 81 99 98;
7 76 63 66 84 94 93 72 92 79 65;
8 76 63 92 69 60 88 79 93 66 73;
9 92 82 77 72 77 81 89 95 80 80;
10 88 89 83 86 69 78 91 64 94 92]

# Définition du modèle
model = Model(GLPK.Optimizer)

# Déclaration des variables
@variable(model, x[i in 1:M, j in 1:N], Bin)  # Variable booléenne xij
@variable(model, d[i in 1:M, j in 1:N], Int)  # Variable booléenne xij

# Fonction objective
@objective(model, Max, sum(t[i, j] * (1 - x[i, j]) for i in 1:M, j in 1:N) + w2 * g * l *sum( (4 * x[i, j] - d[i, j]) for i in 1:M, j in 1:N))

# Contraintes

#initialisation
for i in 1:M
    for j in 1:N
        if i==1 && j==1
            @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i+1,j) (i,j+1)]) - 2 * (1 - x[i, j]))
        end
        if i==1 && j==N
            @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i+1,j) (i,j-1)]) - 2 * (1 - x[i, j]))
        end
        if i==M && j==1
            @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i,j+1)]) - 2 * (1 - x[i, j]))
        end
        if i==M && j==N
            @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i,j-1)]) - 2 * (1 - x[i, j]))
        end
        if i==1 && j!=1 && j!=N
            @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i+1,j) (i,j-1) (i,j+1)]) - 3 * (1 - x[i, j]) )
        end
        if i==M && j!=1 && j!=N
            @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i,j-1) (i,j+1)]) - 3 * (1 - x[i, j]) )
        end
        if j==1 && i!=1 && i!=M
            @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i+1,j)  (i,j+1)]) - 3 * (1 - x[i, j]) )
        end
        if j==N && i!=1 && i!=M
            @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i+1,j) (i,j-1) ]) - 3 * (1 - x[i, j]) )
        end
        if i!=1 && i!=M && j!=1 && j!=N
            @constraint(model, d[i, j] >= sum(x[k, l] for (k,l) in [(i-1,j) (i+1,j) (i,j-1) (i,j+1)]) - 4 * (1 - x[i, j]) )
        end
    end
end

# Résolution du modèle
optimize!(model)

# Récupération des résultats
solution_x = value.(x)

# Affichage des résultats ou autre traitement nécessaire
# Affichage solution
if termination_status(model) == MOI.OPTIMAL
    # Afficher la valeur de l'objectif
    objectiveValue = objective_value(model)
    println("Valeur de l'objectif : ", objectiveValue)

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
