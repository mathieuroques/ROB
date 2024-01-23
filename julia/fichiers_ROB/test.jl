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

# Affichage de la matrice Aij
println("Matrice Aij : ")
for i in 1:M
    for j in 1:N
        println("A[", i, ",", j, "] = ", A[i, j])
    end
end

# Récupération de la taille du vecteur A[1,1]
n = size(A[1,1], 1)