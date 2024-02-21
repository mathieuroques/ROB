using Random
using Distributions

nb_instances = 1

function generate_instance()
    N = rand(2:2:100)
    Nm = floor(Int,N/2)
    Nf = floor(Int,N/2)
    C = 1
    G = rand(N:2N)
    A = rand(2:10)
    T = 50
    init=0.001;
    M = Array{Array{Int64,2},1}(undef, N)
    for i in 1:N
        matrice = zeros(Int64, 2G, 5)
        for j in 1:2G
            matrice[j,1]=i
            matrice[j, 2] = 1
            matrice[j, 3] = 1 + ((j - 1) ÷ 2) % N
            matrice[j, 4] = 2 -  (j % 2)
            matrice[j, 5] = rand(1:A)
        end
        M[i]=matrice
    end

    return N, Nm, Nf, C, G,A, T, init, M
end

instances = [generate_instance() for i in 1:nb_instances]

for i in 1:nb_instances
    println("Instance ", i)
    N = instances[i][1]
    Nm = instances[i][2]
    Nf =  instances[i][3]
    C = instances[i][4]
    G = instances[i][5]
    A = instances[i][6]
    T = instances[i][7]    
    init = instances[i][8]
    M = instances[i][9]

    # Ecrire dans le fichier    
    file_name = "Instance_Projet_3"*"-"*string(N)*"-"*string(G)*"-"*string(A)*".txt"
    file = open(file_name,"w")

    write(file, "N=")
    write(file,string(N))
    write(file, "\n")
    write(file, "Nm=")
    write(file,string(Nm))
    write(file, "\n")
    write(file, "Nf=")
    write(file,string(Nf))
    write(file, "\n")
    write(file, "C=")
    write(file,string(C))
    write(file, "\n")
    write(file, "G=")
    write(file,string(G))
    write(file, "\n")
    write(file, "A=")
    write(file,string(A))
    write(file, "\n")
    write(file, "T=")
    write(file,string(T))
    write(file, "\n")
    write(file, "init=")
    write(file,string(init))
    write(file, "\n")
    write(file, "M=")
    write(file,string(M))
    write(file, "\n")

    close(file)
end

