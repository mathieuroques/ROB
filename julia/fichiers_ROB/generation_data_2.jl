using Random
using Distributions

nb_instances = 1

function generate_instance()
    lambda = 20
    n = rand(20:1:20)
    Amin = floor(Int,rand(0.1*(n^2):1:0.9*(n^2)))
    Amax = floor(Int,Amin*1.1)
    B = floor(Int,rand((n^2)*5:1:(n^2)*9))
    couts = rand(1:10, n, n)

    return lambda, n, Amin, Amax, B, couts
end

instances = [generate_instance() for i in 1:nb_instances]

for i in 1:nb_instances
    println("Instance ", i)
    lambda = instances[i][1]
    n = instances[i][2]
    Amin =  instances[i][3]
    Amax = instances[i][4]
    B = instances[i][5]
    couts = instances[i][6]

    # Ecrire dans le fichier    
    file_name = "Instance_Projet_4"*"-"*string(n)*"-"*string(Amin)*"-"*string(Amax)*string(B)*".txt"
    file = open(file_name,"w")

    write(file, "lambda=")
    write(file,string(lambda))
    write(file, "\n")
    write(file, "n=")
    write(file,string(n))
    write(file, "\n")
    write(file, "Amin=")
    write(file,string(Amin))
    write(file, "\n")
    write(file, "Amax=")
    write(file,string(Amax))
    write(file, "\n")
    write(file, "B=")
    write(file,string(B))
    write(file, "\n")
    write(file, "couts=")
    write(file,string(couts))
    write(file, "\n")

    close(file)
end

