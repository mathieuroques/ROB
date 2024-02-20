using Random
using Distributions

nb_instances = 1

function generate_instance()
    M = rand(500:1000)
    N = M
    w1 = rand(1:5)
    w2 = rand(w1:5w1)
    l =3
    g = 1.26157 
    t = rand(50:100, (N, M))

    return M, N, w1, w2, l, g, t
end

instances = [generate_instance() for i in 1:nb_instances]

for i in 1:nb_instances
    println("Instance ", i)
    M = instances[i][1]
    N = instances[i][2]
    w1 =  instances[i][3]
    w2 = instances[i][4]
    l = instances[i][5]
    g = instances[i][6]
    t = instances[i][7]

    # Ecrire dans le fichier    
    file_name = "Instance_Projet_4"*"-"*string(N)*"-"*string(w1)*"-"*string(w2)*".txt"
    file = open(file_name,"w")

    write(file, "N=")
    write(file,string(N))
    write(file, "\n")
    write(file, "M=")
    write(file,string(M))
    write(file, "\n")
    write(file, "w1=")
    write(file,string(w1))
    write(file, "\n")
    write(file, "w2=")
    write(file,string(w2))
    write(file, "\n")
    write(file, "l=")
    write(file,string(l))
    write(file, "\n")
    write(file, "g=")
    write(file,string(g))
    write(file, "\n")
    write(file, "t=")
    write(file,string(t))
    write(file, "\n")

    close(file)
end

