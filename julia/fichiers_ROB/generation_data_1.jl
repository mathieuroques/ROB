using Random
using Distributions

nb_instances = 1

function generate_instance()
    m = rand(1:1000)
    n = m
    p = rand(1:floor(Int, n/2)+1)
    q = rand(1:floor(Int, n/2)+1)
    alpha = rand(0.1:0.01:0.9, p+q)  
    data = rand(Float64, (rand(2n:12n), 3))  # Les trois premières colonnes
    data[:, 1] .= rand(1:p+q, size(data, 1))
    data[:, 2] .= rand(1:n, size(data, 1))
    data[:, 3] .= rand(1:m, size(data, 1))
    data = hcat(data, rand(size(data, 1)))  # Ajouter une colonne pour les valeurs non entières


    c = rand(1:8, (n, m))

    return m, n, p, q, alpha, data, c
end

instances = [generate_instance() for i in 1:nb_instances]

for i in 1:nb_instances
    println("Instance ", i)
    m = instances[i][1]
    n = instances[i][2]
    p =  instances[i][3]
    q = instances[i][4]
    alpha = instances[i][5]
    data = instances[i][6]
    c = instances[i][7]

    # Ecrire dans le fichier    
    file_name = "Instance"*"-"*string(n)*"-"*string(p)*"-"*string(q)*"-"*string(floor(size(data)[1]/(p+q)))*".txt"
    file = open(file_name,"w")

    write(file, "n=")
    write(file,string(n))
    write(file, "\n")
    write(file, "m=")
    write(file,string(m))
    write(file, "\n")
    write(file, "p=")
    write(file,string(p))
    write(file, "\n")
    write(file, "q=")
    write(file,string(q))
    write(file, "\n")
    write(file, "alpha=")
    write(file,string(alpha))
    write(file, "\n")
    write(file, "data=")
    write(file,string(data))
    write(file, "\n")
    write(file, "c=")
    write(file,string(c))
    write(file, "\n")

    close(file)
end

