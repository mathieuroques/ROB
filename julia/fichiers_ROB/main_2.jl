using JuMP
using CPLEX
using Colors
using Plots

start_time = time()
m = Model(CPLEX.Optimizer)

lambda=20
n=20
Amin=193
Amax=212
B=3089
couts=[4 3 6 10 7 1 1 2 1 3 8 9 7 1 8 3 9 2 6 3; 7 5 4 2 2 3 2 8 5 9 2 2 5 6 6 5 1 5 8 4; 4 1 9 7 6 10 9 4 3 9 6 1 3 5 9 7 9 9 3 3; 7 2 1 10 9 4 5 1 4 9 2 8 8 4 8 6 5 10 6 9; 8 4 7 10 5 8 10 3 10 10 2 6 4 3 4 1 9 4 9 6; 8 10 1 7 4 2 2 3 10 3 2 8 2 3 8 8 1 9 9 1; 4 4 9 8 3 5 10 10 4 3 4 2 9 3 9 3 2 4 1 4; 7 4 8 5 3 10 8 4 9 3 7 6 6 8 3 5 3 5 8 8; 9 7 10 6 4 8 8 7 10 7 6 2 4 5 2 5 2 9 6 1; 1 10 10 5 9 4 10 6 6 3 8 8 4 7 6 9 10 4 8 6; 8 3 9 7 1 7 9 4 7 1 6 10 3 3 9 6 6 9 4 2; 1 4 9 1 2 8 9 1 6 1 2 2 8 10 10 5 9 9 2 8; 5 6 9 2 5 8 10 8 6 10 4 10 2 8 8 2 1 1 6 3; 3 9 5 4 8 7 1 6 8 7 5 4 3 1 2 10 4 2 2 1; 2 1 1 2 4 7 1 9 1 3 10 7 8 7 6 6 4 9 2 4; 1 8 3 4 5 8 4 10 5 4 2 6 10 2 4 10 7 9 3 10; 10 6 8 9 6 10 3 6 10 10 3 8 5 7 8 4 10 4 7 9; 2 5 7 2 4 2 3 8 1 5 4 9 6 9 9 8 7 3 10 9; 9 9 6 5 3 1 2 8 4 3 8 9 4 4 9 6 4 9 2 10; 6 2 9 2 5 4 5 7 6 9 1 9 3 1 9 9 8 8 7 6]


# lambda = 20
# n = 10
# Amin = 30
# Amax = 35
# B = 920
# Amin = 20
# Amax = 21
# B = 520
# Amin = 70
# Amax = 75
# B = 3500

# couts = Matrix{Int64}([
#     7 3 10 10 2 8 6 4 5 5;
#     7 7 10 5 2 8 6 3 9 9;
#     7 3 4 6 3 2 4 9 7 8;
#     6 2 7 6 4 7 5 10 7 8;
#     2 4 3 4 9 6 4 9 8 4;
#     7 5 2 9 8 9 5 6 10 10;
#     5 2 3 7 9 9 4 9 6 3;
#     5 2 9 4 2 8 6 9 3 4;
#     9 6 5 4 5 6 8 9 6 6;
#     8 8 7 7 3 5 8 3 9 9
# ]
# )

couts = 10*couts

@variable(m, x[1:n, 1:n], Bin)
@variable(m, y[1:n,1:n,1:n,1:n], Bin)

@constraint(m, Amin <= sum(x[i,j] for i in 1:n for j in 1:n))
@constraint(m, sum(x[i,j] for i in 1:n for j in 1:n) <= Amax)

@constraint(m, sum(couts[i,j]*x[i,j] for i in 1:n for j in 1:n) <= B)

@constraint(m, [i in 1:n, j in 1:n], sum( y[i,j,k,l] for k in 1:n for l in 1:n if (i!=k || j!=l) ) == x[i,j])
@constraint(m, [i in 1:n, j in 1:n, k in 1:n, l in 1:n], y[i,j,k,l] <= x[k,l])


v_lambda = 10

function distance(y)
    return sum(sqrt((i-k)^2+(j-l)^2)*y[i,j,k,l] for i in 1:n for j in 1:n for k in 1:n for l in 1:n)
end

function somme(x)
    return sum(x[i,j] for i in 1:n for j in 1:n)
end

global compteur=0
while  (v_lambda > 0.001 || v_lambda < - 0.001) 
    global compteur +=1
    @objective(m, Min, sum(sqrt((i-k)^2+(j-l)^2)*y[i,j,k,l] for i in 1:n for j in 1:n for k in 1:n for l in 1:n if ((i!=k) || (j!=l))) - lambda*sum(x[i,j] for i in 1:n for j in 1:n))
    optimize!(m)
    global feasible_solution_found = (primal_status(m) == MOI.FEASIBLE_POINT)
    global isOpt = termination_status(m) == MOI.OPTIMAL
    if feasible_solution_found
        global v_lambda = JuMP.objective_value(m)
        global vy = value.(y)
        global vx = value.(x)

        if v_lambda > 0 
            global lambda = distance(vy)/somme(vx)
            # println("vlambda ", v_lambda)
            # println("Valeur fonction objective : ", JuMP.objective_value(m))
            # println("valeur distance moyenne : ", lambda)
            # println("v_lambda supérieur à 0")
            continue
        else 
            # println("Valeurs x :")
            # for i in 1:size(vx, 1)
            #     for j in 1:size(vx, 1)
            #         print(Int64(vx[i,j]>0.9), " ")
            #     end
            #     println()
            # end
            println("valeur distance moyenne : ")
            global vy = value.(y)
            global vx = value.(x)
            println(sum(sqrt((i-k)^2+(j-l)^2)*vy[i,j,k,l] for i in 1:n for j in 1:n for k in 1:n for l in 1:n)/sum(vx[i,j] for i in 1:n for j in 1:n))

            break
        end
    end
end

println("Nombre de noeuds explorés : ",JuMP.node_count(m))

println("Itérations : ", compteur)

end_time = time()
println("Execution time : ", end_time-start_time)

println("Paramètres")
println("n : ", n)
println("Amin : ", Amin)
println("Amax : ", Amax)
println("B : ", B)