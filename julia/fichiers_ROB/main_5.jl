using JuMP
using CPLEX


NbParcelles=40;

T=10;
SURF=1;
lmax=2;
amax=2;

C=[1,2]; # 1=riz, 2=haricot
Demande=[[1200,0,1200,0,1200,0,1200,0,1200,0],[0,400,0,400,0,400,0,400,0,400]];


gamma = [2,1]

sommets = [
    [2, 0, 0],
    [1, 0, 0],
    [2, 1, 1],
    [2, 2, 1],
    [1, 1, 1],
    [1, 2, 1],
    [2, 1, 2],
    [2, 2, 2],
    [1, 1, 2],
    [1, 2, 2]
]
NbSommets = size(sommets, 1)

arcs = [
    [[2, 0, 0], [2, 0, 0], 0],
    [[2, 0, 0], [2, 1, 1], 120],
    [[2, 0, 0], [2, 1, 2], 90],
    [[1, 0, 0], [2, 0, 0], 0],
    [[1, 0, 0], [1, 1, 1], 72],
    [[1, 0, 0], [1, 1, 2], 54],
    [[2, 1, 1], [1, 0, 0], 0],
    [[2, 1, 1], [2, 2, 2], 65],
    [[2, 2, 1], [1, 0, 0], 0],
    [[1, 1, 1], [1, 0, 0], 0],
    [[1, 1, 1], [1, 2, 2], 39],
    [[1, 2, 1], [1, 0, 0], 0],
    [[2, 1, 2], [1, 0, 0], 0],
    [[2, 1, 2], [2, 2, 1], 90],
    [[2, 2, 2], [1, 0, 0], 0],
    [[1, 1, 2], [1, 0, 0], 0],
    [[1, 1, 2], [1, 2, 1], 54],
    [[1, 2, 2], [1, 0, 0], 0]
]
NbArc = size(arcs, 1)


function solve(NbParcelles::Int, NbArc::Int, T::Int, SURF::Int, lmax::Int, amax::Int, C::Array{Int,1}, Demande::Array{Array{Int,1},1}, gamma::Array{Int,1})        

    time_start = time()
    model = Model(CPLEX.Optimizer)

    @variable(model, x[1:NbParcelles, 1:T, 1:NbArc] >= 0, Bin)
    @objective(model, Min, sum(sum(x[p, 1, i] for i in 1:NbArc if arcs[i][1] == sommets[1]) for p in 1:NbParcelles))

    # Contrainte 0
    for p in 1:NbParcelles
        for arc in 1:NbArc
            if arcs[arc][1] != sommets[1]
                @constraint(model, x[p, 1, arc] == 0)
            end
        end
    end


    # Contrainte 1
    for p in 1:NbParcelles
        @constraint(model, sum(x[p, 1, arc] for arc in 1:NbArc) <= 1)
    end

    # Contrainte 2
    for t in 1:T
        for j in 1:length(C)
                @constraint(model, sum(x[p, t, u]* arcs[u][3] for p in 1:NbParcelles for u in 1:NbArc) >= Demande[j][t])
        end
    end

    # Contrainte 3
    for p in 1:NbParcelles
        for t in 2:T
            for s in sommets
                a = sum(x[p, t-1, arc] for arc in 1:NbArc if arcs[arc][2]==s)
                b = sum(x[p, t, arc] for arc in 1:NbArc if arcs[arc][1]==s)
                @constraint(model,  a == b)
            end
        end
    end


    optimize!(model)
    solution = value.(x)


    NbParcellesSelected = 0
    for p in 1:NbParcelles
        NbParcellesSelected += (sum(solution[p,:,:])>0)
    end

    
    println(termination_status(model))

    # Check if the problem was solved to optimality
    if termination_status(model) == MOI.OPTIMAL
        println("Objective value: ", objective_value(model))
        println("Number of selected parcels: ", NbParcellesSelected)
    else
        println("Optimization problem could not be solved.")
        println("MOI termination status: ", termination_status(model))
    end

    time_end = time()
    execution_time = time_end - time_start
    println("Temps d'execution : ", execution_time)

    return value.(x)

end

solution = solve(NbParcelles, NbArc, T, SURF, lmax, amax, C, Demande, gamma)

print("")
