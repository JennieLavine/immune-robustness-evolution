function compute_fitness_landscape(
    S0, C0,
    optimization_function, fitness_function;
    output_callback = nothing
)
    initial_fitness = fitness_function(S0, C0)

    if !ismissing(output_callback)
        output_callback(0, nothing, nothing, S0, C0, initial_fitness)
    end

    # S matrices
    # dictionary mapping from matrix to fitness value
    fitness_landscape = Dict{BitMatrix, Float64}()

    next_id = 1

    function compute_fitness_landscape_recursive(parent_id, S_parent)
        edges = findall(S_parent)
        if length(edges) <= 16
            return nothing
        end

        for edge in edges
            S = remove_edge(S_parent, edge)
            id = next_id
            next_id += 1

            if haskey(fitness_landscape, S)
#                 println("found repeat matrix")
                continue
            end

            C = optimization_function(S)
            fitness = fitness_function(S, C)

            fitness_landscape[S] = fitness
            if !isnothing(output_callback)
                output_callback(id, parent_id, edge, S, C, fitness)
            end

            compute_fitness_landscape_recursive(id, S)

        end

        nothing
    end

    compute_fitness_landscape_recursive(0, S0)
    (fitness_landscape)
end
