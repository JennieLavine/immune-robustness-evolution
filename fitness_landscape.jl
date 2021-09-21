#!/usr/bin/env julia

# This script generates a receptor-effector matrix, enumerates the neutral
# network, and writes it all out to a SQLite database.

include("julia/robustness-genomesize_play.jl")
include("redundancy_metrics.jl")

using SQLite: DB, Stmt
import SQLite.DBInterface.execute

const N_RECEPTORS = 8
const N_EFFECTORS = 16
const MAX_EDGES_PER_RECEPTOR = 5
const MAX_WEIGHT_PER_RECEPTOR_INIT = 1.0
const WEIGHT_PRIOR_SD_REGRESSION = 100.0
const PROP_CAPACITY = 0.9
const MIN_EDGES = 25
const N_NETWORKS = 1000
const MAX_ROW_SUM = 5
const MAX_COMBO_SIZE = 4

const RECORD_EVERYTHING = false

const N_DISRUPTIONS_HSC = 4

const V0 = 10
const I0 = 100
const k0 = 10e-3
const μ = 1.5
const β = 1.0
const ϕ = 10e3
const α = 0
const tspan = (0.0, 100)

Random.seed!(1234)
s_mats = generate_networks(
    N_RECEPTORS, N_EFFECTORS, N_NETWORKS, MAX_ROW_SUM
    )

high_redun_ind = [ 2  28  94 364 437 445 510 646 802]
low_redun_ind = [19 223 438 679]



function run_sims(mat_ind)
    fileName = string("output_mat", mat_ind, ".sqlite")

    # Generate an initial receptor-effector network as a matrix.
    S0 = s_mats[:,:,mat_ind]

    # Generate weights as a matrix
    C0 = generate_weights_with_max_rowsum(S0, MAX_WEIGHT_PER_RECEPTOR_INIT, PROP_CAPACITY)
    println("Initial weights C0:")
    display(C0)
    print("\n\n")

    # Compute the initial phenotype
    println("Initial phenotype:")
    initial_phenotype = compute_phenotype(C0)
    display(initial_phenotype)
    print("\n\n")

    # # Define an optimization function for the search
    # function optimize(S)
    #     C, w_mean, w_cov = optimize_weights_regression_analytical(
    #         initial_phenotype, S,
    #         WEIGHT_PRIOR_SD_REGRESSION;
    #         apply_prior_per_receptor = true
    #     )
    #     max.(0.0, C)
    # end

    # Define an optimization function for the search
    function optimize(S)
        C, obj, w = optimize_weights_real_whitacre(
            initial_phenotype, S,
            MAX_WEIGHT_PER_RECEPTOR_INIT;
            n_initial_conditions = 100, return_all=true
        )
        max.(0.0, C)
    end


    # Define a fitness function for the search
    function compute_fitness(S, C)
        phenotype = min.(initial_phenotype, compute_phenotype(C))
        k = k0 * sum(phenotype) # Or whatever
        n_disruptions = sum(S .!= S0)

        µ_modified = µ - n_disruptions / (n_disruptions + N_DISRUPTIONS_HSC)
        t, Y, F, G = solve_host_virus(V0, I0, μ_modified, k, β, ϕ, α, tspan)

        G
    end

    # Define an output callback to write things to database
    db = init_database(fileName)
    summary_statement = Stmt(db, "INSERT INTO summary VALUES (?,?,?,?,?,?,?)")
    edges_statement = Stmt(db, "INSERT INTO edges VALUES (?,?,?,?)")

    function write_output(id, parent_id, removed_edge, S, C, fitness)
        @show id, parent_id, removed_edge, fitness, sum(S)

        parent_id_sql = if isnothing(parent_id) missing else parent_id end
        receptor_sql = if isnothing(removed_edge) missing else removed_edge[1] end
        effector_sql = if isnothing(removed_edge) missing else removed_edge[2] end
        imm_pwr_sql = sum(compute_phenotype(C))
        n_edges_sql = sum(S)
        execute(
            summary_statement,
            (id, parent_id_sql, receptor_sql, effector_sql, imm_pwr_sql, fitness, n_edges_sql)
        )

        if id == 0 || RECORD_EVERYTHING
            for edge in findall(S)
                execute(
                    edges_statement,
                    (id, edge[1], edge[2], C[edge])
                )
            end
        end
    end

    compute_fitness_landscape(
        S0, C0, optimize, compute_fitness, MIN_EDGES;
        output_callback = write_output
    )
end

"""
Initializes database tables to hold all the output we want to save.
"""
function init_database(file_name)
    if isfile(file_name)
        error("Output database already exists; delete first")
    end
    db = DB(file_name)

    # Metadata (dimensions)
    execute(db, """
        CREATE TABLE meta (
            n_receptors INTEGER,
            n_effectors INTEGER
        )
    """)
    execute(db, "INSERT INTO meta VALUES (?, ?)", (N_RECEPTORS, N_EFFECTORS))

    # Space-efficient representation of neutral network:
    # * ID of genotype
    # * parent ID (0 = initial genotype)
    # * receptor/effector removed relative to parent
    # * fitness
    execute(db, """
        CREATE TABLE summary (
            id INTEGER PRIMARY KEY,
            parent_id INTEGER,
            removed_receptor INTEGER,
            removed_effector INTEGER,
            immune_power REAL,
            fitness REAL,
            n_edges INTEGER
        )
    """)

    # All receptor-effector edges and fitted weights
    # Aside from initial genotype, all information here is implied by
    # the `summary` table.
    # Save space by setting `RECORD_EVERYTHING = false` above.
    execute(db, """
        CREATE TABLE edges (
            id INTEGER,
            receptor INTEGER,
            effector INTEGER,
            weight REAL
        )
    """)

    db
end



# Fix the seed
#Random.seed!(128)
#run_sims(high_redun_ind[3])



# high_redun_ind = [ 2  28  94 364 437 445 510 646 802]
# low_redun_ind = [19 223 438 679]
#
# for mat in high_redun_ind
#     run_sims(mat)
# end
