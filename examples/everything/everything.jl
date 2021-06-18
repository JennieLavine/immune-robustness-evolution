#!/usr/bin/env julia

# This script generates a receptor-effector matrix, enumerates the neutral
# network, and writes it all out to a SQLite database.

include("../../julia/robustness-genomesize.jl")

using SQLite: DB, Stmt
import SQLite.DBInterface.execute

const N_RECEPTORS = 8
const N_EFFECTORS = 16
const MAX_EDGES_PER_RECEPTOR = 4
const MAX_WEIGHT_PER_RECEPTOR_INIT = 1.0
const WEIGHT_PRIOR_SD_REGRESSION = 100.0

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

function run(neutrality_threshold)
    # Fix the seed
    Random.seed!(1)

    # Generate an initial receptor-effector network as a matrix.
    S0 = generate_matrix_with_max_rowsum(
        N_RECEPTORS, N_EFFECTORS, MAX_EDGES_PER_RECEPTOR
    )
    println("Initial network S0:")
    display(S0)
    print("\n\n")

    # Generate weights as a matrix
    C0 = generate_weights_with_max_rowsum(S0, MAX_WEIGHT_PER_RECEPTOR_INIT)
    println("Initial weights C0:")
    display(C0)
    print("\n\n")

    # Compute the initial phenotype
    println("Initial phenotype:")
    initial_phenotype = compute_phenotype(S0, C0)
    display(initial_phenotype)
    print("\n\n")

    # Define an optimization function for the search
    function optimize(S)
        C, w_mean, w_cov = optimize_weights_regression_analytical(
            initial_phenotype, S,
            WEIGHT_PRIOR_SD_REGRESSION;
            apply_prior_per_receptor = true
        )
        max.(0.0, C)
    end

    # Define a fitness function for the search
    function compute_fitness(S, C)
        phenotype = min.(initial_phenotype, compute_phenotype(S, C))
        k = k0 * sum(phenotype) # Or whatever
        n_disruptions = sum(S .!= S0)

        µ_modified = µ - n_disruptions / (n_disruptions + N_DISRUPTIONS_HSC)
        t, Y, F, G = solve_host_virus(V0, I0, μ_modified, k, β, ϕ, α, tspan)

        log(G)
    end

    # Define an output callback to write things to database
    db = init_database()
    summary_statement = Stmt(db, "INSERT INTO summary VALUES (?,?,?,?,?)")
    edges_statement = Stmt(db, "INSERT INTO edges VALUES (?,?,?,?)")

    function write_output(id, parent_id, removed_edge, S, C, fitness)
        @show id, parent_id, removed_edge, fitness

        parent_id_sql = if isnothing(parent_id) missing else parent_id end
        receptor_sql = if isnothing(removed_edge) missing else removed_edge[1] end
        effector_sql = if isnothing(removed_edge) missing else removed_edge[2] end
        execute(
            summary_statement,
            (id, parent_id_sql, receptor_sql, effector_sql, fitness)
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

    compute_neutral_network(
        S0, C0, neutrality_threshold, optimize, compute_fitness;
        output_callback = write_output
    )
end

"""
Initializes database tables to hold all the output we want to save.
"""
function init_database()
    if isfile("output.sqlite")
        error("Output database already exists; delete first")
    end
    db = DB("output.sqlite")

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
            fitness REAL
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

neutral_network, one_neighborhood = run(0.02)
