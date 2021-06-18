#!/usr/bin/env julia

# This script generates a matrix, applies some disruptions, and optimizes
# weights via the Whitacre-type objective function or via a linear regression.

include("../../julia/robustness-genomesize.jl")

# Fix the seed
Random.seed!(1)

n_receptors = 8
n_effectors = 16
max_edges_per_receptor = 4
max_weight_per_receptor_init = 1.0 # WLOG
max_weight_per_receptor_whitacre = 1.0
weight_prior_sd_regression = 100.0

# Generate an initial receptor-effector network as a matrix.
S0 = generate_matrix_with_max_rowsum(
    n_receptors, n_effectors, max_edges_per_receptor
)
println("Initial network S0:")
display(S0)
print("\n\n")

# Generate weights as a matrix
C0 = generate_weights_with_max_rowsum(S0, max_weight_per_receptor_init)
println("Initial weights C0:")
display(C0)
print("\n\n")

# Compute the initial phenotype
println("Initial phenotype:")
initial_phenotype = compute_phenotype(S0, C0)
display(initial_phenotype)
print("\n\n")

# Add 10 random disruptions
S1 = add_random_disruption(S0)
for i in 1:9
    global S1 = add_random_disruption(S1)
end
println("Disrupted network S1:")
display(S1)
print("\n\n")
println("Missing entry in S1:")
display(findall(S0 .!= S1))
print("\n\n")

# Define a function to run the Whitacre-type objective function on this data
function optimize_whitacre_local(n_initial_conditions)
    # Optimize with the Whitacre-like objective function
    C1, obj_all, w_all = optimize_weights_real_whitacre(
        initial_phenotype, S1,
        max_weight_per_receptor_whitacre; n_initial_conditions = n_initial_conditions, return_all = true
    )

    obj_min, min_index = findmin(obj_all)
    obj_max, max_index = findmax(obj_all)
    w_min = w_all[min_index]

    println("objective function range: ($(obj_min), $(obj_max))")

    println("Value of objective function: $(obj_min)")
    println("Weights:")
    display(w_min)

    # Show new phenotype:
    println("New phenotype:")
    phenotype = compute_phenotype(S1, C1)
    display(phenotype)
    print("\n\n")

    # Show the difference between the phenotypes
    dphenotype = phenotype .- initial_phenotype
    println("Difference in phenotype (new - old):")
    display(dphenotype)
    print("\n\n")

    # Show which receptors are below capacity
    println("Receptors below capacity:")
    display(findall(dphenotype .< 0.0))
    print("\n\n")
end

println("Running Whitacre-like optimization 100 times...")
optimize_whitacre_local(100)

# Define a function to run the linear regression-based optimization
function optimize_regression_local()
    C1, w_mean, w_cov = optimize_weights_regression_analytical(
        initial_phenotype, S1,
        weight_prior_sd_regression;
        apply_prior_per_receptor = true
    )

    println("Weights:")
    display(w_mean)

    # Show new phenotype:
    println("New phenotype:")
    phenotype = compute_phenotype(S1, C1)
    display(phenotype)
    print("\n\n")

    # Show the difference between the phenotypes
    dphenotype = phenotype .- initial_phenotype
    println("Difference in phenotype (new - old):")
    display(dphenotype)
    print("\n\n")

    # Show which receptors are below capacity
    println("Receptors below capacity:")
    display(findall(dphenotype .< 0.0))
    print("\n\n")
end

println("Running linear regression-based optimization...")
optimize_regression_local()
