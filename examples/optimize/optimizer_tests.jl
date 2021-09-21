#Playing with the optimizer

include("../../julia/robustness-genomesize_play.jl")
# Fix the seed
#Random.seed!(10)

n_receptors = 8
n_effectors = 16
max_edges_per_receptor = 5
max_weight_per_receptor_init = 1.0 # WLOG
max_weight_per_receptor_whitacre = 1.0
weight_prior_sd_regression = 100.0
proportion_capacity = 0.5
capacity = proportion_capacity
n_initial_conditions = 100

# Generate an initial receptor-effector network as a matrix.
S0 = generate_matrix_with_max_rowsum(
    n_receptors, n_effectors, max_edges_per_receptor
)


S0 = s_mats[:,:,94]
heatmap(S0, title="Undisrupted immune connections",
    xlabel="Effectors", ylabel="Receptors")

#for capacity in proportion_capacity
    C0 = generate_weights_with_max_rowsum(
        S0, max_weight_per_receptor_init, capacity
    )

heatmap(C0, title="Undisrupted immune power",
        xlabel="Effectors", ylabel="Receptors")


    initial_phenotype = compute_phenotype(C0)
    bar(initial_phenotype, title = "Initial phenotype",
        ylims=(0,3), labels="Effector contribution")
    xlabel!("Effectors")
    ylabel!("Contribution per effector")

    res_obj = []
    res_w = []
    res_C1 = []

    n_disrupts = n_effectors-10
    n_reps = 100

    for n_dels in repeat(1:n_disrupts, inner=n_reps)
        global S1 = add_random_disruption(S0)
        for i in 1:n_dels
            global S1 = add_random_disruption(S1)
        end

        C1, obj_all, w_all = optimize_weights_real_whitacre(
        initial_phenotype, S1,
        max_weight_per_receptor_whitacre; n_initial_conditions = n_initial_conditions, return_all = true
        )

        push!(res_obj, obj_all)
        push!(res_w, w_all)
        push!(res_C1, C1)
    end

    res_phenotype = []
    for i in 1:length(res_C1)
        push!(res_phenotype, sum(res_C1[i], dims=1)[:])
    end

    res_sum = []
    res_mean = []

    for i in 1:length(res_phenotype)
        push!(res_sum, sum(res_phenotype[i]))
    end

    mean_phenotype = []
    sd_phenotype = []
    for i in 0:(n_disrupts-1)
        push!(mean_phenotype, mean(res_sum[((i*n_reps)+1):((i*n_reps)+n_reps)]))
        push!(sd_phenotype, std(res_sum[((i*n_reps)+1):((i*n_reps)+n_reps)]))
        println(i)
    end

    plot(mean_phenotype, lw=3, label="mean phenotype")
    plot!(mean_phenotype - 1.96*sd_phenotype, color="gray")
    plot!(mean_phenotype + 1.96*sd_phenotype, color="gray")
    xlabel!("number of disruptions")
    ylabel!("phenotype (sum of effector power)")
    fig_title = string("test5")
    title!(fig_title)

    savefig(string(fig_title, ".png"))

#end
