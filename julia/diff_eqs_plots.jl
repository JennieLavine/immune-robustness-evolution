###plotting trajectories

include("julia/robustness-genomesize_play.jl")
default(linewidth=3)

const V0 = 10
const I0 = 100
const k0 = 10e-3
const μ = 1.5
const β = 1.0
const ϕ = 10e3
const α = 0
const tspan = (0.0, 100)

t, Y, F, G = solve_host_virus(V0, I0, μ, k0, β, ϕ, α, tspan)

plot(t, Y[:,1],
    label="viral titers",
    ylim=(0,1200),
    title = string("No extra viral genes or immune disruptions \n
    μ=", μ, "; k=", k0)
    ) # plot viral titers
plot!(t, Y[:,2], label = "immunity") # immunity
plot!(
    map(x -> hill(x, 10, 500), collect(1:1000)),
    inset = (1, bbox(0.05, 0.05, 0.5, 0.25, :center, :right)),
    #ticks = nothing,
    subplot = 2,
    bg_inside = nothing,
    label = "fitness function",
    legend = :bottomright,
    color=:purple
)
plot!(twinx(), t, Y[:,3],
    label= "cumulative fitness",
    color=:green, xticks=:none, legend= :top,
    ylim=(0,10),
    yforeground_color_axis = :green,
    yforeground_color_text = :green) #



t2, Y2, F2, G2 = solve_host_virus(V0, I0, μ-0.1, k0, β, ϕ, α, tspan)
plot(t2, Y2[:,1],
    label="viral titers",
    ylim = (0,1200),
    title = string("1 extra viral gene \n
         robustness --> no immune effector disuption \n
        μ=", μ-0.1, "; k=", k0)
        ) # plot viral titers
plot!(t2, Y2[:,2], label = "immunity") # immunity

plot!(twinx(), t2, Y2[:,3],
    label= "cumulative fitness",
    color=:green, xticks=:none, legend= :top,
    ylim=(0,10),
    yforeground_color_axis = :green,
    yforeground_color_text = :green)#




t3, Y3, F3, G3 = solve_host_virus(V0, I0, μ-0.2, k0-(0.1*k0), β, ϕ, α, tspan)
plot(t3, Y3[:,1],
    label="viral titers",
    ylim = (0,1200), # plot viral titers
    title = string("2 extra viral genes \n
         some immune disruption \n
        μ=", μ-0.2, "; k=", round(k0-(0.1*k0), sigdigits=2))
        )
plot!(t3, Y3[:,2], label = "immunity") # immunity
plot!(twinx(), t3, Y3[:,3],
    label= "cumulative fitness",
    color=:green, xticks=:none, legend= :top,
    ylim=(0,10),
    yforeground_color_axis = :green,
    yforeground_color_text = :green) #



t4, Y4, F4, G4 = solve_host_virus(V0, I0, μ-0.3, k0-(0.5*k0), β, ϕ, α, tspan)
    plot(t4, Y4[:,1],
        label="viral titers",
        ylim = (0,1200), # plot viral titers
        title = string("3 extra viral genes \n
        lots of immune disruption \n
            μ=", μ-0.4, "; k=", round(k0-(0.5*k0), sigdigits=2))
            )
    plot!(t4, Y4[:,2], label = "immunity") # immunity
    plot!(twinx(), t4, Y4[:,3],
        label= "cumulative fitness",
        color=:green, xticks=:none, legend= :top,
        ylim=(0,10),
        yforeground_color_axis = :green,
        yforeground_color_text = :green) #


mu_vals = collect(0.9:0.01:1.5)
k_vals = collect(0.006:0.00001:0.012)

fitness_res = Matrix{Float64}(undef, length(k_vals), length(mu_vals))

col_ind = 1
for mu in mu_vals
    row_ind=1
    println(string("mu=",mu))
    for k in k_vals
        println(string("k=", k))
        t, Y, F, G = solve_host_virus(V0, I0, mu, k, β, ϕ, α, tspan)
        fitness_res[row_ind, col_ind] = Y[end,4]
        row_ind += 1
    end
    col_ind += 1
end


heatmap(mu_vals, k_vals, fitness_res,
    xlabel = "μ", ylabel="k",
    title = "Fitness")
