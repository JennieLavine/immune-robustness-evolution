#!/usr/bin/env julia

println("This takes a long time only because of initial compilation...")

include("../../julia/robustness-genomesize.jl")

using Plots

V0 = 10
I0 = 100
μ = 1.5
k = 10e-3
β = 1.0
ϕ = 10e3
α = 0
tspan = (0.0, 100)

t, Y, F, G = solve_host_virus(V0, I0, μ, k, β, ϕ, α, tspan)

# Plot V, I, F, G
VI_plot = plot(t, Y[:,1:2], labels = ["V" "I"])
savefig(VI_plot, "VI.png")
F_plot = plot(t, Y[:,3], labels = "F")
savefig(F_plot, "F.png")
G_plot = plot(t, log.(Y[:,4]), labels = "log(G)")
savefig(G_plot, "G.png")

@show F
@show G
@show log(G)
