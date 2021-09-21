include("run_redundancy.jl")

high_redun_ind = [ 2  28  94 364 437 445 510 646 802]
low_redun_ind = [19 223 438 679]

random_samp = sample(1:1000, 10, replace = false)
##########################
# Explore networks we have generated
###########################


#Name vectors for the receptors (R) and effetors (E)
R_names = map(string, repeat("R", N_RECEPTORS), 1:N_RECEPTORS)
E_names = map(string, repeat("E", N_EFFECTORS), 1:N_EFFECTORS)

###########################
# Plot of the network as bipartite graph with two rows
##########################
function graph_two_row(N)
    net = BipartiteNetwork(N, R_names, E_names)
    I_net = initial(BipartiteInitialLayout, net)
    position!(NestedBipartiteLayout(0.4), I_net, net)
    plot(I_net, net, aspectratio=1)
    scatter!(I_net, net, bipartite=true,
        nodefill=brim(lp(net)...)[2], c=:lighttest,
        )
end


function graph_mod_bip(N)
    net = BipartiteNetwork(N, R_names, E_names)
    I_net = initial(RandomInitialLayout, net)
    for step in 1:(100richness(net))
      position!(ForceDirectedLayout(0.3, 0.3; gravity=0.2), I_net, net)
    end
    plot(I_net, net, aspectratio=1)
    scatter!(I_net, net, bipartite=true,
        nodefill=brim(lp(net)...)[2], c=:lighttest,
        )
end

function hist_bip_deg(N)
    net = BipartiteNetwork(N, R_names, E_names)
    degree_bip_net = EcologicalNetworks.degree(net)
    vals = collect(values(degree_bip_net))
    histogram(vals, bins=range(0, stop = maximum(vals), length = 8))
    annotate!(0.5, 1, string("skewness=\n", round(skewness(vals), sigdigits=2)))
end

function graph_bot_proj(N)
    bot_proj = generate_bottom_projection(N)
    net = UnipartiteNetwork(bot_proj)
    I_net = initial(RandomInitialLayout, net)
    for step in 1:(100richness(net))
      position!(ForceDirectedLayout(0.3, 0.3; gravity=0.2), I_net, net)
    end
    plot(I_net, net, aspectratio=1)
    scatter!(I_net, net, bipartite=true,
        nodefill=brim(lp(net)...)[2], c=:lighttest,
        )
end

function hist_bot_deg(N)
    bot_proj = generate_bottom_projection(N)
    net = UnipartiteNetwork(bot_proj)
    degree_bot_net = EcologicalNetworks.degree(net, dims=1)
    vals = collect(values(degree_bot_net))
    histogram(vals, bins=range(0, stop = maximum(vals), length = 8))
    annotate!(0.5, 1, string("skewness=\n", round(skewness(vals), sigdigits=2)))
end


#
# for s_mat_ind in high_redun_ind
#     s_mat = s_mats[:,:,s_mat_ind]
#     p1 = graph_two_row(s_mat)
#     p2 = graph_mod_bip(s_mat)
#     p3 = hist_bip_deg(s_mat)
#     p4 = graph_bot_proj(s_mat)
#     p5 = hist_bot_deg(s_mat)
#     println(s_mat_ind)
#     l = @layout [a b
#                 c d e]
#     plot(p1, p4, p3, p2, p5, layout = l)
#     fig_title = string("network ",s_mat_ind, " high redundancy")
#     savefig(string(fig_title, ".png"))
# end
#
#
# for s_mat_ind in low_redun_ind
#     s_mat = s_mats[:,:,s_mat_ind]
#     p1 = graph_two_row(s_mat)
#     p2 = graph_mod_bip(s_mat)
#     p3 = hist_bip_deg(s_mat)
#     p4 = graph_bot_proj(s_mat)
#     p5 = hist_bot_deg(s_mat)
#     println(s_mat_ind)
#
#     l = @layout [a b
#                 c d e]
#     plot(p1, p4, p3, p2, p5, layout = l)
#     fig_title = string("network ",s_mat_ind, " low redundancy")
#     savefig(string(fig_title, ".png"))
# end


### another high redundancy network, as detected from the full algorithm with virus
#
# s_mat_ind=341
# s_mat = s_mats[:,:,s_mat_ind]
# p1 = graph_two_row(s_mat)
# p2 = graph_mod_bip(s_mat)
# p3 = hist_bip_deg(s_mat)
# p4 = graph_bot_proj(s_mat)
# p5 = hist_bot_deg(s_mat)
# println(s_mat_ind)
# l = @layout [a b
#             c d e]
# plot(p1, p4, p3, p2, p5, layout = l)
# fig_title = string("network ",s_mat_ind, " high redundancy")
# savefig(string(fig_title, ".png"))




for s_mat_ind in random_samp
    s_mat = s_mats[:,:,s_mat_ind]
    p1 = graph_two_row(s_mat)
    p2 = graph_mod_bip(s_mat)
    p3 = hist_bip_deg(s_mat)
    p4 = graph_bot_proj(s_mat)
    p5 = hist_bot_deg(s_mat)
    println(s_mat_ind)

    l = @layout [a b
                c d e]
    plot(p1, p4, p3, p2, p5, layout = l)
    fig_title = string("network ",s_mat_ind, " random redundancy")
    savefig(string(fig_title, ".png"))
end
