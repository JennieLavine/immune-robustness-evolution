using LightGraphs
using GraphPlot
using EcologicalNetworks
using EcologicalNetworksPlots
using Mangal
using DataFrames

include("julia/robustness-genomesize_play.jl")

function vec_triu_loop(M::AbstractMatrix{T}) where T
    m, n = size(M)
    m == n || throw(error("not square"))
    l = n*(n-1) ÷ 2
    v = zeros(Int16, l)
    k = 0
    @inbounds for i in 2:n
        for j in 1:(i-1)
            v[k + j] = M[j, i]
        end
        k += i-1
    end
    v
end


function generate_networks(n_receptors, n_effectors, n_networks, max_row_sum)
    S_mats = zeros(Bool, (n_receptors, n_effectors, n_networks) )

    mat_n = 1
    ind = 0
    while mat_n <= n_networks

        ind += 1
        if ind > 1000 * n_networks
            println("\n can't find enough non-degenerate networks \n")
            break
        end

        mat_prop = generate_matrix_with_max_rowsum(n_receptors, n_effectors, max_row_sum)

        #if the matrix is degenerate, replace with a new one
        if length(intersect(union(sum(mat_prop, dims=1), sum(mat_prop, dims=2)), 0))>=1
            continue
        end

        S_mats[ : , : , mat_n] = mat_prop

        mat_n += 1
    end
    (S_mats)
end


function generate_bottom_projection(N)

    n_effectors = size(N)[2]
    n_receptors = size(N)[1]

    bottom_projection = zeros(Bool, (n_effectors, n_effectors))

    for i in 1:n_effectors
        rec_cnxn = findall(N[:,i])
        for rec in rec_cnxn
            eff_cnxn = findall(N[rec,:])
            for eff in eff_cnxn
                bottom_projection[i, eff] = 1
            end
        end
    end
    bottom_projection
end

##### Write a function to calculate summary statistics on the generated graphs
#N is a boolean matrix, either representing the bipartite graph as a matrix,
# or a matrix representing the bottom projection of the bipartite graph

function calculate_graph_stats(N)

    n_receptors = size(N)[1]
    n_effectors = size(N)[2]

    R_names = map(string, repeat("R", n_receptors), 1:n_receptors)
    E_names = map(string, repeat("E", n_effectors), 1:n_effectors)

    bot_proj = generate_bottom_projection(N)
    net = UnipartiteNetwork(bot_proj)

    #path length
    path_length = shortest_path(net)

    upper_tri_vec = vec_triu_loop(path_length)
    bot_mean_path_length = mean(upper_tri_vec)
    bot_sd_path_length = std(upper_tri_vec)

    bot_degeneracy = isdegenerate(net)

    #modularity
    B, P = brim(lp(net)...)
    bot_n_modules = maximum(values(P))
    bot_modularity_net= Q(net, P)

    #degree distribution
    degree_net = EcologicalNetworks.degree(net)
    bot_mean_deg = mean(values(degree_net))
    bot_sd_deg = std(values(degree_net))
    bot_skew_deg = skewness(collect(values(degree_net)))

    #density
    bot_dens_net = connectance(net)


    ### Caclulate everything except path length for the bipartite network
    net = BipartiteNetwork(N, R_names, E_names)
    bip_degeneracy = isdegenerate(net)

    #modularity
    B, P = brim(lp(net)...)
    bip_n_modules = maximum(values(P))
    bip_modularity_net= Q(net, P)

    #degree distribution
    degree_net = EcologicalNetworks.degree(net)
    bip_mean_deg = mean(values(degree_net))
    bip_sd_deg = std(values(degree_net))
    bip_skew_deg = skewness(collect(values(degree_net)))
    #density
    bip_dens_net = connectance(net)

    (n_receptors, n_effectors,
        bot_degeneracy, bot_dens_net,
        bot_mean_deg, bot_sd_deg,
        bot_n_modules, bot_modularity_net,
        bot_mean_path_length, bot_sd_path_length,
        bip_degeneracy, bip_dens_net, bip_mean_deg, bip_sd_deg,
        bip_n_modules, bip_modularity_net, bot_skew_deg, bip_skew_deg)
end



function quantify_redundancy(N0, max_combo_size)

    n_receptors = size(N0)[1]
    n_effectors = size(N0)[2]

    n_orphan_receptors =  Dict{BitMatrix, Float64}()
    n_orphan_effectors =  Dict{BitMatrix, Float64}()

    next_id = 1

    function quantify_redundancy_recursive(id, N_parent)
        edges = findall(N_parent)

        if sum(N0)-sum(N_parent) > max_combo_size
            return nothing
        end

        for edge in edges
            N = remove_edge(N_parent, edge)
            id = next_id
            next_id += 1

            if haskey(n_orphan_receptors, N) || haskey(n_orphan_effectors, N)
                continue
            end

            bot_proj = generate_bottom_projection(N)

            if sum(N0) - sum(N) <= max_combo_size
                quantify_redundancy_recursive(id, N)
                n_orphan_effectors[N] = sum(sum(bot_proj, dims=1).==0)
                n_orphan_receptors[N] = sum(sum(bot_proj, dims=2).==0)
                continue
            end

        end

        nothing
    end

    quantify_redundancy_recursive(0, N0)
    (n_orphan_effectors, n_orphan_receptors)
end

#
#
# nets = generate_networks(8, 16, 10, 4)
# net1 = nets[:,:,1]
# (orphE, orphR) = quantify_redundancy(net1, 4)
# histogram(collect(values(orphE)))
#
#
#
#
#
#
#
#
#
#
#
# ##########################
# # Explore networks we have generated
# ###########################
#
#
# #Name vectors for the receptors (R) and effetors (E)
# R_names = map(string, repeat("R", N_RECEPTORS), 1:N_RECEPTORS)
# E_names = map(string, repeat("E", N_EFFECTORS), 1:N_EFFECTORS)
#
# ###########################
# # Plot of the network as bipartite graph with two rows
# ##########################
#
#
# function graph_two_row(N)
#     net = BipartiteNetwork(N, R_names, E_names)
#     I = initial(BipartiteInitialLayout, net)
#     position!(NestedBipartiteLayout(0.4), I, net)
#     plot(I, net, aspectratio=1)
#     scatter!(I, net, bipartite=true,
#         nodefill=brim(lp(net)...)[2], c=:lighttest,
#         )
#     fig_title = map(string, size(N), (" Receptors", " Effectors"))
#     title!(fig_title)
#     savefig(string(fig_title, ".png"))
# end
#
# net = BipartiteNetwork(s_mats[:,:,2], R_names, E_names)
# isdegenerate(net)
# I = initial(BipartiteInitialLayout, net)
# position!(NestedBipartiteLayout(0.4), I, net)
#
# #plot(I, net, aspectratio=1)
# #scatter!(I, net, bipartite=true,
#     #nodefill=salp(net)[2], #c=:cividis,
#     #nodesize = EcologicalNetworks.degree(net),
#     #series_annotations = collect(values(salp(net)[2]))
#     )
#
#
# plot(I, net, aspectratio=1)
# scatter!(I, net, bipartite=true, msw=0.0, nodefill=P, c=:lighttest)
#
#
# #uni_net = convert(UnipartiteNetwork, net)
# I = initial(RandomInitialLayout, net)
# for step in 1:(100richness(net))
#   position!(ForceDirectedLayout(0.3, 0.3; gravity=0.2), I, net)
# end
# plot(I, net, aspectratio=1)
# scatter!(I, net, bipartite=true, msw=0.0, nodefill=P, c=:lighttest)
#
#
# #############################
# # Metrics of the bipartite network
# ##############################
# #degree distribution
# degree_net = EcologicalNetworks.degree(net)
# histogram(collect(values(degree_net)))
#
# #density
# connectance(net)
#
# ##############################
# #Plot the bottom projection
# #############################
# bot_proj_m = bot_proj[:,:,2]
# #sum(bot_proj_m, dims=2)
# bot_proj_m[diagind(bot_proj_m)] .= 0
# bot_net = UnipartiteNetwork(bot_proj_m)
# #simplify!(bot_net)
#
# B, P = brim(lp(bot_net)...)
#
# I_bot = initial(RandomInitialLayout, bot_net)
# for step in 1:(100richness(bot_net))
#   position!(ForceDirectedLayout(0.3, 0.3; gravity=0.2), I_bot, bot_net)
# end
# plot(I_bot, bot_net, aspectratio=1)
# scatter!(I_bot, bot_net, bipartite=false,
#     nodefill=P, c=:Set2)
#
# ##############################
# #Metrics of the bottom projection
# #############################
#
# path_length = shortest_path(bot_net)
#
# function vec_triu_loop(M::AbstractMatrix{T}) where T
#     m, n = size(M)
#     m == n || throw(error("not square"))
#     l = n*(n-1) ÷ 2
#     v = zeros(Int16, l)
#     k = 0
#     @inbounds for i in 2:n
#         for j in 1:(i-1)
#             v[k + j] = M[j, i]
#         end
#         k += i-1
#     end
#     v
# end
#
# upper_tri_vec = vec_triu_loop(path_length)
# mean_path_length = mean(upper_tri_vec)
# std_path_length = std(upper_tri_vec)
#
#
# ################
# #Effect of removing nodes on network structure
# ################
#
#
#
#
#
# #plot the bipartite network
# net_1 = BipartiteNetwork(test_mat, R_names, E_names)
# I = initial(BipartiteInitialLayout, net_1)
# position!(NestedBipartiteLayout(0.4), I, net_1)
# plot(I, net_1, aspectratio=1)
# scatter!(I, net_1, bipartite=true)
#
#
# proj_net = UnipartiteNetwork(bottom_projection)
# test = shortest_path(proj_net, nmax=10)
#
#
#
# proj_I = initial(CircularInitialLayout, proj_1)
# plot(proj_I, proj_1, aspectratio=1)
# scatter!(proj_I, proj_1, bipartite=false)
#
#
# #calculate the density
# density = sum(test_mat)/(n_effectors*n_receptors)
#
#
#
#
# bot_overlap = overlap(net_1; dims=2)
# bot_mat = zeros(Bool, (n_effectors, n_effectors))
# for edge in bot_overlap
#     bot_mat
#
# df = DataFrame(E_names)
#
#
# UnipartiteNetwork(bot_overlap)
#
#
# net_uni_1 = convert(UnipartiteNetwork, net_1)
# I = initial(RandomInitialLayout, net_uni_1)
# plot(I, net_uni_1, aspectratio=1)
# scatter!(I, net_uni_1, bipartite=true)
#
# test = overlap(net_1, dims=1)
#
# initial()
#
#
# test_mat_edges = matrix_to_edges(test_mat)
#
# tops = 1:n_receptors
# bottoms = 1:n_effectors
# edges = test_mat_edges
#
# n_nodes = n_receptors + n_effectors
# G1 = Graph(n_nodes)
#     for edge in edges
#         add_edge!(G1, edge[1], n_receptors+edge[2])
#     end
#
# gplot(G1, layout=spring_layout)
#
#
# function bottom_projection_stats(edges)
#     tops = zeros(Int8, length(edges))
#     for i in 1:length(edges)
#         tops[i] = edges[i][1]
#     end
#
#     top_degrees = countmap(tops)
#
#     for i in 1:length(top_degrees)
#         findall()
#
#
#
#     degrees
# end
