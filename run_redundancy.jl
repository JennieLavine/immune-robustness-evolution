

using SQLite: DB, Stmt
import SQLite.DBInterface.execute
using DelimitedFiles
using DataFrames

include("redundancy_metrics.jl")

Random.seed!(1234)

N_NETWORKS = 1000
N_RECEPTORS = 8
N_EFFECTORS = 16
MAX_ROW_SUM = 5
MAX_COMBO_SIZE = 4




########
#Call functions and generate matrix with results
########

s_mats = generate_networks(
    N_RECEPTORS, N_EFFECTORS, N_NETWORKS, MAX_ROW_SUM
    )

#writedlm("matrices.txt", reshape(s_mats,
#    (N_RECEPTORS * N_NETWORKS, N_EFFECTORS)))



stats_res = zeros(N_NETWORKS, 18)
redundancy_res = zeros(N_NETWORKS, MAX_COMBO_SIZE+1)

for i in 1:N_NETWORKS
    sys=s_mats[:,:,i]
    stats_res[i,:] = collect(calculate_graph_stats(sys))
    println(i)
    #orphE, orphR = quantify_redundancy(sys, MAX_COMBO_SIZE)
    #Evec = zeros(MAX_COMBO_SIZE+1)
    #for j in 1:(MAX_COMBO_SIZE+1)
    #    Evec[j]=sum(collect(values(orphE)).==(j-1))
    #end
    #redundancy_res[i,:] = Evec
end

writedlm("graph_stats.txt", stats_res)
#writedlm("redundancy_stats.txt", redundancy_res)








"""
Initializes database tables to hold all the output we want to save.
"""

function init_database()

    if isfile("redundancy_output.sqlite")
        error("Redundancy output database already exists; delete first")
    end
    db = DB("redundancy_output.sqlite")

    execute(db, """
        CREATE TABLE summary (
            id INTEGER PRIMARY KEY,
            n_receptors INTEGER,
            n_effectors INTEGER,
            bot_degeneracy INTEGER,
            bot_dens_net REAL ,
            bot_mean_deg REAL,
            bot_sd_deg REAL,
            bot_n_modules INTEGER,
            bot_modularity_net REAL,
            bot_mean_path_length REAL,
            bot_sd_path_length REAL,
            bip_degeneracy INTEGER,
            bip_dens_net REAL,
            bip_mean_deg REAL,
            bip_sd_deg REAL,
            bip_n_modules INTEGER,
            bip_modularity_net REAL)
        )
    """)
    db
end
