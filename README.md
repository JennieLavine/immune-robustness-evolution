# Code for project on robustness and genome size

Jennie Lavine<br>
Ed Baskerville (initial example code)


Folder: Julia
(1) robustness-genomesize_play.jl -- main file with the core functions in it.  This is Jennie's edited version of Ed's robustness-genomesize.jl
(2) diff_eqs_plots.jl -- calls (1), runs simulations and plots trajectories, no network model used in these.

(3) redundancy_metrics.jl -- Uses the EcologicalNetworks + related packages in Julia define functions to calculate statistics on a set of randomly generated graphs.  Generates graphs using some of the core function in (1)
(4) run_redunancy.jl -- Runs many of the functions in (3) and saves results to a .txt file.
(5) graph_visualization.jl -- Generates  plots of the graphs and graph statistics computed in (4) including degree distributions, two-row representations of bipartite networks, bottom projections, and more.

(6) fitness_landscape.jl -- Define functions to run simulations and compute the fitness of all possible combinations of interruptions that leave at least MIN_EDGES in the graph.



functional_rltnshps.R -- plots graphs for relationship between 