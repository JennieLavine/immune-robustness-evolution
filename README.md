# Code for project on robustness and genome size

Jennie Lavine<br>
Ed Baskerville (initial example code)


Folder: Julia -- this folder has all the Julia code -- including writing functions, running sims and any plotting done in Julia.
(1) robustness-genomesize_play.jl -- main file with the core functions in it.  This is Jennie's edited version of Ed's robustness-genomesize.jl
(2) diff_eqs_plots.jl -- calls (1), runs simulations and plots trajectories, no network model used in these.

(3) redundancy_metrics.jl -- Uses the EcologicalNetworks + related packages in Julia define functions to calculate statistics on a set of randomly generated graphs.  Generates graphs using some of the core function in (1)
(4) run_redunancy.jl -- Runs many of the functions in (3) and saves results to a .txt file.
(5) graph_visualization.jl -- Generates  plots of the graphs and graph statistics computed in (4) including degree distributions, two-row representations of bipartite networks, bottom projections, and more.

(6) fitness_landscape.jl -- Define functions to run simulations and compute the fitness of all possible combinations of interruptions that leave at least MIN_EDGES in the graph.


Folder: R_RMarkdowns -- this folder contains R scripts and .Rmd files that run analyses on data in folders inside this one, mostly generated from Julia scripts.
(1) functional_rltnshps.R -- plots graphs for relationship between (immune power \& viral reproduction rate) and viral genome size.
(2) graph_stats.Rmd -- calculates statistics on graphs and runs linear models identifying important predictors of robustness / redundancy.
(3) plotting_fitness_landscapes.Rmd -- INCOMPLETE CODE! Writes functions to take output of fitness_landscape.jl and estimate fitness landscapes / deadends /etc.


Folder: R_Ed_original_files -- Contains code that Ed provided for running Julia through R.  Jennie never used these.

Folder: graph_plots -- Contains plots of graphs from initial set of simulations.

Folder: examples -- Ed's code on which Jennie's code is based containing examples of running most of the functions.

Folder: databases -- Downloaded data from some online databases of innate immune molecules

Folder:articles -- PDFs of articles of interest to this project.