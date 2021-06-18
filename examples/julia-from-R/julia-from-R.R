#!/usr/bin/env Rscript

source("../../R/call_julia.R")

# NB: Julia's JSON library and R's jsonlite library have opposite defaults
# for column- vs. row-major storage in JSON nested arrays, hence the calls to
# `t(...)`.

# Call an external Julia script file
cat("Calling Julia file:\n")
X1 <- t(call_julia_file_jsonlite("julia_file.jl"))
print(X1)

# Call Julia code stored in a string
julia_code <- "
using JSON

X = zeros(2, 3)
X[1,2] = true

println(JSON.json(X))
"
cat("Calling Julia string:\n")
X2 <- t(call_julia_string_jsonlite(julia_code))
print(X2)
