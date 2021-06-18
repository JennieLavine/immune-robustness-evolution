library(jsonlite)

# Return the standard output, as a string, from running a Julia script
# located at `julia_filename`.
call_julia_file <- function(julia_filename) {
  system2("julia", c(julia_filename), stdout = TRUE)
}

# Return the standard output from running a Julia script as a object parsed
# via the `fromJSON()` function in the `jsonlite` package.
#
# NB: By default, Julia's `JSON` package generates column-major nested arrays
# from matrices, but R's `jsonlite` package parses them as row-major matrices.
# Therefore, for matrices, the output must be transposed.
call_julia_file_jsonlite <- function(julia_filename) {
  fromJSON(call_julia_file(julia_filename))
}

# Return the standard output, as a string, from running the Julia code in the
# string `julia_code`.
#
# To avoid extraneous REPL output, the code is first written to a temporary
# file, and then `call_julia_file` is called.
call_julia_string <- function(julia_code) {
  julia_filename <- tempfile()
  write(julia_code, julia_filename)
  call_julia_file(julia_filename)
}

# Return the standard output, as a string, from running the Julia code in the
# string `julia_code`.
#
# To avoid extraneous REPL output, the code is first written to a temporary
# file, and then `call_julia_file` is called.
call_julia_string_jsonlite <- function(julia_code) {
  julia_filename <- tempfile()
  write(julia_code, julia_filename)
  call_julia_file_jsonlite(julia_filename)
}

# Return the Unix error code from running a Julia script; do not transfer any
# data.
run_julia_file <- function(julia_filename) {
  system2("julia", c(julia_filename))
}

# Return the Unix error code from running a Julia string; do not transfer any
# data.
run_julia_string <- function(julia_code) {
  julia_filename <- tempfile()
  write(julia_code, julia_filename)
  run_julia_file(julia_filename)
}
