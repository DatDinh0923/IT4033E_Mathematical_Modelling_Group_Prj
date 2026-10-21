"""
Orange / ROADEF 2026 "Keep the flow!" — segment routing, time periods 0 and 1.

Files follow the steps of `project/readme`:
`graph.jl` (Step 2), `rcoef.jl` (Step 3), `data.jl` (Step 4), `model_t0.jl` (Step 5),
`data_t1.jl` (Step 7), `model_t01.jl` (Step 8).
"""
module OrangeSR

using Graphs
using HiGHS
using JSON3
using JuMP
using SparseArrays

include("graph.jl")
include("rcoef.jl")
include("data.jl")
include("model_t0.jl")
include("data_t1.jl")
include("model_t01.jl")

export NetworkData, json_to_graph, metric_matrix, capacity_matrix
export all_pairs_shortest_data, count_shortest_paths, split_coefficients, check_split_coefficients
export InstanceT0, load_instance_t0, num_vertices, arc_loads, mlu, direct_paths
export ModelT0, build_model_t0, solve_t0!, extract_paths, write_srpaths
export PeriodT1, load_period1, routing_distance
export ModelT01, build_model_t01, solve_t01!

end
