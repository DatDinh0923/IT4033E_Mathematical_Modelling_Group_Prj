"""
Orange / ROADEF 2026 "Keep the flow!" — segment routing, time period 0.

Files follow the steps of `project/readme`:
`graph.jl` (Step 2), `rcoef.jl` (Step 3), `data.jl` (Step 4), `model_t0.jl` (Step 5).
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

export NetworkData, json_to_graph, metric_matrix, capacity_matrix
export all_pairs_shortest_data, count_shortest_paths, split_coefficients, check_split_coefficients
export InstanceT0, load_instance_t0, num_vertices, arc_loads, mlu, direct_paths
export ModelT0, build_model_t0, solve_t0!, extract_paths, write_srpaths

end
