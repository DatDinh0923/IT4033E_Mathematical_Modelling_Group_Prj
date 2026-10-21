# test_toy.jl — sanity checks on the toy instance of the challenge PDF
#   julia --project=. test_toy.jl
# Expected with ECMP r: MLU without waypoints 0.4375 (t=0), optimal 0.375,
# r(v0,v5,a35,0) = 0.75 as in Table 2 of the subject.
include(joinpath(@__DIR__, "src", "data.jl")); include(joinpath(@__DIR__, "src", "model.jl"))
dir = joinpath(@__DIR__, "toy")
for meth in (:ecmp, :pathcount)
    I = read_instance(dir, "toy"; rmethod = meth)
    v = I.id_to_vertex
    println("[$meth] r(v0,v5,a35,0) = ", round(r_value(I, 0, v[0], v[5], 6), digits = 4))
    println("[$meth] MLU without waypoints t=0: ", mlu_from_waypoints(I, Dict(1 => Int[], 2 => Int[]), 0))
    for (per, bud) in (([0], false), ([0, 1], true))
        m, x, λ, S = build_model(I; periods = per, budget = bud, timelimit = 60.0)
        res = solve!(m)
        println("[$meth] periods $per: $(res.status), MLU = $(round(res.obj, digits = 4))")
    end
end