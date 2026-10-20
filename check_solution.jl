# Đọc một file nghiệm *-srpaths.json (định dạng của đề) và tính lại MỌI THỨ một cách độc lập:
# MLU từng chu kỳ, link nghẽn nhất, số waypoint, số thay đổi giữa t=0 và t=1, rồi so với bảng.
#   julia --project=. check_solution.jl                      # kiểm tra mọi nghiệm trong results/
#   julia --project=. check_solution.jl setA-09 9            # in chi tiết một nghiệm
include(joinpath(@__DIR__, "src", "orange.jl"))

const DATA = joinpath(@__DIR__, "..", "project", "setA")
const RES = joinpath(@__DIR__, "results")

"srpaths.json -> P[ti][d] = [s, w1, ..., t] (đỉnh Julia); demand không có trong file = đi thẳng."
function read_srpaths(inst::Instance, file, periods)
    P = [baseline_paths(inst) for _ in periods]
    for e in JSON3.read(read(file, String))["srpaths"]
        ti = findfirst(==(Int(e["t"])), periods)
        d = Int(e["d"]) + 1
        P[ti][d] = vcat(inst.src[d], [inst.net.id_to_vertex[Int(w)] for w in e["w"]], inst.dst[d])
    end
    return P
end

function check(name, step; verbose = false, sub = "step$step")
    periods = step == 6 ? [0] : [0, 1]
    inst = load_instance(joinpath(DATA, name))
    pd = [period_data(inst, t) for t in periods]
    alives = [p[1] for p in pd]; sds = [p[2] for p in pd]
    P = read_srpaths(inst, joinpath(RES, "solutions", sub, name * "-srpaths.json"), periods)
    chk = verify_solution(inst, sds, alives, periods, P)      # lỗi nếu nghiệm không hợp lệ
    raw = JSON3.read(read(joinpath(RES, "raw", sub, name * "_s1.json"), String))
    ok = abs(chk.mlu - raw["best_obj"]) <= 1e-6 * max(1, raw["best_obj"])
    @printf("%-8s %-9s MLU tính lại = %.6f | bảng = %.6f  %s | waypoint = %d | thay đổi = %s\n",
            name, sub, chk.mlu, raw["best_obj"], ok ? "KHỚP" : "LỆCH", chk.waypoints, chk.changes)
    if verbose
        for (ti, t) in enumerate(periods)
            λ = arc_loads(inst, sds[ti], t, P[ti])
            λ[.!alives[ti]] .= 0
            top = sortperm(λ; rev = true)[1:3]
            println("  t = $t: MLU = $(round(chk.mlu_per_period[ti], digits = 4)); 3 link tải cao nhất:")
            for a in top
                u, v = inst.net.vertex_to_id[inst.net.arc_from[a]], inst.net.vertex_to_id[inst.net.arc_to[a]]
                @printf("     link id %-5d (v%d -> v%d)  tải = %.4f\n", inst.net.link_id[a], u, v, λ[a])
            end
            nwp = [length(p) - 2 for p in P[ti]]
            println("  t = $t: số demand dùng 0/1/2/3+ waypoint = ",
                    [count(==(k), nwp) for k in 0:2], " / ", count(>=(3), nwp))
        end
    end
    return ok
end

if length(ARGS) >= 2
    check(ARGS[1], parse(Int, ARGS[2]); verbose = true, sub = length(ARGS) >= 3 ? ARGS[3] : "step$(ARGS[2])")
else
    bad = 0
    for sub in ("step6", "step9", "step6_ipm", "step9_ipm"), f in sort(readdir(joinpath(RES, "raw", sub)))
        step = parse(Int, sub[5:5])
        check(f[1:7], step; sub) || (global bad += 1)
    end
    println(bad == 0 ? "TẤT CẢ NGHIỆM KHỚP VỚI BẢNG" : "$bad nghiệm LỆCH")
end
