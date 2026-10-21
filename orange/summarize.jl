# summarize.jl — turn results/<mode>-<rmethod>.csv into the report table
# (averages over the runs of each instance), printed as Markdown and LaTeX.
#   julia --project=. summarize.jl results/p0-pathcount.csv
using Statistics, Printf

file = ARGS[1]
lines = readlines(file)
hdr = split(lines[1], ",")
col(name) = findfirst(==(name), hdr)
rows = [split(l, ",") for l in lines[2:end] if !isempty(strip(l))]

groups = Dict{String,Vector{Vector{SubString{String}}}}()
for r in rows; push!(get!(groups, r[1], []), r); end

num(s) = (x = tryparse(Float64, s); x === nothing ? NaN : x)
nanmean(v) = (w = filter(!isnan, v); isempty(w) ? NaN : mean(w))
fmt(x; d = 4) = isnan(x) ? "–" : @sprintf("%.*f", d, x)

println("| Instance | Vertices | Links | Demands | Runs | Objective (MLU) | CPU (s) | Status | Gap | B&B nodes |")
println("|---|---|---|---|---|---|---|---|---|---|")
tex = String[]
for name in sort(collect(keys(groups)))
    g = groups[name]; r1 = g[1]
    obj  = nanmean(num.(getindex.(g, col("obj"))))
    cpu  = nanmean(num.(getindex.(g, col("cpu"))))
    gap  = nanmean(num.(getindex.(g, col("gap"))))
    node = nanmean(num.(getindex.(g, col("nodes"))))
    sts  = unique(getindex.(g, col("status")))
    st   = join(sts, "/")
    gs   = all(==("Optimal"), sts) ? "–" : @sprintf("%.2f%%", 100gap)
    println("| $name | $(r1[col("V")]) | $(r1[col("A")]) | $(r1[col("D")]) | $(length(g)) | ",
            fmt(obj), " | ", fmt(cpu; d = 1), " | $st | $gs | ", fmt(node; d = 0), " |")
    push!(tex, "$name & $(r1[col("V")]) & $(r1[col("A")]) & $(r1[col("D")]) & $(fmt(obj)) & " *
               "$(fmt(cpu; d = 1)) & $st & $(replace(gs, "%" => "\\%")) & $(fmt(node; d = 0)) \\\\")
end
println("\n% LaTeX rows: Instance & Vertices & Links & Demands & Obj & CPU & Status & Gap & Nodes")
foreach(println, tex)