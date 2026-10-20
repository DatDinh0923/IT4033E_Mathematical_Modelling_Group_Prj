# =============================================================================
#  run_all.jl — chạy TOÀN BỘ project từ đầu đến cuối:
#      julia --project=. run_all.jl
#  1. test toy  2. prepare 20 instance  3. step 6 & step 9 (NSEEDS seed mỗi instance)
#  4. aggregate -> results/table_step6.{csv,md}, results/table_step9.{csv,md}
#
#  GIỚI HẠN RAM 15 GB (= 15·10^9 byte ≈ 13.97 GiB) cho TOÀN BỘ quá trình. Các ngưỡng dưới đây
#  tính theo GiB (2^30 byte). Chạy bằng:
#      ./run.sh        (= systemd-run --scope -p MemoryMax=15000000000 -p MemorySwapMax=0 julia run_all.jl)
#  cgroup của systemd là chốt chặn cứng: tổng RAM không bao giờ vượt 15 GB.
#  Mỗi lần giải là một process Julia riêng (run_one.jl, HiGHS threads = 1) để chạy
#  song song nhiều lần giải trong 15 GB đó. Bộ điều phối giữ RAM an toàn bằng 3 lớp:
#   (a) ước lượng RAM mỗi job (lần đầu: theo kích thước mô hình; sau khi seed đầu tiên
#       của cùng (instance, step) chạy xong: 1.3 × peak RSS đo được),
#       tổng ước lượng của các job đang chạy ≤ MEM_BUDGET_GB;
#   (b) chỉ khởi động job mới khi RAM đang dùng (đo trong cgroup) + ước lượng ≤ MEM_START_GB;
#   (c) watchdog: nếu RAM đang dùng > MEM_KILL_GB thì kill job đang chiếm RAM nhiều
#       nhất và đưa lại vào hàng đợi với ước lượng RAM lớn hơn.
#  Kết quả đã có (results/raw/...json) được bỏ qua -> chạy lại lệnh này để resume.
# =============================================================================
using JSON3, Printf, Dates

const ROOT = @__DIR__
const JULIA = joinpath(Sys.BINDIR, Base.julia_exename())
const RES = joinpath(ROOT, "results")
const INSTANCES = ["setA-" * lpad(i, 2, '0') for i in 1:20]

# ---------------- cấu hình (có thể đổi bằng biến môi trường) ----------------
envint(k, d) = parse(Int, get(ENV, k, string(d)))
envf(k, d) = parse(Float64, get(ENV, k, string(d)))
const NSEEDS        = envint("ORANGE_NSEEDS", 1)        # 1 lần chạy / instance (readme: 10)
const MAX_JOBS      = envint("ORANGE_MAX_JOBS", 12)     # số process giải song song
const MEM_BUDGET_GB = envf("ORANGE_MEM_BUDGET_GB", 12.5) # tổng ước lượng của các job đang chạy
const MEM_START_GB  = envf("ORANGE_MEM_START_GB", 12.0)  # dùng thực tế + ước lượng job mới
const MEM_KILL_GB   = envf("ORANGE_MEM_KILL_GB", 13.3)  # vượt ngưỡng này -> kill job lớn nhất
const DEADLINE_H    = envf("ORANGE_DEADLINE_H", 40.0)   # tổng thời gian tối đa
const TL = Dict(6 => envf("ORANGE_TL6", 1800.0), 9 => envf("ORANGE_TL9", 3600.0))
# mô hình đầy đủ khi số biến x ≤ ngưỡng; lớn hơn -> mô hình rút gọn (K waypoint ứng viên)
const EXACT_MAX_X   = envf("ORANGE_EXACT_MAX_X", 2.5e6)
const RESTRICTED_X  = envf("ORANGE_RESTRICTED_X", 6.0e5)  # số biến x mục tiêu của mô hình rút gọn

# Nếu cgroup chạm 15 GB, OOM killer của kernel phải giết một process giải (RSS lớn), không phải
# bộ điều phối: bộ điều phối có oom_score_adj = -900, process con được đặt lại về 0.
try write("/proc/self/oom_score_adj", "-900") catch end
child(cmd) = `sh -c 'echo 0 > /proc/self/oom_score_adj; exec "$0" "$@"' $(cmd.exec)`

const RAMKILLS = joinpath(RES, "ram_kills.txt")   # (instance, step) từng bị kill vì RAM
const T0 = time()
const LOG = open(joinpath(mkpath(RES), "run_all.log"), "a")
logmsg(s) = (line = "[$(Dates.format(now(), "mm-dd HH:MM:SS")) +$(round((time()-T0)/3600, digits=2))h] $s";
             println(line); println(LOG, line); flush(LOG); flush(stdout))

# cgroup v2 của process này (systemd-run --scope tạo riêng một cgroup có memory.max)
const CGROUP = let l = first(eachline("/proc/self/cgroup"))
    joinpath("/sys/fs/cgroup", lstrip(split(l, ":"; limit = 3)[3], '/'))
end
# Đo RAM ở cấp orange.slice (tổng 15 GB của MỌI process trong slice) nếu có, không chỉ scope này.
const CGMEM = let par = dirname(CGROUP)
    endswith(par, ".slice") && isfile(joinpath(par, "memory.max")) &&
        strip(read(joinpath(par, "memory.max"), String)) != "max" ? par : CGROUP
end
const CG_LIMITED = isfile(joinpath(CGMEM, "memory.max")) && strip(read(joinpath(CGMEM, "memory.max"), String)) != "max"

"RAM đang dùng (GB) của cả nhóm process: memory.current của cgroup trừ page cache (thu hồi được)."
function mem_used_gb(running)
    if CG_LIMITED
        cur = parse(Int, strip(read(joinpath(CGMEM, "memory.current"), String)))
        file = 0
        for l in eachline(joinpath(CGMEM, "memory.stat"))
            startswith(l, "file ") && (file = parse(Int, split(l)[2]); break)
        end
        return (cur - file) / 2^30
    end
    return 0.5 + sum((rss_gb(getpid(j.proc)) for j in running); init = 0.0)
end
function rss_gb(pid)
    f = "/proc/$pid/status"
    isfile(f) || return 0.0
    for l in eachline(f)
        startswith(l, "VmRSS:") && return parse(Int, split(l)[2]) / 2^20
    end
    return 0.0
end

# ---------------- 1. test toy + 2. prepare ----------------
function run_cmd(cmd, logfile)
    ok = success(pipeline(child(cmd); stdout = logfile, stderr = logfile))
    ok || logmsg("LỖI: $(cmd) (xem $logfile)")
    return ok
end

function prepare_all()
    logdir = mkpath(joinpath(RES, "logs"))
    run_cmd(`$JULIA --project=$ROOT $(joinpath(ROOT, "test_toy.jl"))`, joinpath(logdir, "test_toy.log")) ||
        error("test toy thất bại")
    logmsg("test toy: OK")
    todo = [nm for nm in INSTANCES if !isfile(joinpath(RES, "instances", nm * ".json"))]
    # Cận dưới MCF (LP interior point) của instance lớn tốn tới ~8 GB (setA-20) -> chạy tuần tự;
    # instance nhỏ (< 1 GB) chạy 3 process cùng lúc.
    prep(nm) = run_cmd(`$JULIA --project=$ROOT $(joinpath(ROOT, "prepare.jl")) $nm`, joinpath(logdir, "prepare_$nm.log"))
    big = ["setA-$i" for i in 14:20]
    asyncmap(prep, filter(!in(big), todo); ntasks = 3)
    foreach(prep, filter(in(big), todo))
    info = Dict(nm => JSON3.read(read(joinpath(RES, "instances", nm * ".json"), String)) for nm in INSTANCES)
    logmsg("prepare: xong $(length(info)) instance")
    return info
end

# ---------------- 3. jobs ----------------
mutable struct Job
    inst::String
    step::Int
    seed::Int
    mode::String
    K::Int
    est_gb::Float64
    proc::Union{Nothing,Base.Process}
    started::Float64
    attempts::Int
    ram_killed::Bool
end

resfile(j::Job) = joinpath(RES, "raw", "step$(j.step)", "$(j.inst)_s$(j.seed).json")
key(j::Job) = (j.inst, j.step)

"Chọn exact / restricted và K cho (instance, step) dựa trên kích thước mô hình."
function choose_mode(info, step)
    nx = step == 6 ? info["x_vars_step6"] : info["x_vars_step9"]
    nx <= EXACT_MAX_X && return ("exact", 0, nx)
    # K giống nhau cho step 6 và step 9 (tính theo kích thước step 9): mô hình step 6 đúng bằng
    # phần t = 0 của mô hình step 9 -> hai bảng so sánh được. (Lần chạy đầu cho step 6 K lớn
    # hơn — tới 30 — thì LP gốc của setA-06/07/08/13 không giải xong trong 1800 s.)
    D = info["demands"]
    per_d = RESTRICTED_X / (2D)                          # số segment / demand / chu kỳ
    K = clamp(floor(Int, sqrt(per_d)) - 1, 3, 30)        # (K+2) đỉnh -> ~ (K+1)^2 segment
    return ("restricted", K, D * (K + 1)^2 * (step == 6 ? 1 : 2))
end

# Ước lượng RAM ban đầu (GB). Số đo của lần chạy đầu (OOM lúc 13:30): mô hình rút gọn
# ~5.8e5 biến x dùng 2.6–7.1 GB (tăng theo n và |D|, không chỉ theo số biến x) -> ước lượng
# bảo thủ 6.5 GB; mô hình exact: setA-05 step 9 (1.98e6 biến x) đo 4.3 GB khi chưa giải xong
# LP gốc, cộng thêm bộ nhớ của IPM -> 1.5 + 4 KB/biến x.
initial_est_gb(mode, nx) = mode == "exact" ? 1.5 + 4.0e-6 * nx : 6.5

function main()
    info = prepare_all()
    jobs = Job[]
    for seed in 1:NSEEDS, step in (9, 6), nm in reverse(INSTANCES)   # instance lớn trước (LPT)
        mode, K, nx = choose_mode(info[nm], step)
        push!(jobs, Job(nm, step, seed, mode, K, initial_est_gb(mode, nx), nothing, 0.0, 0, false))
    end
    for step in (6, 9), nm in INSTANCES
        mode, K, nx = choose_mode(info[nm], step)
        logmsg(@sprintf("plan step%d %s: %s K=%d  ~%.2e biến x  ước lượng RAM %.1f GB", step, nm, mode, K, nx, initial_est_gb(mode, nx)))
    end
    filter!(j -> !isfile(resfile(j)), jobs)
    # các lần bị kill vì RAM ở những lần chạy trước -> thu nhỏ mô hình tương ứng
    if isfile(RAMKILLS)
        for l in eachline(RAMKILLS), j in jobs
            isempty(strip(l)) && continue
            nm, st = split(l)
            (j.inst == nm && j.step == parse(Int, st)) || continue
            j.mode == "exact" ? (j.mode = "restricted"; j.K = 20) : (j.K = max(3, round(Int, 0.6 * j.K)))
            logmsg("đã từng bị kill vì RAM -> $(j.inst) step$(j.step): $(j.mode) K=$(j.K)")
        end
    end
    # job còn đang chạy từ lần điều phối trước (process mồ côi): không chạy trùng; nó tự ghi kết quả
    alive = read(`ps -eo args`, String)
    for j in copy(jobs)
        if occursin("run_one.jl $(j.inst) $(j.step) $(j.seed) ", alive)
            logmsg("ĐANG CHẠY từ trước, bỏ qua: $(j.inst) step$(j.step) s$(j.seed)")
            filter!(!=(j), jobs)
        end
    end
    logmsg("còn $(length(jobs)) job; MAX_JOBS=$MAX_JOBS MEM_BUDGET=$(MEM_BUDGET_GB)GB NSEEDS=$NSEEDS cgroup đo RAM=$CGMEM giới hạn=$CG_LIMITED")

    measured = Dict{Tuple{String,Int},Float64}()   # peak RSS đo được (GB)
    overhead = Dict{Tuple{String,Int},Float64}()   # thời gian ngoài optimize! (s)
    for nm in INSTANCES, step in (6, 9), s in 1:NSEEDS   # nạp số đo từ các run đã có (khi resume)
        f = joinpath(RES, "raw", "step$step", "$(nm)_s$s.json")
        isfile(f) || continue
        r = JSON3.read(read(f, String))
        measured[(nm, step)] = max(get(measured, (nm, step), 0.0), r["maxrss_gb"])
    end
    running = Job[]
    deadline = T0 + DEADLINE_H * 3600
    logdir = mkpath(joinpath(RES, "logs", "runs"))
    probing = Set{Tuple{String,Int}}()
    last_report = 0.0

    while !isempty(jobs) || !isempty(running)
        # --- thu job đã xong ---
        for j in copy(running)
            process_exited(j.proc) || continue
            filter!(!=(j), running)
            delete!(probing, key(j))
            if isfile(resfile(j))
                r = JSON3.read(read(resfile(j), String))
                measured[key(j)] = max(get(measured, key(j), 0.0), r["maxrss_gb"])
                overhead[key(j)] = (time() - j.started) - r["wall_time"]
                logmsg(@sprintf("XONG  %s step%d s%d [%s] %s obj=%.5f gap=%.4f nodes=%d cpu=%.0fs rss=%.1fGB",
                                j.inst, j.step, j.seed, j.mode, r["status"], r["best_obj"], r["gap"],
                                r["nodes"], r["cpu_time"], r["maxrss_gb"]))
            else
                j.attempts += 1
                logmsg("THẤT BẠI $(j.inst) step$(j.step) s$(j.seed) (lần $(j.attempts), exit=$(j.proc.exitcode))")
                if j.attempts < 3
                    if j.ram_killed
                        # bị kill vì RAM: thu nhỏ mô hình (exact -> rút gọn K = 20; rút gọn -> K × 0.6)
                        oldK = j.K
                        j.mode == "exact" ? (j.mode = "restricted"; j.K = 20) : (j.K = max(3, round(Int, 0.6 * j.K)))
                        logmsg("  -> chạy lại với mô hình nhỏ hơn: $(j.mode) K=$(j.K) (trước: K=$oldK)")
                        j.ram_killed = false
                    end
                    j.est_gb = min(1.5 * j.est_gb, MEM_BUDGET_GB); j.proc = nothing
                    push!(jobs, j)
                end
            end
        end
        # --- watchdog RAM ---
        if mem_used_gb(running) > MEM_KILL_GB && !isempty(running)
            victim = argmax(j -> rss_gb(getpid(j.proc)), running)
            logmsg(@sprintf("RAM cao (%.1f GB đang dùng): kill %s step%d s%d (rss %.1f GB)", mem_used_gb(running),
                            victim.inst, victim.step, victim.seed, rss_gb(getpid(victim.proc))))
            victim.est_gb = max(victim.est_gb, 1.5 * rss_gb(getpid(victim.proc)))
            victim.ram_killed = true
            open(io -> println(io, victim.inst, " ", victim.step), RAMKILLS, "a")   # nhớ qua các lần resume
            kill(victim.proc); wait(victim.proc)
        end
        # --- khởi động job mới ---
        used = sum((j.est_gb for j in running); init = 0.0)
        for j in copy(jobs)
            length(running) >= MAX_JOBS && break
            k = key(j)
            # seed khác của cùng (instance, step) chờ lần đầu đo xong RAM
            (k in probing) && continue
            est = haskey(measured, k) ? max(1.3 * measured[k], 1.0) :
                  # step 6 cùng instance, cùng kiểu mô hình nhưng chỉ một chu kỳ -> không lớn hơn step 9
                  (j.step == 6 && haskey(measured, (j.inst, 9)) && j.mode == "restricted") ?
                      max(1.1 * measured[(j.inst, 9)], 1.0) : j.est_gb
            j.attempts > 0 && (est = max(est, j.est_gb))
            expected = TL[j.step] + get(overhead, k, 600.0)
            if time() + expected > deadline
                filter!(!=(j), jobs)
                logmsg("BỎ (hết thời gian) $(j.inst) step$(j.step) s$(j.seed)")
                continue
            end
            (used + est <= MEM_BUDGET_GB && mem_used_gb(running) + est <= MEM_START_GB) || continue
            j.est_gb = est
            lf = joinpath(logdir, "$(j.inst)_step$(j.step)_s$(j.seed).log")
            cmd = `$JULIA --project=$ROOT $(joinpath(ROOT, "run_one.jl")) $(j.inst) $(j.step) $(j.seed) $(j.mode) $(TL[j.step]) $(max(j.K, 1))`
            j.proc = run(pipeline(child(cmd); stdout = lf, stderr = lf); wait = false)
            j.started = time()
            haskey(measured, k) || push!(probing, k)
            push!(running, j); filter!(!=(j), jobs)
            used += est
            logmsg(@sprintf("CHẠY  %s step%d s%d [%s K=%d] ước lượng %.1f GB (đang chạy %d, RAM dùng %.1f GB)",
                            j.inst, j.step, j.seed, j.mode, j.K, est, length(running), mem_used_gb(running)))
        end
        # --- job treo (quá 1.5 × time limit + 1 giờ) ---
        for j in running
            if time() - j.started > 1.5 * TL[j.step] + 3600
                logmsg("TREO: kill $(j.inst) step$(j.step) s$(j.seed)"); kill(j.proc)
            end
        end
        if time() - last_report > 900
            last_report = time()
            logmsg(@sprintf("trạng thái: %d đang chạy, %d chờ, RAM ước lượng %.1f GB, RAM dùng %.1f GB",
                            length(running), length(jobs), used, mem_used_gb(running)))
        end
        sleep(2)
    end
    logmsg("tất cả job đã xong — tổng hợp kết quả")
    run_cmd(`$JULIA --project=$ROOT $(joinpath(ROOT, "aggregate.jl"))`, joinpath(RES, "logs", "aggregate.log"))
    logmsg("HOÀN TẤT")
end

main()
