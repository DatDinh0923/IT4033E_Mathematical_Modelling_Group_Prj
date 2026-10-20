# Project Orange — Flow A→Z, giải thích output, đánh giá kết quả và quy trình suy luận

> File này đi kèm `GIAI_THICH.md` (kiến thức nền + mô hình toán). Ở đây tập trung vào 4 câu hỏi:
>
> 1. **Code chạy thế nào từ đầu tới cuối?** (mục 1 — lần theo một instance thật, setA-01, với số liệu thật)
> 2. **Mỗi file output nghĩa là gì, đọc thế nào?** (mục 2)
> 3. **Kết quả đạt được đến đâu, tốt hay chưa, vì sao?** (mục 3)
> 4. **Đã suy luận thế nào để ra các quyết định và kết luận?** (mục 4)
>
> Mọi con số trong file đều lấy từ `orange/results/` và đã được kiểm tra lại bằng
> `julia --project=. check_solution.jl` (kết quả: **44/44 nghiệm khớp với bảng**).

---

## Mục lục

1. [Flow A→Z: lần theo setA-01 từ file JSON tới con số trong bảng](#1-flow-az)
2. [Giải thích output đầu ra](#2-giải-thích-output-đầu-ra)
3. [Đánh giá kết quả đạt được](#3-đánh-giá-kết-quả-đạt-được)
4. [Quy trình suy luận](#4-quy-trình-suy-luận)
5. [Tóm tắt một trang](#5-tóm-tắt-một-trang)

---

## 1. Flow A→Z

### 1.0 Bức tranh tổng thể

```
                    ┌──────────────────────── ./run.sh ─────────────────────────┐
                    │ cgroup 15 GB  →  run_all.jl (bộ điều phối)                │
                    └───────────────────────────────────────────────────────────┘
  INPUT                   GIAI ĐOẠN                                        OUTPUT
  ─────                   ─────────                                        ──────
  toy/*.json      ──►  [A] test_toy.jl: 39 kiểm tra trên ví dụ của đề  ──► logs/test_toy.log (dừng nếu sai)
  setA-XX-*.json  ──►  [B] prepare.jl ×20: đọc dữ liệu, tính r, baseline,
                           cận dưới MCF, kích thước mô hình             ──► instances/setA-XX.json
                       [C] run_all: lập 40 job, chọn exact/restricted,
                           xếp lịch theo RAM                            ──► run_all.log
                       [D] run_one.jl ×40 (mỗi job một process):
                           D1 đọc JSON → đồ thị G(0), G(1)
                           D2 tính r(i,j,a,t)
                           D3 baseline (không waypoint)
                           D4 [exact] pha 1: mô hình rút gọn → warm start
                           D5 build MIP (1)(2)(3)[(4)]
                           D6 HiGHS optimize!                           ──► highs_logs/stepX/…log
                           D7 decode x → waypoint, kiểm tra độc lập     ──► solutions/stepX/…-srpaths.json
                           D8 ghi kết quả                               ──► raw/stepX/setA-XX_s1.json
                       [E] aggregate.jl: gom 40 file raw                ──► table_step6/9.{md,csv},
                                                                            all_runs_step6/9.csv
                       [F] check_solution.jl (tùy chọn): đọc lại
                           srpaths, tính lại MLU, so với bảng           ──► "TẤT CẢ NGHIỆM KHỚP VỚI BẢNG"
```

Dưới đây lần theo **setA-01, step 9** (t ∈ {0,1} + budget) — instance nhỏ nhất, đủ để thấy mọi bước.

### 1.1 [Input] Ba file JSON

| File | Nội dung (setA-01) | Dùng làm gì |
|---|---|---|
| `setA-01-net.json` | 20 node (`id` 0…19), 80 link có `from`, `to`, `metric`, `capacity` | đồ thị G, độ dài ω(a), capacity c(a) |
| `setA-01-tm.json` | `num_time_slots = 2`; 40 demand, mỗi demand `s`, `t`, `v = [ν(d,0), ν(d,1)]` (vd. demand 0: s=7, t=8, v=[0.69, 46.9]) | lưu lượng ν(d,t) |
| `setA-01-scenario.json` | `max_segments = 6`, `budget = [{t:1, value:51}]`, `interventions = [{t:1, links:[22]}]` | maxSeg, κ(1) = 51, q(1) = {link 22} |

Lưu ý: lưu lượng ở t = 1 **khác hẳn** t = 0 (demand 0: 0,69 → 46,9), nên t = 1 không chỉ khác ở link bị tắt.

### 1.2 [A] Test trên instance toy (`test_toy.jl`)

Trước khi đụng vào setA, code chạy toàn bộ pipeline trên ví dụ của đề (7 node, 22 link, 2 demand)
và so với con số tính tay. Ví dụ những gì được kiểm tra:

- r theo `calculating_r.pdf` cho demand (v0, v5): cung (v1,v3) = 1/3, (v1,v4) = 2/3, (v3,v5) = 2/3
  (có 3 đường ngắn nhất, 1 đi qua v3, 2 đi qua v4).
- Hàm tính tải khớp đúng Table 2 và Fig. 4 của đề khi dùng ECMP (MLU 0,4375 → 0,375).
- Mô hình (cả 2 cách viết ràng buộc tải × 2 cách viết budget) giải ra tối ưu 0,375 < baseline 0,4167,
  nghiệm: demand (v0,v5) đi qua waypoint v2 (step 6) / v4 (step 9).
- Budget κ = 0 ⇒ routing t = 1 phải trùng t = 0.

Sai bất kỳ điều nào → `run_all.jl` dừng ngay ("test toy thất bại").

### 1.3 [B] Prepare (`prepare.jl setA-01`) → `results/instances/setA-01.json`

Mục đích: biết instance "to cỡ nào" và có các mốc so sánh **trước khi** giải.

| Bước | Hàm (file:dòng) | Kết quả setA-01 |
|---|---|---|
| đọc 3 JSON | `load_instance` (`src/orange.jl:230`) → `json_to_graph` (`:34`) | n = 20, 80 cung, 40 demand |
| G(0), G(1) | `graph_at` (`:75`): bỏ link trong `interventions` | G(1) có 79 cung (bỏ link 22) |
| d, σ | `all_pairs_shortest_data` (`:94`, chép từ `calculating_r.pdf`): Dijkstra `allpaths=true` từ mọi đỉnh | ma trận 20×20 khoảng cách và số đường ngắn nhất |
| r | `compute_r` (`:163`): r = σ_iu·σ_vj/σ_ij nếu d_iu+ω_a+d_vj = d_ij | 380 cặp (i,j), tổng 935 phần tử khác 0 ("ΣFG") |
| kiểm tra r | `check_r_conservation` (`:193`): r(i,j,·) phải là luồng đơn vị i→j | sai số 0 (với metric không nguyên — setA-08, 14 — Graphs.jl đếm thiếu, code tự chuyển sang Dijkstra có tolerance) |
| baseline | `arc_loads` (`:273`) với mọi demand đi thẳng s→t | MLU t=0 = 1,0000007; t=1 = 0,566 |
| cận dưới MCF | `mcf_lower_bound` (`:654`): LP luồng đa hàng hóa, chia luồng tùy ý | t=0: 0,5604; t=1: 0,2736 ⇒ step 9: 0,5604 |
| kích thước | số biến x = (#demand có volume > 0) × #cặp × #chu kỳ | step 6: 14 820; step 9: 30 400 |

### 1.4 [C] Lập kế hoạch (`run_all.jl`, `choose_mode` dòng 127)

- Số biến x ≤ 2,5 triệu ⇒ **exact** (mô hình đầy đủ). setA-01 step 9: 30 400 ⇒ exact.
- Lớn hơn ⇒ **restricted** với K waypoint ứng viên/demand, K = ⌊√(600 000 / 2|D|)⌋ − 1 ∈ [3, 30],
  dùng chung cho step 6 và 9 (vd. setA-20: |D| = 6000 ⇒ K = 6).
- Mỗi job có ước lượng RAM; job chỉ được khởi động khi tổng RAM (đo trong cgroup) còn đủ. Ghi vào
  `run_all.log`: `CHẠY setA-01 step9 s1 [exact K=0] ước lượng 1.6 GB (đang chạy 3, RAM dùng 4.8 GB)`.

### 1.5 [D] Một lần giải (`run_one.jl setA-01 9 1 exact 3600`)

**D0 — khởi động JIT**: chạy toy 2 lần (exact + restricted) để Julia biên dịch xong mọi hàm; nhờ đó
CPU time đo được là của *giải*, không lẫn thời gian biên dịch.

**D1–D3 — dữ liệu và baseline**: như prepare (mục 1.3), mất 0,0007 s (`time_data`).

**D4 — pha 1 (chỉ exact): tìm warm start.** Giải nhanh mô hình rút gọn K = 20 (≤ 10 % time limit =
360 s), thu được routing có MLU **0,9294** sau 6,1 s CPU (`warm_cpu_time`). Routing này khả thi cho
mô hình đầy đủ (mọi segment của nó đều có trong mô hình đầy đủ).

**D5 — build mô hình** (`build_model`, `src/orange.jl:369`), 1,6 s (`time_build`). Mô hình setA-01
step 9 có **46 401 biến**, **17 817 ràng buộc** — đếm cụ thể:

| Thành phần | Số lượng | Giải thích |
|---|---|---|
| x^{dt}_{ij} (nhị phân) | 30 400 | 40 demand × 380 cặp × 2 chu kỳ |
| U | 1 | MLU |
| y^t_{ij} (liên tục) | 760 | luồng tổng trên segment (i,j) — 380 × 2 (dạng gộp của ràng buộc tải) |
| z (liên tục) | 15 200 | z ≥ x¹ − x⁰ cho budget — 40 × 380 |
| b_d (liên tục) | 40 | số segment thay đổi của từng demand |
| **Tổng biến** | **46 401** | |
| (1) bảo toàn luồng | 1 600 | 40 demand × 20 đỉnh × 2 chu kỳ |
| (2) số segment ≤ 6 | 80 | 40 × 2 |
| định nghĩa y | 760 | y = Σ_d ν x |
| (3) tải ≤ U·c(a) | 136 | chỉ cho cung nằm trên ít nhất một forwarding graph (12 trong 80 cung của setA-01 không bao giờ nằm trên đường ngắn nhất nào → bỏ) |
| z ≥ x¹ − x⁰ | 15 200 | |
| b_d = Σ(2z − x¹ + x⁰) | 40 | budget tách theo demand |
| Σ b_d ≤ 51 | 1 | budget (4) |
| **Tổng ràng buộc** | **17 817** | |

**D6 — giải** (`solve!`, `:540`). Log HiGHS `highs_logs/step9/setA-01_s1.log` (rút gọn):

```
MIP has 17817 rows; 46401 cols; 215280 nonzeros; 30400 integer variables
Presolve reductions: rows 17681(-136); columns 46265(-136)       ← presolve bỏ 136 hàng/cột thừa
MIP start solution is feasible, objective value is 0.929383899821 ← nạp warm start từ pha 1
   Nodes |  BestBound       BestSol   Gap  | ...  Time
       0 |  0               0.9293838  100% |       1.4s   ← chưa giải LP: cận dưới = 0
       1 |  0.9293838998    0.9293838  0.00%|       8.5s   ← LP gốc (IPM) cho bound = UB ⇒ xong
Status  Optimal      Nodes 1      Timing 8.52
```

Đọc: cận trên (BestSol) có ngay từ warm start; cận dưới (BestBound) từ LP relaxation ở nút gốc bằng
đúng cận trên ⇒ gap = 0 ⇒ **chứng minh tối ưu** mà không cần rẽ nhánh (Nodes = 1). CPU 8,57 s.

**D7 — decode và kiểm tra độc lập.** `decode_paths` (`:580`) đọc các x = 1, đi BFS từ s tới t ⇒
segment path cho từng demand. `verify_solution` (`:614`) kiểm tra *không dùng mô hình*: đường đơn,
≤ 6 segment, mọi segment tồn tại trong G(t), budget; và tính lại MLU bằng `arc_loads`:
t = 0: **0,9294**, t = 1: **0,6385**, số thay đổi = **4 ≤ 51**. Nếu MLU tính lại khác objective, job
báo lỗi và **không** ghi kết quả.

**D8 — ghi output**: `raw/step9/setA-01_s1.json` (mục 2.3) và `solutions/step9/setA-01-srpaths.json`
(mục 2.4).

### 1.6 [E] Tổng hợp (`aggregate.jl`)

Đọc 40 file raw + 20 file instances → một dòng/instance: `table_step9.md`:

```
| setA-01 | 20 | 80 | 40 | exact | 1 | Optimal×1 | 0.9294 | 0.9294 | 0.00% | 8.6 | 1 | 1.0000 | 0.9294 | 7.1% | 0.5604 | 39.7% | 4.0 / 51 |
```

### 1.7 [F] Kiểm tra lại từ file nghiệm (`check_solution.jl`)

Script độc lập đọc **file nghiệm theo định dạng của đề** (không đọc mô hình), dựng lại đồ thị, r, tính
lại MLU và so với bảng. `julia --project=. check_solution.jl setA-01 9`:

```
setA-01  step9  MLU tính lại = 0.929384 | bảng = 0.929384  KHỚP | waypoint = 114 | thay đổi = [4]
  t = 0: MLU = 0.9294; 3 link tải cao nhất:
     link id 29    (v7 -> v14)  tải = 0.9294
     link id 32    (v16 -> v3)  tải = 0.9065
     link id 42    (v1 -> v4)   tải = 0.7428
  t = 0: số demand dùng 0/1/2/3+ waypoint = [9, 12, 13] / 6
  t = 1: MLU = 0.6385; ...
```

⇒ setA-01: SR giảm MLU từ 1,0000 xuống 0,9294 (−7,1 %); nút thắt là link 29 (v7→v14) ở t = 0;
31/40 demand dùng ít nhất một waypoint.

### 1.8 Khác biệt khi instance lớn (restricted), ví dụ setA-20

- Không có pha 1; warm start = baseline (mọi demand đi thẳng).
- `candidate_segments(:restricted, K=6)`: mỗi demand chỉ được dùng 6 waypoint "ít vòng nhất" ⇒
  893 911 biến thay vì ~1,9 tỉ.
- LP gốc giải bằng dual simplex; ở setA-20 nghiệm LP đã nguyên ⇒ "Optimal" sau 62 s, MLU 0,9913.
- "Optimal" này là tối ưu **trong không gian K = 6**, không phải của bài gốc (mục 3.2).

---

## 2. Giải thích output đầu ra

### 2.1 Cây thư mục `orange/results/`

```
results/
├── table_step6.md / .csv        ★ BẢNG KẾT QUẢ CHÍNH step 6 (một dòng / instance)
├── table_step9.md / .csv        ★ BẢNG KẾT QUẢ CHÍNH step 9
├── all_runs_step6.csv / 9.csv     mọi lần chạy, mọi cột (nếu chạy 10 seed sẽ có 10 dòng/instance)
├── instances/setA-XX.json         thông tin instance từ prepare (mục 2.2)
├── raw/step6|step9/setA-XX_s1.json   kết quả chi tiết từng lần giải (mục 2.3) — nguồn của mọi bảng
├── raw/step6_ipm|step9_ipm/       thử nghiệm bổ sung IPM cho setA-06/08 (không vào bảng chính)
├── solutions/step6|step9/setA-XX-srpaths.json   NGHIỆM: waypoint của từng demand (mục 2.4)
├── highs_logs/step6|step9/setA-XX_s1.log        log HiGHS của đúng lần giải trong bảng (mục 2.5)
├── highs_logs_full/…              log gốc: HiGHS ghi NỐI TIẾP qua các lần chạy lại (xem 2.5)
├── run_all.log / run_all.out      nhật ký bộ điều phối (mục 2.6)
├── ram_kills.txt                  (instance, step) từng bị kill vì RAM → lần sau dùng mô hình nhỏ hơn
├── logs/                          log từng process: test_toy, prepare_*, runs/*, aggregate, ipm_extra
└── run_ipm_extra.sh               script của thử nghiệm bổ sung IPM
```

### 2.2 `instances/setA-XX.json` (từ prepare)

| Trường | Nghĩa |
|---|---|
| `n`, `arcs`, `demands` | số đỉnh, số link, số demand |
| `zero_volume_t0` | số demand có ν(d,0) = 0 (bị bỏ khỏi mô hình step 6 — không ảnh hưởng tải) |
| `maxseg`, `budget_t1`, `down_links_t1` | maxSeg, κ(1), id các link bị tắt ở t = 1 |
| `pairs_t0/t1` | số cặp (i,j) có đường đi ở G(t) — số segment khả dĩ |
| `sum_fg_t0/t1` | tổng số phần tử khác 0 của r (= Σ kích thước forwarding graph) |
| `r_conservation_err_t0/t1` | sai số bảo toàn luồng của r (phải ≈ 0) |
| `r_err_graphsjl_t0`, `tolerant_dijkstra` | sai số nếu chỉ dùng Graphs.jl; `true` = đã phải dùng Dijkstra có tolerance |
| `baseline_mlu_t0/t1/t01` | MLU khi không dùng waypoint (t01 = max hai chu kỳ = baseline của step 9) |
| `mcf_lb_t0/t1/t01` | cận dưới MCF (null nếu LP không giải xong) |
| `x_vars_step6/9` | số biến x của mô hình đầy đủ |
| `nnz_load_direct_t0`, `nnz_load_aggregated_t0` | số hệ số ràng buộc tải nếu viết thẳng / dạng gộp |
| `time_r_s`, `time_mcf_s`, `maxrss_gb` | thời gian tính r, giải MCF, RAM đỉnh của prepare |

### 2.3 `raw/stepX/setA-XX_s1.json` (một lần giải)

| Nhóm | Trường | Nghĩa |
|---|---|---|
| định danh | `instance`, `step`, `seed`, `time_limit` | |
| mô hình | `mode` (`exact`/`restricted`), `K` (0 = exact), `lp_solver` (`ipm`/`simplex`) | cấu hình đã dùng |
| kích thước | `n`, `arcs`, `demands`, `demands_in_model`, `nvars`, `ncons`, `maxseg`, `budget` (−1 ở step 6) | |
| solver | `status` (`Optimal` / `Not Optimal` / `Infeasible` / `No solution`), `termination` (mã gốc của JuMP: `OPTIMAL`, `TIME_LIMIT`) | |
| | `best_obj` | MLU của nghiệm tốt nhất (cận trên) |
| | `best_bound` | cận dưới của **mô hình đã giải** (null nếu chưa có) |
| | `gap` | 1 − bound/obj (readme); bound null ⇒ 1 |
| | `nodes` | số node branch-and-bound (0 = chưa xong nút gốc; 1 = xong ở nút gốc) |
| | `cpu_time` | **CPU time của `optimize!`** (con số readme yêu cầu); `wall_time` = thời gian thực |
| ngoài optimize | `time_data` (tính r), `time_build` (build + warm start), `warm_obj`, `warm_cpu_time` (pha 1 của exact) | |
| chất lượng | `baseline_mlu`, `improvement` = 1 − best_obj/baseline | |
| kiểm tra | `mlu_check` (MLU tính lại độc lập), `mlu_per_period`, `changes` (dist(P⁰,P¹) thực tế), `waypoints` (tổng số waypoint) | |
| tài nguyên | `maxrss_gb` | RAM đỉnh của process |

### 2.4 `solutions/stepX/setA-XX-srpaths.json` (nghiệm, đúng định dạng `toy-srpaths.json` của đề)

```json
{"srpaths": [ {"d": 0, "t": 0, "w": [11, 10, 5]},
              {"d": 1, "t": 0, "w": [9]}, ... ]}
```

- `d` = chỉ số demand (0-based, theo thứ tự trong `tm.json`), `t` = chu kỳ, `w` = danh sách waypoint
  theo thứ tự (node id gốc trong `net.json`).
- Demand **không xuất hiện** ở chu kỳ t = đi thẳng s → t (không waypoint).
- Ví dụ trên (setA-01 step 9): demand 0 (s = 7, t = 8) ở t = 0 đi v7 → v11 → v10 → v5 → v8: 3 waypoint,
  4 segment ≤ 6. Giữa hai waypoint, gói tin đi theo các đường ngắn nhất (chia luồng theo r).

### 2.5 `highs_logs/…log` — đọc log HiGHS

```
MIP has R rows; C cols; NZ nonzeros; B integer variables      kích thước mô hình
Presolve reductions: rows R'(−a); columns C'(−b)               sau rút gọn
MIP start solution is feasible, objective value is X           warm start được chấp nhận
   Src Proc. InQueue | Leaves Expl. | BestBound BestSol Gap | Cuts InLp Confl. | LpIters Time
```

| Cột | Nghĩa |
|---|---|
| `Src` | nguồn của nghiệm mới: `H` heuristic, `L` sub-MIP, `C` central rounding, `T` node, trống = dòng tiến độ |
| `Proc.` / `InQueue` | số node đã xử lý / đang chờ |
| `BestBound` / `BestSol` / `Gap` | cận dưới / cận trên / gap hiện tại |
| `Cuts`, `InLp` | số lát cắt đã sinh / đang nằm trong LP |
| `LpIters`, `Time` | tổng vòng lặp simplex, thời gian |

Dòng cuối có `Status` (`Optimal` / `Time limit reached`), `Nodes`, `Timing … (Presolve) … (Solve)`.

Dấu hiệu cần chú ý khi đọc log:
- `BestBound` = 0 suốt và `Proc.` = 0 tới hết giờ ⇒ **LP gốc chưa giải xong** (setA-06, 08).
- `BestSol` không đổi so với dòng đầu ⇒ không tìm được nghiệm tốt hơn warm start.

> **Lưu ý**: HiGHS ghi **nối tiếp** vào file log khi một job được chạy lại. Thư mục `highs_logs/` đã
> được lọc để mỗi file chỉ còn đúng khối của lần giải nằm trong bảng (khớp theo số rows/cols và Status).
> Bản gốc đầy đủ ở `highs_logs_full/` (vd. `step9/setA-20_s1.log` có 2 khối: khối 1 là lần simplex
> được giữ, khối 2 là lần thử IPM đã bị hủy).

### 2.6 `run_all.log` — nhật ký bộ điều phối

| Dòng | Nghĩa |
|---|---|
| `test toy: OK`, `prepare: xong 20 instance` | giai đoạn A, B |
| `plan step9 setA-20: restricted K=6 ~5.88e+05 biến x ước lượng RAM 6.5 GB` | kế hoạch cho từng (instance, step) |
| `CHẠY setA-20 step9 s1 [restricted K=6] ước lượng … (đang chạy N, RAM dùng X GB)` | khởi động một job |
| `XONG setA-20 step9 s1 [restricted] Optimal obj=0.99131 gap=0.0000 nodes=1 cpu=62s rss=7.1GB` | job xong |
| `RAM cao (13.8 GB đang dùng): kill setA-04 step9 s1 (rss 9.0 GB)` | watchdog giết job lớn nhất |
| `THẤT BẠI … (lần 1, exit=0)` + `-> chạy lại với mô hình nhỏ hơn: restricted K=20` | job hỏng/bị kill và được chạy lại |
| `đã từng bị kill vì RAM -> setA-07 step9: restricted K=11` | đọc `ram_kills.txt` khi resume |
| `ĐANG CHẠY từ trước, bỏ qua: …` | job của lần điều phối trước vẫn còn chạy, không chạy trùng |
| `trạng thái: 2 đang chạy, 13 chờ, RAM ước lượng …, RAM dùng …` | báo cáo mỗi 15 phút |
| `=== RESTART n ===` (trong `run_all.out`) | các lần khởi động lại bộ điều phối (mục 4.2) |
| `HOÀN TẤT` | đã tổng hợp xong |

### 2.7 Các cột của bảng `table_step6/9.md`

| Cột | Nghĩa | Cách hiểu |
|---|---|---|
| Mô hình | `exact` / `restricted(K=…)` | quyết định ý nghĩa của Status và Gap |
| Status | theo HiGHS | `Optimal` của restricted = tối ưu trong không gian rút gọn |
| Best obj | MLU tốt nhất tìm được | càng nhỏ càng tốt; < 1 = hết nghẽn |
| Best bound, Gap | cận dưới và gap readme | chỉ là cận của **mô hình đã giải** |
| CPU time | của `optimize!` | |
| Nodes | node B&B | 1 = xong ở nút gốc; 0 = chưa xong nút gốc |
| MLU baseline | không waypoint | mốc so sánh |
| Warm start | MLU sau pha 1 (chỉ exact) | |
| Giảm MLU | 1 − Best obj / baseline | lợi ích của Segment Routing |
| Cận dưới MCF | cận dưới **hợp lệ cho bài gốc** | nếu = Best obj ⇒ tối ưu bài gốc |
| Gap theo MCF | 1 − MCF / Best obj | giới hạn trên của gap thật (thường bi quan) |
| Thay đổi / κ (step 9) | dist(P⁰,P¹) / budget | luôn ≤ κ |

---

## 3. Đánh giá kết quả đạt được

### 3.1 So với yêu cầu của readme

| Yêu cầu | Trạng thái | Ghi chú |
|---|---|---|
| Step 2: graph G(0) từ JSON (script Pluto) | ✅ | `json_to_graph` chuyển từ notebook v2; bỏ lỗi `multigraph=true` (không có cung song song thật) |
| Step 3: r(i,j,a,0) theo `calculating_r.pdf` | ✅ | đúng công thức; thêm tolerance cho metric không nguyên; kiểm tra bảo toàn luồng = 0 |
| Step 4: Data section (r, c, ν, maxSeg) | ✅ | `load_instance`, `period_data` |
| Step 5: mô hình t = 0, không budget | ✅ | test toy đúng tối ưu |
| Step 6: chạy setA, 1800 s, bảng | ✅ | 20/20 instance; bảng có đủ cột readme yêu cầu |
| Step 7: G(1), r(·,1) | ✅ | bỏ link trong `interventions` |
| Step 8: thêm budget (4) | ✅ | tuyến tính hóa chính xác |
| Step 9: chạy lại, 3600 s | ✅ | 20/20 instance |
| "Average on ten runs" | ⚠️ 1 lần | theo yêu cầu của bạn; `ORANGE_NSEEDS=10 ./run.sh` để chạy 10 lần |
| Mô hình đúng trên mọi instance | ⚠️ | exact chỉ cho 01–05; 06–20 dùng mô hình rút gọn (bắt buộc vì 15 GB RAM) |

### 3.2 Phân loại 40 kết quả theo độ tin cậy

```
                        ┌─ exact + Optimal ─────────────────────────► TỐI ƯU BÀI GỐC (chứng minh bởi HiGHS)
Kết quả một lần giải ───┼─ Best obj = cận MCF ─────────────────────► TỐI ƯU BÀI GỐC (kẹp giữa 2 cận)
                        ├─ step 9: Best obj = tối ưu exact step 6 ──► TỐI ƯU BÀI GỐC (step 6 là nới lỏng của step 9)
                        ├─ restricted + Optimal, obj > cận MCF ─────► tối ưu trong không gian rút gọn;
                        │                                              là nghiệm tốt, chưa biết cách tối ưu bao xa
                        └─ Not Optimal ─────────────────────────────► nghiệm khả thi + gap
```

| Nhóm | Step 6 | Step 9 |
|---|---|---|
| **A. Tối ưu cho bài gốc** (chứng minh được) | 01, 02, 03, 07, 15, 17, 18 (7) | 01, 03, 05, 07, 09, 11, 15, 17, 18 (9) |
| **B. Tối ưu trong không gian rút gọn** | 09, 10, 11, 12, 13, 14, 16, 19, 20 (9) | 04, 10, 12, 13, 14, 16, 19, 20 (8) |
| **C. Có cải thiện nhưng chưa chứng minh tối ưu** | 05 (gap 12,8 %) | 02 (30,9 %) |
| **D. Không cải thiện được baseline** | 04, 06, 08 (solver) | 06, 08 (solver) |

(07 và 18 ở step 6, 18 ở step 9 nằm ở nhóm A nhưng "giảm MLU 0 %": baseline **đã là tối ưu**.)

### 3.3 Đánh giá từng instance

| Inst | Step 6 | Step 9 | Nhận xét |
|---|---|---|---|
| 01 | 0,9294 tối ưu (−7 %) | 0,9294 tối ưu (−7 %) | t = 0 là nút thắt; cận MCF 0,56 rất lỏng |
| 02 | 0,5486 tối ưu (−45 %) | 0,9031 (−18 %), gap 31 % | step 9 khó hơn nhiều: t = 1 baseline 1,098 (> 1: nghẽn) |
| 03 | 0,9435 tối ưu (−6 %) | 0,9435 (−6 %) ⇒ **tối ưu** dù HiGHS báo gap 33 % | bỏ chu kỳ t = 1 của step 9 được đúng bài step 6 ⇒ opt(step 9) ≥ opt(step 6) = 0,9435 đã chứng minh; step 9 đạt 0,9435 ⇒ tối ưu |
| 04 | 0,6918 = baseline, gap 34 % | 0,7421 (−26 %) | step 6: solver không tìm được nghiệm (step 9 cho thấy t = 0 đạt được 0,643) |
| 05 | 0,1324 (−60 %), gap 13 % | 0,2050 (−80 %) = cận MCF ⇒ tối ưu | giảm mạnh nhất; κ = 1 nhưng dùng 0 thay đổi |
| 06 | 1,0 — LP gốc không xong | 1,0 — LP gốc không xong | thất bại của solver; IPM cho bound 0,13 / 0,25 nhưng không có nghiệm |
| 07 | 0,9080 = baseline = MCF | 0,9080 (−9 %) = MCF | baseline t = 0 đã tối ưu; step 9 tối ưu |
| 08 | 1,0 — LP gốc không xong | 1,0 — LP gốc không xong | IPM: step 6 **tối ưu 0,3189 (−68 %)** |
| 09 | 0,5054 (−17 %) | 0,8497 = MCF ⇒ tối ưu (−15 %) | |
| 10 | 0,4239 (−58 %) | 0,4239 (−58 %) | cận MCF 0,07 quá lỏng để kết luận |
| 11 | 0,6446 (−36 %), 203 node | 0,7858 = MCF ⇒ tối ưu (−21 %) | instance hiếm hoi phải rẽ nhánh |
| 12 | 1,0001 = baseline | 1,0001 = baseline | **mô hình rút gọn thiếu waypoint v161**; có nghiệm 0,9898 (mục 3.4) |
| 13 | 0,6197 (−38 %) | 0,6497 (−35 %) | budget dùng 9/12 — gần chặn |
| 14 | 0,4397 (−56 %) | 0,5331 (−47 %) | cách cận MCF chỉ 3,6 % / 4,8 % |
| 15 | 0,8252 = MCF ⇒ tối ưu (−17 %) | 0,8987 = MCF ⇒ tối ưu (−10 %) | |
| 16 | 0,8730 (−10 %) | 0,8811 (−12 %) | K = 6 nhỏ (4800 demand) |
| 17 | 0,4242 = MCF ⇒ tối ưu (−58 %) | 0,4242 = MCF ⇒ tối ưu (−58 %) | |
| 18 | 0,6334 = baseline = MCF | 1,0000 = baseline = MCF | nguồn v138 chỉ có 1 link ra ⇒ không thể tốt hơn |
| 19 | 0,7377 (−26 %) | 0,7377 (−26 %) | |
| 20 | 0,9913 (−1 %) | 0,9913 (−1 %) | K = 6 nhỏ; cận MCF 0,648 ⇒ có thể còn cải thiện |

### 3.4 Vì sao có những kết quả đó — giải thích theo nguyên nhân

**(1) Vì sao SR giảm được MLU (trung bình ~21 %, tới 80 %).** Baseline đi theo đường ngắn nhất ⇒
nhiều demand dồn vào cùng vài link "trung tâm". Waypoint buộc demand đi vòng qua một đỉnh khác ⇒ dùng
các đường song song đang trống. Instance có nhiều đường song song và lưu lượng tập trung (05, 10, 14,
17) được lợi nhiều nhất. Dữ liệu được sinh sao cho baseline ≈ 1,0000 (vừa chạm 100 %), nên mọi giá trị
< 1 nghĩa là SR gỡ được nghẽn.

**(2) Vì sao có instance mà baseline đã tối ưu (07, 18).** Nút thắt là **giới hạn vật lý**: ở setA-18,
nguồn v138 chỉ có một link đi ra (v138→v90); mọi demand từ v138 bắt buộc qua link đó dù chọn waypoint
nào. Cận MCF (chia luồng tùy ý) cũng bằng baseline ⇒ không phương pháp định tuyến nào làm tốt hơn.

**(3) Vì sao mô hình rút gọn có thể "Optimal" mà vẫn chưa tốt (12).** Ứng viên waypoint được chọn theo
"độ vòng" d(s,w)+d(w,t)−d(s,t). Ở setA-12, link còn trống đi vào v175 là (v161→v175) có metric 0,97
(gấp ~27 lần link thường) ⇒ v161 có độ vòng rất lớn, xếp hạng 161/198, bị loại khỏi K = 26 ứng viên.
Kiểm tra tay: đưa 2 demand qua v161 cho MLU 0,9898 < 1,0001. Đây là **giới hạn của heuristic chọn
ứng viên**, không phải của mô hình toán.

**(4) Vì sao solver thất bại ở 06, 08 (và 04 step 6).** Mô hình rút gọn với K lớn (23–30) cho mỗi
demand ~600–960 segment; LP relaxation rất lớn và **suy biến** (nhiều nghiệm tối ưu tương đương), dual
simplex lặp mãi không hội tụ ⇒ hết giờ vẫn ở nút gốc, không có bound, heuristic không có thông tin ⇒
giữ nguyên warm start = baseline. Bằng chứng: chỉ đổi sang IPM thì setA-08 step 6 giải tối ưu trong
172 s.

**(5) Vì sao exact cỡ vừa (02, 03 step 9; 04, 05) có gap lớn.** 100 000 – 3 000 000 biến nhị phân:
LP gốc chỉ vừa giải xong (IPM) là gần hết giờ ⇒ Nodes = 0, không kịp rẽ nhánh. Gap theo HiGHS lớn nhưng
nghiệm thực tế có thể đã tối ưu (05 step 9: gap 100 % nhưng obj = cận MCF).

**(6) Vì sao cận MCF đôi khi rất lỏng (01, 10, 13, 16, 19).** MCF cho chia luồng của *mỗi demand* tùy
ý trên mọi đường; SR chỉ chia theo ECMP/đường ngắn nhất của tối đa 5 waypoint. Với nhiều demand nhỏ
và mạng nhiều đường (|D| lớn), MCF trải luồng mỏng khắp mạng ⇒ cận rất thấp (0,04–0,07), không phản
ánh khả năng thật của SR. Vì vậy "gap theo MCF" chỉ là giới hạn trên (bi quan) của gap thật.

**(7) Vì sao step 9 ≥ step 6 và budget ít khi chặn.** Mô hình step 6 chính là phần t = 0 của step 9
(cùng K) ⇒ obj(step 9) ≥ obj(step 6) ở mọi instance. Chu kỳ t = 1 thường là nút thắt (link bị tắt +
lưu lượng khác). Budget κ thường lớn so với nhu cầu: một bộ waypoint tốt cho cả hai chu kỳ thường tồn
tại (nhiều instance dùng 0 thay đổi); chỉ setA-13 dùng 9/12.

### 3.5 Kết luận đánh giá

- **Đạt**: pipeline đầy đủ step 2–9, chạy hết 20 instance cho cả hai step trong 15 GB, mọi nghiệm hợp
  lệ và được kiểm tra độc lập; SR giảm MLU trung bình ~21 %; 16/40 kết quả được chứng minh tối ưu cho
  bài gốc.
- **Chưa đạt / hạn chế**: mô hình đầy đủ không mở rộng quá ~100 đỉnh trong 15 GB; mô hình rút gọn là
  heuristic (setA-12 cho thấy nó có thể bỏ sót nghiệm tốt); 2 instance (06, 08) thất bại vì thuật toán
  LP; chỉ 1 lần chạy nên không có độ biến thiên.
- **Hướng cải thiện có bằng chứng**: (i) LP gốc chạy simplex, nếu chậm thì chuyển IPM (08: 1,0 → 0,3189);
  (ii) chọn ứng viên waypoint có tính đến link còn trống (12: 1,0001 → ≤ 0,9898); (iii) warm start
  step 9 bằng nghiệm step 6; (iv) column generation thay cho mô hình đầy đủ cho mạng lớn.

---

## 4. Quy trình suy luận

### 4.1 Suy luận khi thiết kế (trước khi chạy)

| Câu hỏi | Suy luận | Quyết định |
|---|---|---|
| Objective lex-min của đề có giải được không? | cần |A|·|T| lần giải lặp, mỗi lần là một MIP lớn | chỉ min MLU (thành phần đầu của lex-min) |
| Có cần ràng buộc chống chu trình rời (subtour)? | r ≥ 0, ν ≥ 0 ⇒ bỏ chu trình chỉ giảm tải và số segment ⇒ luôn có nghiệm tối ưu không chu trình | không cần; decode bằng BFS bỏ chu trình |
| Segment (i,j) không có đường đi? | r = 0 mọi cung ⇒ "chở luồng miễn phí" — sai | cấm (fix x = 0) |
| Mô hình có vừa 15 GB? | x = |D|·n²·|T|; setA-20 step 9 ≈ 1,9 tỉ biến | exact chỉ khi ≤ 2,5 triệu; còn lại rút gọn |
| Ràng buộc tải quá dày? | |D|·Σ|FG| hệ số (setA-08: 23 triệu) | đặt y = Σ_d ν x (tương đương), còn 4,6 triệu |
| Làm sao biết nghiệm rút gọn tốt hay không? | bound của mô hình con không phải bound bài gốc | thêm cận dưới MCF độc lập |
| Làm sao tin con số? | objective của solver có thể sai do mô hình sai | test toy với số tính tay + tính lại MLU độc lập mỗi lần giải |

### 4.2 Suy luận khi chạy: quan sát → giả thuyết → thử nghiệm → quyết định

Đây là phần quan trọng nhất để hiểu vì sao code có dạng như hiện tại. Mỗi thay đổi đều xuất phát từ
một quan sát cụ thể và được kiểm chứng bằng thử nghiệm trước khi áp dụng.

**① Presolve treo.**
- *Quan sát*: setA-20 step 9 (time limit 60 s) chạy > 20 phút, log dừng ở "Presolving model".
- *Giả thuyết 1*: nhiều hàng budget → viết 1 hàng/cặp thay vì 2. *Thử*: vẫn treo ⇒ sai.
- *Giả thuyết 2*: "probing" (log gợi ý). *Thử*: tắt probing, setA-02 vẫn presolve 206 s ⇒ chưa đủ.
- *Giả thuyết 3*: **một hàng budget chứa mọi biến** (hàng dày) làm presolve chậm. *Thử*: tách
  b_d theo demand ⇒ presolve 5 s, setA-13 Optimal trong 38 s. ⇒ **Áp dụng**; bật lại probing (đo thấy
  không còn ảnh hưởng).

**② Hết RAM làm chết cả lượt chạy.**
- *Quan sát*: 13:30 mọi job nhận SIGTERM. Journal: "oom-kill" và scope "Failed".
- *Suy luận*: cgroup đã chặn đúng 15 GB, nhưng systemd mặc định dừng **cả scope** khi một process bị
  OOM; và ước lượng RAM mô hình rút gọn (2,7 GB) sai — đo thực tế 2,6–7,1 GB.
- *Quyết định*: `OOMPolicy=continue`; ước lượng bảo thủ; watchdog ở cấp slice; bị kill vì RAM thì chạy
  lại với mô hình nhỏ hơn và ghi nhớ vào `ram_kills.txt`.

**③ Mô hình đầy đủ không có tiến triển.**
- *Quan sát*: setA-03 step 9 hết 1 giờ: Nodes 0, bound 0, 752 481 vòng simplex.
- *Thử*: `mip_lp_solver=ipm` ⇒ LP gốc xong 36 s, bound 0,6337; nhưng UB vẫn = baseline.
- *Suy luận*: heuristic HiGHS không tìm được nghiệm trên mô hình lớn ⇒ cần nghiệm tốt cho sẵn. Mô hình
  rút gọn giải rất nhanh và nghiệm của nó khả thi cho mô hình đầy đủ.
- *Quyết định*: exact = pha 1 (rút gọn K=20) + pha 2 (IPM). Kết quả: setA-02 step 6 từ "LP gốc 191 s"
  thành Optimal trong 3 s.

**④ Mô hình rút gọn với K lớn bị kẹt; IPM không phải lúc nào cũng tốt.**
- *Quan sát*: setA-13 step 6 (K=16) cho 1,0 trong khi step 9 (K=11, khó hơn) cho 0,6496 — mâu thuẫn.
- *Suy luận*: step 6 được cho K lớn hơn (ngân sách biến chia cho 1 chu kỳ) ⇒ LP lớn hơn ⇒ simplex kẹt.
- *Thử IPM*: setA-08/13 step 6 tối ưu; **nhưng** setA-20 step 9 tệ hơn (simplex 62 s Optimal, IPM 20 phút
  chưa có nghiệm) — vì nghiệm đỉnh của simplex thường đã nguyên, còn IPM cho nghiệm trong miền.
- *Quyết định*: cùng K cho step 6 & 9 (hai bảng so sánh được, mô hình step 6 = phần t=0 của step 9),
  simplex cho mô hình rút gọn, IPM cho exact; IPM cho 06/08 báo cáo riêng như thử nghiệm bổ sung.

**⑤ Kết quả "không cải thiện" ở setA-12 — mô hình sai hay dữ liệu?**
- *Quan sát*: Optimal 1,0001 = baseline, trong khi cận MCF 0,58.
- *Điều tra*: tìm link nghẽn (v150→v175), các demand đi qua, link thay thế (v161→v175, metric 0,97).
- *Thử tay*: đưa 2 demand qua v161 ⇒ 0,9898 — tốt hơn. v161 xếp hạng 161/198 theo độ vòng.
- *Kết luận*: giới hạn của cách chọn ứng viên, không phải lỗi mô hình.

**⑥ Kiểm tra output có đúng là của lần chạy trong bảng.**
- *Quan sát*: log HiGHS bị ghi nối tiếp qua các lần chạy lại; khối cuối của `setA-20` step 9 là lần IPM
  đã bị hủy.
- *Quyết định*: lọc log theo (rows, cols, Status) khớp với file raw; viết `check_solution.jl` để tính
  lại MLU từ file nghiệm ⇒ 44/44 khớp.

### 4.3 Suy luận khi đọc kết quả — cây quyết định

```
Đọc một dòng của bảng:
1. Mô hình = exact?
   ├─ có, Status Optimal ─────────────► tối ưu bài gốc. Xong.
   └─ không / Not Optimal ──► 2
2. Best obj = Cận dưới MCF (gap theo MCF ≈ 0)?
   ├─ có ─────────────────────────────► tối ưu bài gốc (mọi routing đều ≥ cận MCF).
   └─ không ──► 2b
2b. (step 9) Best obj = MLU tối ưu ĐÃ CHỨNG MINH của step 6 (exact) cùng instance?
   ├─ có ─────────────────────────────► tối ưu bài gốc: bỏ chu kỳ t = 1 và budget của step 9 thì được
   │                                    đúng bài step 6, nên opt(step 9) ≥ opt(step 6) (vd. setA-03)
   └─ không ──► 3
3. Best obj < baseline?
   ├─ có ──► SR có hiệu quả; độ tốt chưa rõ:
   │          - restricted Optimal: tối ưu trong K ứng viên; tăng K / đổi cách chọn có thể tốt hơn
   │          - Not Optimal: xem Gap (HiGHS) — cận trên của khoảng cách tới tối ưu của mô hình đó
   └─ không (= baseline) ──► 4
4. Vì sao không cải thiện?
   ├─ Nodes = 0 & bound ≈ 0 ──────────► solver chưa giải xong LP gốc (06, 08) → thử IPM / mô hình nhỏ hơn
   ├─ baseline = cận MCF ─────────────► giới hạn vật lý của mạng (07, 18) → không thể tốt hơn
   └─ Optimal (restricted) nhưng MCF thấp ► ứng viên thiếu (12) → xem link nghẽn bằng check_solution.jl
```

---

## 5. Tóm tắt một trang

- **Input**: 3 JSON/instance (mạng, lưu lượng 2 chu kỳ, kịch bản tắt link + budget).
- **Xử lý**: đồ thị G(t) → đếm đường ngắn nhất → hệ số chia luồng r → MIP min MLU với biến x (segment
  nào được dùng) → HiGHS → decode ra waypoint → kiểm tra độc lập.
- **Output chính**: `results/table_step6.md`, `results/table_step9.md` (bảng), `results/solutions/…`
  (waypoint của từng demand), `results/raw/…` (mọi số đo).
- **Kết quả**: MLU giảm trung bình ~21 % (tới 80 %); 16/20 (step 6) và 15/20 (step 9) "Optimal" theo
  HiGHS; được chứng minh tối ưu cho bài gốc: 7 (step 6) và 9 (step 9); 2 instance thất bại do thuật toán LP.
- **Giới hạn**: 15 GB RAM ⇒ mô hình đầy đủ chỉ cho ≤ 100 đỉnh, còn lại là heuristic rút gọn; 1 lần chạy.
- **Kiểm chứng**: `julia --project=. test_toy.jl` (39/39) và `julia --project=. check_solution.jl`
  (44/44 nghiệm khớp với bảng).
