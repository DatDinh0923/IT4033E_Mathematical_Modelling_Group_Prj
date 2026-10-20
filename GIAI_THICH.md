# Project Orange (T-ASR) — Giải thích flow chạy, kết quả và kiến thức nền

> Xem thêm **`KET_QUA_VA_FLOW.md`**: flow A→Z lần theo một instance thật, giải thích từng file output,
> đánh giá chi tiết từng instance và quy trình suy luận.

> Code: thư mục `orange/` (Julia 1.13 + JuMP + HiGHS). Kế hoạch gốc: `project/PLAN.md`.
> Dữ liệu: 20 instance `project/setA/setA-XX-{net,tm,scenario}.json`.
> Điều kiện chạy: **toàn bộ quá trình bị giới hạn cứng 15 GB RAM** (cgroup của systemd) và
> **mỗi instance chỉ chạy 1 lần** (1 seed HiGHS) thay vì trung bình 10 lần như readme.

---

## Mục lục

1. [Bài toán nói gì (tóm tắt đề)](#1-bài-toán-nói-gì)
2. [Kiến thức nền](#2-kiến-thức-nền)
3. [Mô hình toán đã cài đặt](#3-mô-hình-toán-đã-cài-đặt)
4. [Cấu trúc code và flow chạy](#4-cấu-trúc-code-và-flow-chạy)
5. [Giới hạn 15 GB RAM được đảm bảo như thế nào](#5-giới-hạn-15-gb-ram)
6. [Kết quả](#6-kết-quả)
7. [Đánh giá và giải thích vì sao có kết quả đó](#7-đánh-giá-và-giải-thích)
8. [Những điểm kỹ thuật quan trọng đã phát hiện khi làm](#8-những-điểm-kỹ-thuật-quan-trọng)
9. [Cách chạy lại](#9-cách-chạy-lại)

---

## 1. Bài toán nói gì

Orange (nhà mạng) có một mạng IP/MPLS: các **router** (đỉnh) nối với nhau bằng các **link có hướng**
(cung). Mỗi link có:

- **metric ω(a)** — "độ dài" dùng để tính đường đi ngắn nhất (giao thức định tuyến IGP như OSPF/IS-IS dùng nó);
- **capacity c(a)** — băng thông tối đa.

Có một tập **demand** d = (s, t) với **lưu lượng ν(d, t)** thay đổi theo chu kỳ thời gian t
(t = 0: mạng bình thường; t = 1: có **intervention** — một số link bị tắt để bảo trì).

Mặc định, gói tin đi theo **đường ngắn nhất** s → t. Nếu có nhiều đường ngắn nhất bằng nhau thì
luồng được chia ra (ECMP). Vấn đề: đường ngắn nhất thường dồn nhiều luồng vào cùng vài link →
nghẽn. **Segment Routing (SR)** cho phép chèn thêm các **waypoint**: gói tin đi ngắn nhất từ s tới
w₁, rồi ngắn nhất từ w₁ tới w₂, …, rồi tới t. Như vậy ta "bẻ" luồng sang đường khác mà không phải
đổi metric của mạng.

**Mục tiêu**: chọn segment path (dãy waypoint) cho mọi demand, ở mọi chu kỳ, sao cho mạng ít nghẽn
nhất. Đề dùng lex-min của mọi tải; thành phần đầu tiên (quan trọng nhất) là
**MLU — Maximum Link Utilization** = tải lớn nhất trên mọi link. Ràng buộc:

1. mỗi demand có đúng một segment path từ s tới t;
2. số segment ≤ maxSeg (ở setA: maxSeg = 6, tức tối đa 5 waypoint);
3. tải mỗi link = (tổng lưu lượng qua nó) / capacity;
4. **budget**: số segment thay đổi giữa t−1 và t ≤ κ(t) (mỗi thay đổi tốn công vận hành).

Readme của thầy chia việc thành các bước:

| Step | Việc | Ở đâu trong code |
|---|---|---|
| 2 | `net.json` → đồ thị G(0) | `json_to_graph`, `graph_at` |
| 3 | tính r(i,j,a,0) theo `calculating_r.pdf` | `all_pairs_shortest_data`, `compute_r` |
| 4 | Data section: r, c(a), ν(d,0), maxSeg | `load_instance`, `period_data` |
| 5 | mô hình chỉ t = 0, không có budget (4) | `build_model(…, [0], …)` |
| 6 | chạy setA, time limit 1800 s | `run_one.jl … 6 …` |
| 7–8 | G(1) (bỏ link bị intervention), r(·,1), thêm budget (4) | `graph_at(net, down)`, `build_model(…, [0,1], …)` |
| 9 | chạy lại, time limit 3600 s | `run_one.jl … 9 …` |

---

## 2. Kiến thức nền

### 2.1 Đường đi ngắn nhất, đếm số đường ngắn nhất

Với đồ thị có trọng số dương, **Dijkstra** từ một nguồn s cho khoảng cách d(s, v) tới mọi đỉnh.
Thêm một nhãn σ(v) = số đường ngắn nhất s → v, cập nhật cùng lúc với d(v):

- tìm được đường **ngắn hơn hẳn** (ℓ < d(v)): d(v) ← ℓ, σ(v) ← σ(u);
- tìm được đường **bằng** (ℓ = d(v)): σ(v) ← σ(v) + σ(u);
- dài hơn: bỏ qua.

Chạy Dijkstra từ mọi đỉnh → ma trận d[i,j] và σ[i,j]. Độ phức tạp O(n·(m + n log n)).

### 2.2 Forwarding graph và hệ số chia luồng r

**Forwarding graph FG(i, j)** = mọi cung nằm trên ít nhất một đường ngắn nhất i → j.
Cung a = (u, v) thuộc FG(i, j) ⇔ d[i,u] + ω(a) + d[v,j] = d[i,j].

**Hệ số chia luồng r(i, j, a, t)** = tỉ lệ luồng đi từ i tới j chạy qua cung a. Theo
`calculating_r.pdf`:

```
r(i,j,a,t) = σ_iu · σ_vj / σ_ij    nếu d_iu + ω_a + d_vj = d_ij,   ngược lại 0
```

Lý do: một đường ngắn nhất i→j đi qua a = (u,v) = (một đường ngắn nhất i→u) + a + (một đường ngắn
nhất v→j), nên số đường như vậy là σ_iu·σ_vj. Chia cho tổng số đường σ_ij ra tỉ lệ đường đi qua a
(giả thiết: luồng chia **đều theo từng đường** ngắn nhất).

> **Lưu ý**: đề Orange dùng ECMP chia đều **tại từng nút** (Figure 1, Table 2 của đề), khác với
> cách chia đều theo từng đường của `calculating_r.pdf`. Ví dụ toy, demand (v0, v5): cung (v1,v3) có
> r = 1/2 theo ECMP nhưng 1/3 theo path-count. Code làm **đúng theo `calculating_r.pdf`** (yêu cầu
> của thầy); `test_toy.jl` có thêm hàm ECMP chỉ để kiểm chứng hàm tính tải khớp đúng con số của đề
> (MLU 0.4375 → 0.375).

**Tính chất kiểm tra được**: với mỗi cặp (i,j), r(i,j,·) là một luồng đơn vị i → j
(tổng ra − tổng vào = 1 tại i, −1 tại j, 0 nơi khác). Code kiểm tra điều này cho mọi cặp, mọi
instance (`check_r_conservation`), sai số = 0.

### 2.3 Tải và MLU

Tải của cung a ở chu kỳ t:

```
λ(a,t) = Σ_d Σ_(i,j) r(i,j,a,t) · ν(d,t) · [segment (i,j) thuộc path của d] / c(a)
```

MLU = max_a λ(a,t). MLU > 1 nghĩa là link bị quá tải (nghẽn).

### 2.4 Quy hoạch nguyên hỗn hợp (MIP) và Branch-and-Bound

Mô hình có biến nhị phân x ∈ {0,1} và biến liên tục U → **MILP**. Solver (HiGHS) giải bằng
**branch-and-bound**:

1. **Presolve**: rút gọn mô hình (bỏ biến/ràng buộc thừa, siết bound).
2. **LP relaxation** ở nút gốc: nới x ∈ [0,1], giải LP → **cận dưới** (best bound, LB).
3. **Cutting planes**: thêm lát cắt để siết LP.
4. **Heuristic**: tìm nghiệm nguyên tốt → **cận trên** (incumbent, UB).
5. **Rẽ nhánh**: chọn một x phân số, tách thành x = 0 và x = 1, lặp lại. Mỗi bài con là một **node**.

Dừng khi **gap** = 1 − LB/UB ≤ ngưỡng (ở đây 1e-4) → **Optimal**, hoặc khi hết time limit →
**Not Optimal** (có nghiệm nhưng chưa chứng minh được tối ưu). Gap cho biết nghiệm hiện có cách
tối ưu nhiều nhất bao nhiêu phần trăm.

### 2.5 Tuyến tính hóa giá trị tuyệt đối

Budget dùng |x¹ − x⁰|, không tuyến tính. Cách sách giáo khoa: thêm z ≥ 0 với z ≥ x¹ − x⁰,
z ≥ x⁰ − x¹ và Σ z ≤ κ. Vì z chỉ bị chặn dưới, luôn chọn được z = |Δ| khi và chỉ khi Σ|Δ| ≤ κ,
nên tương đương chính xác. Code dùng một biến thể chỉ cần **một** ràng buộc / cặp
(|Δ| = 2Δ⁺ − Δ), xem mục 3.

### 2.6 Cận dưới bằng Multi-Commodity Flow (MCF)

Nếu bỏ hết luật đường ngắn nhất/ECMP/SR, cho phép mỗi demand chia luồng **tùy ý** trên mọi đường,
ta được bài toán luồng đa hàng hóa — là một **LP** (dễ giải). Mọi segment routing đều là một luồng
hợp lệ của bài MCF, nên **MLU tối ưu của MCF ≤ MLU tối ưu của T-ASR**: đây là một cận dưới hợp lệ,
độc lập với MIP. Dùng nó để đánh giá chất lượng nghiệm cho các instance mà MIP không cho bound tốt
(cột "Gap theo MCF"). Cận này thường lỏng (SR bị ràng buộc bởi ECMP nên không thể chia luồng tùy ý).

---

## 3. Mô hình toán đã cài đặt

**Tập**: V (router), A (link), D (demand), T = {0} (step 6) hoặc {0,1} (step 9),
P_t = {(i,j): i ≠ j, j tới được từ i trong G(t)} (segment hợp lệ).

**Dữ liệu (Data section — step 4)**:

| Ký hiệu | Ý nghĩa | Lấy từ |
|---|---|---|
| c(a) | capacity | `net.json` → `links[].capacity` |
| ω(a) | metric | `net.json` → `links[].metric` |
| r(i,j,a,t) | hệ số chia luồng | tính từ metric trên G(t) (mục 2.2) |
| ν(d,t) | lưu lượng | `tm.json` → `demands[].v[t]` |
| maxSeg | số segment tối đa | `scenario.json` → `max_segments` |
| κ(1) | budget | `scenario.json` → `budget` |
| q(1) | link bị tắt | `scenario.json` → `interventions` |

**Biến**: x^{dt}_{ij} ∈ {0,1} — segment path của d ở chu kỳ t có dùng segment (i,j);
U ≥ 0 — MLU.

**Mục tiêu**: min U.

**Ràng buộc**:

```
(1) Σ_j x^{dt}_{ij} − Σ_j x^{dt}_{ji} = 1 (i = s_d), −1 (i = t_d), 0 (khác)     ∀d, i, t
(2) Σ_{ij} x^{dt}_{ij} ≤ maxSeg                                                  ∀d, t
(3) Σ_d Σ_{ij} r(i,j,a,t) ν(d,t) x^{dt}_{ij} ≤ U · c(a)                          ∀t, ∀a ∈ A \ q(t)
(4) Σ_d Σ_{ij} |x^{d1}_{ij} − x^{d0}_{ij}| ≤ κ(1)                               (chỉ step 9)
    x^{dt}_{ij} = 0 nếu (i,j) ∉ P_t
```

**Giải thích từng ràng buộc**:

- **Vì sao min U**: lex-min của đề có thành phần đầu là MLU; giải trọn lex-min cần giải lặp |A|·|T|
  lần → không khả thi. MLU là tiêu chí chuẩn của Traffic Engineering ("the most commonly used" — đề).
  Đây là **xấp xỉ mức 1** của lex-min. Ở step 9, cùng một U cho cả hai chu kỳ ⇒ U = max(MLU₀, MLU₁).
- **(1) bảo toàn luồng**: coi mọi cặp (i,j) là một "cung ảo" (segment). Các x = 1 tạo thành một
  đường từ s tới t trên đồ thị đầy đủ các segment; đỉnh trung gian chính là waypoint.
- **(2)**: số cung ảo trên đường = số segment = số waypoint + 1.
- **(3)**: luồng của d chạy trên segment (i,j) được trải lên các cung thật theo r. Vì đúng với mọi
  a nên U ≥ max tải; minimize ⇒ U = MLU. Link bị tắt ở t = 1 không có ràng buộc (và r của G(1) không
  bao giờ dùng nó).
- **(4) budget**: đếm segment xuất hiện ở đúng một trong hai routing (dist(P⁰,P¹) của đề).
- **Chỉ cho phép (i,j) ∈ P_t**: nếu j không tới được từ i thì r = 0 với mọi a ⇒ segment đó "chở
  luồng miễn phí" — sai thực tế; phải cấm (quan trọng ở t = 1 khi mất link).
- **Chu trình rời (subtour)**: (1) không cấm chu trình tách rời đường s→t. Nhưng r ≥ 0, ν ≥ 0 nên bỏ
  chu trình chỉ làm giảm tải và số segment ⇒ luôn có nghiệm tối ưu không có chu trình. Khi decode,
  code đi BFS từ s theo các x = 1 và bỏ chu trình thừa; MLU tính lại không tăng.
- **Demand volume 0**: không ảnh hưởng tải → gán đi thẳng s→t, bỏ khỏi mô hình (rút gọn chính xác).

**Ba cách viết lại tương đương để mô hình nhỏ hơn** (đều đã kiểm chứng trên toy cho cùng tối ưu):

1. **Ràng buộc tải dạng gộp theo segment**: đặt y^t_{ij} = Σ_d ν(d,t) x^{dt}_{ij} (luồng tổng trên
   segment (i,j)), rồi (3) thành Σ_{ij} r(i,j,a,t) y^t_{ij} ≤ U c(a). Chỉ là đặt tên cho một biểu
   thức nên tương đương tuyệt đối, nhưng số phần tử khác 0 giảm từ |D|·Σ|FG| xuống |D|·|P| + Σ|FG|
   (setA-08: 23,2 triệu → 4,6 triệu).
2. **Budget một dòng / cặp**: |Δ| = 2·Δ⁺ − Δ nên dùng z ≥ Δ, z ≥ 0 và Σ(2z − Δ) ≤ κ.
3. **Budget tách theo demand**: b_d = Σ_{ij}(2z − Δ), Σ_d b_d ≤ κ. Lý do xem mục 8.2.

**Kích thước**: số biến x = |D|·n(n−1)·|T|. setA-20 (n = 400, |D| = 6000) ở step 9 cần
**1,9 tỉ biến nhị phân** — không thể build trong 15 GB (ước lượng ~2 KB/biến ⇒ hàng TB RAM).

### 3.1 Mô hình rút gọn cho instance lớn (`restricted`)

Khi số biến x của mô hình đầy đủ vượt 2,5 triệu (setA-06 trở đi), code dùng **mô hình con**: mỗi
demand chỉ được chọn waypoint trong K đỉnh "ít vòng nhất" — các w có d(s,w) + d(w,t) − d(s,t) nhỏ
nhất; segment là mọi cặp trong {s} ∪ W_d ∪ {t} (bỏ cung đi vào s, đi ra khỏi t, vì một đường đơn
s→t không bao giờ dùng chúng). K chọn sao cho mô hình step 9 có ≈ 600 000 biến x
(K = ⌊√(600 000 / 2|D|)⌋ − 1, kẹp trong [3, 30]) và **dùng cùng K cho step 6 và step 9** — nhờ đó
mô hình step 6 đúng bằng phần t = 0 của mô hình step 9, hai bảng so sánh được với nhau.

- Mọi nghiệm của mô hình con **khả thi cho bài gốc** ⇒ objective là **cận trên** (một heuristic).
- "Optimal" của mô hình con chỉ là tối ưu **trong không gian bị giới hạn**, *không* phải tối ưu của
  bài gốc; best bound của nó **không** phải cận dưới của bài gốc. Vì vậy cần thêm cận dưới MCF.
- Nếu objective = cận dưới MCF thì nghiệm **chắc chắn tối ưu cho bài gốc** (kẹp giữa hai cận).

### 3.2 Cách giải mô hình đầy đủ (`exact`, setA-01…05)

Hai pha trong cùng một lần chạy:

1. **Pha 1 (warm start)**: giải mô hình rút gọn K = 20 với ≤ 10 % time limit. Nghiệm của nó khả
   thi cho mô hình đầy đủ (mục 3.1).
2. **Pha 2**: build mô hình đầy đủ, nạp nghiệm pha 1 làm MIP start, giải với LP gốc bằng
   **interior point (IPM)** thay cho dual simplex.

Lý do (đều đo được, xem mục 8): với mô hình đầy đủ, heuristic của HiGHS hầu như không tìm được nghiệm
tốt hơn baseline, và dual simplex không giải xong LP gốc trong 1 giờ (setA-03 step 9). CPU time trong
bảng là CPU time của `optimize!` ở **pha 2** (đúng yêu cầu readme); CPU của pha 1 ghi riêng
(`warm_cpu_time` trong `results/raw/…json`, ≤ 360 s).

---

## 4. Cấu trúc code và flow chạy

```
orange/
  Project.toml, Manifest.toml   môi trường Julia (JSON3, Graphs, JuMP, HiGHS)
  src/orange.jl                 TOÀN BỘ logic (đọc dữ liệu, r, mô hình, giải, decode, kiểm tra)
  test_toy.jl                   kiểm tra trên instance toy của đề (39 assert)
  prepare.jl                    phân tích 1 instance: r, kích thước mô hình, baseline, cận dưới MCF
  run_one.jl                    1 lần giải = (instance, step 6|9, seed)
  run_all.jl                    điều phối: test → prepare → các lần giải song song → tổng hợp
  aggregate.jl                  gom kết quả thành bảng
  check_solution.jl             đọc file nghiệm, tính lại MLU độc lập, so với bảng
  run.sh                        chạy run_all.jl trong cgroup giới hạn 15 GB
  toy/                          instance toy chép từ đề
  results/
    instances/setA-XX.json      thông tin instance (prepare): kích thước, r, baseline, cận MCF
    raw/step6|step9/*.json      kết quả từng lần giải (mọi con số trong bảng lấy từ đây)
    raw/step6_ipm|step9_ipm/    thử nghiệm bổ sung (mục 7.5)
    highs_logs/…                log HiGHS của đúng lần giải trong bảng (đã lọc)
    highs_logs_full/…           log gốc (HiGHS ghi nối tiếp qua các lần chạy lại)
    solutions/…-srpaths.json    nghiệm (waypoint của từng demand) theo định dạng của đề
    table_step6.{md,csv}        BẢNG KẾT QUẢ STEP 6
    table_step9.{md,csv}        BẢNG KẾT QUẢ STEP 9
    all_runs_step*.csv          mọi cột của mọi lần chạy
    run_all.log / run_all.out   nhật ký điều phối (cả các lần khởi động lại, mục 8)
```

### 4.1 Flow tổng thể (`./run.sh`)

```
run.sh ── systemd-run (cgroup MemoryMax = 15·10⁹ byte, không swap, OOMPolicy=continue)
  └─ run_all.jl
      1. test_toy.jl                       → dừng nếu sai
      2. prepare.jl × 20 instance          → results/instances/setA-XX.json
         (instance nhỏ 3 process song song; setA-14…20 tuần tự vì cận MCF tốn tới 7,8 GB)
      3. tạo 40 job = 20 instance × {step 9, step 6} × 1 seed
         chọn exact / restricted (+K) theo kích thước; xếp instance lớn trước
         vòng lặp mỗi 2 s: khởi động job khi RAM còn đủ → theo dõi RAM → thu kết quả
           └─ run_one.jl <inst> <step> <seed> <mode> <time_limit> [K]   (process riêng)
                → results/raw/step<6|9>/<inst>_s1.json
                → results/highs_logs/step<6|9>/<inst>_s1.log
                → results/solutions/step<6|9>/<inst>-srpaths.json
         job đã có file kết quả thì bỏ qua → chạy lại lệnh = resume
      4. aggregate.jl                      → results/table_step6.{md,csv}, table_step9.{md,csv}
```

### 4.2 Bên trong một lần giải (`run_one.jl`)

```
(chạy toy 2 lần trước để loại thời gian biên dịch JIT khỏi phép đo)
load_instance(prefix)                          đọc 3 file JSON (step 4)
  └─ json_to_graph(net.json)                   Network: đỉnh, cung, metric, capacity (step 2)
period_data(inst, t) cho t ∈ periods
  ├─ graph_at(net, q(t))                       G(t) = G bỏ các link bị intervention (step 7)
  ├─ compute_r → all_pairs_shortest_data       Dijkstra allpaths (Graphs.jl) → d, σ → r (step 3)
  └─ check_r_conservation                      r không bảo toàn luồng (sai số float khi metric
                                               không nguyên) → tính lại bằng Dijkstra có tolerance
baseline = MLU khi mọi demand đi thẳng s→t (không waypoint)
[exact] pha 1: mô hình rút gọn K=20, ≤10% time limit → nghiệm warm start
candidate_segments(:exact | :restricted, K)
build_model                                    (1)(2)(3) [+ (4)]   (step 5 / step 8)
set_warm_start!                                nghiệm pha 1 (exact) hoặc baseline (restricted)
solve!                                         HiGHS, threads=1, seed=1, time limit;
                                               LP gốc: IPM (exact) / dual simplex (restricted)
                                               đo CPU time CHỈ của optimize!
decode_paths                                   x → danh sách waypoint (BFS s→t)
verify_solution                                kiểm tra ĐỘC LẬP: đường đơn, ≤ maxSeg, segment tồn tại
                                               ở G(t), budget, tính lại MLU bằng arc_loads
                                               → phải khớp objective, sai thì báo lỗi
```

### 4.3 Thiết lập solver và môi trường

| Mục | Giá trị |
|---|---|
| Máy | AMD EPYC 7742 (64 nhân), Linux 6.8 — chỉ được dùng ≤ 15 GB RAM |
| Phần mềm | Julia 1.13.2, JuMP 1.32.1, HiGHS.jl 1.26.0, Graphs.jl 1.15.0, JSON3 1.14.3 |
| Solver | HiGHS qua `direct_model` (không giữ bản copy → tiết kiệm RAM), `threads = 1`, `random_seed = 1`, `mip_rel_gap = 1e-4` |
| Time limit | 1800 s (step 6), 3600 s (step 9) |
| Số lần chạy | **1** / instance / step (readme yêu cầu trung bình 10 lần) |
| CPU time | CPU time của process trong khoảng `optimize!` (đo bằng `clock()`) |

---

## 5. Giới hạn 15 GB RAM

Máy chạy có nhiều RAM hơn, nên giới hạn 15 GB được **ép cứng** chứ không chỉ "cố gắng":

1. **cgroup v2 của systemd** — `run.sh` chạy `run_all.jl` bằng
   `systemd-run --scope --slice=orange.slice -p MemoryMax=15000000000 -p MemorySwapMax=0 -p OOMPolicy=continue`.
   Mọi process con nằm trong cgroup đó; tổng RAM **không thể** vượt giới hạn (kernel OOM-kill).
   Slice `~/.config/systemd/user/orange.slice` có cùng giới hạn nên cả các lần chạy thử, hiệu chỉnh và
   thử nghiệm bổ sung đều nằm chung trong 15 GB.
2. **Bộ điều phối (run_all.jl)** đọc RAM thực dùng của cả slice (`memory.current` − page cache) mỗi 2 s:
   - ước lượng RAM: mô hình exact 1,5 GB + 4 KB/biến x; mô hình rút gọn 6,5 GB; step 6 rút gọn dùng
     1,1 × RAM đo được của step 9 cùng instance;
   - chỉ khởi động job mới khi tổng ước lượng ≤ 12,5 GiB **và** RAM thực + ước lượng ≤ 12 GiB;
   - watchdog: RAM thực > 13,3 GiB → kill job lớn nhất, ghi vào `results/ram_kills.txt`, chạy lại với
     mô hình nhỏ hơn (exact → rút gọn K = 20; rút gọn → K × 0,6);
   - bộ điều phối có `oom_score_adj = −900` để nếu kernel phải OOM-kill thì nó giết một job giải.
3. **Giảm RAM của chính mô hình**: `direct_model`, ràng buộc tải dạng gộp (mục 3), r lưu thưa dạng
   CSR, mô hình rút gọn, cận MCF giải từng chu kỳ riêng, prepare của instance lớn chạy tuần tự.

> **Trung thực về con số**: trong lần chạy này giới hạn cứng ban đầu là `MemoryMax=15G`, mà systemd
> hiểu là 15 GiB = 16,1·10⁹ byte. RAM chạm đúng giới hạn đó **một lần** (13:30, sự cố OOM ở mục 8.5)
> và một lần watchdog đo được 14,1 GiB (≈ 15,1·10⁹ byte) ngay trước khi kill job. Sau đó `run.sh` và
> slice đã được đổi sang đúng 15·10⁹ byte và hạ ngưỡng watchdog, để lần chạy lại tuân thủ chặt "15 GB".
> RAM lớn nhất của **một** lần giải hoàn tất: 7,5 GB (setA-05 step 9, mô hình đầy đủ 3 triệu biến).

---

## 6. Kết quả

Bảng đầy đủ: `results/table_step6.md`, `results/table_step9.md` (và bản CSV). Các cột:

- **Mô hình**: `exact` = mô hình đầy đủ; `restricted(K)` = mô hình rút gọn với K waypoint ứng viên.
- **Status / Best obj / Best bound / Gap**: theo HiGHS, gap = 1 − LB/UB (readme). Với mô hình rút gọn,
  bound và gap là của **mô hình con**.
- **CPU time**: của `optimize!`; **Nodes**: số node branch-and-bound.
- **MLU baseline**: không waypoint. **Giảm MLU** = 1 − Best obj / baseline.
- **Cận dưới MCF** (hợp lệ cho bài gốc) và **Gap theo MCF** = 1 − LB_MCF / Best obj — là gap *thật*
  tối đa so với bài gốc (thường rất bi quan vì cận MCF lỏng).
- Step 9: **Thay đổi / κ** = dist(P⁰, P¹) thực tế / budget.

### 6.1 Step 6 — chỉ t = 0, không budget, time limit 1800 s

| Inst | \|V\| | \|A\| | \|D\| | Mô hình | Status | Best obj | Best bound | Gap | CPU (s) | Nodes | Baseline | Giảm MLU | Cận MCF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 01 | 20 | 80 | 40 | exact | Optimal | 0.9294 | 0.9294 | 0 % | 0.8 | 1 | 1.0000 | 7.1 % | 0.5604 |
| 02 | 30 | 150 | 45 | exact | Optimal | 0.5486 | 0.5486 | 0 % | 3.1 | 1 | 1.0000 | 45.1 % | 0.5388 |
| 03 | 50 | 250 | 20 | exact | Optimal | 0.9435 | 0.9435 | 0 % | 404.6 | 1 | 1.0000 | 5.6 % | 0.6337 |
| 04 | 50 | 250 | 200 | exact | Not Optimal | 0.6918 | 0.4591 | 33.6 % | 1799.2 | 0 | 0.6918 | 0 % | 0.3817 |
| 05 | 100 | 396 | 100 | exact | Not Optimal | 0.1324 | 0.1155 | 12.8 % | 3050.1* | 0 | 0.3333 | 60.3 % | 0.1155 |
| 06 | 100 | 500 | 500 | restr. K=23 | Not Optimal | 1.0000 | 0 | 100 % | 1799.8 | 0 | 1.0000 | 0 % | 0.0986 |
| 07 | 100 | 500 | 800 | restr. K=11 | Optimal | 0.9080 | 0.9080 | 0 % | 3.4 | 1 | 0.9080 | 0 % | **0.9080** |
| 08 | 150 | 654 | 200 | restr. K=30 | Not Optimal | 1.0000 | 0.0098 | 99.0 % | 1799.5 | 0 | 1.0000 | 0 % | 0.3189 |
| 09 | 150 | 750 | 200 | restr. K=30 | Optimal | 0.5054 | 0.5054 | 0 % | 27.9 | 1 | 0.6111 | 17.3 % | 0.4254 |
| 10 | 150 | 966 | 1000 | restr. K=16 | Optimal | 0.4239 | 0.4239 | 0 % | 238.9 | 1 | 1.0000 | 57.6 % | 0.0683 |
| 11 | 200 | 1000 | 400 | restr. K=26 | Optimal | 0.6446 | 0.6446 | 0 % | 130.9 | 203 | 1.0000 | 35.5 % | 0.6131 |
| 12 | 200 | 898 | 400 | restr. K=26 | Optimal | 1.0001 | 1.0001 | 0 % | 10.7 | 1 | 1.0001 | 0 % | 0.5789 |
| 13 | 200 | 1000 | 2000 | restr. K=11 | Optimal | 0.6197 | 0.6197 | 0 % | 58.8 | 1 | 1.0000 | 38.0 % | 0.0389 |
| 14 | 250 | 1108 | 600 | restr. K=21 | Optimal | 0.4397 | 0.4397 | 0 % | 40.2 | 1 | 1.0000 | 56.0 % | 0.4239 |
| 15 | 250 | 1250 | 600 | restr. K=21 | Optimal | 0.8252 | 0.8252 | 0 % | 10.9 | 1 | 1.0000 | 17.5 % | **0.8252** |
| 16 | 250 | 1452 | 4800 | restr. K=6 | Optimal | 0.8730 | 0.8730 | 0 % | 19.2 | 1 | 0.9672 | 9.7 % | 0.0386 |
| 17 | 300 | 1270 | 2000 | restr. K=11 | Optimal | 0.4242 | 0.4242 | 0 % | 35.8 | 1 | 1.0001 | 57.6 % | **0.4242** |
| 18 | 300 | 1500 | 2000 | restr. K=11 | Optimal | 0.6334 | 0.6334 | 0 % | 10.3 | 1 | 0.6334 | 0 % | **0.6334** |
| 19 | 300 | 1998 | 6000 | restr. K=6 | Optimal | 0.7377 | 0.7377 | 0 % | 25.8 | 1 | 1.0000 | 26.2 % | 0.0468 |
| 20 | 400 | 2000 | 6000 | restr. K=6 | Optimal | 0.9913 | 0.9913 | 0 % | 12.9 | 1 | 1.0000 | 0.9 % | 0.6482 |

\* HiGHS vượt time limit (xem mục 7.4). **In đậm**: Best obj = cận MCF ⇒ tối ưu cho bài gốc.

### 6.2 Step 9 — t ∈ {0, 1}, có budget (4), time limit 3600 s

| Inst | Mô hình | Status | Best obj | Best bound | Gap | CPU (s) | Nodes | Baseline | Giảm MLU | Cận MCF | MLU t=0 / t=1 | Thay đổi / κ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 01 | exact | Optimal | 0.9294 | 0.9294 | 0 % | 8.6 | 1 | 1.0000 | 7.1 % | 0.5604 | 0.929 / 0.639 | 4 / 51 |
| 02 | exact | Not Optimal | 0.9031 | 0.6243 | 30.9 % | 3594.9 | 0 | 1.0981 | 17.8 % | 0.6243 | 0.718 / 0.903 | 18 / 63 |
| 03 | exact | Not Optimal | 0.9435 | 0.6337 | 32.8 % | 3597.3 | 0 | 1.0000 | 5.6 % | 0.6337 | 0.944 / 0.848 | 22 / 53 |
| 04 | restr. K=20† | Optimal | 0.7421 | 0.7421 | 0 % | 50.0 | 1 | 1.0000 | 25.8 % | 0.4832 | 0.643 / 0.742 | 0 / 44 |
| 05 | exact | Not Optimal | 0.2050 | 0 | 100 % | 3632.4* | 0 | 1.0000 | 79.5 % | **0.2050** | 0.167 / 0.205 | 0 / 1 |
| 06 | restr. K=23 | Not Optimal | 1.0000 | 0 | 100 % | 3598.6 | 0 | 1.0000 | 0 % | 0.0986 | 1.000 / 0.775 | 0 / 13 |
| 07 | restr. K=11† | Optimal | 0.9080 | 0.9080 | 0 % | 126.8 | 1 | 1.0000 | 9.2 % | **0.9080** | 0.908 / 0.858 | 14 / 90 |
| 08 | restr. K=30 | Not Optimal | 1.0000 | 0.0147 | 98.5 % | 3596.5 | 0 | 1.0000 | 0 % | 0.3189 | 1.000 / 0.702 | 0 / 13 |
| 09 | restr. K=30 | Optimal | 0.8497 | 0.8497 | 0 % | 61.6 | 1 | 1.0000 | 15.0 % | **0.8497** | 0.709 / 0.850 | 0 / 18 |
| 10 | restr. K=16 | Optimal | 0.4239 | 0.4239 | 0 % | 628.6 | 1 | 1.0000 | 57.6 % | 0.0683 | 0.424 / 0.424 | 0 / 1 |
| 11 | restr. K=26 | Optimal | 0.7858 | 0.7858 | 0 % | 89.0 | 1 | 1.0000 | 21.4 % | **0.7858** | 0.656 / 0.786 | 5 / 89 |
| 12 | restr. K=26 | Optimal | 1.0001 | 1.0001 | 0 % | 70.4 | 1 | 1.0001 | 0 % | 0.8799 | 1.000 / 0.880 | 0 / 13 |
| 13 | restr. K=11 | Optimal | 0.6496 | 0.6496 | 0 % | 482.5 | 1 | 1.0000 | 35.0 % | 0.0408 | 0.645 / 0.650 | 9 / 12 |
| 14 | restr. K=21 | Optimal | 0.5331 | 0.5331 | 0 % | 315.3 | 1 | 1.0000 | 46.7 % | 0.5074 | 0.455 / 0.533 | 0 / 13 |
| 15 | restr. K=21 | Optimal | 0.8987 | 0.8987 | 0 % | 71.1 | 1 | 1.0000 | 10.1 % | **0.8987** | 0.880 / 0.899 | 0 / 54 |
| 16 | restr. K=6 | Optimal | 0.8811 | 0.8811 | 0 % | 196.0 | 1 | 1.0000 | 11.9 % | 0.0441 | 0.873 / 0.881 | 4 / 13 |
| 17 | restr. K=11 | Optimal | 0.4242 | 0.4242 | 0 % | 328.9 | 1 | 1.0001 | 57.6 % | **0.4242** | 0.424 / 0.366 | 0 / 1 |
| 18 | restr. K=11 | Optimal | 1.0000 | 1.0000 | 0 % | 46.6 | 1 | 1.0000 | 0 % | **1.0000** | 0.633 / 1.000 | 3 / 89 |
| 19 | restr. K=6 | Optimal | 0.7377 | 0.7377 | 0 % | 249.6 | 1 | 1.0000 | 26.2 % | 0.0468 | 0.738 / 0.733 | 0 / 13 |
| 20 | restr. K=6 | Optimal | 0.9913 | 0.9913 | 0 % | 62.4 | 1 | 1.0000 | 0.9 % | 0.6482 | 0.991 / 0.862 | 28 / 90 |

† Bị watchdog kill vì RAM ở lần đầu (04: exact 9 GB; 07: K = 18, 12 GB) → chạy lại với mô hình nhỏ hơn.
\* Vượt time limit. Kích thước (|V|, |A|, |D|) như bảng 6.1.

### 6.3 Tóm tắt

| | Step 6 | Step 9 |
|---|---|---|
| Status "Optimal" (theo HiGHS) | 16 / 20 | 15 / 20 |
| Tối ưu cho **bài gốc** (exact Optimal; hoặc obj = cận MCF; hoặc step 9 obj = tối ưu exact của step 6) | 7 (01, 02, 03, 07, 15, 17, 18) | 9 (01, 03, 05, 07, 09, 11, 15, 17, 18) |
| Giảm được MLU so với baseline | 14 / 20 | 16 / 20 |
| Giảm MLU trung bình (20 instance) | 21,7 % | 21,4 % |
| Không cải thiện (MLU = baseline) | 04, 06, 07, 08, 12, 18 | 06, 08, 12, 18 |
| Tổng CPU time của `optimize!` | ≈ 2,6 giờ | ≈ 5,8 giờ |

Mọi nghiệm đều qua kiểm tra độc lập (mục 8.4): MLU tính lại khớp objective, số segment ≤ 6,
không dùng link bị tắt ở t = 1, số thay đổi ≤ κ.

---

## 7. Đánh giá và giải thích

### 7.1 Segment Routing giảm MLU đáng kể — nhưng không đều

Trung bình MLU giảm ~21 %, nhiều instance giảm 45–80 % (02, 05, 10, 14, 17). Cơ chế: baseline định
tuyến mọi demand theo đường ngắn nhất, nên các link "trung tâm" gánh quá nhiều; thêm 1–2 waypoint
đẩy một phần demand sang các đường song song còn trống. Số waypoint dùng nhiều ở instance có nhiều
demand (setA-20 step 9: 9 906 waypoint cho 12 000 cặp (demand, chu kỳ)).

Baseline của hầu hết instance có MLU = 1,0000 chính xác (hoặc 1,0001) — dữ liệu được sinh bằng cách
**co giãn lưu lượng sao cho baseline vừa chạm 100 %**, nên MLU < 1 nghĩa là SR đã gỡ được nghẽn.

### 7.2 Vì sao một số instance KHÔNG giảm được MLU — 3 nguyên nhân khác nhau

**(a) Baseline đã tối ưu — giới hạn vật lý của mạng** (07 step 6, 18 cả hai step).
Ở đây Best obj = cận MCF, tức ngay cả khi chia luồng tùy ý cũng không làm tốt hơn. Ví dụ setA-18:
link nghẽn nhất là (v138 → v90) và **v138 chỉ có đúng một link đi ra** (outdeg = 1); mọi demand xuất
phát từ v138 bắt buộc đi qua link đó, waypoint nào cũng vô ích. Ở step 9, link bị tắt làm MLU t = 1
lên 1,000 — cũng là giới hạn vật lý (cận MCF = 1,000).

**(b) Mô hình rút gọn bỏ sót waypoint cần thiết** (setA-12). Link nghẽn là (v150 → v175), một trong
hai link đi vào v175; link kia (v161 → v175) có metric 0,97 — gấp ~27 lần các link xung quanh (0,036).
Bốn demand đổ về v175 phải chia sang link (v161 → v175), tức phải dùng **waypoint v161**. Kiểm tra
tay: đưa demand 140 và 351 qua v161 cho MLU = **0,9898 < 1,0001**. Nhưng vì metric lớn, độ vòng qua
v161 xếp hạng 161/198 → không nằm trong K = 26 ứng viên "ít vòng nhất". Mô hình con báo "Optimal"
một cách đúng *trong không gian của nó*, nhưng không gian đó thiếu nghiệm tốt. Bài học: tiêu chí
chọn ứng viên theo độ vòng không nhìn thấy link nghẽn; nên ưu tiên các đỉnh kề link còn trống.

**(c) Solver không giải xong LP gốc trong time limit** (06 và 08, cả hai step; 04 step 6).
Log HiGHS: hết 1800/3600 s vẫn ở node 0, chưa có LP bound → không có thông tin để rẽ nhánh hay chạy
heuristic, nên giữ nguyên warm start (baseline). Với 04 step 6 (exact): pha 1 chỉ có 180 s nên mô
hình rút gọn K = 20 chưa giải xong, mô hình đầy đủ (475 000 biến) có bound 0,459 nhưng không tìm
được nghiệm tốt hơn — trong khi ở step 9, cùng mô hình rút gọn K = 20 giải 50 s cho MLU t = 0 là
0,643 < 0,6918, chứng tỏ có nghiệm tốt hơn baseline. Đây là giới hạn của **thuật toán LP (dual simplex)
trên LP suy biến nặng**, không phải của mô hình: thử nghiệm bổ sung ở mục 7.5 với IPM giải setA-08
step 6 tối ưu (= cận MCF 0,3189) trong 172 s.

### 7.3 Mô hình đầy đủ (exact) và quy mô

- Giải tối ưu được khi mô hình nhỏ (01–03 step 6, 01 step 9: ≤ 50 đỉnh, ≤ 50 000 biến x).
- Từ ~100 000 biến x trở lên (02/03 step 9, 04, 05), LP gốc rất lớn và suy biến: IPM giải được LP gốc
  (cho bound) nhưng B&B không đi được node nào trong time limit (Nodes = 0). Gap 13–34 %.
- Gap theo HiGHS ở exact **lớn nhưng không có nghĩa nghiệm tệ**: setA-05 step 9 báo gap 100 % (LP
  gốc chưa xong nên bound = 0) nhưng objective 0,2050 **bằng đúng cận MCF** ⇒ thực ra đã tối ưu.
  Tương tự setA-03 step 9 báo gap 33 % nhưng obj 0,9435 = tối ưu đã chứng minh của step 6 — mà step 6
  là nới lỏng của step 9 (bỏ chu kỳ t = 1) ⇒ step 9 cũng tối ưu.
  Ngược lại setA-01 có tối ưu chứng minh được 0,9294 trong khi cận MCF chỉ 0,5604: **cận MCF có thể
  rất lỏng** (SR không chia luồng tùy ý được), nên "gap theo MCF" chỉ là giới hạn trên của gap thật.
- Kích thước tăng theo |D|·n²: setA-20 cần 1,9 tỉ biến cho mô hình đầy đủ step 9 — không thể với
  15 GB. Mô hình rút gọn (≈ 0,6–0,9 triệu biến) giải được mọi instance lớn trong vài phút, phần lớn
  đạt "Optimal" **ngay tại node gốc** (Nodes = 1): LP relaxation của bài toán này rất chặt — nghiệm
  đỉnh của simplex thường đã nguyên.
- Pha warm start (mô hình rút gọn K = 20) là yếu tố quyết định cho exact: ở 01, 02 step 6 và 05,
  nghiệm cuối chính là nghiệm pha 1 (pha 2 chỉ chứng minh tối ưu hoặc thu hẹp bound); 03 cải thiện
  nhẹ ở pha 2 (0,9444 → 0,9435); 02 step 9 pha 2 cải thiện mạnh (1,0981 → 0,9031) vì pha 1 không
  cải thiện được gì (n = 30 nên K = 20 gần bằng mô hình đầy đủ, LP không xong trong 360 s).

### 7.4 Thời gian

- Mô hình rút gọn: 3–630 s; mô hình exact cỡ vừa chạm time limit.
- **HiGHS không phải lúc nào cũng dừng đúng time limit**: setA-05 step 6 dùng 3050 s CPU cho time
  limit 1800 s (IPM chỉ kiểm tra thời gian giữa các vòng lặp lớn, và một vòng IPM trên LP 1 triệu
  biến rất lâu). Bảng báo đúng CPU đo được, không cắt.
- Thời gian ngoài `optimize!` (tính r, build mô hình) nhỏ: r ≤ 4 s, build ≤ 95 s.

### 7.5 Thử nghiệm bổ sung: IPM cho setA-06 và setA-08

Cùng mô hình rút gọn, cùng K, chỉ đổi thuật toán LP gốc sang IPM (`ORANGE_LP=ipm`,
`results/raw/step6_ipm`, `step9_ipm`):

| Instance / step | Simplex (bảng chính) | IPM | CPU IPM (s) |
|---|---|---|---|
| setA-06 step 6 | 1,0000, không có bound (gap 100 %) | 1,0000, bound 0,1303 (gap 87 %) | 1799 |
| setA-08 step 6 | 1,0000, bound 0,0098 (gap 99 %) | **Optimal 0,3189** (= cận MCF ⇒ tối ưu bài gốc, giảm 68 %) | 172 |
| setA-06 step 9 | 1,0000, không có bound (gap 100 %) | 1,0000, bound 0,2535 (gap 75 %) | 3600 |
| setA-08 step 9 | 1,0000, bound 0,0147 (gap 99 %) | 1,0000, bound 0,3189 (gap 68 %) | 3770 |

IPM luôn giải xong LP gốc (có bound tốt hơn hẳn), nhưng chỉ ở setA-08 step 6 heuristic mới tìm ra
nghiệm nguyên tốt hơn baseline. Ở setA-08 step 9, chỉ cần dùng lại routing tối ưu của step 6 cho cả
hai chu kỳ (0 thay đổi) là đã có nghiệm MLU = max(0,3189, MLU t=1 của routing đó) — gợi ý một warm
start tốt hơn cho step 9: **khởi động từ nghiệm step 6** thay vì từ baseline.

Simplex và IPM bổ sung cho nhau: simplex cho nghiệm đỉnh (thường đã nguyên → có ngay incumbent tối
ưu), IPM giải LP suy biến lớn nhanh hơn nhiều nhưng nghiệm nằm trong miền (không nguyên). Ví dụ
ngược: setA-20 step 9 với IPM có bound 0,9913 sau 164 s nhưng sau 20 phút vẫn chưa tìm được
nghiệm < 1,0, trong khi simplex cho Optimal 0,9913 sau 62 s. Một chiến lược tốt hơn (chưa làm):
chạy simplex, nếu sau một khoảng thời gian chưa có LP bound thì chuyển sang IPM (hoặc chạy song song).

### 7.6 Step 9 so với step 6

- Vì mô hình step 6 = phần t = 0 của step 9 (cùng K), luôn có obj(step 9) ≥ obj(step 6) — đúng ở mọi
  instance. Chu kỳ quyết định MLU thường là **t = 1** (link bị tắt + ma trận lưu lượng khác).
- Budget ít khi chặn: phần lớn số thay đổi thực tế ≪ κ, nhiều instance dùng 0 thay đổi (cùng một bộ
  waypoint tốt cho cả hai chu kỳ). Ngoại lệ đáng chú ý: setA-13 dùng 9/12, setA-05 κ = 1 (gần như
  bắt buộc giữ nguyên routing).
- Ở 01, 07, 10, 17, 19, 20 MLU step 9 = step 6: chu kỳ t = 0 là nút thắt, intervention ở t = 1 không làm
  tình hình tệ hơn khi đã có waypoint.

### 7.7 Hạn chế của kết quả

- **1 lần chạy** / instance thay vì trung bình 10 lần (theo yêu cầu của bạn): không có độ lệch chuẩn;
  với các instance Optimal ở node gốc thì seed hầu như không ảnh hưởng, với instance chạm time limit
  thì kết quả có thể dao động.
- Objective là MLU (mức 1 của lex-min), không phải lex-min đầy đủ của đề.
- r theo `calculating_r.pdf` (chia theo số đường), không phải ECMP chia tại nút như đề Orange.
- Mô hình rút gọn là heuristic; "Optimal" của nó không phải tối ưu của bài gốc.

---

## 8. Những điểm kỹ thuật quan trọng

### 8.1 Sai số số thực khi đếm đường ngắn nhất

Ở setA-08 và setA-14, metric không nguyên (vd. `0.013814041391015053`). Graphs.jl so sánh
`ℓ == d(v)` tuyệt đối, nên hai đường dài bằng nhau về toán học có thể lệch 1 ulp → **đếm thiếu
đường** → r không còn là luồng đơn vị (sai số bảo toàn luồng = 1,0). Code phát hiện bằng
`check_r_conservation` và tính lại bằng Dijkstra tự viết so sánh có tolerance (rtol = 1e-9), đúng
thuật toán Case 1/2/3 của `calculating_r.pdf`. Sau đó sai số = 0 ở mọi instance.

### 8.2 Hàng budget dày làm presolve của HiGHS bị treo

Viết budget thành **một** ràng buộc Σ z ≤ κ chứa mọi biến tạo ra một hàng có hàng trăm nghìn hệ số.
Thử nghiệm (time limit 120 s):

| Mô hình | Budget một hàng dày | Budget tách theo demand |
|---|---|---|
| setA-02 step 9 (exact) | presolve 204 s, không có bound | presolve 5 s, bound 0,624 sau 55 s |
| setA-13 step 9 (rút gọn) | presolve > 20 phút (bỏ qua time limit) | **Optimal trong 38 s** |

Tách thành b_d = Σ_{ij}|Δ| cho từng demand rồi Σ_d b_d ≤ κ là **tương đương chính xác** nhưng mỗi
hàng chỉ vài chục phần tử. (Đã thử tắt "probing" như log HiGHS gợi ý — không đủ; nguyên nhân là
hàng dày.)

### 8.3 Dual simplex không giải xong LP gốc của mô hình đầy đủ

setA-03 step 9 (exact, 150 000 cột): dual simplex chạy 752 481 vòng lặp trong 3600 s vẫn chưa xong
LP gốc → gap 100 %. Cùng mô hình với `mip_lp_solver = "ipm"`: LP gốc xong trong 36 s, bound 0,6337.
Nhưng cận trên vẫn kẹt ở baseline (heuristic không tìm được nghiệm) → thêm pha warm start (mục 3.2):
UB 0,944 ngay từ đầu.

### 8.4 Kiểm tra độc lập nghiệm

Không tin objective của solver: mọi nghiệm được decode thành segment path và tính lại MLU bằng hàm
`arc_loads` (không dùng mô hình), kiểm tra số segment, segment tồn tại ở G(t), budget. Run báo lỗi
nếu MLU tính lại lớn hơn objective. Warm start "không waypoint" (luôn khả thi, budget dùng 0) đảm bảo
solver luôn có incumbent và objective không bao giờ tệ hơn baseline.

### 8.5 Nhật ký các lần khởi động lại (vì sao `run_all.out` có nhiều đoạn RESTART)

| Lúc | Sự kiện | Xử lý |
|---|---|---|
| 13:30 | cgroup chạm giới hạn → kernel OOM-kill 1 job; systemd mặc định (`OOMPolicy=stop`) **dừng cả scope** → mọi job chết | thêm `OOMPolicy=continue`; ước lượng RAM bảo thủ hơn (mô hình rút gọn dùng 2,6–7,1 GB chứ không phải 2,7 GB) |
| 13:54 | exact setA-03 step 9 gap 100 % (mục 8.3) | thêm pha warm start + IPM cho exact; chạy lại 01–03 |
| 14:38 | setA-07 step 9 (K=18) phình 12,2 GB → watchdog kill; lần chạy lại bị chặn vì ước lượng 18 GB > ngân sách | khi bị kill vì RAM: thu nhỏ K (×0,6), lưu vào `ram_kills.txt`; đo RAM ở cấp slice |
| 18:44 | thấy step 6 dùng K lớn hơn step 9 (tới 30) → LP gốc của 06/07/08/13 không xong; thử IPM cho mô hình rút gọn thì tốt ở 08/13 nhưng tệ ở setA-20 | chốt: **cùng K cho step 6 & 9**, simplex cho rút gọn, IPM cho exact; chạy lại step 6 rút gọn |

Kết quả cuối trong bảng đều được chạy với **cùng một cấu hình** (mục 3.1, 3.2, 4.3), trừ: setA-05
step 9 có pha warm start giải bằng IPM (chạy trước khi chốt cấu hình; kết quả = cận MCF nên đã tối ưu,
không chạy lại vì tốn 1,7 giờ và 7,5 GB).

---

## 9. Cách chạy lại

```bash
cd orange
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # lần đầu
julia --project=. test_toy.jl                          # kiểm tra (39 assert)
./run.sh                                               # toàn bộ, giới hạn 15 GB, 1 seed (~8 giờ)
ORANGE_NSEEDS=10 ./run.sh                              # như readme: 10 seed
julia --project=. aggregate.jl                         # chỉ tổng hợp lại bảng

# một lần giải đơn lẻ:
julia --project=. run_one.jl setA-09 9 1 restricted 3600 30
julia --project=. run_one.jl setA-02 6 1 exact 1800
```

Chạy lại `./run.sh` sẽ **resume**: job đã có `results/raw/…json` được bỏ qua. Muốn chạy lại từ đầu
thì xóa `results/raw`, `results/ram_kills.txt`.
