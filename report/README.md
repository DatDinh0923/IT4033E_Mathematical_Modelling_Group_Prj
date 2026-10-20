# Report (LaTeX)

- `main.tex` — report (≤ 10 trang, tiếng Anh). Tên tác giả để placeholder ở `\author{...}`.
- `make_tables.py` — sinh `tables/*.tex` (bảng + dữ liệu biểu đồ) trực tiếp từ `../results/`.
- `tables/` — bảng đã sinh; không sửa tay, chạy lại script nếu kết quả thay đổi.

Build:

```bash
cd orange/report
python3 make_tables.py          # chỉ cần khi results/ thay đổi
latexmk -pdf main.tex           # -> main.pdf
```
