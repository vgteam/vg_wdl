version 1.0

task plotAardvarkComparison {
    input {
        File baseline_summary
        File candidate_summary
        String baseline_name
        String candidate_name
    }

    command <<<
        set -euo pipefail
        cat > plot.py <<'EOF'
import csv, sys

def load(path):
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["comparison"] == "GT" and r["region_label"] == "ALL" and r["filter"] == "ALL":
                out[r["variant_type"]] = r
    return out

def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

base, cand = load(sys.argv[1]), load(sys.argv[2])
base_name, cand_name, out = sys.argv[3], sys.argv[4], sys.argv[5]

types = [t for t in ["ALL", "Snv", "Insertion", "Deletion"] if t in base or t in cand]
metrics = [("metric_recall", "Recall"), ("metric_precision", "Precision"), ("metric_f1", "F1")]
colors = {"base": "#4C72B0", "cand": "#DD8452"}

W, H = 1100, 460
pad_l, pad_t, pad_b = 60, 70, 90
panel_w = (W - pad_l - 20) / len(metrics)
plot_h = H - pad_t - pad_b
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="Helvetica,Arial,sans-serif" font-size="12">',
       f'<rect width="{W}" height="{H}" fill="white"/>',
       f'<text x="{W/2}" y="24" text-anchor="middle" font-size="16" font-weight="bold">Aardvark GT metrics: baseline vs candidate</text>']
# legend
for i, (k, label) in enumerate([("base", "baseline: " + base_name), ("cand", "candidate: " + cand_name)]):
    x = pad_l + i * 420
    svg.append(f'<rect x="{x}" y="38" width="12" height="12" fill="{colors[k]}"/>')
    svg.append(f'<text x="{x+18}" y="49">{label}</text>')

for m_i, (col, title) in enumerate(metrics):
    x0 = pad_l + m_i * panel_w
    vals = [num(d.get(t, {}).get(col)) for d in (base, cand) for t in types]
    vals = [v for v in vals if v is not None]
    # zoom the y axis when all values are high, so small differences are visible
    lo = 0.0
    if vals and min(vals) > 0.9:
        lo = max(0.0, round(min(vals) - 0.01, 2))
    span = 1.0 - lo
    y = lambda v: pad_t + plot_h * (1 - (v - lo) / span)
    svg.append(f'<text x="{x0 + panel_w/2}" y="{pad_t - 8}" text-anchor="middle" font-weight="bold">{title}</text>')
    for k in range(5):
        v = lo + span * k / 4
        svg.append(f'<line x1="{x0}" x2="{x0 + panel_w - 20}" y1="{y(v):.1f}" y2="{y(v):.1f}" stroke="#ddd"/>')
        svg.append(f'<text x="{x0 - 4}" y="{y(v)+4:.1f}" text-anchor="end" fill="#555">{v:.3f}</text>')
    group_w = (panel_w - 30) / max(1, len(types))
    bar_w = group_w * 0.35
    for t_i, t in enumerate(types):
        gx = x0 + 5 + t_i * group_w
        for b_i, (k, d) in enumerate([("base", base), ("cand", cand)]):
            v = num(d.get(t, {}).get(col))
            bx = gx + group_w * 0.12 + b_i * bar_w
            if v is None:
                svg.append(f'<text x="{bx + bar_w/2}" y="{pad_t + plot_h - 4}" text-anchor="middle" fill="#999">n/a</text>')
                continue
            top = y(max(v, lo))
            svg.append(f'<rect x="{bx:.1f}" y="{top:.1f}" width="{bar_w:.1f}" height="{pad_t + plot_h - top:.1f}" fill="{colors[k]}"/>')
            svg.append(f'<text x="{bx + bar_w/2:.1f}" y="{top - 3:.1f}" text-anchor="middle" font-size="9">{v:.4f}</text>')
        svg.append(f'<text x="{gx + group_w*0.12 + bar_w:.1f}" y="{pad_t + plot_h + 16}" text-anchor="middle">{t}</text>')
    if lo > 0:
        svg.append(f'<text x="{x0 + panel_w/2 - 10}" y="{pad_t + plot_h + 34}" text-anchor="middle" fill="#999" font-size="10">y axis starts at {lo:.2f}</text>')

svg.append(f'<text x="{pad_l}" y="{H - 20}" fill="#555" font-size="11">Rows: comparison=GT, region=ALL, filter=ALL. Values from each run\'s Aardvark summary.tsv.</text>')
svg.append("</svg>")
open(out, "w").write("\n".join(svg))
print("wrote", out)
EOF
        python3 plot.py ~{baseline_summary} ~{candidate_summary} "~{baseline_name}" "~{candidate_name}" graph_comparison.svg
    >>>

    output {
        File plot = "graph_comparison.svg"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 5 SSD"
    }
}