import openpyxl, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np

XLSX = "/root/.claude/uploads/00ec2b76-52d9-5cc4-ad12-c189fdae984f/99969136-Table_Step1p5_L1TMenu.xlsx"
wb = openpyxl.load_workbook(XLSX, data_only=True)
s15, s1 = wb.worksheets[0], wb.worksheets[1]   # 'Step 1.5 July26 Workshop', 'Step 1 July26 Workshop - In Progress'
assert s15["E79"].value == "Total Menu" and s1["C65"].value == "Total:"

# label: (Step1 140PU col, Step1 200PU col, Step1.5 140PU col, Step1.5 200PU col, Run3 2024 col [% already]); None = not in sheet
samples = {
 "VBF H→ττ":  ("F","G","S","T","U"),
 "VBF H→bb":  (None,"H",None,"V","W"),
 "VBF H→inv": (None,"I",None,"X","Y"),
 "ggH→γγ":    ("J","K","Z","AA","AB"),
 "HH→bbττ":   (None,"L",None,"AC","AD"),
 "HH→bbbb":   (None,"M",None,"AE","AF"),
}
def val(ws, c, r, scale):
    v = ws[f"{c}{r}"].value if c else None
    return np.nan if v is None else float(v)*scale

series = [  # label, sheet, row, column index in tuple, scale, colour
 ("Step 1, 140 PU",   s1, 65, 0, 100, "#d3d8de"),
 ("Step 1.5, 140 PU", s15,79, 2, 100, "#a9cbe8"),
 ("Step 1, 200 PU (rate %.0f kHz)"%float(s1["E65"].value),   s1, 65, 1, 100, "#9aa5b1"),
 ("Step 1.5, 200 PU (rate %.0f kHz)"%float(s15["P79"].value), s15,79, 3, 100, "#1f6fb2"),
 ("Run 3 2024 menu (Run 3 samples)", s15, 79, 4, 1, "#e08a1e"),
]
fig, ax = plt.subplots(figsize=(12,4.8))
x = np.arange(len(samples)); w = 0.16
for i,(lab,ws,row,k,sc,col) in enumerate(series):
    v = np.array([val(ws,cols[k],row,sc) for cols in samples.values()])
    b = ax.bar(x+(i-2)*w, np.nan_to_num(v), w, label=lab, color=col)
    # "n/a" only where the sheet should have a value (200 PU / Run 3); 140 PU simply doesn't exist for those samples
    ax.bar_label(b, labels=[("" if "140" in lab else "n/a") if np.isnan(t) else "%.0f"%t for t in v], fontsize=8, padding=2)
ax.set_xticks(x); ax.set_xticklabels(samples.keys())
ax.set_ylabel("Total L1 menu efficiency [%]"); ax.set_ylim(0,125)
ax.spines[["top","right"]].set_visible(False)
ax.legend(frameon=False,fontsize=8,ncol=2,loc="upper right")
ax.set_title("L1 menu signal efficiency: Run 3 2024 vs Phase-2 Step 1 / Step 1.5 (170pre3, July26 workshop)",fontsize=11)
fig.tight_layout(); fig.savefig("plots/menu_efficiency_step1_vs_step1p5.png",dpi=150)
