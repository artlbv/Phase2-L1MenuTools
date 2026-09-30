import openpyxl, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np

XLSX = "/root/.claude/uploads/00ec2b76-52d9-5cc4-ad12-c189fdae984f/99969136-Table_Step1p5_L1TMenu.xlsx"
wb = openpyxl.load_workbook(XLSX, data_only=True)
s15, s1 = wb.worksheets[0], wb.worksheets[1]   # 'Step 1.5 July26 Workshop', 'Step 1 July26 Workshop - In Progress'
assert s15["E79"].value == "Total Menu" and s1["C65"].value == "Total:"

# (label, Step1 col, Step1.5 col), 200 PU unless noted
# (label, Step1 col, Step1.5 col, Run3 2024 col in Step 1.5 tab [% already]); Run 3 is None if not in the sheet
samples200 = [("VBF H→ττ","G","T","U"),("VBF H→bb","H","V","W"),("VBF H→inv","I","X","Y"),
              ("ggH→γγ","K","AA","AB"),("HH→bbττ","L","AC","AD"),("HH→bbbb","M","AE","AF")]
samples140 = [("VBF H→ττ","F","S",None),("ggH→γγ","J","Z",None)]
g = lambda ws,c,r: float(ws[f"{c}{r}"].value)*100

fig, axs = plt.subplots(1,2,figsize=(11,4.6),gridspec_kw={"width_ratios":[6,2]},sharey=True)
for ax,samples,title in [(axs[0],samples200,"Spring24 200 PU"),(axs[1],samples140,"Spring24 140 PU")]:
    x=np.arange(len(samples)); w=0.27
    v1=[g(s1,a,65) for _,a,_,_ in samples]; v15=[g(s15,b,79) for _,_,b,_ in samples]
    vr3=[float(s15[f"{c}79"].value) if c and s15[f"{c}79"].value is not None else np.nan for *_,c in samples]
    bars=[(-w,v1,"Step 1 (Phase-2, 200PU rate: %.0f kHz)"%float(s1["E65"].value),"#9aa5b1"),
          (0,v15,"Step 1.5 (Phase-2, 200PU rate: %.0f kHz)"%float(s15["P79"].value),"#1f6fb2")]
    if ax is axs[0]: bars.append((w,vr3,"Run 3 2024 menu (Run 3 samples)","#e08a1e"))
    for off,v,lab,col in bars:
        b=ax.bar(x+off,np.nan_to_num(v),w,label=lab if ax is axs[0] else None,color=col)
        ax.bar_label(b,labels=["n/a" if np.isnan(t) else "%.0f"%t for t in v],fontsize=8,padding=2)
    ax.set_xticks(x); ax.set_xticklabels([s[0] for s in samples],rotation=25,ha="right")
    ax.set_title(title,fontsize=10); ax.spines[["top","right"]].set_visible(False)
axs[0].set_ylabel("Total L1 menu efficiency [%]"); axs[0].set_ylim(0,105)
axs[0].legend(frameon=False,fontsize=8,loc="upper left",bbox_to_anchor=(0,1.0))
fig.suptitle("L1 menu signal efficiency: Run 3 2024 vs Phase-2 Step 1 / Step 1.5 (170pre3, July26 workshop)",fontsize=11)
fig.tight_layout(); fig.savefig("plots/menu_efficiency_step1_vs_step1p5.png",dpi=150)
for lab,a,b,c in samples200+samples140: print(lab, round(g(s1,a,65),1), round(g(s15,b,79),1), s15[f"{c}79"].value if c else None)
