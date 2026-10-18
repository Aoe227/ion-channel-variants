import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")      # write files, never open a blocking window
import matplotlib.pyplot as plt
from pathlib import Path
Path("figures").mkdir(exist_ok=True)

BLUE, CORAL, INK, MUTED, GRID = "#185FA5", "#D85A30", "#2C2C2A", "#73726C", "#E5E4E0"
plt.rcParams.update({"font.size":10,"axes.edgecolor":MUTED,"axes.labelcolor":INK,
    "text.color":INK,"xtick.color":MUTED,"ytick.color":MUTED,"axes.spines.top":False,
    "axes.spines.right":False,"figure.facecolor":"white","axes.facecolor":"white"})

# Figure 1 - the leakage demonstration
res = pd.read_csv("data/results.csv")
order = ["gene identity only","T1 conservation","T2 +classic","T3 +ensembles"]
lab   = ["gene identity\nonly","T1\nconservation","T2\n+classic","T3\n+ensembles"]
gw = [res.query("model==@m and split=='gene-wise (correct)'").auc.iloc[0] for m in order]
vw = [res.query("model==@m and split=='variant-wise (leaky)'").auc.iloc[0] for m in order]
ye, h = np.arange(4), 0.34
fig, ax = plt.subplots(figsize=(7.8,4.3))
b1 = ax.barh(ye+h/2+0.012, vw, h, color=CORAL, label="variant-wise split (leaky)")
b2 = ax.barh(ye-h/2-0.012, gw, h, color=BLUE,  label="gene-wise split (correct)")
for i,(a,b) in enumerate(zip(gw,vw)):
    ax.text(b+.008, i+h/2+0.012, f"{b:.3f}", va="center", fontsize=9, color=MUTED)
    ax.text(a+.008, i-h/2-0.012, f"{a:.3f}", va="center", fontsize=9, color=MUTED)
ax.axvline(0.5, color=MUTED, lw=1, ls="--")
ax.set_yticks(ye, lab); ax.set_xlim(0.45,1.05); ax.set_ylim(3.72,-0.92)
ax.annotate("chance", xy=(0.504,-0.86), fontsize=8.5, color=MUTED, va="top")
ax.set_xlabel("ROC-AUC"); ax.grid(axis="x", lw=.4, color=GRID); ax.set_axisbelow(True)
ax.legend(handles=[b2,b1], frameon=False, fontsize=9, ncol=2, loc="lower left",
          bbox_to_anchor=(0,1.01,1,.1))
fig.suptitle("Leakage, measured: a model given only gene identity scores 0.855\n"
             "under the naive split and exactly chance under the correct one",
             x=0.013, ha="left", fontsize=11.5, y=0.99)
fig.tight_layout(rect=[0,0,1,0.88]); fig.savefig("figures/leakage.png", dpi=160)

# Figure 2 - what the honest model uses
imp = pd.read_csv("data/importances.csv").query("tier=='T2'").sort_values("importance")
fig, ax = plt.subplots(figsize=(7,3.8))
ax.barh(imp.feature, imp.importance, 0.62, color=BLUE)
for i,v in enumerate(imp.importance): ax.text(v+.006, i, f"{v:.2f}", va="center", fontsize=9, color=MUTED)
ax.set_xlim(0, imp.importance.max()*1.18)
ax.set_xlabel("share of permutation importance (AUC drop on held-out genes)")
ax.grid(axis="x", lw=.4, color=GRID); ax.set_axisbelow(True)
ax.set_title("What the honest model relies on", loc="left", pad=10)
fig.tight_layout(); fig.savefig("figures/importance.png", dpi=160)

# Figure 3 - VUS predictions
v = pd.read_csv("data/vus_predictions.csv")
fig, ax = plt.subplots(figsize=(7,3.6))
ax.hist(v.p_path, bins=40, color=BLUE, edgecolor="white", linewidth=.5)
ax.axvline(0.5, color=MUTED, lw=1, ls="--")
ax.set_xlabel("predicted probability pathogenic"); ax.set_ylabel("variants")
ax.grid(axis="y", lw=.4, color=GRID); ax.set_axisbelow(True)
ax.set_title(f"{len(v):,} variants of uncertain significance, scored by the T2 model",
             loc="left", pad=10)
fig.tight_layout(); fig.savefig("figures/vus_distribution.png", dpi=160)
print("saved 3 figures")
