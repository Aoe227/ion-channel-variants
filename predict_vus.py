import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline

lab = pd.read_csv("data/features_labelled.csv")
vus = pd.read_csv("data/features_vus.csv")

# T2 deliberately, not T3: T3's importance is 82% AlphaMissense + REVEL, so scoring
# VUS with it would just be re-reporting AlphaMissense under another name.
COLS = ["gerp","phylop","phastcons","siphy","sift","polyphen2","provean","mutassessor"]

pipe = make_pipeline(SimpleImputer(strategy="median"),
                     RandomForestClassifier(n_estimators=400, random_state=42, n_jobs=-1))
pipe.fit(lab[COLS].values, lab.label.values)
vus["p_path"] = pipe.predict_proba(vus[COLS].values)[:, 1]
p = vus.p_path.values

print(f"{len(vus):,} VUS scored\n")
for lo, hi, name in [(0,.1,"very likely benign"), (.1,.3,"likely benign"), (.3,.7,"uncertain"),
                     (.7,.9,"likely pathogenic"), (.9,1.01,"very likely pathogenic")]:
    n = ((p >= lo) & (p < hi)).sum()
    print(f"  {name:<24}{n:>7,}  ({100*n/len(p):4.1f}%)")
print(f"\nmedian {np.median(p):.3f}  mean {p.mean():.3f}   "
      f"(labelled set is {lab.label.mean():.1%} pathogenic)")

t = vus.groupby("GeneSymbol").p_path.agg(n="size", mean_p="mean").query("n >= 100")
print("\nhighest mean VUS pathogenicity:"); print(t.nlargest(5,"mean_p").round(3).to_string())
print("\nlowest:");                          print(t.nsmallest(4,"mean_p").round(3).to_string())

vus[["VariationID","GeneSymbol","Name","p_path"]].to_csv("data/vus_predictions.csv", index=False)
print("\nsaved data/vus_predictions.csv")
