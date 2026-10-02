import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

d = pd.read_csv("data/features_labelled.csv")
y, groups = d.label.values, d.GeneSymbol.values

CONS    = ["gerp","phylop","phastcons","siphy"]          # evolution only
CLASSIC = ["sift","polyphen2","provean","mutassessor"]   # classic predictors
ENSEM   = ["cadd","revel","metasvm","alphamissense"]     # modern ensembles
TIERS = {"T1 conservation": CONS,
         "T2 +classic":     CONS + CLASSIC,
         "T3 +ensembles":   CONS + CLASSIC + ENSEM}

gkf = GroupKFold(n_splits=5)
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
SPLITS = {"gene-wise (correct)": (gkf, {"groups": groups}),
          "variant-wise (leaky)": (skf, {})}

print(f"{len(d)} variants | {d.GeneSymbol.nunique()} genes | {y.mean():.3f} pathogenic\n")

def run(X, splitter, split_args):
    aucs, bals = [], []
    for tr, te in splitter.split(X, y, **split_args):
        # imputer inside the pipeline: medians learned per fold, never from test data
        pipe = make_pipeline(SimpleImputer(strategy="median"),
                             RandomForestClassifier(n_estimators=400, random_state=42, n_jobs=-1))
        pipe.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], pipe.predict_proba(X[te])[:, 1]))
        bals.append(balanced_accuracy_score(y[te], pipe.predict(X[te])))
    return np.mean(aucs), np.std(aucs), np.mean(bals), np.std(bals)

rows = []
print(f"{'model':<20}{'split':<24}{'ROC-AUC':>16}{'bal acc':>16}")
print("-" * 76)
for tname, cols in TIERS.items():
    for sname, (sp, args) in SPLITS.items():
        a, asd, b, bsd = run(d[cols].values, sp, args)
        rows.append(dict(model=tname, split=sname, auc=a, auc_sd=asd, bal=b, bal_sd=bsd))
        print(f"{tname:<20}{sname:<24}{a:>10.3f}\u00b1{asd:.3f}{b:>10.3f}\u00b1{bsd:.3f}")
    print()

# The probe: a model given ONLY gene identity - nothing about the variant at all.
# Under a correct split it must score chance. Under a leaky split it does not.
Xg = pd.get_dummies(d.GeneSymbol).values.astype(float)
print("gene-identity-only probe (no variant information whatsoever)")
for sname, (sp, args) in SPLITS.items():
    a, asd, b, bsd = run(Xg, sp, args)
    rows.append(dict(model="gene identity only", split=sname, auc=a, auc_sd=asd, bal=b, bal_sd=bsd))
    print(f"{'gene id only':<20}{sname:<24}{a:>10.3f}\u00b1{asd:.3f}{b:>10.3f}\u00b1{bsd:.3f}")

print("\nsingle-score baselines (same gene-wise folds):")
for s in ["revel","alphamissense","cadd","phylop"]:
    aucs = []
    for tr, te in gkf.split(d[[s]].values, y, groups=groups):
        m = ~np.isnan(d[s].values[te])
        if m.sum() > 10 and len(np.unique(y[te][m])) > 1:
            aucs.append(roc_auc_score(y[te][m], d[s].values[te][m]))
    rows.append(dict(model=f"{s} alone", split="gene-wise (correct)",
                     auc=np.mean(aucs), auc_sd=np.std(aucs), bal=np.nan, bal_sd=np.nan))
    print(f"  {s:<16}{np.mean(aucs):.3f} \u00b1 {np.std(aucs):.3f}")

dm, bals = DummyClassifier(strategy="most_frequent"), []
for tr, te in gkf.split(d[CONS].values, y, groups=groups):
    dm.fit(d[CONS].values[tr], y[tr])
    bals.append(balanced_accuracy_score(y[te], dm.predict(d[CONS].values[te])))
print(f"  {'majority class':<16}bal acc {np.mean(bals):.3f}")

pd.DataFrame(rows).to_csv("data/results.csv", index=False)
print("\nsaved data/results.csv")
