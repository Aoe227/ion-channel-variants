import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import GroupKFold
from sklearn.inspection import permutation_importance

d = pd.read_csv("data/features_labelled.csv")
y, g = d.label.values, d.GeneSymbol.values
CONS    = ["gerp","phylop","phastcons","siphy"]
CLASSIC = ["sift","polyphen2","provean","mutassessor"]
ENSEM   = ["cadd","revel","metasvm","alphamissense"]
TIERS = {"T1": CONS, "T2": CONS+CLASSIC, "T3": CONS+CLASSIC+ENSEM}

rows = []
for tier, cols in TIERS.items():
    X = d[cols].values
    imp = np.zeros(len(cols))
    # Permutation importance measured on HELD-OUT genes: how much AUC is lost when
    # a feature is shuffled. Unlike impurity importance this is computed out-of-sample.
    for tr, te in GroupKFold(n_splits=5).split(X, y, groups=g):
        pipe = make_pipeline(SimpleImputer(strategy="median"),
                             RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=-1)).fit(X[tr], y[tr])
        imp += permutation_importance(pipe, X[te], y[te], n_repeats=5,
                                      random_state=42, scoring="roc_auc", n_jobs=-1).importances_mean
    imp = np.clip(imp / 5, 0, None)
    imp = imp / imp.sum()
    print(f"\n{tier} permutation importance (AUC drop on held-out genes, normalised):")
    for i in np.argsort(imp)[::-1]:
        print(f"   {cols[i]:<15}{imp[i]:.3f}")
        rows.append(dict(tier=tier, feature=cols[i], importance=imp[i]))

pd.DataFrame(rows).to_csv("data/importances.csv", index=False)
print("\nsaved data/importances.csv")
