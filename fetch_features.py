import time, requests, pandas as pd

# record path -> (column name, is a HIGHER score more damaging?)
FEATURES = {
    "dbnsfp.gerp.91_mammals.score":             ("gerp",          True),
    "dbnsfp.phylop.100way_vertebrate.score":    ("phylop",        True),
    "dbnsfp.phastcons.100way_vertebrate.score": ("phastcons",     True),
    "dbnsfp.siphy_29way.logodds_score":         ("siphy",         True),
    "dbnsfp.sift.score":                        ("sift",          False),
    "dbnsfp.polyphen2.hdiv.score":              ("polyphen2",     True),
    "dbnsfp.provean.score":                     ("provean",       False),
    "dbnsfp.mutationassessor.score":            ("mutassessor",   True),
    "dbnsfp.cadd.phred":                        ("cadd",          True),
    "dbnsfp.revel.score":                       ("revel",         True),
    "dbnsfp.metasvm.score":                     ("metasvm",       True),
    "dbnsfp.alphamissense.score":               ("alphamissense", True),
}
FIELDS = ",".join(FEATURES)
URL = "https://myvariant.info/v1/query"
BATCH = 1000

def dig(rec, path):
    d = rec
    for key in path.split("."):
        if not isinstance(d, dict):
            return None
        d = d.get(key)
        if d is None:
            return None
    return d

def collapse(value, higher_is_damaging):
    """One variant can sit in several transcripts, so dbNSFP returns a list of
    scores. Keep the most damaging transcript."""
    if value is None:
        return None
    if isinstance(value, list):
        nums = [v for v in value if isinstance(v, (int, float))]
        if not nums:
            return None
        return max(nums) if higher_is_damaging else min(nums)
    return value if isinstance(value, (int, float)) else None

def fetch(ids, retries=3):
    for attempt in range(retries):
        try:
            r = requests.post(URL, timeout=180, data={
                "q": ",".join(ids), "scopes": "_id",
                "fields": FIELDS, "assembly": "hg38"})
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if attempt == retries - 1:
                raise
            print(f"    retry {attempt+1} after {type(e).__name__}")
            time.sleep(5 * (attempt + 1))

def annotate(df, name):
    df = df.copy()
    df["_id"] = ("chr" + df.Chromosome.astype(str) + ":g." + df.PositionVCF.astype(str)
                 + df.ReferenceAlleleVCF.astype(str) + ">" + df.AlternateAlleleVCF.astype(str))
    ids = df._id.tolist()
    rows, t0 = {}, time.time()
    for i in range(0, len(ids), BATCH):
        for rec in fetch(ids[i:i + BATCH]):
            if rec.get("notfound"):
                continue
            rows[rec["query"]] = {col: collapse(dig(rec, path), hi)
                                  for path, (col, hi) in FEATURES.items()}
        print(f"  {name}: {min(i+BATCH, len(ids))}/{len(ids)}", end="\r")
    print(f"  {name}: {len(ids)}/{len(ids)} in {time.time()-t0:.0f}s      ")

    out = df.join(pd.DataFrame.from_dict(rows, orient="index"), on="_id")
    print(f"  coverage (of {len(out)}):")
    for col, _ in FEATURES.values():
        print(f"    {col:<14} {out[col].notna().sum():>6}  ({100*out[col].notna().mean():.0f}%)")
    return out

if __name__ == "__main__":
    lab = annotate(pd.read_csv("data/labelled.csv"), "labelled")
    lab.to_csv("data/features_labelled.csv", index=False)
    vus = annotate(pd.read_csv("data/vus.csv"), "vus")
    vus.to_csv("data/features_vus.csv", index=False)
    print(f"\nsaved features_labelled.csv {lab.shape} and features_vus.csv {vus.shape}")
