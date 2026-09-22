import re
import pandas as pd
from pathlib import Path

PANEL = {
    # cardiac
    "SCN5A","KCNQ1","KCNH2","KCNE1","KCNE2","KCNJ2","CACNA1C","RYR2","CASQ2",
    "SCN1B","SCN4B","HCN4","KCNJ5","TRDN",
    # neuro
    "SCN1A","SCN2A","SCN3A","SCN8A","SCN9A","SCN11A","KCNQ2","KCNQ3","KCNA1",
    "KCNA2","KCNB1","KCNT1","KCNMA1","CACNA1A","CACNA1H","HCN1","GRIN2A","GRIN2B",
    # muscle
    "SCN4A","CACNA1S","CLCN1","RYR1",
    # other
    "CFTR","KCNJ11","ABCC8","TRPV4",
}

COLS = ["GeneSymbol","ClinicalSignificance","ReviewStatus","Type","Name","Assembly",
        "Chromosome","PositionVCF","ReferenceAlleleVCF","AlternateAlleleVCF","VariationID"]

# Only labels backed by stated assertion criteria. Drops "no assertion criteria
# provided" and, deliberately, "criteria provided, conflicting classifications".
GOOD_REVIEW = {
    "criteria provided, single submitter",
    "criteria provided, multiple submitters, no conflicts",
    "reviewed by expert panel",
    "practice guideline",
}

PATHOGENIC = {"Pathogenic","Likely pathogenic","Pathogenic/Likely pathogenic"}
BENIGN     = {"Benign","Likely benign","Benign/Likely benign"}

MIN_PER_GENE = 30          # a gene needs enough variants to work as a CV group
MISSENSE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")

# ---- read in chunks: the file is ~4 GB uncompressed, so never load it whole ----
kept = []
for chunk in pd.read_csv("data/variant_summary.txt.gz", sep="\t", usecols=COLS,
                         chunksize=500_000, low_memory=False, dtype=str):
    kept.append(chunk[
        (chunk.Assembly == "GRCh38") &          # each variant is listed twice: 37 AND 38
        (chunk.GeneSymbol.isin(PANEL)) &
        (chunk.Type == "single nucleotide variant")
    ])
df = pd.concat(kept, ignore_index=True)
print(f"panel SNVs (GRCh38):        {len(df):,}")

df = df[df.ReviewStatus.isin(GOOD_REVIEW)]
print(f"with confident review status:{len(df):,}")

# missense = protein change of the form p.AaaNNNBbb, excluding stop-gain
aa = df.Name.fillna("").str.extract(MISSENSE)
df = df.assign(aa_ref=aa[0], aa_pos=aa[1], aa_alt=aa[2])
df = df[df.aa_ref.notna() & (df.aa_alt != "Ter") & (df.aa_ref != df.aa_alt)]
print(f"missense:                   {len(df):,}")

sig = df.ClinicalSignificance.str.split(";").str[0].str.strip()
df = df.assign(sig=sig)

vus = df[sig == "Uncertain significance"].copy()
lab = df[sig.isin(PATHOGENIC | BENIGN)].copy()
lab["label"] = lab.sig.isin(PATHOGENIC).astype(int)     # 1 = pathogenic, 0 = benign

counts = lab.GeneSymbol.value_counts()
lab = lab[lab.GeneSymbol.isin(counts[counts >= MIN_PER_GENE].index)]

Path("data").mkdir(exist_ok=True)
lab.to_csv("data/labelled.csv", index=False)
vus.to_csv("data/vus.csv", index=False)

print(f"\nlabelled: {len(lab):,}  "
      f"(pathogenic {lab.label.sum():,} / benign {(1-lab.label).sum():,})")
print(f"genes (CV groups): {lab.GeneSymbol.nunique()}")
print(f"VUS held out: {len(vus):,}")
print(f"majority-class baseline: {max(lab.label.mean(), 1-lab.label.mean()):.3f}")

rates = lab.groupby("GeneSymbol").label.agg(n="size", frac_path="mean")
print("\nper-gene pathogenic rate (the reason for gene-wise splitting):")
print(rates.sort_values("frac_path").to_string(float_format=lambda x: f"{x:.2f}"))