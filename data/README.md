# Data Sources

Raw data files are **not** tracked in git (see `.gitignore`).
Run `python src/data_download.py` to fetch everything automatically.

## Primary Datasets

### 1. Tox21 (MoleculeNet mirror)

| Field | Value |
|-------|-------|
| **URL** | `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz` |
| **Records** | ~8,000 compounds |
| **Labels** | 12 nuclear-receptor and stress-response pathway assay outcomes (binary) |
| **Format** | CSV (SMILES + 12 label columns) |
| **Use** | *In vitro* HTS readouts — the core "digital twin" feature set |

### 2. CATMoS — Acute Oral Toxicity (LD50)

| Field | Value |
|-------|-------|
| **Source** | NICEATM / ICE — Collaborative Acute Toxicity Modeling Suite |
| **URL** | `https://ice.ntp.niehs.nih.gov/` (search "CATMoS") |
| **Records** | ~11,000 compounds with rat oral LD50 values |
| **Format** | CSV/XLSX |
| **Use** | *In vivo* ground-truth labels — binary toxic/non-toxic based on GHS categories |

### 3. EPA CompTox — DSSTox Substance Identifiers

| Field | Value |
|-------|-------|
| **URL** | `https://comptox.epa.gov/dashboard/` |
| **Use** | DTXSID cross-referencing to link Tox21 assay IDs to CATMoS LD50 records |

## Derived Files (in `data/processed/`)

| File | Description |
|------|-------------|
| `tox21_clean.csv` | Tox21 with valid SMILES, de-duplicated |
| `ld50_labels.csv` | CATMoS LD50 values mapped to binary labels |
| `merged_features.csv` | Combined: Morgan fingerprints + assay readouts + LD50 label |
| `train.csv` / `val.csv` / `test.csv` | Scaffold-aware splits |

## Reproducibility

- **Access date:** *(fill in when you download)*
- **Python environment:** see `requirements.txt`
- **Download script:** `python src/data_download.py`
