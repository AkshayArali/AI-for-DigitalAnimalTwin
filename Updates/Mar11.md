## Week of Mar 11 — Progress Update

**1. Data acquisition and labeling**

- Implemented `src/data_download.py` to fully automate dataset downloads:
  - Tox21 (*in vitro* HTS): 7,831 compounds, 12 toxicity-related assays.
  - NCI Acute Toxicity Database: 80,081 compounds, 59 LD50 endpoints.
- Extracted **rat oral LD50** endpoint and converted values from `-log(mol/kg)` to `mg/kg` using RDKit molecular weights.
- Created **binary toxicity label** using the GHS 2000 mg/kg cutoff:
  - 10,189 compounds with valid rat oral LD50.
  - 63.1% labeled toxic (LD50 ≤ 2000 mg/kg), 36.9% non-toxic.
- Saved processed labels to `data/processed/ld50_labels.csv`.

**2. Feature engineering (“digital twin” feature space)**

- Implemented `src/feature_pipeline.py`:
  - Generated **Morgan fingerprints** (ECFP4, radius 2, 2048 bits) from SMILES.
  - Computed basic **physicochemical descriptors**: MW, LogP, TPSA, HBD, HBA, rotatable bonds.
  - Integrated all 12 **Tox21 assay readouts** as in-vitro biological features.
- Built a combined feature matrix:
  - `data/processed/merged_features.csv` with 7,823 molecules × 2,067 features.
- Performed **Murcko scaffold–aware splitting** to avoid structural leakage:
  - Train: 6,259 compounds
  - Validation: 783 compounds
  - Test: 781 compounds
  - Saved to `train.csv`, `val.csv`, `test.csv` in `data/processed/`.

**3. Environment, structure, and reproducibility**

- Created a dedicated **virtual environment** (`venv/`) and pinned dependencies in `requirements.txt`:
  - RDKit, pandas, numpy, scikit-learn, requests, tqdm, etc.
- Expanded `.gitignore` to exclude:
  - `data/raw/`, `data/processed/`, virtual environments, and notebook checkpoints.
- Added `data/README.md` documenting:
  - Data sources (Tox21 + NCI LD50) with URLs.
  - Processed file descriptions and reproducibility notes.
- Created `notebooks/01_data_exploration.ipynb` for:
  - Inspecting Tox21 assay sparsity and activation rates.
  - Reviewing LD50 label distribution (once labels are loaded).
  - Validating SMILES parsing.

**4. Modeling plan (next steps)**

- Baseline model: **XGBoost** using structure-only features (fingerprints + descriptors) to predict binary LD50 toxicity.
- Digital Twin model: **XGBoost and/or Graph Neural Network** using:
  - Molecular structure (fingerprints or graph embedding) **plus** Tox21 assay readouts.
- Evaluation focus:
  - AUPRC and sensitivity (recall) on the test set.
  - A “confidence threshold vs. tests-avoided” curve to estimate how many animal tests can be skipped at high sensitivity (e.g., ≥95%).