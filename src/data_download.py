"""
Download and prepare raw datasets for the Digital Animal Twin project.

Datasets:
    1. Tox21 (MoleculeNet mirror) — 12 in-vitro HTS assay labels + SMILES
    2. NCI Acute Toxicity DB — rat oral LD50 (in-vivo ground truth)

Both are downloaded automatically — no manual steps required.

Usage:
    python src/data_download.py
"""

from __future__ import annotations

import gzip
import logging
from datetime import date
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = PROJECT_ROOT / "data" / "raw"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

TOX21_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz"
NCI_LD50_URL = "https://cactus.nci.nih.gov/download/acute-toxicity-db/LD50_dataset.csv"

TOX21_ASSAYS = [
    "NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase",
    "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma",
    "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53",
]

RAT_ORAL_COL = "rat_oral_LD50_(?log(mol/kg))"


def _download(url: str, dest: Path) -> None:
    """Stream-download a file with progress bar."""
    if dest.exists():
        log.info("Already cached: %s", dest.name)
        return
    log.info("Downloading %s …", url)
    resp = requests.get(url, stream=True, timeout=120)
    resp.raise_for_status()
    total = int(resp.headers.get("content-length", 0))
    with open(dest, "wb") as f, tqdm(total=total, unit="B", unit_scale=True) as bar:
        for chunk in resp.iter_content(chunk_size=8192):
            f.write(chunk)
            bar.update(len(chunk))
    log.info("Saved → %s", dest)


# ── Tox21 ────────────────────────────────────────────────────────────────────

def download_tox21() -> pd.DataFrame:
    """Fetch Tox21 CSV from the MoleculeNet S3 mirror and return a DataFrame."""
    gz_path = RAW_DIR / "tox21.csv.gz"
    _download(TOX21_URL, gz_path)

    with gzip.open(gz_path, "rt") as f:
        df = pd.read_csv(f)

    log.info("Tox21 raw shape: %s", df.shape)
    return df


def clean_tox21(df: pd.DataFrame) -> pd.DataFrame:
    """Basic cleaning: drop rows without SMILES, rename columns."""
    smiles_col = [c for c in df.columns if "smiles" in c.lower()][0]
    df = df.rename(columns={smiles_col: "smiles"})
    df = df.dropna(subset=["smiles"])
    df = df.drop_duplicates(subset=["smiles"])
    log.info("Tox21 after dedup: %d compounds, %d assay columns",
             len(df), len([c for c in df.columns if c in TOX21_ASSAYS]))
    return df


# ── NCI Acute Toxicity DB (rat oral LD50) ────────────────────────────────────

def download_ld50() -> pd.DataFrame:
    """
    Download the NCI Acute Toxicity Database (80k compounds, 59 endpoints)
    and extract the rat oral LD50 subset (~10k compounds).

    Source: https://cactus.nci.nih.gov/download/acute-toxicity-db/
    Reference: Zhu et al., Chem. Res. Toxicol. 2009
    """
    csv_path = RAW_DIR / "ld50_nci.csv"
    _download(NCI_LD50_URL, csv_path)
    df = pd.read_csv(csv_path, low_memory=False)
    log.info("NCI LD50 raw shape: %s (59 endpoint columns)", df.shape)
    return df


def clean_ld50(df: pd.DataFrame, threshold_mgkg: float = 2000.0) -> pd.DataFrame:
    """
    Extract rat oral LD50 values, convert from -log(mol/kg) to mg/kg,
    and assign a binary toxic label using the GHS 2000 mg/kg cutoff.

    Conversion:
        LD50 (mol/kg)  = 10^(-value)
        LD50 (mg/kg)   = LD50 (mol/kg) * MW (g/mol) * 1000
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    sub = df[["SMILES", RAT_ORAL_COL]].dropna(subset=[RAT_ORAL_COL]).copy()
    sub = sub.rename(columns={"SMILES": "smiles", RAT_ORAL_COL: "neg_log_molkg"})
    sub = sub.dropna(subset=["smiles"])
    sub = sub.drop_duplicates(subset=["smiles"])
    log.info("Rat oral LD50 records with SMILES: %d", len(sub))

    log.info("Computing molecular weights for unit conversion …")
    mws = []
    for smi in tqdm(sub["smiles"], desc="MW calc"):
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        mws.append(Descriptors.MolWt(mol) if mol else np.nan)
    sub["mw"] = mws

    before = len(sub)
    sub = sub.dropna(subset=["mw"])
    log.info("Valid MW: %d / %d", len(sub), before)

    # -log(mol/kg) → mol/kg → mg/kg
    sub["ld50_molkg"] = 10.0 ** (-sub["neg_log_molkg"])
    sub["ld50_mgkg"] = sub["ld50_molkg"] * sub["mw"] * 1000

    # GHS binary label: Cat 1-4 (<=2000 mg/kg) → toxic=1
    sub["toxic"] = (sub["ld50_mgkg"] <= threshold_mgkg).astype(int)

    n_toxic = sub["toxic"].sum()
    log.info("LD50 labeling (threshold=%d mg/kg): %d toxic / %d total (%.1f%%)",
             threshold_mgkg, n_toxic, len(sub), 100 * n_toxic / len(sub))

    out = sub[["smiles", "mw", "ld50_mgkg", "toxic"]].reset_index(drop=True)
    return out


# ── Orchestration ────────────────────────────────────────────────────────────

def main():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # ── Tox21 (in-vitro HTS) ──
    tox21_raw = download_tox21()
    tox21 = clean_tox21(tox21_raw)
    tox21.to_csv(PROCESSED_DIR / "tox21_clean.csv", index=False)
    log.info("Wrote tox21_clean.csv (%d rows)", len(tox21))

    # ── Rat oral LD50 (in-vivo ground truth) ──
    ld50_raw = download_ld50()
    ld50 = clean_ld50(ld50_raw)
    ld50.to_csv(PROCESSED_DIR / "ld50_labels.csv", index=False)
    log.info("Wrote ld50_labels.csv (%d rows)", len(ld50))

    log.info("Access date: %s", date.today().isoformat())
    log.info("Done. Check data/processed/ for cleaned files.")


if __name__ == "__main__":
    main()
