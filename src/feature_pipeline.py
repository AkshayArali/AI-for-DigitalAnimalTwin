"""
Feature engineering pipeline for the Digital Animal Twin project.

Generates:
    1. Morgan fingerprints (ECFP4, 2048-bit) from SMILES — chemical structure
    2. In-vitro assay feature matrix from Tox21 — biological context
    3. Merged feature set for the "Digital Twin" model

Usage:
    python src/feature_pipeline.py
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from tqdm import tqdm

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

TOX21_ASSAYS = [
    "NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase",
    "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma",
    "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53",
]

FP_RADIUS = 2       # ECFP4
FP_NBITS = 2048


# ── SMILES → RDKit Mol ──────────────────────────────────────────────────────

def smiles_to_mol(smiles: str):
    """Parse SMILES; return None for invalid strings."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    return Chem.MolFromSmiles(smiles.strip())


def validate_smiles(df: pd.DataFrame, smiles_col: str = "smiles") -> pd.DataFrame:
    """Drop rows with un-parseable SMILES and add a canonical SMILES column."""
    log.info("Validating SMILES for %d rows …", len(df))
    mols = [smiles_to_mol(s) for s in tqdm(df[smiles_col], desc="Parsing SMILES")]
    df = df.copy()
    df["mol"] = mols
    before = len(df)
    df = df.dropna(subset=["mol"])
    df["canon_smiles"] = df["mol"].apply(lambda m: Chem.MolToSmiles(m))
    log.info("Valid SMILES: %d / %d (dropped %d)", len(df), before, before - len(df))
    return df


# ── Morgan Fingerprints ─────────────────────────────────────────────────────

def compute_morgan_fingerprints(
    df: pd.DataFrame,
    radius: int = FP_RADIUS,
    n_bits: int = FP_NBITS,
) -> np.ndarray:
    """
    Compute Morgan (ECFP) circular fingerprints.

    Returns an (N, n_bits) numpy array of uint8.
    """
    log.info("Computing Morgan FP (radius=%d, bits=%d) for %d molecules …",
             radius, n_bits, len(df))
    fps = []
    for mol in tqdm(df["mol"], desc="Morgan FP"):
        arr = np.zeros(n_bits, dtype=np.uint8)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        fp_array = np.zeros(n_bits, dtype=np.uint8)
        Chem.DataStructs.ConvertToNumpyArray(fp, fp_array)
        fps.append(fp_array)
    return np.stack(fps)


def compute_descriptors(df: pd.DataFrame) -> pd.DataFrame:
    """Compute a handful of physicochemical descriptors (MW, LogP, TPSA, HBD, HBA)."""
    log.info("Computing physicochemical descriptors …")
    records = []
    for mol in tqdm(df["mol"], desc="Descriptors"):
        records.append({
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        })
    return pd.DataFrame(records)


# ── Assay Feature Matrix ────────────────────────────────────────────────────

def build_assay_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract the 12 Tox21 assay columns as the in-vitro feature block.

    Missing assay values are left as NaN (handled downstream by the model
    or imputed with median / -1 sentinel).
    """
    present = [c for c in TOX21_ASSAYS if c in df.columns]
    log.info("Found %d / %d Tox21 assay columns", len(present), len(TOX21_ASSAYS))
    assay_df = df[present].copy()
    pct_missing = assay_df.isna().mean().mean() * 100
    log.info("Assay matrix sparsity: %.1f%% missing", pct_missing)
    return assay_df


# ── Merge & Save ─────────────────────────────────────────────────────────────

def build_feature_matrix(tox21_path: Optional[Path] = None) -> pd.DataFrame:
    """
    End-to-end pipeline: load cleaned Tox21 → fingerprints + assays → save.

    Returns the merged feature DataFrame.
    """
    if tox21_path is None:
        tox21_path = PROCESSED_DIR / "tox21_clean.csv"

    df = pd.read_csv(tox21_path)
    df = validate_smiles(df)

    # Fingerprints
    fp_matrix = compute_morgan_fingerprints(df)
    fp_cols = [f"fp_{i}" for i in range(fp_matrix.shape[1])]
    fp_df = pd.DataFrame(fp_matrix, columns=fp_cols, index=df.index)

    # Physicochemical descriptors
    desc_df = compute_descriptors(df)
    desc_df.index = df.index

    # Assay readouts
    assay_df = build_assay_features(df)
    assay_df.index = df.index

    # Combine everything
    meta = df[["canon_smiles"]].copy()
    meta.index = df.index

    # Keep assay labels if present (for baseline/digital-twin comparison)
    label_cols = [c for c in TOX21_ASSAYS if c in df.columns]

    merged = pd.concat([meta, fp_df, desc_df, assay_df.add_prefix("assay_")], axis=1)

    out_path = PROCESSED_DIR / "merged_features.csv"
    merged.to_csv(out_path, index=False)
    log.info("Wrote merged feature matrix → %s  (%d rows × %d cols)",
             out_path.name, *merged.shape)

    return merged


# ── Scaffold-aware train/val/test split ──────────────────────────────────────

def scaffold_split(
    df: pd.DataFrame,
    train_frac: float = 0.8,
    val_frac: float = 0.1,
    seed: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Split by Murcko scaffold so structurally similar compounds stay together.
    Falls back to random split if scaffolds can't be computed.
    """
    from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles

    log.info("Computing Murcko scaffolds for splitting …")
    scaffolds: dict[str, list[int]] = {}
    for idx, smi in tqdm(df["canon_smiles"].items(), total=len(df), desc="Scaffolds"):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            scaf = "NONE"
        else:
            try:
                scaf = MurckoScaffoldSmiles(mol=mol, includeChirality=False)
            except Exception:
                scaf = "NONE"
        scaffolds.setdefault(scaf, []).append(idx)

    rng = np.random.RandomState(seed)
    scaffold_keys = sorted(scaffolds.keys(), key=lambda k: len(scaffolds[k]), reverse=True)

    train_idxs, val_idxs, test_idxs = [], [], []
    n = len(df)
    for key in scaffold_keys:
        idxs = scaffolds[key]
        if len(train_idxs) / n < train_frac:
            train_idxs.extend(idxs)
        elif len(val_idxs) / n < val_frac:
            val_idxs.extend(idxs)
        else:
            test_idxs.extend(idxs)

    train = df.loc[train_idxs]
    val = df.loc[val_idxs]
    test = df.loc[test_idxs]
    log.info("Split sizes — train: %d, val: %d, test: %d", len(train), len(val), len(test))
    return train, val, test


def main():
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    merged = build_feature_matrix()

    train, val, test = scaffold_split(merged)
    train.to_csv(PROCESSED_DIR / "train.csv", index=False)
    val.to_csv(PROCESSED_DIR / "val.csv", index=False)
    test.to_csv(PROCESSED_DIR / "test.csv", index=False)
    log.info("Saved train/val/test splits to data/processed/")


if __name__ == "__main__":
    main()
