#!/usr/bin/env python3
"""Hybrid raw-first + legacy-rich cancer drug-response pipeline."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import subprocess
import warnings
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET

import h5py
import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.decomposition import TruncatedSVD
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import average_precision_score, mean_absolute_error, mean_squared_error, r2_score, roc_auc_score
from sklearn.model_selection import GroupKFold, KFold
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


ROOT = Path(__file__).resolve().parents[1]

KEY_COLUMNS = {
    "pair_id",
    "sample_id",
    "cell_line_name",
    "model_id",
    "SANGER_MODEL_ID",
    "canonical_drug_id",
    "DRUG_ID",
    "drug_name",
    "TCGA_DESC",
    "cohort",
    "label_regression",
    "label_binary",
    "label_main",
    "label_aux",
    "label_main_type",
    "label_aux_type",
    "binary_threshold",
    "LN_IC50",
    "AUC",
    "Z_SCORE",
    "canonical_smiles",
    "drug__canonical_smiles",
    "drug__target_list",
    "target",
    "target_pathway",
    "putative_target",
    "pathway_name",
}

SMILES_DERIVED_PREFIXES = (
    "drug_morgan_",
    "smiles_svd_",
    "drug__descriptor_",
)

GENERIC_TARGET_TOKENS = {
    "ABL",
    "ACTIVITY",
    "AGONIST",
    "ANTAGONIST",
    "BROMODOMAIN",
    "CELL",
    "CYTOSKELETON",
    "DESTABILISER",
    "EPIGENETIC",
    "GENOME",
    "HISTONE",
    "INHIBITOR",
    "INHIBITORS",
    "KINASE",
    "MICROTUBULE",
    "OTHER",
    "PATHWAY",
    "PROTEIN",
    "SIGNALING",
    "STABILISER",
    "TARGET",
    "UNKNOWN",
}

TARGET_ALIASES = {
    "PI3KALPHA": "PIK3CA",
    "PI3KBETA": "PIK3CB",
    "PI3KDELTA": "PIK3CD",
    "PI3KGAMMA": "PIK3CG",
    "MEK1": "MAP2K1",
    "MEK2": "MAP2K2",
    "ERK1": "MAPK3",
    "ERK2": "MAPK1",
    "BCLXL": "BCL2L1",
    "VEGFR": "KDR",
}


@dataclass
class RunPaths:
    root: Path
    raw_cache: Path
    run_dir: Path
    step1: Path
    step2: Path
    step3: Path
    step4: Path
    reports: Path


class ResidualBlock(nn.Module):
    def __init__(self, dim: int, dropout: float) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(dim, dim),
        )
        self.norm = nn.LayerNorm(dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.norm(x + self.net(x))


class ResidualMLP(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int = 128, blocks: int = 3, dropout: float = 0.1) -> None:
        super().__init__()
        layers: list[nn.Module] = [nn.Linear(input_dim, hidden_dim), nn.ReLU(), nn.Dropout(dropout)]
        layers.extend(ResidualBlock(hidden_dim, dropout) for _ in range(blocks))
        layers.append(nn.Linear(hidden_dim, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class FlatMLP(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int = 192, dropout: float = 0.15) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--stage",
        default="all",
        choices=["all", "inventory", "download", "build", "train", "report", "upload"],
    )
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--upload-s3", action="store_true")
    parser.add_argument("--variants", default="")
    parser.add_argument("--models", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    paths = make_paths(args.root, config)
    ensure_dirs(paths)
    stages = expand_stage(args.stage)

    if "inventory" in stages:
        write_source_inventory(config, paths)
    if "download" in stages and not args.skip_download:
        download_sources(config, paths)
    if "build" in stages:
        build_all_inputs(config, paths)
    if "train" in stages:
        train_variants(config, paths, args)
    if "report" in stages:
        write_report(config, paths)
    if args.upload_s3 or "upload" in stages:
        upload_outputs(config, paths)

    print(json.dumps(build_index(config, paths), indent=2, ensure_ascii=False))


def load_config(path: Path) -> dict[str, Any]:
    config = json.loads(path.read_text(encoding="utf-8"))
    if "extends" not in config:
        return config
    parent_path = path.parent / str(config.pop("extends"))
    parent = load_config(parent_path)
    return deep_update(parent, config)


def make_paths(root: Path, config: dict[str, Any]) -> RunPaths:
    cancer = str(config["cancer_id"])
    run_id = str(config["run_id"])
    run_dir = root / "data/processed/runs" / run_id
    return RunPaths(
        root=root,
        raw_cache=root / "data/raw_cache" / f"{cancer}_raw",
        run_dir=run_dir,
        step1=run_dir / "step1_raw_inventory",
        step2=run_dir / "step2_intermediate",
        step3=run_dir / "step3_features",
        step4=run_dir / "step4_model_inputs",
        reports=root / "outputs/reports" / run_id,
    )


def ensure_dirs(paths: RunPaths) -> None:
    for value in paths.__dict__.values():
        if isinstance(value, Path):
            value.mkdir(parents=True, exist_ok=True)


def expand_stage(stage: str) -> list[str]:
    if stage == "all":
        return ["inventory", "download", "build", "train", "report"]
    if stage == "build":
        return ["inventory", "build"]
    return [stage]


def write_source_inventory(config: dict[str, Any], paths: RunPaths) -> None:
    inventory_path = paths.step1 / "raw_s3_inventory.txt"
    if config.get("s3_raw_prefix") and not inventory_path.exists():
        cmd = ["aws", "s3", "ls", ensure_s3_slash(config["s3_raw_prefix"]), "--recursive", "--summarize"]
        with inventory_path.open("w", encoding="utf-8") as fh:
            subprocess.run(cmd, check=True, stdout=fh)
    manifest = []
    for key, rel, s3_uri in source_file_specs(config):
        local = paths.raw_cache / rel
        manifest.append(
            {
                "source_key": key,
                "s3_uri": s3_uri,
                "local_path": str(local),
                "exists_local": local.exists(),
                "local_size_bytes": local.stat().st_size if local.exists() else None,
            }
        )
    write_json(paths.step1 / "raw_source_manifest.json", {"sources": manifest})


def download_sources(config: dict[str, Any], paths: RunPaths) -> None:
    for _, rel, s3_uri in source_file_specs(config):
        local = paths.raw_cache / rel
        if local.exists():
            continue
        local.parent.mkdir(parents=True, exist_ok=True)
        run(["aws", "s3", "cp", s3_uri, str(local), "--only-show-errors"])


def build_all_inputs(config: dict[str, Any], paths: RunPaths) -> None:
    labels, cells, drugs = build_intermediate_tables(config, paths)
    sample_crispr = build_sample_crispr(config, paths, cells)
    drug_lincs = build_drug_lincs(config, paths, drugs, cells)
    drug_features = build_drug_features(config, paths, drugs, drug_lincs)
    train = build_train_table(config, paths, labels, cells, sample_crispr, drug_features)
    train.to_parquet(paths.step3 / "train_table_full.parquet", index=False)
    write_input_qc(config, paths, labels, cells, drugs, sample_crispr, drug_features, train)
    build_variants(config, paths, train)


def build_intermediate_tables(
    config: dict[str, Any],
    paths: RunPaths,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    labels_raw = read_source(config, paths, "gdsc_labels")
    depmap = read_source(config, paths, "depmap_model")

    if config.get("source_tier") in {"original_core", "raw_only"}:
        gdsc = labels_raw.copy()
        gdsc_drugs = gdsc[["DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME"]].drop_duplicates("DRUG_ID").copy()
    else:
        smiles = read_source(config, paths, "drug_smiles_catalog")
        gdsc_cells = read_source(config, paths, "gdsc_cell_annotations")
        gdsc_drugs = read_source(config, paths, "gdsc_drug_annotations")
        gdsc = labels_raw.merge(gdsc_cells, on="SANGER_MODEL_ID", how="left")
        gdsc = gdsc.merge(gdsc_drugs, on="DRUG_ID", how="left")
    gdsc = gdsc[gdsc["TCGA_DESC"].isin(config["tcga_codes"])].copy()
    gdsc = gdsc.rename(columns=canonical_gdsc_column_map(gdsc.columns))
    gdsc["cell_line_name"] = gdsc["cell_line_name"].fillna(gdsc["SANGER_MODEL_ID"]).astype(str)
    gdsc["drug_name"] = gdsc["drug_name"].fillna("").astype(str)

    cells = build_cell_metadata(config, gdsc, depmap)
    labels = build_labels(config, gdsc, cells)
    if config.get("source_tier") == "raw_only":
        smiles = build_raw_smiles_catalog(config, paths, labels)
    elif config.get("source_tier") == "original_core":
        smiles = read_source(config, paths, "drug_smiles_catalog")
    drugs = build_drug_master(labels, gdsc_drugs, smiles)
    target_mapping = build_drug_target_mapping(drugs)

    labels.to_parquet(paths.step2 / "labels.parquet", index=False)
    cells.to_parquet(paths.step2 / "cell_metadata.parquet", index=False)
    drugs.to_parquet(paths.step2 / "drug_master.parquet", index=False)
    target_mapping.to_parquet(paths.step2 / "drug_target_mapping.parquet", index=False)
    return labels, cells, drugs


def build_cell_metadata(config: dict[str, Any], gdsc: pd.DataFrame, depmap: pd.DataFrame) -> pd.DataFrame:
    base = gdsc[["SANGER_MODEL_ID", "cell_line_name", "TCGA_DESC"]].drop_duplicates().copy()
    dep = depmap.copy()
    dep["SangerModelID"] = dep["SangerModelID"].fillna("").astype(str)
    keep = [
        "ModelID",
        "CellLineName",
        "StrippedCellLineName",
        "OncotreeCode",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "SangerModelID",
        "COSMICID",
    ]
    for column in keep:
        if column not in dep.columns:
            dep[column] = ""
    cells = base.merge(dep[keep], left_on="SANGER_MODEL_ID", right_on="SangerModelID", how="left")
    missing = cells["ModelID"].isna() | cells["ModelID"].astype(str).eq("")
    if missing.any():
        name_map = dep.dropna(subset=["CellLineName"]).copy()
        name_map["cell_line_norm"] = name_map["CellLineName"].map(norm_key)
        fallback = cells.loc[missing, ["cell_line_name"]].copy()
        fallback["cell_line_norm"] = fallback["cell_line_name"].map(norm_key)
        fallback = fallback.merge(name_map[keep + ["cell_line_norm"]], on="cell_line_norm", how="left", suffixes=("", "_fb"))
        for col in keep:
            cells.loc[missing, col] = fallback[col].to_numpy()
    cells["sample_id"] = cells["cell_line_name"].astype(str)
    cells["model_id"] = cells["ModelID"].fillna("").astype(str)
    cells["is_depmap_mapped"] = cells["model_id"].ne("").astype(int)
    cells["cohort"] = str(config["cancer_id"])
    cells = cells.rename(
        columns={
            "OncotreeCode": "depmap_oncotree_code",
            "OncotreePrimaryDisease": "depmap_primary_disease",
            "OncotreeSubtype": "depmap_oncotree_subtype",
        }
    )
    keep_out = [
        "sample_id",
        "cell_line_name",
        "SANGER_MODEL_ID",
        "model_id",
        "is_depmap_mapped",
        "TCGA_DESC",
        "cohort",
        "depmap_oncotree_code",
        "depmap_primary_disease",
        "depmap_oncotree_subtype",
    ]
    for column in keep_out:
        if column not in cells.columns:
            cells[column] = ""
    return cells[keep_out].sort_values("sample_id").reset_index(drop=True)


def build_labels(config: dict[str, Any], gdsc: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    labels = gdsc.copy()
    label_col = config.get("label", {}).get("regression_column", "LN_IC50")
    labels["sample_id"] = labels["cell_line_name"].astype(str)
    labels["canonical_drug_id"] = labels["DRUG_ID"].astype(int).astype(str)
    labels["pair_id"] = labels["sample_id"] + "__" + labels["canonical_drug_id"]
    labels["label_regression"] = pd.to_numeric(labels[label_col], errors="coerce")
    labels = labels.dropna(subset=["label_regression"]).copy()
    threshold = float(labels["label_regression"].quantile(float(config.get("label", {}).get("binary_quantile", 0.3))))
    labels["label_binary"] = (labels["label_regression"] <= threshold).astype(int)
    labels["label_main"] = labels["label_regression"]
    labels["label_aux"] = labels["label_binary"]
    labels["label_main_type"] = "regression"
    labels["label_aux_type"] = "binary"
    labels["binary_threshold"] = threshold
    labels["gdsc_version"] = labels.get("DATASET", "GDSC2")
    labels["cohort"] = str(config["cancer_id"])
    labels = labels.merge(cells[["sample_id", "model_id", "is_depmap_mapped"]], on="sample_id", how="left")
    keep = [
        "pair_id",
        "sample_id",
        "cell_line_name",
        "SANGER_MODEL_ID",
        "model_id",
        "is_depmap_mapped",
        "canonical_drug_id",
        "DRUG_ID",
        "drug_name",
        "TCGA_DESC",
        "cohort",
        "gdsc_version",
        "putative_target",
        "pathway_name",
        "LN_IC50",
        "AUC",
        "Z_SCORE",
        "label_regression",
        "label_binary",
        "label_main",
        "label_aux",
        "label_main_type",
        "label_aux_type",
        "binary_threshold",
    ]
    for column in keep:
        if column not in labels.columns:
            labels[column] = ""
    return labels[keep].reset_index(drop=True)


def build_drug_master(labels: pd.DataFrame, gdsc_drugs: pd.DataFrame, smiles: pd.DataFrame) -> pd.DataFrame:
    gdsc = gdsc_drugs.rename(
        columns={
            "DRUG_NAME": "drug_name",
            "PUTATIVE_TARGET_NORMALIZED": "putative_target",
            "PUTATIVE_TARGET": "putative_target",
            "PATHWAY_NAME_NORMALIZED": "pathway_name",
            "PATHWAY_NAME": "pathway_name",
        }
    ).copy()
    label_drugs = labels[["DRUG_ID", "drug_name", "putative_target", "pathway_name"]].drop_duplicates("DRUG_ID")
    master = label_drugs.merge(gdsc[["DRUG_ID", "drug_name", "putative_target", "pathway_name"]], on="DRUG_ID", how="left", suffixes=("", "_gdsc"))
    for column in ["drug_name", "putative_target", "pathway_name"]:
        master[column] = master[column].fillna(master[f"{column}_gdsc"]).fillna("").astype(str)
    smi = smiles.rename(columns={"DRUG_NAME": "drug_name_catalog"}).copy()
    smi["DRUG_ID"] = pd.to_numeric(smi["DRUG_ID"], errors="coerce").astype("Int64")
    smi = smi.dropna(subset=["DRUG_ID"]).copy()
    smi["DRUG_ID"] = smi["DRUG_ID"].astype(int)
    master = master.merge(smi[["DRUG_ID", "canonical_smiles", "match_source", "has_smiles"]], on="DRUG_ID", how="left")
    master["canonical_drug_id"] = master["DRUG_ID"].astype(int).astype(str)
    master["canonical_smiles"] = master["canonical_smiles"].fillna("").astype(str)
    master["has_smiles"] = pd.to_numeric(master["has_smiles"], errors="coerce").fillna(0).astype(int)
    master["has_smiles"] = master["canonical_smiles"].str.strip().ne("").astype(int)
    master["match_source"] = master["match_source"].fillna("unmatched").astype(str)
    master["target"] = master["putative_target"].fillna("").astype(str)
    master["target_pathway"] = master["pathway_name"].fillna("").astype(str)
    master["drug_name_norm"] = master["drug_name"].map(norm_key)
    master["drug__target_list"] = master["target"].map(lambda value: ",".join(extract_targets(value)))
    keep = [
        "canonical_drug_id",
        "DRUG_ID",
        "drug_name",
        "drug_name_norm",
        "canonical_smiles",
        "has_smiles",
        "match_source",
        "target",
        "target_pathway",
        "drug__target_list",
    ]
    return master[keep].sort_values("canonical_drug_id").reset_index(drop=True)


def build_drug_target_mapping(drugs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in drugs.iterrows():
        for target in extract_targets(row.get("target", "")):
            rows.append(
                {
                    "canonical_drug_id": str(row["canonical_drug_id"]),
                    "drug_name": row["drug_name"],
                    "target_gene_symbol": target,
                    "source": "gdsc_putative_target",
                }
            )
    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def build_raw_smiles_catalog(config: dict[str, Any], paths: RunPaths, labels: pd.DataFrame) -> pd.DataFrame:
    drugs = labels[["DRUG_ID", "drug_name"]].drop_duplicates("DRUG_ID").copy()
    candidates: dict[int, list[tuple[str, str]]] = {}
    for _, row in drugs.iterrows():
        drug_id = int(row["DRUG_ID"])
        candidates.setdefault(drug_id, [])
        candidates[drug_id].append((str(row["drug_name"]), "gdsc_drug_name"))

    if source_key_exists(config, "gdsc_screened_compounds"):
        screened = read_source(config, paths, "gdsc_screened_compounds")
        for _, row in screened.iterrows():
            drug_id_raw = pd.to_numeric(row.get("DRUG_ID"), errors="coerce")
            if pd.isna(drug_id_raw):
                continue
            drug_id = int(drug_id_raw)
            if drug_id not in candidates:
                continue
            candidates[drug_id].append((str(row.get("DRUG_NAME", "")), "gdsc_screened_name"))
            for synonym in split_synonyms(row.get("SYNONYMS", "")):
                candidates[drug_id].append((synonym, "gdsc_screened_synonym"))

    wanted_names = {norm_key(name) for names in candidates.values() for name, _ in names if norm_key(name)}
    lincs_index = build_lincs_smiles_name_index(config, paths, wanted_names)
    drugbank_index = build_drugbank_smiles_name_index(config, paths, wanted_names)
    source_indexes = [
        ("lincs_pert_info", lincs_index),
        ("drugbank_xml", drugbank_index),
    ]

    rows = []
    for _, row in drugs.iterrows():
        drug_id = int(row["DRUG_ID"])
        canonical_smiles = ""
        match_source = "unmatched"
        matched_name = ""
        for candidate_name, candidate_source in candidates.get(drug_id, []):
            key = norm_key(candidate_name)
            if not key:
                continue
            for source_name, index in source_indexes:
                smiles = index.get(key, "")
                if smiles:
                    canonical_smiles = smiles
                    match_source = f"{source_name}:{candidate_source}"
                    matched_name = candidate_name
                    break
            if canonical_smiles:
                break
        rows.append(
            {
                "DRUG_ID": drug_id,
                "DRUG_NAME": str(row["drug_name"]),
                "canonical_smiles": canonical_smiles,
                "has_smiles": int(bool(canonical_smiles)),
                "match_source": match_source,
                "matched_name": matched_name,
            }
        )
    out = pd.DataFrame(rows).sort_values("DRUG_ID").reset_index(drop=True)
    out.to_parquet(paths.step2 / "raw_smiles_catalog.parquet", index=False)
    return out


def build_lincs_smiles_name_index(config: dict[str, Any], paths: RunPaths, wanted_names: set[str]) -> dict[str, str]:
    if not source_key_exists(config, "lincs_pert_info"):
        return {}
    pert = read_source(config, paths, "lincs_pert_info")
    if "pert_type" in pert.columns:
        pert = pert[pert["pert_type"].astype(str).eq("trt_cp")].copy()
    smiles_col = first_existing_column(pert, ["canonical_smiles", "smiles"])
    if smiles_col is None:
        return {}
    name_cols = [col for col in ["pert_iname", "cmap_name", "pert_name"] if col in pert.columns]
    alias_cols = [col for col in ["compound_aliases", "pert_alias", "alias"] if col in pert.columns]
    out: dict[str, str] = {}
    for _, row in pert.iterrows():
        smiles = canonicalize_smiles(row.get(smiles_col, ""))
        if not smiles:
            continue
        for name in iter_lincs_names(row, name_cols, alias_cols):
            key = norm_key(name)
            if key in wanted_names and key not in out:
                out[key] = smiles
    return out


def build_drugbank_smiles_name_index(config: dict[str, Any], paths: RunPaths, wanted_names: set[str]) -> dict[str, str]:
    if not source_key_exists(config, "drugbank_xml_zip"):
        return {}
    zip_path = source_local_path(config, paths, "drugbank_xml_zip")
    if not zip_path.exists():
        raise FileNotFoundError(f"Missing DrugBank XML zip: {zip_path}")
    out: dict[str, str] = {}
    with zipfile.ZipFile(zip_path) as zf:
        xml_names = [name for name in zf.namelist() if name.lower().endswith(".xml")]
        if not xml_names:
            return {}
        with zf.open(xml_names[0]) as fh:
            for _, elem in ET.iterparse(fh, events=("end",)):
                if local_xml_name(elem.tag) != "drug":
                    continue
                names = extract_drugbank_names(elem)
                keys = [norm_key(name) for name in names]
                matched_keys = [key for key in keys if key in wanted_names]
                if matched_keys:
                    smiles = canonicalize_smiles(extract_drugbank_smiles(elem))
                    if smiles:
                        for key in matched_keys:
                            out.setdefault(key, smiles)
                elem.clear()
    return out


def extract_drugbank_names(elem: ET.Element) -> list[str]:
    names: list[str] = []
    primary = direct_child_text(elem, "name")
    if primary:
        names.append(primary)
    for child in elem:
        if local_xml_name(child.tag) != "synonyms":
            continue
        for synonym in child:
            if local_xml_name(synonym.tag) == "synonym" and synonym.text:
                names.append(synonym.text)
    return names


def extract_drugbank_smiles(elem: ET.Element) -> str:
    for container_name in ["calculated-properties", "experimental-properties"]:
        for container in elem:
            if local_xml_name(container.tag) != container_name:
                continue
            for prop in container:
                if local_xml_name(prop.tag) != "property":
                    continue
                kind = ""
                value = ""
                for child in prop:
                    child_name = local_xml_name(child.tag)
                    if child_name == "kind":
                        kind = str(child.text or "")
                    elif child_name == "value":
                        value = str(child.text or "")
                if kind.strip().lower() == "smiles" and value.strip():
                    return value
    return ""


def direct_child_text(elem: ET.Element, name: str) -> str:
    for child in elem:
        if local_xml_name(child.tag) == name and child.text:
            return str(child.text)
    return ""


def local_xml_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def build_sample_crispr(config: dict[str, Any], paths: RunPaths, cells: pd.DataFrame) -> pd.DataFrame:
    crispr = read_source(config, paths, "depmap_crispr_gene_effect")
    if "ModelID" not in crispr.columns and "Unnamed: 0" in crispr.columns:
        crispr = crispr.rename(columns={"Unnamed: 0": "ModelID"})
    crispr["ModelID"] = crispr["ModelID"].astype(str)
    model_ids = cells["model_id"].astype(str).tolist()
    sub = crispr[crispr["ModelID"].isin(set(model_ids))].copy()
    gene_cols = [c for c in sub.columns if c != "ModelID"]
    selected_raw = select_top_variance_columns(sub, gene_cols, int(config["feature_limits"].get("full_crispr", 3500)))
    rename = {col: f"sample__crispr__{parse_gene_symbol(col)}" for col in selected_raw}
    sub = sub[["ModelID"] + selected_raw].rename(columns=rename)
    sub = sub.rename(columns={"ModelID": "model_id"})
    out = cells[["sample_id", "cell_line_name", "model_id"]].merge(sub, on="model_id", how="left")
    feature_cols = [c for c in out.columns if c.startswith("sample__crispr__")]
    out["sample__has_crispr_profile"] = out[feature_cols].notna().any(axis=1).astype(int)
    out[feature_cols] = out[feature_cols].fillna(0.0).astype(np.float32)
    out.to_parquet(paths.step2 / "sample_crispr.parquet", index=False)
    return out


def build_drug_lincs(
    config: dict[str, Any],
    paths: RunPaths,
    drugs: pd.DataFrame | None = None,
    cells: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if config.get("source_tier") == "raw_only":
        if drugs is None:
            raise ValueError("raw_only LINCS construction requires the drug master table")
        return build_raw_lincs_features(config, paths, drugs, cells)
    raw = read_source(config, paths, "lincs_pancancer").copy()
    raw["canonical_drug_id"] = raw["canonical_drug_id"].astype(str)
    lincs_cols = [c for c in raw.columns if c != "canonical_drug_id" and pd.api.types.is_numeric_dtype(raw[c])]
    selected = select_top_variance_columns(raw, lincs_cols, int(config["feature_limits"].get("legacy_lincs", 768)))
    out = raw[["canonical_drug_id"] + selected].copy()
    rename = {c: "drug__lincs__" + clean_feature_token(c.replace("drug__lincs__", "")) for c in selected}
    out = out.rename(columns=rename)
    feature_cols = [c for c in out.columns if c.startswith("drug__lincs__")]
    out["drug__has_lincs_signature"] = 1
    out = add_lincs_summary(out, feature_cols)
    out.to_parquet(paths.step2 / "drug_lincs.parquet", index=False)
    return out


def build_raw_lincs_features(
    config: dict[str, Any],
    paths: RunPaths,
    drugs: pd.DataFrame,
    cells: pd.DataFrame | None = None,
) -> pd.DataFrame:
    lincs_cfg = config.get("lincs", {})
    if lincs_cfg.get("mode", "raw") == "none":
        write_json(
            paths.step2 / "raw_lincs_build_summary.json",
            {
                "mode": "none",
                "matched_drugs": 0,
                "matched_perturbagens": 0,
                "matched_signatures": 0,
                "selected_gene_features": 0,
                "cell_filter": "none",
                "cell_column": "",
                "allowed_cell_ids": [],
                "source_total_cell_line_count_after_drug_filter": 0,
                "included_total_cell_line_count": 0,
                "available_cancer_overlap_cell_ids": [],
                "available_cancer_overlap_cell_line_count": 0,
                "available_representative_cell_ids": [],
                "available_representative_cell_line_count": 0,
            },
        )
        return empty_lincs_features(paths, drugs)

    if not all(source_key_exists(config, key) for key in ["lincs_pert_info", "lincs_sig_info", "lincs_gene_info", "lincs_gctx"]):
        return empty_lincs_features(paths, drugs)

    pert = read_source(config, paths, "lincs_pert_info")
    sig = read_source(config, paths, "lincs_sig_info")
    gene = read_source(config, paths, "lincs_gene_info")
    gctx_path = source_local_path(config, paths, "lincs_gctx")
    if not gctx_path.exists():
        raise FileNotFoundError(f"Missing LINCS GCTX source: {gctx_path}")

    drug_to_pert = map_drugs_to_lincs_perturbagens(drugs, pert)
    matched_pert_ids = sorted({pert_id for ids in drug_to_pert.values() for pert_id in ids})
    if not matched_pert_ids:
        return empty_lincs_features(paths, drugs)

    sig = sig[sig["pert_id"].astype(str).isin(set(matched_pert_ids))].copy()
    if "pert_type" in sig.columns:
        sig = sig[sig["pert_type"].astype(str).eq("trt_cp")].copy()
    sig, cell_filter_info = filter_lincs_signatures_by_cell(config, sig, cells)
    if sig.empty:
        write_json(
            paths.step2 / "raw_lincs_build_summary.json",
            {
                "mode": "raw",
                "matched_drugs": 0,
                "matched_perturbagens": int(len(matched_pert_ids)),
                "matched_signatures": 0,
                "selected_gene_features": 0,
                "source": str(gctx_path),
                **cell_filter_info,
            },
        )
        return empty_lincs_features(paths, drugs)

    matrix, row_ids, col_ids = read_gctx_for_signatures(config, paths, gctx_path, sig["sig_id"].astype(str).tolist(), gene)
    if matrix.size == 0:
        return empty_lincs_features(paths, drugs)

    sig_to_pos = {sig_id: idx for idx, sig_id in enumerate(col_ids)}
    pert_to_positions: dict[str, list[int]] = {}
    for _, row in sig.iterrows():
        sig_id = str(row["sig_id"])
        if sig_id in sig_to_pos:
            pert_to_positions.setdefault(str(row["pert_id"]), []).append(sig_to_pos[sig_id])

    rows = []
    feature_names = ["drug__lincs__" + clean_feature_token(symbol) for symbol in row_ids]
    for _, drug_row in drugs.iterrows():
        drug_id = str(drug_row["canonical_drug_id"])
        positions = sorted({pos for pert_id in drug_to_pert.get(drug_id, []) for pos in pert_to_positions.get(pert_id, [])})
        if positions:
            values = matrix[:, positions].mean(axis=1).astype(np.float32)
            has_signature = 1
        else:
            values = np.zeros((len(feature_names),), dtype=np.float32)
            has_signature = 0
        row = {"canonical_drug_id": drug_id, "drug__has_lincs_signature": has_signature}
        row.update({name: float(value) for name, value in zip(feature_names, values)})
        rows.append(row)

    out = pd.DataFrame(rows)
    lincs_cols = [c for c in out.columns if c.startswith("drug__lincs__")]
    keep = select_top_variance_columns(out, lincs_cols, int(config["feature_limits"].get("legacy_lincs", 768)))
    out = out[["canonical_drug_id", "drug__has_lincs_signature"] + keep].copy()
    out = add_lincs_summary(out, keep)
    out.to_parquet(paths.step2 / "drug_lincs.parquet", index=False)
    write_json(
        paths.step2 / "raw_lincs_build_summary.json",
        {
            "mode": "raw",
            "matched_drugs": int(out["drug__has_lincs_signature"].sum()),
            "matched_perturbagens": int(len(matched_pert_ids)),
            "matched_signatures": int(len(col_ids)),
            "selected_gene_features": int(len(keep)),
            "source": str(gctx_path),
            **cell_filter_info,
        },
    )
    return out


def filter_lincs_signatures_by_cell(
    config: dict[str, Any],
    sig: pd.DataFrame,
    cells: pd.DataFrame | None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    lincs_cfg = config.get("lincs", {})
    policy_cfg = config.get("lincs_policy", {})
    policy = str(lincs_cfg.get("cell_filter", "all"))
    cell_col = first_existing_column(sig, ["cell_id", "cell_iname"])
    representative_norm = {norm_key(value) for value in policy_cfg.get("representative_cell_lines", [])}
    available_cancer_overlap: list[str] = []
    available_representative: list[str] = []
    source_total_cell_lines = int(sig[cell_col].astype(str).nunique()) if cell_col else 0
    if cells is not None and cell_col is not None:
        disease_names = {norm_key(value): str(value) for value in cells["cell_line_name"].astype(str)}
        lincs_names = {norm_key(value): str(value) for value in sig[cell_col].astype(str).unique()}
        available_cancer_overlap = sorted(lincs_names[key] for key in disease_names.keys() & lincs_names.keys())
        available_representative = sorted(
            value for value in available_cancer_overlap if norm_key(value) in representative_norm
        )
    if policy in {"all", "none", ""}:
        return sig, {
            "cell_filter": "all",
            "cell_column": cell_col or "",
            "allowed_cell_ids": [],
            "source_total_cell_line_count_after_drug_filter": source_total_cell_lines,
            "included_total_cell_line_count": source_total_cell_lines,
            "available_cancer_overlap_cell_ids": available_cancer_overlap,
            "available_cancer_overlap_cell_line_count": int(len(available_cancer_overlap)),
            "available_representative_cell_ids": available_representative,
            "available_representative_cell_line_count": int(len(available_representative)),
            "signatures_after_cell_filter": int(len(sig)),
            "perturbagens_after_cell_filter": int(sig["pert_id"].astype(str).nunique()) if "pert_id" in sig.columns else 0,
        }
    if cell_col is None:
        raise ValueError("LINCS cell filtering requested, but sig_info has no cell_id/cell_iname column")
    if policy == "explicit":
        allowed = [str(value) for value in lincs_cfg.get("allowed_cell_ids", [])]
    elif policy == "cancer_overlap":
        if cells is None:
            raise ValueError("cancer_overlap LINCS cell filtering requires cell metadata")
        paad_names = {norm_key(value): str(value) for value in cells["cell_line_name"].astype(str)}
        lincs_names = {norm_key(value): str(value) for value in sig[cell_col].astype(str).unique()}
        allowed = sorted(lincs_names[key] for key in paad_names.keys() & lincs_names.keys())
    else:
        raise ValueError(f"Unknown LINCS cell filter policy: {policy}")
    allowed_set = set(allowed)
    filtered = sig[sig[cell_col].astype(str).isin(allowed_set)].copy()
    return filtered, {
        "cell_filter": policy,
        "cell_column": cell_col,
        "allowed_cell_ids": allowed,
        "source_total_cell_line_count_after_drug_filter": source_total_cell_lines,
        "included_total_cell_line_count": int(filtered[cell_col].astype(str).nunique()),
        "available_cancer_overlap_cell_ids": available_cancer_overlap,
        "available_cancer_overlap_cell_line_count": int(len(available_cancer_overlap)),
        "available_representative_cell_ids": available_representative,
        "available_representative_cell_line_count": int(len(available_representative)),
        "signatures_after_cell_filter": int(len(filtered)),
        "perturbagens_after_cell_filter": int(filtered["pert_id"].astype(str).nunique()) if "pert_id" in filtered.columns else 0,
    }


def empty_lincs_features(paths: RunPaths, drugs: pd.DataFrame) -> pd.DataFrame:
    out = drugs[["canonical_drug_id"]].copy()
    out["drug__has_lincs_signature"] = 0
    out = add_lincs_summary(out, [])
    out.to_parquet(paths.step2 / "drug_lincs.parquet", index=False)
    return out


def map_drugs_to_lincs_perturbagens(drugs: pd.DataFrame, pert: pd.DataFrame) -> dict[str, list[str]]:
    pert = pert.copy()
    if "pert_type" in pert.columns:
        pert = pert[pert["pert_type"].astype(str).eq("trt_cp")].copy()
    smiles_col = first_existing_column(pert, ["canonical_smiles", "smiles"])
    name_cols = [col for col in ["pert_iname", "cmap_name", "pert_name"] if col in pert.columns]
    alias_cols = [col for col in ["compound_aliases", "pert_alias", "alias"] if col in pert.columns]
    pert["smiles_canonical_raw"] = pert[smiles_col].map(canonicalize_smiles) if smiles_col else ""
    by_smiles: dict[str, set[str]] = {}
    by_name: dict[str, set[str]] = {}
    for _, row in pert.iterrows():
        pert_id = str(row.get("pert_id", ""))
        smiles = str(row.get("smiles_canonical_raw", ""))
        if smiles:
            by_smiles.setdefault(smiles, set()).add(pert_id)
        for name in iter_lincs_names(row, name_cols, alias_cols):
            name_key = norm_key(name)
            if name_key:
                by_name.setdefault(name_key, set()).add(pert_id)

    out: dict[str, list[str]] = {}
    for _, row in drugs.iterrows():
        drug_id = str(row["canonical_drug_id"])
        smiles = canonicalize_smiles(row.get("canonical_smiles", ""))
        name_key = norm_key(row.get("drug_name", ""))
        matches = set()
        if smiles:
            matches.update(by_smiles.get(smiles, set()))
        if name_key:
            matches.update(by_name.get(name_key, set()))
        out[drug_id] = sorted(matches)
    return out


def read_gctx_for_signatures(
    config: dict[str, Any],
    paths: RunPaths,
    gctx_path: Path,
    sig_ids: list[str],
    gene: pd.DataFrame,
) -> tuple[np.ndarray, list[str], list[str]]:
    h5_path = ensure_hdf5_gctx(config, paths, gctx_path)
    with h5py.File(h5_path, "r") as h5:
        row_ids = decode_h5_values(h5["0/META/ROW/id"][()])
        col_ids = decode_h5_values(h5["0/META/COL/id"][()])
        matrix_ds = h5["0/DATA/0/matrix"]
        gene_ids = select_lincs_gene_ids(config, gene, row_ids)
        row_lookup = {value: idx for idx, value in enumerate(row_ids)}
        col_lookup = {value: idx for idx, value in enumerate(col_ids)}
        row_idx = sorted(row_lookup[gene_id] for gene_id in gene_ids if gene_id in row_lookup)
        col_idx = sorted(col_lookup[sig_id] for sig_id in set(sig_ids) if sig_id in col_lookup)
        if not row_idx or not col_idx:
            return np.zeros((0, 0), dtype=np.float32), [], []
        if matrix_ds.shape == (len(row_ids), len(col_ids)):
            sub = np.asarray(matrix_ds[row_idx, :], dtype=np.float32)[:, col_idx]
        elif matrix_ds.shape == (len(col_ids), len(row_ids)):
            sub = np.asarray(matrix_ds[col_idx, :], dtype=np.float32)[:, row_idx].T
        else:
            raise ValueError(f"Unexpected GCTX matrix shape: {matrix_ds.shape}")
        selected_row_ids = [row_ids[i] for i in row_idx]
        selected_col_ids = [col_ids[i] for i in col_idx]
    gene_symbol_map = build_lincs_gene_symbol_map(gene)
    row_symbols = [gene_symbol_map.get(row_id, row_id) for row_id in selected_row_ids]
    return np.nan_to_num(sub, nan=0.0, posinf=0.0, neginf=0.0), row_symbols, selected_col_ids


def ensure_hdf5_gctx(config: dict[str, Any], paths: RunPaths, path: Path) -> Path:
    if h5py.is_hdf5(str(path)):
        return path
    if path.suffix != ".gz":
        raise ValueError(f"LINCS source is not HDF5 and is not gzip-compressed: {path}")
    out = paths.raw_cache / "raw_only_work" / path.name[:-3]
    if out.exists() and h5py.is_hdf5(str(out)):
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "rb") as src, out.open("wb") as dst:
        while True:
            chunk = src.read(16 * 1024 * 1024)
            if not chunk:
                break
            dst.write(chunk)
    if not h5py.is_hdf5(str(out)):
        raise ValueError(f"Decompressed LINCS file is not HDF5: {out}")
    return out


def decode_h5_values(values: Any) -> list[str]:
    out: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8"))
        else:
            out.append(str(value))
    return out


def select_lincs_gene_ids(config: dict[str, Any], gene: pd.DataFrame, row_ids: list[str]) -> list[str]:
    gene = gene.copy()
    gene_id_col = first_existing_column(gene, ["pr_gene_id", "gene_id"])
    if gene_id_col is None:
        return []
    gene[gene_id_col] = gene[gene_id_col].astype(str)
    if "pr_is_lm" in gene.columns:
        gene = gene[pd.to_numeric(gene["pr_is_lm"], errors="coerce").fillna(0).astype(int).eq(1)].copy()
    elif "feature_space" in gene.columns:
        gene = gene[gene["feature_space"].astype(str).str.lower().eq("landmark")].copy()
    row_id_set = set(row_ids)
    gene = gene[gene[gene_id_col].isin(row_id_set)].copy()
    limit = int(config["feature_limits"].get("raw_lincs_genes", 978))
    return gene[gene_id_col].astype(str).head(limit).tolist()


def build_lincs_gene_symbol_map(gene: pd.DataFrame) -> dict[str, str]:
    out: dict[str, str] = {}
    gene_id_col = first_existing_column(gene, ["pr_gene_id", "gene_id"])
    symbol_col = first_existing_column(gene, ["pr_gene_symbol", "gene_symbol"])
    if gene_id_col is None or symbol_col is None:
        return out
    for _, row in gene.iterrows():
        out[str(row[gene_id_col])] = str(row[symbol_col])
    return out


def build_drug_features(
    config: dict[str, Any],
    paths: RunPaths,
    drugs: pd.DataFrame,
    drug_lincs: pd.DataFrame,
) -> pd.DataFrame:
    rdkit_df = build_rdkit_features(drugs, n_bits=2048)
    smiles_svd = build_smiles_svd_features(drugs, int(config["feature_limits"].get("smiles_svd", 64)))
    out = drugs.merge(rdkit_df, on="canonical_drug_id", how="left")
    out = out.merge(smiles_svd, on="canonical_drug_id", how="left")
    out = out.merge(drug_lincs, on="canonical_drug_id", how="left")
    lincs_cols = [c for c in out.columns if c.startswith("drug__lincs__") or c.startswith("drug__lincs_summary_")]
    out[lincs_cols] = out[lincs_cols].fillna(0.0)
    out["drug__has_lincs_signature"] = pd.to_numeric(out.get("drug__has_lincs_signature", 0), errors="coerce").fillna(0).astype(int)
    out["drug__target_count"] = out["drug__target_list"].map(lambda s: len([x for x in str(s).split(",") if x]))
    out["drug__has_target_mapping"] = out["drug__target_count"].gt(0).astype(int)
    out.to_parquet(paths.step2 / "drug_features.parquet", index=False)
    return out


def build_train_table(
    config: dict[str, Any],
    paths: RunPaths,
    labels: pd.DataFrame,
    cells: pd.DataFrame,
    sample_crispr: pd.DataFrame,
    drug_features: pd.DataFrame,
) -> pd.DataFrame:
    train = labels.merge(cells.drop(columns=["cell_line_name", "TCGA_DESC"], errors="ignore"), on=["sample_id", "SANGER_MODEL_ID", "model_id", "is_depmap_mapped", "cohort"], how="left")
    train = train.merge(sample_crispr.drop(columns=["cell_line_name", "model_id"], errors="ignore"), on="sample_id", how="left")
    train = train.merge(drug_features, on=["canonical_drug_id", "DRUG_ID"], how="left", suffixes=("", "_drug"))
    train["drug_name"] = train["drug_name"].fillna(train.get("drug_name_drug", "")).astype(str)
    train["drug__canonical_smiles"] = train["canonical_smiles"].fillna("").astype(str)
    train["drug__has_smiles"] = train["drug__canonical_smiles"].str.strip().ne("").astype(int)
    train = add_context_features(config, train)
    train = add_pair_target_features(train)
    numeric_cols = train.select_dtypes(include=[np.number]).columns
    train[numeric_cols] = train[numeric_cols].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return train.reset_index(drop=True)


def add_context_features(config: dict[str, Any], train: pd.DataFrame) -> pd.DataFrame:
    out = train.copy()
    out[f"sample__is_{str(config['cancer_id']).lower()}"] = 1
    for value in sorted(out["cell_line_name"].fillna("").astype(str).unique()):
        out[f"sample__cell_line__{clean_feature_token(value)}"] = (out["cell_line_name"].astype(str) == value).astype(int)
    for col, prefix in [
        ("pathway_name", "ctxcat__pathway_name_normalized__"),
        ("target_pathway", "ctxcat__target_pathway__"),
        ("depmap_oncotree_code", "ctxcat__depmap_oncotree_code__"),
    ]:
        if col not in out.columns:
            continue
        values = sorted(v for v in out[col].fillna("").astype(str).unique() if v)
        for value in values:
            out[prefix + clean_feature_token(value)] = (out[col].astype(str) == value).astype(int)
    return out


def add_pair_target_features(train: pd.DataFrame) -> pd.DataFrame:
    out = train.copy()
    crispr_cols = {c.replace("sample__crispr__", "").upper(): c for c in out.columns if c.startswith("sample__crispr__")}
    rows = []
    for _, row in out[["drug__target_list"]].iterrows():
        targets = [x for x in str(row["drug__target_list"]).split(",") if x]
        matched = [crispr_cols[t] for t in targets if t in crispr_cols]
        rows.append((targets, matched))
    counts = np.array([len(t[0]) for t in rows], dtype=np.float32)
    matched_counts = np.array([len(t[1]) for t in rows], dtype=np.float32)
    out["pair__target_count"] = counts
    out["pair__matched_target_count"] = matched_counts
    out["pair__target_match_fraction"] = np.divide(matched_counts, np.maximum(counts, 1.0))
    means = np.zeros(len(out), dtype=np.float32)
    abs_means = np.zeros(len(out), dtype=np.float32)
    max_abs = np.zeros(len(out), dtype=np.float32)
    for idx, (_, matched) in enumerate(rows):
        if not matched:
            continue
        values = pd.to_numeric(out.loc[idx, matched], errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)
        means[idx] = float(values.mean())
        abs_means[idx] = float(np.abs(values).mean())
        max_abs[idx] = float(np.abs(values).max())
    out["pair__mean_target_dependency"] = means
    out["pair__mean_abs_target_dependency"] = abs_means
    out["pair__max_abs_target_dependency"] = max_abs
    out["pair__top50_target_hits"] = (out["pair__mean_abs_target_dependency"] > 0.5).astype(int)
    out["pair__top200_target_hits"] = (out["pair__mean_abs_target_dependency"] > 0.2).astype(int)
    return out


def build_variants(config: dict[str, Any], paths: RunPaths, train: pd.DataFrame) -> None:
    for variant in config["variants"]:
        variant_dir = paths.step4 / variant
        variant_dir.mkdir(parents=True, exist_ok=True)
        if variant == "legacy_rich_valid_smiles_only":
            mask = pd.to_numeric(train["drug_has_valid_smiles"], errors="coerce").fillna(0).gt(0)
            variant_train = train[mask].copy()
            mode = "legacy"
        elif variant == "legacy_rich_all_drugs_zero_smiles":
            variant_train = zero_invalid_smiles_features(train.copy())
            mode = "legacy"
        elif variant == "colonstyle_compact_baseline":
            variant_train = zero_invalid_smiles_features(train.copy())
            mode = "compact"
        else:
            raise ValueError(f"Unknown variant: {variant}")
        feature_names = select_feature_names(variant_train, config, mode)
        save_variant(variant_dir, variant, variant_train, feature_names)


def zero_invalid_smiles_features(train: pd.DataFrame) -> pd.DataFrame:
    valid = pd.to_numeric(train["drug_has_valid_smiles"], errors="coerce").fillna(0).gt(0)
    cols = smiles_derived_columns(train)
    train.loc[~valid, cols] = 0.0
    return train


def smiles_derived_columns(df: pd.DataFrame) -> list[str]:
    cols = [c for c in df.columns if c.startswith(SMILES_DERIVED_PREFIXES)]
    for col in ["drug__smiles_length", "drug__has_smiles"]:
        if col in df.columns:
            cols.append(col)
    return cols


def select_feature_names(train: pd.DataFrame, config: dict[str, Any], mode: str) -> list[str]:
    limits = config["feature_limits"]
    numeric = [
        c
        for c in train.select_dtypes(include=[np.number]).columns
        if c not in KEY_COLUMNS and not c.startswith("label_") and c not in {"binary_threshold", "LN_IC50", "AUC", "Z_SCORE"}
    ]
    numeric = [c for c in numeric if float(train[c].var(ddof=0)) > 1e-12]

    crispr = [c for c in numeric if c.startswith("sample__crispr__")]
    morgan = [c for c in numeric if c.startswith("drug_morgan_")]
    lincs_gene = [c for c in numeric if c.startswith("drug__lincs__")]
    lincs_summary = [c for c in numeric if c.startswith("drug__lincs_summary_") or c in {"drug__has_lincs_signature"}]
    smiles = [c for c in numeric if c.startswith("smiles_svd_")]
    descriptors = [c for c in numeric if c.startswith("drug__descriptor_") or c in {"drug__smiles_length", "drug__has_smiles", "drug_has_valid_smiles"}]
    target = [c for c in numeric if c.startswith("drug__target") or c.startswith("pair__target") or c.startswith("pair__mean") or c.startswith("pair__max") or c.startswith("pair__top")]
    context = [c for c in numeric if c.startswith("sample__cell_line__") or c.startswith("sample__is_") or c.startswith("ctxcat__") or c in {"is_depmap_mapped", "sample__has_crispr_profile"}]

    if mode == "compact":
        crispr_keep = top_by_variance(train, crispr, int(limits.get("compact_crispr", 5000)))
        morgan_keep = top_by_variance(train, morgan, int(limits.get("compact_morgan", 1024)))
        feature_names = context + descriptors + smiles + lincs_summary + crispr_keep + morgan_keep
    else:
        crispr_keep = top_by_variance(train, crispr, int(limits.get("legacy_crispr", 3000)))
        morgan_keep = top_by_variance(train, morgan, int(limits.get("legacy_morgan", 1100)))
        lincs_keep = top_by_variance(train, lincs_gene, int(limits.get("legacy_lincs", 768)))
        other = [
            c
            for c in numeric
            if c not in set(crispr + morgan + lincs_gene + lincs_summary + smiles + descriptors + target + context)
        ]
        feature_names = other + context + descriptors + smiles + target + lincs_summary + lincs_keep + crispr_keep + morgan_keep
    deduped = []
    seen = set()
    for col in feature_names:
        if col not in seen and col in train.columns:
            deduped.append(col)
            seen.add(col)
    return deduped


def save_variant(variant_dir: Path, variant: str, train: pd.DataFrame, feature_names: list[str]) -> None:
    X = train[feature_names].replace([np.inf, -np.inf], np.nan).fillna(0.0).astype(np.float32)
    y = train["label_regression"].to_numpy(dtype=np.float32)
    train_out = train.copy()
    train_out[feature_names] = X
    train_out.to_parquet(variant_dir / "train_table.parquet", index=False)
    np.save(variant_dir / "X.npy", X.to_numpy(dtype=np.float32))
    np.save(variant_dir / "y.npy", y)
    write_json(variant_dir / "feature_names.json", feature_names)
    valid = pd.to_numeric(train.get("drug_has_valid_smiles", 0), errors="coerce").fillna(0).gt(0)
    write_json(
        variant_dir / "input_summary.json",
        {
            "variant": variant,
            "shape": [int(train.shape[0]), int(train.shape[1])],
            "feature_count": int(len(feature_names)),
            "drugs": int(train["canonical_drug_id"].astype(str).nunique()),
            "cell_lines": int(train["cell_line_name"].astype(str).nunique()),
            "valid_smiles_drugs": int(train.loc[valid, "canonical_drug_id"].astype(str).nunique()),
            "valid_smiles_pairs": int(valid.sum()),
            "invalid_smiles_pairs_kept": int((~valid).sum()),
            "label_mean": float(train["label_regression"].mean()),
            "label_std": float(train["label_regression"].std()),
            "feature_blocks": summarize_feature_blocks(feature_names),
        },
    )


def train_variants(config: dict[str, Any], paths: RunPaths, args: argparse.Namespace) -> None:
    variants = [v.strip() for v in (args.variants or ",".join(config["variants"])).split(",") if v.strip()]
    models = [m.strip() for m in (args.models or ",".join(config["benchmark"]["models"])).split(",") if m.strip()]
    for variant in variants:
        variant_dir = paths.step4 / variant
        if not (variant_dir / "train_table.parquet").exists():
            raise FileNotFoundError(f"Missing variant input: {variant_dir}")
        train_one_variant(config, variant_dir, variant, models)


def train_one_variant(config: dict[str, Any], variant_dir: Path, variant: str, models: list[str]) -> None:
    train = pd.read_parquet(variant_dir / "train_table.parquet")
    feature_names = json.loads((variant_dir / "feature_names.json").read_text(encoding="utf-8"))
    X = train[feature_names].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=np.float32)
    y = train["label_regression"].to_numpy(dtype=np.float32)
    folds = int(config["benchmark"].get("folds", 4))
    seed = int(config["benchmark"].get("random_state", 42))
    drug_groups = train["canonical_drug_id"].astype(str).to_numpy()
    scaffold_groups = build_scaffold_groups(train, variant_dir)
    splits = {
        "random_4fold": make_splits(X, y, None, folds, seed, "random"),
        "drug_group_4fold": make_splits(X, y, drug_groups, folds, seed, "group"),
        "scaffold_group_4fold": make_splits(X, y, scaffold_groups, folds, seed, "group"),
    }
    results: list[dict[str, Any]] = []
    oof_store: dict[tuple[str, str], np.ndarray] = {}
    step_dir = variant_dir / "step5_benchmark"
    step_dir.mkdir(parents=True, exist_ok=True)
    for split_name, split_indices in splits.items():
        for model in models:
            print(f"[{variant}] Training {model} / {split_name}", flush=True)
            metrics, oof = fit_model_cv(model, X, y, split_indices, config)
            if metrics.get("status") == "skipped":
                results.append({"variant": variant, "split": split_name, "model": model, **metrics})
                continue
            oof_store[(split_name, model)] = oof
            save_oof(step_dir, train, y, oof, split_name, model)
            results.append({"variant": variant, "split": split_name, "model": model, **metrics})
        ensemble = build_ensemble(step_dir, train, y, split_name, variant, results, oof_store)
        if ensemble:
            results.append(ensemble)
    all_df = pd.DataFrame(results)
    completed = all_df[all_df["status"].eq("completed")].copy()
    if not completed.empty:
        completed = completed.sort_values(["split", "spearman"], ascending=[True, False])
    all_df.to_csv(step_dir / "benchmark_summary_all.csv", index=False)
    completed.to_csv(step_dir / "benchmark_summary_completed.csv", index=False)
    write_json(
        step_dir / "metrics_summary.json",
        {
            "variant": variant,
            "n_rows": int(X.shape[0]),
            "n_features": int(X.shape[1]),
            "models": models,
            "results": results,
        },
    )


def fit_model_cv(
    model: str,
    X: np.ndarray,
    y: np.ndarray,
    split_indices: list[tuple[np.ndarray, np.ndarray]],
    config: dict[str, Any],
) -> tuple[dict[str, Any], np.ndarray]:
    oof = np.zeros(len(y), dtype=np.float32)
    folds = []
    for fold_idx, (tr, va) in enumerate(split_indices, start=1):
        try:
            pred = fit_one_model(model, X[tr], y[tr], X[va], int(config["benchmark"].get("random_state", 42)) + fold_idx, config)
        except ModuleNotFoundError as exc:
            return {"status": "skipped", "reason": str(exc)}, oof
        oof[va] = pred.astype(np.float32)
        fold_metrics = regression_metric_bundle(y[va], pred)
        folds.append({"fold": fold_idx, "n_train": int(len(tr)), "n_valid": int(len(va)), **fold_metrics})
    fold_spearman = [f["spearman"] for f in folds if np.isfinite(f["spearman"])]
    metrics = regression_metric_bundle(y, oof)
    return (
        {
            "status": "completed",
            **metrics,
            "fold_spearman_mean": float(np.mean(fold_spearman)) if fold_spearman else float("nan"),
            "fold_spearman_std": float(np.std(fold_spearman)) if fold_spearman else float("nan"),
            "folds": folds,
        },
        oof,
    )


def fit_one_model(
    model_name: str,
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_valid: np.ndarray,
    seed: int,
    config: dict[str, Any],
) -> np.ndarray:
    n_estimators = int(config["benchmark"].get("tree_estimators", 120))
    if model_name == "lightgbm":
        import lightgbm as lgb

        model = lgb.LGBMRegressor(
            objective="regression",
            n_estimators=n_estimators,
            learning_rate=0.05,
            num_leaves=31,
            subsample=0.85,
            colsample_bytree=0.85,
            random_state=seed,
            n_jobs=-1,
            verbosity=-1,
        )
        model.fit(X_train, y_train)
        return model.predict(X_valid).astype(np.float32)
    if model_name == "xgboost":
        from xgboost import XGBRegressor

        model = XGBRegressor(
            objective="reg:squarederror",
            tree_method="hist",
            n_estimators=n_estimators,
            learning_rate=0.05,
            max_depth=5,
            subsample=0.85,
            colsample_bytree=0.85,
            reg_lambda=1.0,
            random_state=seed,
            n_jobs=4,
        )
        model.fit(X_train, y_train)
        return model.predict(X_valid).astype(np.float32)
    if model_name == "extratrees":
        model = ExtraTreesRegressor(
            n_estimators=n_estimators,
            random_state=seed,
            n_jobs=-1,
            max_features="sqrt",
            min_samples_leaf=1,
        )
        model.fit(X_train, y_train)
        return model.predict(X_valid).astype(np.float32)
    if model_name == "randomforest":
        model = RandomForestRegressor(
            n_estimators=max(80, n_estimators // 2),
            random_state=seed,
            n_jobs=-1,
            max_features="sqrt",
            min_samples_leaf=1,
        )
        model.fit(X_train, y_train)
        return model.predict(X_valid).astype(np.float32)
    if model_name == "residual_mlp":
        return fit_mlp(X_train, y_train, X_valid, seed, int(config["benchmark"].get("mlp_epochs", 10)), residual=True)
    if model_name == "flat_mlp":
        return fit_mlp(X_train, y_train, X_valid, seed, int(config["benchmark"].get("mlp_epochs", 10)), residual=False)
    raise ValueError(f"Unsupported model: {model_name}")


def fit_mlp(X_train: np.ndarray, y_train: np.ndarray, X_valid: np.ndarray, seed: int, epochs: int, residual: bool) -> np.ndarray:
    torch.manual_seed(seed)
    np.random.seed(seed)
    X_train_s, X_valid_s = scale_by_train(X_train, X_valid)
    y_mean = float(y_train.mean())
    y_std = float(y_train.std()) or 1.0
    y_train_s = ((y_train - y_mean) / y_std).astype(np.float32)
    dataset = TensorDataset(torch.from_numpy(X_train_s), torch.from_numpy(y_train_s).unsqueeze(1))
    loader = DataLoader(dataset, batch_size=512, shuffle=True)
    model: nn.Module = ResidualMLP(X_train.shape[1]) if residual else FlatMLP(X_train.shape[1])
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    loss_fn = nn.MSELoss()
    model.train()
    for _ in range(epochs):
        for bx, by in loader:
            opt.zero_grad()
            loss = loss_fn(model(bx), by)
            loss.backward()
            opt.step()
    model.eval()
    with torch.no_grad():
        pred = model(torch.from_numpy(X_valid_s)).squeeze(1).cpu().numpy()
    return (pred * y_std + y_mean).astype(np.float32)


def write_report(config: dict[str, Any], paths: RunPaths) -> None:
    completed_frames = []
    sections = [
        f"# {config['cancer_id']} Hybrid Pipeline Report",
        "",
        f"- Cancer: {config.get('english_name')} ({config.get('korean_name')})",
        f"- Run ID: `{config['run_id']}`",
        f"- Raw S3: `{config['s3_raw_prefix']}`",
        f"- Upload S3: `{config['upload_prefix']}`",
        f"- Fold policy: {config['benchmark'].get('folds', 4)}-fold random/drug/scaffold",
        "",
    ]
    qc_path = paths.step2 / "input_qc.json"
    if qc_path.exists():
        qc = json.loads(qc_path.read_text(encoding="utf-8"))
        qc = augment_qc_with_lincs_policy(config, paths, qc)
        write_json(qc_path, qc)
        sections.extend(
            [
                "## Input QC",
                "",
                f"- Labels: {qc['labels']['rows']} rows, {qc['labels']['cell_lines']} cell lines, {qc['labels']['drugs']} drugs",
                f"- Valid SMILES drugs: {qc['drugs']['valid_smiles_drugs']} / {qc['drugs']['drugs']}",
                f"- DepMap mapped cell lines: {qc['cells']['depmap_mapped']} / {qc['cells']['cell_lines']}",
                f"- CRISPR features: {qc['features']['sample_crispr_features']}",
                f"- LINCS mapped drugs: {qc['features']['lincs_mapped_drugs']}",
                "",
            ]
        )
        if "lincs_policy" in qc:
            sections.extend(render_lincs_policy(qc["lincs_policy"]))
    for variant in config["variants"]:
        variant_dir = paths.step4 / variant
        summary_path = variant_dir / "input_summary.json"
        if not summary_path.exists():
            continue
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        sections.extend(
            [
                f"## {variant}",
                "",
                f"- Input shape: {summary['shape']}",
                f"- Feature count: {summary['feature_count']}",
                f"- Drugs: {summary['drugs']}",
                f"- Cell lines: {summary['cell_lines']}",
                f"- Invalid SMILES pairs kept: {summary['invalid_smiles_pairs_kept']}",
                "",
            ]
        )
        completed_path = variant_dir / "step5_benchmark/benchmark_summary_completed.csv"
        if completed_path.exists():
            completed = pd.read_csv(completed_path)
            completed_frames.append(completed)
            for split in sorted(completed["split"].dropna().unique()):
                top = completed[completed["split"].eq(split)].sort_values("spearman", ascending=False).head(5)
                sections.append(f"### {split}")
                for _, row in top.iterrows():
                    sections.append(f"- {row['model']}: Spearman {float(row['spearman']):.4f}, RMSE {float(row['rmse']):.4f}")
                sections.append("")
    if completed_frames:
        combined = pd.concat(completed_frames, ignore_index=True)
        combined = enrich_completed_metrics(combined, paths)
        combined.to_csv(paths.reports / "benchmark_completed_all_variants.csv", index=False)
        best = combined.sort_values(["split", "spearman"], ascending=[True, False]).groupby("split").head(8)
        best.to_csv(paths.reports / "benchmark_best_by_split.csv", index=False)
        extended_best = combined.sort_values(["split", "spearman"], ascending=[True, False]).groupby("split").head(1)
        extended_best.to_csv(paths.reports / "benchmark_extended_best_by_split.csv", index=False)
        policy = build_policy_comparison(combined)
        if not policy.empty:
            policy.to_csv(paths.reports / "variant_policy_comparison.csv", index=False)
            sections.extend(render_policy_comparison(policy))
        recommendation = build_variant_recommendation(combined)
        write_json(paths.reports / "variant_recommendation.json", recommendation)
        sections.extend(render_variant_recommendation(recommendation))
        sections.extend(render_extended_metrics_table(extended_best))
    (paths.reports / "summary.md").write_text("\n".join(sections).rstrip() + "\n", encoding="utf-8")
    write_json(paths.reports / "run_index.json", build_index(config, paths))


def upload_outputs(config: dict[str, Any], paths: RunPaths) -> None:
    prefix = ensure_s3_slash(config["upload_prefix"])
    run(["aws", "s3", "sync", str(paths.run_dir), prefix + "run/", "--only-show-errors"])
    run(["aws", "s3", "sync", str(paths.reports), prefix + "reports/", "--only-show-errors"])
    run(["aws", "s3", "sync", str(paths.root / "configs"), prefix + "configs/", "--only-show-errors"])
    run(["aws", "s3", "sync", str(paths.root / "docs"), prefix + "docs/", "--only-show-errors"])
    run(["aws", "s3", "cp", str(paths.root / "scripts/run_hybrid_pipeline.py"), prefix + "code/scripts/run_hybrid_pipeline.py", "--only-show-errors"])


def read_source(config: dict[str, Any], paths: RunPaths, key: str) -> pd.DataFrame:
    path = source_local_path(config, paths, key)
    if not path.exists():
        raise FileNotFoundError(f"Missing source {key}: {path}")
    if path.suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if path.suffix == ".csv":
        return pd.read_csv(path)
    if path.suffix == ".gz" and path.name.endswith(".txt.gz"):
        return pd.read_csv(path, sep="\t", compression="gzip")
    if path.suffix == ".txt":
        return pd.read_csv(path, sep="\t")
    if path.suffix == ".json":
        return pd.read_json(path)
    return pd.read_parquet(path)


def build_rdkit_features(drugs: pd.DataFrame, n_bits: int) -> pd.DataFrame:
    ids = drugs["canonical_drug_id"].astype(str).tolist()
    morgan = np.zeros((len(drugs), n_bits), dtype=np.float32)
    rows = []
    for idx, row in drugs.reset_index(drop=True).iterrows():
        smiles = str(row.get("canonical_smiles", "") or "")
        mol = mol_from_smiles(smiles)
        base = {
            "canonical_drug_id": str(row["canonical_drug_id"]),
            "drug_has_valid_smiles": int(mol is not None),
            "drug__smiles_length": float(len(smiles)),
            "drug__descriptor_mol_wt": 0.0,
            "drug__descriptor_logp": 0.0,
            "drug__descriptor_tpsa": 0.0,
            "drug__descriptor_hbd": 0.0,
            "drug__descriptor_hba": 0.0,
            "drug__descriptor_rot_bonds": 0.0,
            "drug__descriptor_ring_count": 0.0,
            "drug__descriptor_heavy_atom_count": 0.0,
            "drug__descriptor_fraction_csp3": 0.0,
        }
        if mol is not None:
            base.update(
                {
                    "drug__descriptor_mol_wt": float(Descriptors.MolWt(mol)),
                    "drug__descriptor_logp": float(Crippen.MolLogP(mol)),
                    "drug__descriptor_tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
                    "drug__descriptor_hbd": float(Lipinski.NumHDonors(mol)),
                    "drug__descriptor_hba": float(Lipinski.NumHAcceptors(mol)),
                    "drug__descriptor_rot_bonds": float(Lipinski.NumRotatableBonds(mol)),
                    "drug__descriptor_ring_count": float(rdMolDescriptors.CalcNumRings(mol)),
                    "drug__descriptor_heavy_atom_count": float(mol.GetNumHeavyAtoms()),
                    "drug__descriptor_fraction_csp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
                }
            )
            bitvect = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
            arr = np.zeros((n_bits,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(bitvect, arr)
            morgan[idx] = arr
        rows.append(base)
    desc = pd.DataFrame(rows)
    morgan_df = pd.DataFrame(morgan, columns=[f"drug_morgan_{i:04d}" for i in range(n_bits)])
    morgan_df.insert(0, "canonical_drug_id", ids)
    return desc.merge(morgan_df, on="canonical_drug_id", how="left")


def build_smiles_svd_features(drugs: pd.DataFrame, requested_dim: int) -> pd.DataFrame:
    base = drugs[["canonical_drug_id", "canonical_smiles"]].copy()
    base["canonical_smiles"] = base["canonical_smiles"].fillna("").astype(str)
    texts = base["canonical_smiles"].where(base["canonical_smiles"].str.strip().ne(""), "__EMPTY__")
    vectorizer = TfidfVectorizer(analyzer="char", ngram_range=(2, 4), min_df=1)
    tfidf = vectorizer.fit_transform(texts)
    dim = max(1, min(requested_dim, tfidf.shape[1] - 1, len(base) - 1))
    svd = TruncatedSVD(n_components=dim, random_state=42)
    vals = svd.fit_transform(tfidf).astype(np.float32)
    out = pd.DataFrame({"canonical_drug_id": base["canonical_drug_id"].astype(str)})
    for i in range(requested_dim):
        out[f"smiles_svd_{i:02d}"] = vals[:, i] if i < vals.shape[1] else 0.0
    return out


def add_lincs_summary(df: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    if not feature_cols:
        for name in ["mean", "std", "abs_mean", "max_abs", "nonzero_frac"]:
            out[f"drug__lincs_summary_{name}"] = 0.0
        return out
    arr = out[feature_cols].fillna(0.0).to_numpy(dtype=np.float32)
    out["drug__lincs_summary_mean"] = arr.mean(axis=1)
    out["drug__lincs_summary_std"] = arr.std(axis=1)
    out["drug__lincs_summary_abs_mean"] = np.abs(arr).mean(axis=1)
    out["drug__lincs_summary_max_abs"] = np.abs(arr).max(axis=1)
    out["drug__lincs_summary_nonzero_frac"] = (arr != 0).mean(axis=1)
    return out


def make_splits(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray | None,
    n_splits: int,
    random_state: int,
    mode: str,
) -> list[tuple[np.ndarray, np.ndarray]]:
    if mode == "group":
        if groups is None:
            raise ValueError("groups required")
        n = min(n_splits, len(np.unique(groups)))
        return list(GroupKFold(n_splits=n).split(X, y, groups=groups))
    n = min(n_splits, len(y))
    return list(KFold(n_splits=n, shuffle=True, random_state=random_state).split(X, y))


def build_scaffold_groups(train: pd.DataFrame, variant_dir: Path) -> np.ndarray:
    rows = []
    for _, row in train[["canonical_drug_id", "drug__canonical_smiles"]].drop_duplicates("canonical_drug_id").iterrows():
        drug_id = str(row["canonical_drug_id"])
        scaffold = compute_scaffold(str(row.get("drug__canonical_smiles", "") or ""))
        rows.append({"canonical_drug_id": drug_id, "scaffold": scaffold or f"NO_SCAFFOLD_{drug_id}"})
    scaffolds = pd.DataFrame(rows)
    scaffolds.to_parquet(variant_dir / "drug_scaffolds.parquet", index=False)
    lookup = scaffolds.set_index("canonical_drug_id")["scaffold"].to_dict()
    return train["canonical_drug_id"].astype(str).map(lookup).fillna("NO_SCAFFOLD_UNKNOWN").to_numpy()


def compute_scaffold(smiles: str) -> str:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol), isomericSmiles=False)
    except Exception:
        return ""


def save_oof(step_dir: Path, train: pd.DataFrame, y: np.ndarray, pred: np.ndarray, split: str, model: str) -> None:
    cols = [c for c in ["pair_id", "sample_id", "canonical_drug_id", "drug_name"] if c in train.columns]
    out = train[cols].copy()
    out["y_true"] = y
    out["y_pred"] = pred
    out["split"] = split
    out["model"] = model
    out.to_parquet(step_dir / f"oof_{split}_{model}.parquet", index=False)


def build_ensemble(
    step_dir: Path,
    train: pd.DataFrame,
    y: np.ndarray,
    split: str,
    variant: str,
    results: list[dict[str, Any]],
    oof_store: dict[tuple[str, str], np.ndarray],
) -> dict[str, Any] | None:
    candidates = [
        r
        for r in results
        if r.get("variant") == variant and r.get("split") == split and r.get("status") == "completed" and r.get("model") != "weighted_top3_ensemble"
    ]
    candidates = sorted(candidates, key=lambda r: float(r.get("spearman", -np.inf)), reverse=True)[:3]
    if not candidates:
        return None
    weights = np.array([max(float(r["spearman"]), 0.0) for r in candidates], dtype=np.float32)
    if float(weights.sum()) <= 1e-8:
        weights = np.ones(len(candidates), dtype=np.float32)
    weights = weights / weights.sum()
    pred = np.zeros(len(y), dtype=np.float32)
    members = []
    for result, weight in zip(candidates, weights):
        model = str(result["model"])
        pred += weight * oof_store[(split, model)]
        members.append({"model": model, "spearman": float(result["spearman"]), "weight": float(weight)})
    save_oof(step_dir, train, y, pred, split, "weighted_top3_ensemble")
    metrics = regression_metric_bundle(y, pred)
    return {
        "variant": variant,
        "split": split,
        "model": "weighted_top3_ensemble",
        "status": "completed",
        **metrics,
        "fold_spearman_mean": np.nan,
        "fold_spearman_std": np.nan,
        "folds": [],
        "members": members,
    }


def write_input_qc(
    config: dict[str, Any],
    paths: RunPaths,
    labels: pd.DataFrame,
    cells: pd.DataFrame,
    drugs: pd.DataFrame,
    sample_crispr: pd.DataFrame,
    drug_features: pd.DataFrame,
    train: pd.DataFrame,
) -> None:
    valid = pd.to_numeric(drug_features.get("drug_has_valid_smiles", 0), errors="coerce").fillna(0).gt(0)
    lincs = pd.to_numeric(drug_features.get("drug__has_lincs_signature", 0), errors="coerce").fillna(0).gt(0)
    qc = {
        "run_id": config["run_id"],
        "cancer_id": config["cancer_id"],
        "labels": {
            "rows": int(labels.shape[0]),
            "cell_lines": int(labels["cell_line_name"].nunique()),
            "drugs": int(labels["canonical_drug_id"].nunique()),
            "label_mean": float(labels["label_regression"].mean()),
            "label_std": float(labels["label_regression"].std()),
            "label_min": float(labels["label_regression"].min()),
            "label_max": float(labels["label_regression"].max()),
        },
        "cells": {
            "cell_lines": int(cells.shape[0]),
            "depmap_mapped": int(pd.to_numeric(cells["is_depmap_mapped"], errors="coerce").fillna(0).sum()),
        },
        "drugs": {
            "drugs": int(drugs.shape[0]),
            "valid_smiles_drugs": int(valid.sum()),
            "lincs_mapped_drugs": int(lincs.sum()),
        },
        "features": {
            "sample_crispr_features": int(sum(c.startswith("sample__crispr__") for c in sample_crispr.columns)),
            "lincs_mapped_drugs": int(lincs.sum()),
            "train_shape": [int(train.shape[0]), int(train.shape[1])],
        },
        "warnings": build_qc_warnings(labels, cells, valid, lincs),
    }
    qc = augment_qc_with_lincs_policy(config, paths, qc)
    write_json(paths.step2 / "input_qc.json", qc)


def load_raw_lincs_summary(paths: RunPaths) -> dict[str, Any] | None:
    summary_path = paths.step2 / "raw_lincs_build_summary.json"
    if not summary_path.exists():
        return None
    return json.loads(summary_path.read_text(encoding="utf-8"))


def build_lincs_qc_details(raw_summary: dict[str, Any]) -> dict[str, Any]:
    return {
        "mode": str(raw_summary.get("mode", "")),
        "source": str(raw_summary.get("source", "")),
        "cell_filter": str(raw_summary.get("cell_filter", "")),
        "cell_column": str(raw_summary.get("cell_column", "")),
        "matched_drugs": int(raw_summary.get("matched_drugs", 0)),
        "matched_perturbagens": int(raw_summary.get("matched_perturbagens", 0)),
        "matched_signatures": int(raw_summary.get("matched_signatures", 0)),
        "included_total_cell_line_count": int(raw_summary.get("included_total_cell_line_count", 0)),
        "available_cancer_overlap_cell_ids": [str(value) for value in raw_summary.get("available_cancer_overlap_cell_ids", [])],
        "available_cancer_overlap_cell_line_count": int(raw_summary.get("available_cancer_overlap_cell_line_count", 0)),
        "available_representative_cell_ids": [str(value) for value in raw_summary.get("available_representative_cell_ids", [])],
        "available_representative_cell_line_count": int(raw_summary.get("available_representative_cell_line_count", 0)),
    }


def build_lincs_policy_decision(config: dict[str, Any], raw_summary: dict[str, Any]) -> dict[str, Any]:
    policy_cfg = config.get("lincs_policy", {})
    representative_cells = [str(value) for value in policy_cfg.get("representative_cell_lines", [])]
    representative_norm = {norm_key(value) for value in representative_cells}
    available_overlap = [str(value) for value in raw_summary.get("available_cancer_overlap_cell_ids", [])]
    available_representative = [value for value in available_overlap if norm_key(value) in representative_norm]
    matched_drugs = int(raw_summary.get("matched_drugs", 0))
    matched_signatures = int(raw_summary.get("matched_signatures", 0))
    available_overlap_count = int(raw_summary.get("available_cancer_overlap_cell_line_count", len(available_overlap)))
    single_cell_min_drugs = int(policy_cfg.get("single_cell_main_min_mapped_drugs", 80))
    single_cell_min_signatures = int(policy_cfg.get("single_cell_main_min_signatures", 5000))
    filter_policy = str(raw_summary.get("cell_filter", ""))
    mode = str(raw_summary.get("mode", ""))

    strategy = "unclassified"
    main_run_eligible = False
    reason = "No LINCS policy rule matched."

    if mode == "none":
        strategy = "no_lincs"
        reason = "LINCS mode is disabled for this run."
    elif filter_policy == "all":
        strategy = "all_cell_main"
        main_run_eligible = True
        if available_overlap_count == 0:
            reason = "No directly overlapping disease LINCS cell line was available after drug matching, so all-cell drug-level LINCS is used as the fallback main setting."
        elif available_overlap_count == 1:
            reason = f"Only one directly overlapping disease LINCS cell line ({available_overlap[0]}) was available after drug matching, so all-cell LINCS remains the main setting and disease-restricted LINCS should be treated as sensitivity/supporting evidence."
        else:
            reason = f"This run keeps all-cell LINCS although {available_overlap_count} overlapping disease LINCS cell lines were available ({', '.join(available_overlap)}). In that situation, a disease-restricted LINCS main run should usually be preferred."
            main_run_eligible = False
    elif available_overlap_count >= 2:
        strategy = "cancer_specific_main"
        main_run_eligible = True
        reason = f"{available_overlap_count} overlapping disease LINCS cell lines were available ({', '.join(available_overlap)}), so cancer-specific LINCS can be used as the main setting."
    elif available_overlap_count == 1 and available_representative and matched_drugs >= single_cell_min_drugs and matched_signatures >= single_cell_min_signatures:
        strategy = "single_representative_main"
        main_run_eligible = True
        reason = (
            f"Only one overlapping disease LINCS cell line ({available_overlap[0]}) was available, but it is declared as a representative disease cell line and coverage is adequate "
            f"({matched_drugs} matched drugs, {matched_signatures} matched signatures), so single-cell cancer-specific LINCS can be used as the main setting."
        )
    elif available_overlap_count == 1:
        strategy = "cancer_specific_supporting_only"
        reason = (
            f"Only one overlapping disease LINCS cell line ({available_overlap[0]}) was available and it does not meet the representative single-cell main criteria "
            f"({single_cell_min_drugs} matched drugs and {single_cell_min_signatures} signatures by default), so this disease-restricted LINCS run should be treated as sensitivity/supporting evidence."
        )
    else:
        strategy = "cancer_specific_unavailable"
        reason = "No overlapping disease LINCS cell line was available, so a cancer-specific LINCS main run is unavailable."

    return {
        "strategy": strategy,
        "main_run_eligible": bool(main_run_eligible),
        "selection_reason": reason,
        "cell_filter": filter_policy,
        "representative_cell_lines_declared": representative_cells,
        "available_cancer_overlap_cell_ids": available_overlap,
        "available_cancer_overlap_cell_line_count": available_overlap_count,
        "available_representative_cell_ids": available_representative,
        "available_representative_cell_line_count": int(len(available_representative)),
        "matched_drugs": matched_drugs,
        "matched_signatures": matched_signatures,
        "single_cell_main_min_mapped_drugs": single_cell_min_drugs,
        "single_cell_main_min_signatures": single_cell_min_signatures,
    }


def augment_qc_with_lincs_policy(config: dict[str, Any], paths: RunPaths, qc: dict[str, Any]) -> dict[str, Any]:
    raw_summary = load_raw_lincs_summary(paths)
    if raw_summary is None:
        return qc
    out = dict(qc)
    out["lincs"] = build_lincs_qc_details(raw_summary)
    out["lincs_policy"] = build_lincs_policy_decision(config, raw_summary)
    write_json(paths.step2 / "lincs_policy_decision.json", out["lincs_policy"])
    return out


def build_qc_warnings(labels: pd.DataFrame, cells: pd.DataFrame, valid: pd.Series, lincs: pd.Series) -> list[str]:
    warnings_out = []
    if labels["cell_line_name"].nunique() < 10:
        warnings_out.append("cell_line_count_lt_10")
    if float(valid.mean()) < 0.7:
        warnings_out.append("valid_smiles_drug_rate_lt_70pct")
    if int(lincs.sum()) < 20:
        warnings_out.append("lincs_mapped_drugs_lt_20")
    if int(pd.to_numeric(cells["is_depmap_mapped"], errors="coerce").fillna(0).sum()) < cells.shape[0]:
        warnings_out.append("some_cell_lines_missing_depmap_mapping")
    return warnings_out


def render_lincs_policy(decision: dict[str, Any]) -> list[str]:
    overlap = decision.get("available_cancer_overlap_cell_ids", [])
    representatives = decision.get("representative_cell_lines_declared", [])
    return [
        "## LINCS Policy",
        "",
        f"- Strategy: `{decision.get('strategy', 'unknown')}`",
        f"- Main-run eligible: {'yes' if decision.get('main_run_eligible') else 'no'}",
        f"- Disease-overlap LINCS cells: {', '.join(overlap) if overlap else 'none'}",
        f"- Representative cells declared: {', '.join(representatives) if representatives else 'none'}",
        f"- Matched drugs: {int(decision.get('matched_drugs', 0))}",
        f"- Matched signatures: {int(decision.get('matched_signatures', 0))}",
        f"- Reason: {decision.get('selection_reason', '')}",
        "",
    ]


def enrich_completed_metrics(combined: pd.DataFrame, paths: RunPaths) -> pd.DataFrame:
    """Backfill extended metrics from saved OOF predictions.

    Older benchmark CSVs only contain Spearman and RMSE. The OOF files are the
    source of truth, so reports can be regenerated without retraining.
    """

    out = combined.copy()
    threshold_cache: dict[str, float] = {}
    metric_rows: list[dict[str, Any]] = []
    for _, row in out.iterrows():
        if row.get("status") != "completed":
            metric_rows.append({})
            continue
        variant = str(row["variant"])
        split = str(row["split"])
        model = str(row["model"])
        oof_path = paths.step4 / variant / "step5_benchmark" / f"oof_{split}_{model}.parquet"
        if not oof_path.exists():
            metric_rows.append({})
            continue
        oof = pd.read_parquet(oof_path)
        y_true = pd.to_numeric(oof["y_true"], errors="coerce").to_numpy(dtype=np.float32)
        y_pred = pd.to_numeric(oof["y_pred"], errors="coerce").to_numpy(dtype=np.float32)
        metrics = regression_metric_bundle(y_true, y_pred)
        metrics["n_eval"] = int(len(oof))

        if variant not in threshold_cache:
            train_table = pd.read_parquet(
                paths.step4 / variant / "train_table.parquet",
                columns=["binary_threshold"],
            )
            thresholds = pd.to_numeric(train_table["binary_threshold"], errors="coerce").dropna().unique()
            threshold_cache[variant] = float(thresholds[0]) if len(thresholds) else float("nan")
        threshold = threshold_cache[variant]
        if math.isnan(threshold):
            y_binary = np.zeros(len(y_true), dtype=int)
        else:
            y_binary = (y_true <= threshold).astype(int)
        metrics.update(binary_metric_bundle(y_binary, -y_pred))
        metric_rows.append(metrics)

    metrics_df = pd.DataFrame(metric_rows)
    for column in metrics_df.columns:
        out[column] = metrics_df[column]
    return out


def build_policy_comparison(combined: pd.DataFrame) -> pd.DataFrame:
    target_variants = ["legacy_rich_valid_smiles_only", "legacy_rich_all_drugs_zero_smiles", "colonstyle_compact_baseline"]
    rows = []
    for split in sorted(combined["split"].dropna().unique()):
        row = {"split": split}
        for variant in target_variants:
            sub = combined[(combined["variant"].eq(variant)) & (combined["split"].eq(split))]
            if sub.empty:
                continue
            best = sub.sort_values("spearman", ascending=False).iloc[0]
            row[f"{variant}_best_model"] = best["model"]
            row[f"{variant}_spearman"] = float(best["spearman"])
            row[f"{variant}_rmse"] = float(best["rmse"])
        if "legacy_rich_valid_smiles_only_spearman" in row and "legacy_rich_all_drugs_zero_smiles_spearman" in row:
            row["valid_minus_zero_delta"] = row["legacy_rich_valid_smiles_only_spearman"] - row["legacy_rich_all_drugs_zero_smiles_spearman"]
        rows.append(row)
    return pd.DataFrame(rows)


def render_policy_comparison(policy: pd.DataFrame) -> list[str]:
    lines = ["## Variant Policy Comparison", ""]
    lines.append("| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |")
    lines.append("|---|---:|---:|---:|")
    for _, row in policy.iterrows():
        valid = row.get("legacy_rich_valid_smiles_only_spearman", np.nan)
        zero = row.get("legacy_rich_all_drugs_zero_smiles_spearman", np.nan)
        delta = row.get("valid_minus_zero_delta", np.nan)
        lines.append(f"| {row['split']} | {valid:.4f} | {zero:.4f} | {delta:+.4f} |")
    lines.append("")
    return lines


def build_variant_recommendation(combined: pd.DataFrame) -> dict[str, Any]:
    completed = combined[combined["status"].eq("completed")].copy()
    best_rows = []
    for (variant, split), sdf in completed.groupby(["variant", "split"]):
        best = sdf.sort_values("spearman", ascending=False).iloc[0]
        best_rows.append(
            {
                "variant": variant,
                "split": split,
                "model": best["model"],
                "spearman": float(best["spearman"]),
                "rmse": float(best["rmse"]),
            }
        )
    best = pd.DataFrame(best_rows)
    generalization = best[best["split"].isin(["drug_group_4fold", "scaffold_group_4fold"])]
    scores = (
        generalization.groupby("variant")["spearman"].mean().sort_values(ascending=False).reset_index(name="mean_drug_scaffold_spearman")
    )
    recommended = scores.iloc[0].to_dict() if not scores.empty else {}
    split_winners = (
        best.sort_values(["split", "spearman"], ascending=[True, False])
        .groupby("split")
        .head(1)
        .to_dict(orient="records")
    )
    return {
        "recommended_variant": recommended.get("variant"),
        "recommendation_basis": "highest mean Spearman across drug_group_4fold and scaffold_group_4fold",
        "mean_generalization_scores": scores.to_dict(orient="records"),
        "split_winners": split_winners,
    }


def render_variant_recommendation(recommendation: dict[str, Any]) -> list[str]:
    lines = ["## Recommended Input Policy", ""]
    variant = recommendation.get("recommended_variant")
    if not variant:
        lines.append("- Recommendation unavailable because benchmark results are incomplete.")
        lines.append("")
        return lines
    lines.append(f"- Recommended variant: `{variant}`")
    lines.append(f"- Basis: {recommendation.get('recommendation_basis')}")
    lines.append("")
    lines.append("| Variant | Mean drug/scaffold Spearman |")
    lines.append("|---|---:|")
    for row in recommendation.get("mean_generalization_scores", []):
        lines.append(f"| {row['variant']} | {float(row['mean_drug_scaffold_spearman']):.4f} |")
    lines.append("")
    return lines


def render_extended_metrics_table(best: pd.DataFrame) -> list[str]:
    lines = ["## Extended Metrics For Split Winners", ""]
    lines.append(
        "| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |"
    )
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for _, row in best.iterrows():
        lines.append(
            "| {split} | {variant} | {model} | {spearman:.4f} | {pearson:.4f} | {kendall:.4f} | "
            "{rmse:.4f} | {mae:.4f} | {r2:.4f} | {auroc:.4f} | {auprc:.4f} |".format(
                split=row["split"],
                variant=row["variant"],
                model=row["model"],
                spearman=float(row.get("spearman", np.nan)),
                pearson=float(row.get("pearson", np.nan)),
                kendall=float(row.get("kendall", np.nan)),
                rmse=float(row.get("rmse", np.nan)),
                mae=float(row.get("mae", np.nan)),
                r2=float(row.get("r2", np.nan)),
                auroc=float(row.get("sensitivity_auroc", np.nan)),
                auprc=float(row.get("sensitivity_auprc", np.nan)),
            )
        )
    lines.append("")
    lines.append("Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.")
    lines.append("")
    return lines


def build_index(config: dict[str, Any], paths: RunPaths) -> dict[str, Any]:
    return {
        "run_id": config["run_id"],
        "cancer_id": config["cancer_id"],
        "root": str(paths.root),
        "run_dir": str(paths.run_dir),
        "reports": str(paths.reports),
        "summary": str(paths.reports / "summary.md"),
        "upload_prefix": config["upload_prefix"],
    }


def summarize_feature_blocks(feature_names: list[str]) -> dict[str, int]:
    return {
        "crispr": sum(c.startswith("sample__crispr__") for c in feature_names),
        "morgan": sum(c.startswith("drug_morgan_") for c in feature_names),
        "smiles_svd": sum(c.startswith("smiles_svd_") for c in feature_names),
        "lincs_gene": sum(c.startswith("drug__lincs__") for c in feature_names),
        "lincs_summary": sum(c.startswith("drug__lincs_summary_") for c in feature_names),
        "target_pair": sum(c.startswith("pair__") or c.startswith("drug__target") for c in feature_names),
        "context": sum((c.startswith("sample__") and not c.startswith("sample__crispr__")) or c.startswith("ctxcat__") for c in feature_names),
    }


def top_by_variance(df: pd.DataFrame, cols: list[str], max_count: int) -> list[str]:
    if not cols:
        return []
    variances = df[cols].var(axis=0, ddof=0).sort_values(ascending=False)
    return variances.head(max_count).index.tolist()


def select_top_variance_columns(df: pd.DataFrame, cols: list[str], max_count: int) -> list[str]:
    numeric = [c for c in cols if pd.api.types.is_numeric_dtype(df[c])]
    if not numeric:
        return []
    variances = df[numeric].var(axis=0, ddof=0).sort_values(ascending=False)
    variances = variances[variances > 1e-12]
    return variances.head(max_count).index.tolist()


def metric_spearman(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 3 or np.std(y_true) <= 1e-12 or np.std(y_pred) <= 1e-12:
        return float("nan")
    value = pd.Series(y_true).corr(pd.Series(y_pred), method="spearman")
    return float(value) if value is not None else float("nan")


def metric_pearson(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 3 or np.std(y_true) <= 1e-12 or np.std(y_pred) <= 1e-12:
        return float("nan")
    value = pd.Series(y_true).corr(pd.Series(y_pred), method="pearson")
    return float(value) if value is not None else float("nan")


def metric_kendall(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 3 or np.std(y_true) <= 1e-12 or np.std(y_pred) <= 1e-12:
        return float("nan")
    value = pd.Series(y_true).corr(pd.Series(y_pred), method="kendall")
    return float(value) if value is not None else float("nan")


def metric_rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(math.sqrt(mean_squared_error(y_true, y_pred)))


def regression_metric_bundle(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    return {
        "spearman": metric_spearman(y_true, y_pred),
        "pearson": metric_pearson(y_true, y_pred),
        "kendall": metric_kendall(y_true, y_pred),
        "rmse": metric_rmse(y_true, y_pred),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "r2": float(r2_score(y_true, y_pred)),
    }


def binary_metric_bundle(y_binary: np.ndarray, score: np.ndarray) -> dict[str, float]:
    if len(np.unique(y_binary)) < 2:
        return {"sensitivity_auroc": float("nan"), "sensitivity_auprc": float("nan")}
    return {
        "sensitivity_auroc": float(roc_auc_score(y_binary, score)),
        "sensitivity_auprc": float(average_precision_score(y_binary, score)),
    }


def scale_by_train(X_train: np.ndarray, X_valid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = X_train.mean(axis=0, keepdims=True)
    std = X_train.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    return ((X_train - mean) / std).astype(np.float32), ((X_valid - mean) / std).astype(np.float32)


def mol_from_smiles(smiles: str) -> Chem.Mol | None:
    value = str(smiles or "").strip()
    if not value:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return Chem.MolFromSmiles(value)
        except Exception:
            return None


def extract_targets(value: Any) -> list[str]:
    text = str(value or "").upper()
    for src, dst in TARGET_ALIASES.items():
        text = text.replace(src, dst)
    tokens = re.split(r"[^A-Z0-9]+", text)
    out = []
    for token in tokens:
        if len(token) < 2 or token in GENERIC_TARGET_TOKENS:
            continue
        if re.fullmatch(r"[A-Z][A-Z0-9]{1,8}", token):
            out.append(TARGET_ALIASES.get(token, token))
    return sorted(set(out))


def parse_gene_symbol(column: str) -> str:
    symbol = str(column).split(" (", 1)[0].strip().upper()
    return clean_feature_token(symbol).upper()


def clean_feature_token(value: Any) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip().lower()).strip("_")
    return token or "unknown"


def canonical_gdsc_column_map(columns: Any) -> dict[str, str]:
    cols = set(columns)
    mapping = {
        "CELL_LINE_NAME": "cell_line_name",
        "DRUG_NAME": "drug_name",
    }
    if "PUTATIVE_TARGET_NORMALIZED" in cols:
        mapping["PUTATIVE_TARGET_NORMALIZED"] = "putative_target"
    elif "PUTATIVE_TARGET" in cols:
        mapping["PUTATIVE_TARGET"] = "putative_target"
    if "PATHWAY_NAME_NORMALIZED" in cols:
        mapping["PATHWAY_NAME_NORMALIZED"] = "pathway_name"
    elif "PATHWAY_NAME" in cols:
        mapping["PATHWAY_NAME"] = "pathway_name"
    return mapping


def norm_key(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def split_synonyms(value: Any) -> list[str]:
    text = str(value or "")
    if not text or text.lower() == "nan":
        return []
    return [item.strip() for item in re.split(r"[,;|]", text) if item.strip()]


def first_existing_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for column in candidates:
        if column in df.columns:
            return column
    return None


def iter_lincs_names(row: pd.Series, name_cols: list[str], alias_cols: list[str]) -> list[str]:
    names: list[str] = []
    for column in name_cols:
        value = row.get(column, "")
        if pd.notna(value) and str(value).strip():
            names.append(str(value).strip())
    for column in alias_cols:
        names.extend(split_synonyms(row.get(column, "")))
    return names


def canonicalize_smiles(value: Any) -> str:
    text = str(value or "").strip()
    if not text or text.lower() in {"nan", "none", "-666", "restricted"}:
        return ""
    mol = mol_from_smiles(text)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def ensure_s3_slash(prefix: str) -> str:
    return prefix if prefix.endswith("/") else prefix + "/"


def source_key_exists(config: dict[str, Any], key: str) -> bool:
    return key in config.get("source_files", {})


def source_file_specs(config: dict[str, Any]) -> list[tuple[str, str, str]]:
    specs: list[tuple[str, str, str]] = []
    for key, spec in config.get("source_files", {}).items():
        if isinstance(spec, str):
            rel = spec
            s3_uri = ensure_s3_slash(config.get("s3_raw_prefix", "")) + rel
        else:
            rel = str(spec.get("local") or spec.get("path") or Path(str(spec["s3_uri"])).name)
            s3_uri = str(spec.get("s3_uri") or ensure_s3_slash(config.get("s3_raw_prefix", "")) + str(spec["path"]))
        specs.append((key, rel, s3_uri))
    return specs


def source_local_path(config: dict[str, Any], paths: RunPaths, key: str) -> Path:
    for source_key, rel, _ in source_file_specs(config):
        if source_key == key:
            return paths.raw_cache / rel
    raise KeyError(f"Unknown source file key: {key}")


def deep_update(base: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_update(out[key], value)
        else:
            out[key] = value
    return out


def run(cmd: list[str]) -> None:
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
