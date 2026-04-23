# Hybrid Cancer Drug Repurposing Pipeline

이 레포는 팀별 암종 파이프라인을 하나의 공통 구조로 정리하기 위한 작업 공간입니다.

핵심 방향은 다음과 같습니다.

- **Raw source-first 품질 관리**는 colon-style처럼 가져갑니다.
- **최종 모델 입력셋**은 성능이 좋았던 legacy-rich feature 구성을 기본으로 가져갑니다.
- 암종마다 `valid SMILES only`, `all drugs + zero SMILES`, `compact colon-style baseline`을 같이 만들어 같은 4-fold split으로 비교합니다.
- 최종 약물 선정 이후의 ADMET, KG, ClinicalTrials, DrugBank 설명은 모델 ranking을 바꾸는 단계가 아니라 해석/보조 설명 단계로 분리합니다.

## Quick Start

PAAD 췌장암 예시는 아래처럼 실행합니다.

```bash
python3 scripts/run_hybrid_pipeline.py \
  --config configs/paad.json \
  --stage all \
  --upload-s3
```

`configs/paad.json`은 `PAAD_raw/`에 이미 놓여 있던 staged parquet/CSV를 사용하는 빠른 재현용 설정입니다. 진짜 원천 파일에서 core input을 다시 검증하려면 아래 설정을 사용합니다.

```bash
python3 scripts/run_hybrid_pipeline.py \
  --config configs/paad_original_core.json \
  --stage all \
  --upload-s3
```

중간부터 다시 실행할 때는 `--stage build`, `--stage train`, `--stage report`, `--stage upload` 중 하나를 사용합니다.

## Pipeline Layout

```text
configs/
  paad.json
  paad_original_core.json
docs/
  hybrid_pipeline_design.md
  paad_original_source_audit.md
scripts/
  run_hybrid_pipeline.py
data/
  raw_cache/                 # S3 raw subset cache, git ignored
  processed/runs/{run_id}/   # model inputs and benchmarks, git ignored
outputs/
  reports/{run_id}/          # compact reports, committed when useful
```

## Standard Variants

| Variant | Role |
|---|---|
| `legacy_rich_valid_smiles_only` | 기본 성능형 입력셋. invalid/no-SMILES 약물 제거 |
| `legacy_rich_all_drugs_zero_smiles` | no-SMILES 약물을 유지하고 SMILES-derived feature를 0 처리 |
| `colonstyle_compact_baseline` | raw-first 품질 검증용 compact control |

## Standard Evaluation

모든 variant는 같은 기준으로 평가합니다.

- `random_4fold`
- `drug_group_4fold`
- `scaffold_group_4fold`
- `weighted_top3_ensemble`

Random CV는 가까운 분포 내 예측력을, drug/scaffold split은 새로운 약물/새로운 화학 구조에 대한 일반화 가능성을 보기 위한 기준입니다.
