# OncoRegimenTester

White-box synthetic OMOP CDM simulator for oncology regimen detection method testing.

## Overview

OncoRegimenTester generates structurally believable oncology treatment histories as OMOP-like tables. It creates realistic drug exposure patterns from configurable regimen templates -- patterns that are intentionally messy (timing jitter, supportive drugs, modifications, gaps) to mimic what real claims or EHR data looks like.

The key insight: in real data we only have drug exposures. We never directly observe regimens. Regimen detection is the hard part. This simulator lets you generate drug exposures from known regimen templates (ground truth), then test whether your detection algorithm can reconstruct the regimens from those exposures alone.

### Use cases

- Testing regimen detection algorithms against known ground truth
- Benchmarking line-of-therapy identification methods
- Developing methods to reconstruct regimens from raw drug exposures
- Creating shareable synthetic datasets for oncology method research

## Installation

```r
# Install from GitHub
remotes::install_github("ablack3033/OncoRegimenTester")

# Or clone and install locally
# git clone https://github.com/ablack3033/OncoRegimenTester.git
# devtools::install("OncoRegimenTester")

# For development (load without installing)
devtools::load_all()
```

### Dependencies

**Required**: `data.table`

**Optional** (for CDMConnector integration): `CDMConnector`, `DBI`, `duckdb`

## Quick start

### Minimal example

```r
library(OncoRegimenTester)

# Run with defaults: 100 patients, mCRC cohort, no noise
result <- runSimulation()

# Examine the generated drug exposures
head(result$therapy$drugExposures)

# Compare ground truth regimens to detected regimens
result$therapy$groundTruthEpisodes  # what we generated
result$detected$episodes             # what the detector found
```

### With configuration

```r
config <- simulatorConfig(
  paramDir = system.file("example_params", "mCRC", package = "OncoRegimenTester"),
  outputDir = "my_output",
  nPatients = 500,
  seed = 42,
  noise = noiseConfig(
    probabilityJitterSd = 0.1,
    durationJitterDays = 3,
    cycleJitterDays = 1,
    missingnessRate = 0.05
  )
)

result <- runSimulation(config)
```

### Step-by-step pipeline

If you want more control over individual steps:

```r
# 1. Load and validate parameters
pack <- loadParameterPack("inst/example_params/mCRC")

# 2. Generate patients
config <- simulatorConfig(nPatients = 200, seed = 42)
patients <- generatePatients(pack, config)

# 3. Generate therapy (line plans + drug exposures)
therapy <- generateTherapy(patients, pack, config)

# 4. Generate supporting events (conditions, procedures, measurements)
supporting <- generateSupportingEvents(patients, therapy$linePlans, pack, config)

# 5. Run reference detector
detected <- detectRegimens(therapy$drugExposures, pack)

# 6. Render OMOP tables to CSV
renderOmop(patients, therapy, supporting, detected, config)

# 7. Validate output vs parameters
validation <- validateOutput(patients, therapy, pack)
cat(paste(validation$report, collapse = "\n"))
```

### CDMConnector integration

Write synthetic data directly into a CDMConnector CDM reference:

```r
library(CDMConnector)
library(DBI)
library(duckdb)

# Create an in-memory DuckDB CDM
con <- dbConnect(duckdb())
cdm <- cdmFromCon(con, cdmSchema = "main", writeSchema = "main")

# Generate and insert synthetic data
result <- runSimulation(simulatorConfig(nPatients = 100, omopOutput = FALSE))
writeToCdm(cdm, result$patients, result$therapy,
           result$supporting, result$detected, result$config)

# Now use the CDM with OHDSI tools
cdm$drug_exposure
```

## CLI usage

```bash
Rscript inst/scripts/run_simulator.R \
  --param_dir inst/example_params/mCRC \
  --output_dir output \
  --n_patients 200 \
  --seed 42 \
  --noise_prob 0.1 \
  --noise_dur 3 \
  --noise_cycle 1 \
  --noise_miss 0.05
```

## Output tables

### OMOP-like tables (lower_snake_case column names)

| File | Description |
|------|-------------|
| `person.csv` | Synthetic patient demographics |
| `observation_period.csv` | Observation windows |
| `drug_exposure.csv` | Individual drug exposure records (the primary output) |
| `condition_occurrence.csv` | Optional condition events |
| `procedure_occurrence.csv` | Optional procedure events |
| `measurement.csv` | Optional measurement events |

### Ground truth (what the simulator intended)

| File | Description |
|------|-------------|
| `ground_truth_regimen_episodes.csv` | The regimens that were generated |
| `ground_truth_lines.csv` | The lines of therapy that were planned |

### Derived (from the reference detector)

| File | Description |
|------|-------------|
| `derived_regimen_episodes.csv` | Regimens recovered from drug exposures |
| `derived_lines.csv` | Lines recovered from drug exposures |

The comparison between ground truth and derived tables is the core evaluation loop.

## Parameter pack

The simulator is driven entirely by interpretable CSV files. See `inst/example_params/mCRC/` for a complete example.

| File | Purpose |
|------|---------|
| `manifest.csv` | Metadata about the parameter pack |
| `cohort_population.csv` | Demographic mix (age, sex, observation length) |
| `treatment_profile_mix.csv` | Latent treatment profile prevalence |
| `drug_feature_catalog.csv` | Drug definitions (antineoplastic vs supportive, route) |
| `regimen_catalog.csv` | Regimen templates |
| `regimen_drug_map.csv` | Which drugs belong to which regimens |
| `profile_regimen_start.csv` | Starting regimen probabilities by profile and line |
| `regimen_transition.csv` | Transition probabilities between regimens |
| `regimen_duration.csv` | Duration distributions per regimen |
| `drug_schedule.csv` | Cycle timing (e.g., oxaliplatin q14d day 1) |
| `drug_drop_add_rules.csv` | Stochastic drug modifications |
| `supportive_drug_rules.csv` | Supportive drugs that should be ignored by detection |
| `condition_signal.csv` | Optional condition events |
| `procedure_signal.csv` | Optional procedure events |
| `measurement_signal.csv` | Optional measurement events |
| `measurement_value_dist.csv` | Measurement value distributions |

## Noise options

Generation-time noise adds realism. This is NOT differential privacy -- it is simulation variability.

| Parameter | What it does |
|-----------|--------------|
| `probabilityJitterSd` | Perturbs sampling probabilities on logit scale |
| `durationJitterDays` | Varies regimen durations within bins |
| `cycleJitterDays` | Shifts individual administration days within cycles |
| `missingnessRate` | Randomly drops exposure records |

## Testing

```r
devtools::test()
```

## Design principles

- **White-box only**: All parameters are human-readable CSVs. No hidden model weights.
- **Deterministic**: Same seed produces same output.
- **Drug-centered**: Drug exposures are the primary signal. Other domains support plausibility.
- **Modular**: Each step is inspectable and replaceable.
- **Honest about scope**: The reference detector is for self-validation only. Production regimen detection belongs in a separate package.

## Code style

See `STYLE_GUIDE.md` for conventions used in this package.
