# OncoRegimenTester

White-box synthetic OMOP CDM simulator for oncology regimen detection method testing.

## What It Does

Generates fake but structurally believable oncology treatment histories in OMOP-like tables. The simulator creates drug exposure patterns from configurable regimen templates, producing data suitable for:

- Testing regimen detection algorithms
- Benchmarking line-of-therapy identification
- Developing methods to reconstruct regimens from raw drug exposures
- Creating shareable synthetic datasets

**Key design principle**: Drugs drive the realism. We generate from known regimen templates (the ground truth), but the output drug exposures are intentionally messy -- with timing jitter, supportive drugs, modifications, and gaps. A separate regimen detection package then attempts to reconstruct regimens from those exposures alone.

## Installation

```r
# Install from source
devtools::install(".")

# Or install dependencies and load
install.packages("data.table")
devtools::load_all()
```

## Quick Start

```r
library(OncoRegimenTester)

# Run with defaults (100 patients, mCRC cohort)
result <- run_simulation()

# Or configure explicitly
config <- simulator_config(
  param_dir = system.file("example_params", "mCRC", package = "OncoRegimenTester"),
  output_dir = "my_output",
  n_patients = 500,
  seed = 42,
  noise = noise_config(
    probability_jitter_sd = 0.1,
    duration_jitter_days = 3,
    cycle_jitter_days = 1,
    missingness_rate = 0.05
  )
)

result <- run_simulation(config)
```

## CLI Usage

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

## Output Tables

### OMOP-like tables
- `person.csv` -- synthetic patient demographics
- `observation_period.csv` -- observation windows
- `drug_exposure.csv` -- individual drug exposure records (the primary output)
- `condition_occurrence.csv` -- optional condition events
- `procedure_occurrence.csv` -- optional procedure events
- `measurement.csv` -- optional measurement events

### Ground truth (what the simulator intended)
- `ground_truth_regimen_episodes.csv` -- the regimens that were generated
- `ground_truth_lines.csv` -- the lines of therapy that were planned

### Detected (from the reference detector)
- `derived_regimen_episodes.csv` -- regimens recovered from drug exposures
- `derived_lines.csv` -- lines recovered from drug exposures

The comparison between ground truth and derived tables is the core evaluation loop for testing regimen detection methods.

## Parameter Pack

The simulator is driven by interpretable CSV parameter files. See `inst/example_params/mCRC/` for a complete example.

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

## Architecture

```
R/
  config.R              -- simulator_config() and noise_config()
  noise.R               -- jitter functions for realism
  parameter_loader.R    -- CSV loading + validation
  patient_generator.R   -- demographic + profile sampling
  therapy_generator.R   -- lines, regimens, drug exposures
  supporting_event_generator.R -- conditions, procedures, measurements
  regimen_detector.R    -- minimal reference detector (NOT production)
  omop_renderer.R       -- OMOP-like CSV output
  validation.R          -- output vs parameter comparison
  run_simulation.R      -- top-level pipeline
```

## Noise Options

Generation-time noise adds realism (this is NOT differential privacy):

- **probability_jitter_sd**: Perturbs sampling probabilities on logit scale
- **duration_jitter_days**: Varies regimen durations within bins
- **cycle_jitter_days**: Shifts individual admin days within cycles
- **missingness_rate**: Randomly drops exposure records

## Testing

```r
devtools::test()
```

## Design Principles

- **White-box only**: All parameters are human-readable CSVs
- **Deterministic**: Same seed produces same output
- **Drug-centered**: Drug exposures are the primary signal
- **Modular**: Each step is inspectable and replaceable
- **Honest about scope**: The reference detector is for self-validation only; production regimen detection belongs in a separate package
