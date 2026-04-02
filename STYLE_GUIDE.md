# Code Style Guide

## Naming conventions

### R functions and variables: `camelCase`

All R function names and local variable names use camelCase.

```r
# Functions
loadParameterPack()
generatePatients()
simulatorConfig()
noiseConfig()

# Variables
paramDir <- "path/to/params"
nPatients <- 100L
calDays <- as.integer(endDate - startDate)
chosenDemo <- pop[demoIdx]
```

### Database table names and column names: `lower_snake_case`

All database/OMOP table names and their column names use lower_snake_case. This matches the OMOP CDM convention and SQL conventions.

```r
# Table names
"person"
"drug_exposure"
"observation_period"
"condition_occurrence"
"ground_truth_regimen_episodes"

# Column names in data.tables that represent database rows
data.table(
  person_id = 1L,
  drug_exposure_start_date = as.Date("2020-01-01"),
  regimen_id = "FOLFOX",
  line_number = 1L,
  is_supportive = FALSE
)
```

### File names: `lower_snake_case`

All file names (R source files, CSV parameter files, output files) use lower_snake_case.

```
R/
  config.R
  parameter_loader.R
  therapy_generator.R

inst/example_params/mCRC/
  cohort_population.csv
  drug_feature_catalog.csv
  regimen_drug_map.csv

output/
  drug_exposure.csv
  ground_truth_lines.csv
```

### CSV column names: `lower_snake_case`

Parameter CSV files use lower_snake_case for all column headers. This keeps them consistent with database conventions and makes them easy to load directly into databases.

```csv
cohort_id,age_group,sex,obs_length_bin,proportion
mCRC,50_64,M,365_730,0.20
```

## Summary table

| Thing | Convention | Example |
|-------|-----------|---------|
| R function names | camelCase | `generatePatients()` |
| R variable names | camelCase | `nPatients` |
| R list element names (config) | camelCase | `config$paramDir` |
| R list element names (data) | lower_snake_case | `pack$regimenCatalog` (the R list slot is camelCase, but the data.table inside has lower_snake_case columns) |
| data.table column names | lower_snake_case | `dt$drug_feature_id` |
| Database table names | lower_snake_case | `drug_exposure` |
| CSV column headers | lower_snake_case | `cohort_id,regimen_id` |
| R source file names | lower_snake_case | `parameter_loader.R` |
| CSV file names | lower_snake_case | `regimen_catalog.csv` |
| Output file names | lower_snake_case | `ground_truth_lines.csv` |

## Other conventions

- **Indentation**: 2 spaces (no tabs)
- **Line length**: aim for under 100 characters
- **Roxygen**: use `#'` for exported function documentation
- **Comments**: use `#` for inline comments explaining non-obvious design decisions
- **data.table**: preferred over base data.frame for performance. Use `:=` for in-place modification, `.(...)` for grouping
- **Error handling**: fail loudly with `stop(..., call. = FALSE)` on malformed input. Use `message()` for progress updates, not `cat()` or `print()`
- **Type coercion**: be explicit. Use `as.integer()`, `as.numeric()`, `as.character()` rather than relying on implicit coercion
- **Seeds**: use `set.seed()` for reproducibility. Offset seeds for different generation phases to avoid correlation
