# Examples

Each folder contains 3 R scripts:

- `clean.R` reads in raw data and outputs two cleaned datasets (`survey_df` and `poststrat_df`)
- `analysis.R` specifies and fits the model using `brms`, and extracts the implicit weights using `mrpaw`
- `diagnostics.R` uses the implicit weights to perform diagnostic checks

The raw data used in `clean.R`, and the output data from `clean.R` and `analysis.R` are stored in [this Google Drive folder](https://drive.google.com/drive/folders/1SKrSM-JWz9MJDjdQEtxR6pxEJ1Q3c4CW).
