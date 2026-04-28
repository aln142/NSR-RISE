# HEAL EEG Biomarkers

This repository contains an R Markdown workflow for biomarker analyses related to acquired epilepsy after acute provoked neonatal seizures.

## Files

- `heal_eeg_biomarkers_github.Rmd`: sanitized R Markdown analysis template.
- `data/`: local input data directory. Raw or identifiable data should not be committed.
- `outputs/`: generated figures and tables.

## Setup

Install required R packages listed in the setup chunk of the R Markdown file.

Create a local `.Renviron` file containing your REDCap credentials:

```bash
REDCAP_URI="https://redcap.ucsf.edu/api/"
REDCAP_TOKEN="your_redcap_token_here"
```

Do not commit `.Renviron`, REDCap tokens, raw data, or identifiable participant information to GitHub.

## Running the analysis

Open `heal_eeg_biomarkers_github.Rmd` in RStudio and knit to HTML. Several code chunks are set to `eval=FALSE` by default because they require local data files or REDCap access. Change these to `eval=TRUE` after configuring local paths and credentials.

## Data privacy

This public template does not include raw data. De-identified datasets should be stored locally or shared through approved institutional mechanisms.
