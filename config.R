## ============================================================================
## config.R  ——  the only file that normally needs editing
## ============================================================================
## Data root: the directory that holds the raw/intermediate inputs described in
## README.md  (section "Data layout"). Outputs are written under ./out, ./tables
## and ./figures in the working directory.
##
## Usage:
##   git clone <repo> && cd <repo>
##   export MANUSCRIPT_DATA=/path/to/data          (Windows: set MANUSCRIPT_DATA=D:\path)
##   Rscript 02_machine_learning/10_ml130.R
## ============================================================================

DATA_ROOT <- Sys.getenv("MANUSCRIPT_DATA", unset = file.path(getwd(), "data"))
WORK      <- Sys.getenv("MANUSCRIPT_WORK", unset = getwd())

## The machine-learning scripts live in 02_machine_learning/ and are sourced by
## scripts that call source(file.path(SCRIPT_DIR, "...")).
SCRIPT_DIR <- file.path(WORK, "02_machine_learning")

## The transcriptome scripts use ROOT as a synonym for the data root.
ROOT <- DATA_ROOT

for (d in c("out", "tables", "figures")) {
  dir.create(file.path(WORK, d), showWarnings = FALSE, recursive = TRUE)
}

message("DATA_ROOT  = ", DATA_ROOT)
message("WORK       = ", WORK)
message("SCRIPT_DIR = ", SCRIPT_DIR)
