# Remove cache if you prefer to recalculate all the alignments
# unlink("tmp", recursive=T, force=T)

source("02-detect_offtargets.R")
detect_offtargets()

source("03-detect_rdc.R")
detect_rdc()

source("04-rdc_published_overlap.R")
rdc_published_overlap()

source("05-APH_concentration.R")
APH_concentration()

source("06-promoter_enhancer_deletion.R")
promoter_enhancer_deletion()

source("07-rdc_pileup.R")
rdc_pileup()

source("08-multiomic_examples.R")
multiomics_examples()

source("09-replication_fork_length.R")
replication_fork_length()