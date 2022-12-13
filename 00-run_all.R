# Remove cache if you prefer to recalculate all the alignments
# unlink("tmp", recursive=T, force=T)

source("01-replication_fork_nn.R")
tzNN_prepare_training_data()
tzNN_train()
tzNN_evaluate()
gc()

source("02-detect_offtargets.R")
detect_offtargets()
gc()

source("03-detect_rdc.R")
detect_rdc()
gc()

source("04-rdc_published_overlap.R")
rdc_published_overlap()
gc()

source("05-APH_concentration.R")
APH_concentration()
gc()

source("06-promoter_enhancer_deletion.R")
promoter_enhancer_deletion()
gc()

source("07-multiomic_examples.R")
multiomics_examples()
gc()

source("08-rdc_pileup.R")
rdc_pileup()
gc()

source("09-replication_fork_length.R")
replication_fork_length()
gc()
