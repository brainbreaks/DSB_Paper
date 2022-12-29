library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(grid)
library(xlsx)
devtools::load_all("breaktools/")
source("00-utils.R")

export_htgts_to_ncbi = function()
{
  dir.create("reports/00-upload_ncbi", recursive=T, showWarnings=F)

  #
  # Load samples data
  #
  htgts_samples_df = tlx_read_samples("data/htgts_samples.tsv", "data/TLX_paper") %>%
    dplyr::filter(!grepl("Wei|Tena", experiment))

  #
  # Find al the raw files
  #
  htgts_raw_df = htgts_samples_df %>%
    dplyr::group_by(run) %>%
    dplyr::do((function(z){
      zz<<-z
      run_path = paste0("/omics/groups/OE0574/internal/peggy/HTGTS/", z$run[1])
      metadata_path = Sys.glob(paste0(run_path, "/", z$run[1], "*metadata*"))
      if(any(file.exists(metadata_path))) {
        metadata_path = metadata_path[order(nchar(basename(metadata_path)))][1]
        metadata_df = z %>%
          dplyr::left_join(readr::read_tsv(metadata_path), by=c("sample"="Library")) %>%
          dplyr::mutate(htgts_parameters=paste0("[", tidyr::replace_na(MID, ""), "]", Chr, ":", Start, "-", End, ":", Strand)) %>%
          dplyr::mutate(htgts_barcode=tidyr::replace_na(MID, "")) %>%
          dplyr::mutate(metadata_exists=T) %>%
          dplyr::mutate(raw_path=paste0(run_path, "/preprocess/", sample, "_", run, "_R1.fq.gz"), raw_path_exists=file.exists(raw_path)) %>%
          dplyr::mutate(tlx_path=paste0(run_path, "/results/", sample, "_", run, "/", sample, "_", run, "_result.tlx"), tlx_path_exists=file.exists(tlx_path)) %>%
          dplyr::select_at(c(colnames(z), "htgts_parameters", "htgts_barcode", "metadata_exists", "raw_path", "raw_path_exists", "tlx_path", "tlx_path_exists"))
      } else {
        metadata_df = z %>% dplyr::bind_cols(metadata_exists=F, raw_path=NA_character_, raw_path_exists=F, tlx_path=NA_character_, tlx_path_exists=F)
      }

      metadata_df
    })(.)) %>%
    dplyr::ungroup()

  geo_htgts_samples_df = htgts_samples_df %>%
    dplyr::left_join(htgts_raw_df %>% dplyr::select(sample, raw_path, htgts_barcode), by="sample") %>%
    dplyr::mutate(cellparental=toupper(gsub(".*(NXP[0-9]+).*", "\\1", description, ignore.case=T))) %>%
    dplyr::mutate(cellculture=gsub("/", "-", ifelse(grepl("ESC-NPC;NXP[0-9]+;([0-9/-]+).*", description, ignore.case=T), gsub("ESC-NPC;NXP[0-9]+;([0-9/-]+).*", "\\1", description, ignore.case=T), ""))) %>%
    dplyr::group_by(organism, celltype, bait_name, treatment, cellculture) %>%
    dplyr::mutate(replicate=paste0("replicate_", 1:dplyr::n())) %>%
    data.frame() %>%
    dplyr::mutate(
      study_geo=tidyr::replace_na(study_geo, "this study"),
      `library name`=paste0("HTGTS Sample ", 1:dplyr::n()),
      sample_id=gsub("_+", "_", gsub(" ", "-", paste0("HTGTS_", bait_name, "_", gsub("[ .]", "", gsub(" uM [0-9]+h$", "", treatment)), "_", cellculture, "_", replicate, "_", sample))),
      `title`=paste0(gsub("; *", "; ", paste(description, "LAM-HTGTS", replicate, sep="; ")), "; Bait ", bait_name, "; Barcode ", htgts_barcode),
      `genotype`=gsub(",$", "", paste("Xrcc4-/-::p53-/-", cellparental, gsub(" \\(.*", "", gsub("\\b(APH|DMSO).*", "", group)), sep=",")),
      `organism`=organism,
      `cell line`="ESC-NPC",
      `cell type`="Neural progenitor cells",
      `treatment`=treatment,
      `molecule`="genomic DNA",
      `single or paired-end`="paired-end",
      `instrument model`=instrument_model,
      description="LAM-HTGTS",
      raw_extention=gsub("[^.]+(\\..*)$", "\\1", basename(raw_path)),
      raw_local_path1=raw_path,
      raw_local_path2=gsub("_R1", "_R2", raw_local_path1),
      raw_file1=paste0(sample_id, "_R1", raw_extention),
      raw_file2=paste0(sample_id, "_R2", raw_extention),
      tlx_local_path=path,
      tlx_clean_local_path=paste0("reports/00-upload_ncbi/TLX/", sample_id, ".tlx"),
      tlx_file=paste0(sample_id, ".tlx")
    )

  #
  # Remove Psyx from TLX files
  #
  geo_htgts_samples_df %>%
    dplyr::group_by(tlx_local_path) %>%
    dplyr::summarize(return=system(paste0("sed '/phiX/d' ", tlx_local_path, "  > ", tlx_clean_local_path)))

  #
  # Do a sanity check
  #
  if(T)
  {
    error = ""
    if(any(!complete.cases(geo_htgts_samples_df))) {
      current_error = paste0("NA values present in sample ids: ", paste(geo_htgts_samples_df$sample_id[!complete.cases(geo_htgts_samples_df)], collapse=", "))
      error = paste0(error, "\n", current_error)
    }

    if(any(duplicated(geo_htgts_samples_df$sample_id))) {
      current_error = paste0("Duplicate sample ids: ", paste(geo_htgts_samples_df$sample_id[duplicated(geo_htgts_samples_df$sample_id)], collapse=", "))
      error = paste0(error, "\n", current_error)
    }

    raw_local_files = c(geo_htgts_samples_df$raw_local_path1, geo_htgts_samples_df$raw_local_path2)
    if(any(duplicated(raw_local_files))) {
      current_error = paste0("Duplicate raw_local_path1 or raw_local_path2: ", paste(raw_local_files[duplicated(raw_local_files)], collapse=", "))
      error = paste0(error, "\n", unique(current_error))
    }

    file_sizes_df = dplyr::bind_rows(file.info(c(geo_htgts_samples_df$raw_local_path1, geo_htgts_samples_df$raw_local_path2, geo_htgts_samples_df$tlx_clean_local_path))) %>%
      tibble::rownames_to_column("file_path") %>%
      dplyr::mutate(filetype=gsub(".*\\.(.*)$", "\\1", as.character(file_path))) %>%
      dplyr::mutate(size_validation=dplyr::case_when(
        filetype=="tlx" & size<500*1e3 ~ "Small TLX", # 500Kb for TLX files
        filetype=="gz" & size<5*1e6 ~ "Small FASTQ", # 5Mb for FASTQ files
        filetype!="gz" & filetype!="tlx" ~ "Unknown filetype",
        T~""))
    if(any(file_sizes_df$size_validation!="")) {
      current_error = paste0("Error validating file size of(raw_local_path1, raw_local_path2, or tlx_path): ", paste(file_sizes_df$file_path[file_sizes_df$size_validation!=""], collapse=", "))
      error = paste0(error, "\n", current_error)
    }

    if(any(duplicated(geo_htgts_samples_df$raw_local_path2))) {
      current_error = paste0("Duplicate raw_local_path2: ", paste(geo_htgts_samples_df$raw_local_path2[duplicated(geo_htgts_samples_df$raw_local_path2)], collapse=", "))
      error = paste0(error, "\n", current_error)
    }

    if(any(duplicated(geo_htgts_samples_df$tlx_clean_local_path))) {
      current_error = paste0("Duplicate tlx_clean_local_path: ", paste(geo_htgts_samples_df$tlx_clean_local_path[duplicated(geo_htgts_samples_df$tlx_clean_local_path)], collapse=", "))
      error = paste0(error, "\n", current_error)
    }

    if(error!="") {
      stop(error)
    }
  }

  geo_htgts_samples_export_df = geo_htgts_samples_df %>%
    dplyr::select(`library name`, `title`, `organism`, `cell line`, `cell type`, `genotype`, `treatment`, `molecule`, `single or paired-end`, `instrument model`, `description`, tlx_file, raw_file1, raw_file2) %>%
    magrittr::set_names(c("library name", "title", "organism", "cell line", "cell type", "genotype", "treatment", "molecule", "single or paired-end", "instrument model", "description", "processed data file", "raw file", "raw file"))

  geo_paired_export_df = geo_htgts_samples_df %>%
    dplyr::select(`file name 1`=raw_file1, `file name 2`=raw_file2)

  geo_md5_df=dplyr::bind_rows(
    geo_htgts_samples_df %>% dplyr::select(local_path=raw_local_path1, file=raw_file1),
    geo_htgts_samples_df %>% dplyr::select(local_path=raw_local_path2, file=raw_file2),
    geo_htgts_samples_df %>% dplyr::select(local_path=tlx_clean_local_path, file=tlx_file)) %>%
    dplyr::mutate(md5=tools::md5sum(local_path)) %>%
    dplyr::select(`file name`=file, `file checksum`=md5)

  #
  # Create script to copy all ftp files
  #
  cmd_ftpcopy = c(
    "set log:file/xfer reports/00-upload_ncbi/upload.log",
    "set xfer:log 1",
    "set xfer:eta-period 60",
    "set xfer:rate-period 60",
    "set ftp:proxy http://www-int2.inet.dkfz-heidelberg.de:80",
    'open -u "geoftp,Haf5Fatoryen" ftp://ftp-private.ncbi.nlm.nih.gov',
    "bin",
    # paste("put -e -O", "uploads/ch164594_Y54fyeWG", geo_htgts_samples_df$tlx_clean_local_path),
    "",
    paste("put -e -O", "uploads/ch164594_Y54fyeWG", c(geo_htgts_samples_df$raw_local_path1, geo_htgts_samples_df$raw_local_path2), "-o", basename(c(geo_htgts_samples_df$raw_file1, geo_htgts_samples_df$raw_file2))),
    "",
    "bye"
  )
  writeLines(cmd_ftpcopy, con="reports/00-upload_ncbi/upload.ftp")
  writeLines("#!/bin/bash\nlftp -f reports/00-upload_ncbi/upload.ftp", con="reports/00-upload_ncbi/upload.sh")

  #
  # Write EXCEL file describing samples uploaded to GEO
  #
  geo_excel_path = "reports/00-upload_ncbi/geo_htgts.xlsx"
  xlsx::write.xlsx(geo_htgts_samples_export_df, geo_excel_path, sheetName="2-1. Metadata Template",  col.names=T, row.names=F, append=F)
  xlsx::write.xlsx(geo_paired_export_df, geo_excel_path, sheetName="2-2. PAIRED-END EXPERIMENTS",  col.names=T, row.names=F, append=T)
  xlsx::write.xlsx(geo_md5_df, geo_excel_path, sheetName="3. MD5 Checksums",  col.names=T, row.names=F, append=T)

  #
  # LAM-HTGTS SRA data
  #
  sra_htgts_samples_export_df = geo_htgts_samples_df %>%
    dplyr::mutate(
      library_strategy="OTHER",
      library_selection="OTHER",
      library_source="GENOMIC",
      library_layout="PAIRED",
      platform="ILLUMINA",
      instrument_model=gsub(" ?V\\d.*", "", `instrument model`),
      dev_stage="Neural progenitor cells",
      tissue="Embrionic stem cells",
      design_description=description,
      cell_line=paste0(cellparental, ifelse(cellculture!="", " ", ""), cellculture),
      raw_path=raw_local_path1, raw_path1=raw_local_path2,
      filename=paste0(sample_id, "_R1.fastq.gz"), filename1=paste0(sample_id, "_R2.fastq.gz") ) %>%
    dplyr::select(sample_name=sample_id, sample_title=title, organism, strain=genotype, cell_line, bait_name, dev_stage, tissue, treatment, bait_name, replicate, dplyr::matches("library_"), platform, instrument_model, design_description, raw_path, raw_path1, filename, filename1)

  #
  # GRO-seq SRA data
  #
  groseq_fastq_paths_df = data.frame(raw_path=Sys.glob("/omics/groups/OE0574/internal/peggy/groseq/*/*/*R1.fastq.gz")) %>%
    dplyr::mutate(ilse_number=gsub("-LR-.*", "", basename(raw_path)))
  groseq_samples_df = readr::read_tsv("data/groseq_samples.tsv", na="<N/A>") %>%
    dplyr::filter(published=="Y") %>%
    dplyr::left_join(groseq_fastq_paths_df, by="ilse_number")
  sra_groseq_samples_export_df = groseq_samples_df %>%
    dplyr::mutate(
      sample_name=paste0(sample_id, "_", sample),
      library_strategy="OTHER",
      library_selection="OTHER",
      library_source="TRANSCRIPTOMIC",
      library_layout="SINGLE",
      bait_name="no bait",
      cell_line=paste0(cell_parental, ifelse(cellculture!="", " ", ""), cellculture),
      replicate=gsub(".*(replicate_\\d)$", "\\1", sample_id),
      sample_title=paste0(gsub(".*,(ESC-NPC)", "\\1", genotype), " (", cellculture, "); ", replicate),
      dev_stage="Neural progenitor cells", tissue="Embrionic stem cells", design_description="GRO-seq",
      raw_path1="",
      filename=paste0(sample_name, ".fastq.gz"), filename1=""
    ) %>%
    dplyr::select(sample_name, sample_title, organism, strain=genotype, cell_line, bait_name, dev_stage, tissue, treatment, replicate, dplyr::matches("library_"), platform, instrument_model, design_description, raw_path, raw_path1, filename, filename1)

  #
  # DRIP-seq SRA data
  #
  dripseq_fastq_paths_df = data.frame(raw_path=Sys.glob("/omics/groups/OE0574/internal/peggy/dripseq/ILSE24885_120422/0.fastq/*R1.fastq.gz")) %>%
    dplyr::mutate(ilse_number=gsub("-LR-.*", "", basename(raw_path)))
  dripseq_samples_df = readr::read_tsv("data/dripseq_samples.tsv", na="<N/A>") %>%
    dplyr::filter(published=="Y") %>%
    dplyr::left_join(dripseq_fastq_paths_df, by="ilse_number")
  sra_dripseq_samples_export_df = dripseq_samples_df %>%
    dplyr::mutate(
      sample_name=paste0(sample_id, "_", sample),
      bait_name="no bait",
      cell_line=cell_parental,
      library_strategy="OTHER",
      library_selection="Restriction Digest",
      library_source="GENOMIC",
      library_layout="SINGLE",
      sample_title=paste0(gsub(".*,(ESC-NPC)", "\\1", genotype), "; ", treatment, "; replicate_", gsub(".*(\\d)$", "\\1", sample_id)),
      dev_stage="Neural progenitor cells", tissue="Embrionic stem cells", design_description="DRIP-seq",
      raw_path1="",
      filename=paste0(sample_name, ".fastq.gz"), filename1="",
    ) %>%
    dplyr::select(sample_name, sample_title, organism=organism, strain=genotype, cell_line, bait_name, dev_stage, tissue, treatment, replicate, dplyr::matches("library_"), platform, instrument_model, design_description, raw_path, raw_path1, filename, filename1)


  #
  # Create directory with links to SRA files
  #
  dir.create("reports/00-upload_ncbi/SRA", showWarnings=F, recursive=T)
  sra_files_df = dplyr::bind_rows(
    sra_htgts_samples_export_df %>% dplyr::mutate(pair_suffix="_R1") %>% dplyr::select(sample_name, pair_suffix, raw_path, filename),
    sra_htgts_samples_export_df %>% dplyr::mutate(pair_suffix="_R2") %>% dplyr::select(sample_name, pair_suffix, raw_path=raw_path1, filename=filename1),
    sra_groseq_samples_export_df %>% dplyr::mutate(pair_suffix="") %>% dplyr::select(sample_name, pair_suffix, raw_path, filename),
    sra_dripseq_samples_export_df %>% dplyr::mutate(pair_suffix="") %>% dplyr::select(sample_name, pair_suffix, raw_path, filename)) %>%
    dplyr::mutate(link=paste0("reports/00-upload_ncbi/SRA/", sample_name, pair_suffix, ".fastq.gz")) %>%
    dplyr::group_by(sample_name, pair_suffix, link, raw_path) %>%
    dplyr::summarize(return=R.utils::createLink(link=link, target=raw_path, overwrite=T), file.info(link)) %>%
    dplyr::ungroup()

  #
  # Copy fastq.gz files to SRA
  #
  dir.create("reports/00-upload_ncbi/SRA_COPY", showWarnings=F, recursive=T)
  sra_files_df = dplyr::bind_rows(
    sra_htgts_samples_export_df %>% dplyr::mutate(pair_suffix="_R1") %>% dplyr::select(sample_name, pair_suffix, raw_path, filename),
    sra_htgts_samples_export_df %>% dplyr::mutate(pair_suffix="_R2") %>% dplyr::select(sample_name, pair_suffix, raw_path=raw_path1, filename=filename1),
    sra_groseq_samples_export_df %>% dplyr::mutate(pair_suffix="") %>% dplyr::select(sample_name, pair_suffix, raw_path, filename),
    sra_dripseq_samples_export_df %>% dplyr::mutate(pair_suffix="") %>% dplyr::select(sample_name, pair_suffix, raw_path, filename)) %>%
    dplyr::mutate(link=paste0("reports/00-upload_ncbi/SRA_COPY/", sample_name, pair_suffix, ".fastq.gz")) %>%
    dplyr::group_by(sample_name, pair_suffix, link, raw_path) %>%
    dplyr::summarize(return=R.utils::copyFile(srcPathname=raw_path, destPathname=link, overwrite=T), file.info(link)) %>%
    # dplyr::summarize(return=R.utils::createLink(link=link, target=raw_path, overwrite=T), file.info(link)) %>%
    dplyr::ungroup()


  sra_small_samples = sra_files_df %>% dplyr::filter(size<=5e6)
  if(nrow(sra_small_samples)>0) {
    stop("Some samples are too small (<5Mb): ", paste(sra_small_samples$sample_name, collapse="\n"))
  }

  #
  # Write TSV file describing samples uploaded to SRA
  #
  sra_samples_path = "reports/00-upload_ncbi/SRA_samples.tsv"
  dplyr::bind_rows(sra_htgts_samples_export_df, sra_groseq_samples_export_df, sra_dripseq_samples_export_df) %>%
    dplyr::mutate(sex='not collected') %>%
    dplyr::select(sample_name, sample_title, organism, strain, cell_line, bait_name, dev_stage, tissue, sex, treatment, replicate) %>%
    readr::write_tsv(file=sra_samples_path)

  #
  # Write TSV file files metadata uploaded to SRA
  #
  sra_metadata_path = "reports/00-upload_ncbi/SRA_metadata.tsv"
  dplyr::bind_rows(sra_htgts_samples_export_df, sra_groseq_samples_export_df, sra_dripseq_samples_export_df) %>%
    dplyr::mutate(library_ID=sample_name, title=paste0(design_description, " of ", organism, ": ", sample_title), filetype="fastq") %>%
    dplyr::select(sample_name, library_ID, title, library_strategy, library_source, library_selection, library_layout, platform, instrument_model, design_description, filetype, filename, filename1) %>%
    readr::write_tsv(file=sra_metadata_path)

  # xlsx::write.xlsx(sra_samples_export_df, sra_excel_path, sheetName="Model.organism.animal.1.0",  col.names=T, row.names=F, append=F)
}