library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(grid)
library(xlsx)
devtools::load_all("breaktools/")
source("00-utils.R")


export_files_to_ncbi = function()
{
  #
  # Load samples data
  #
  samples_df = tlx_read_paper_samples("data/htgts_samples.tsv", "data") %>%
    dplyr::filter(celltype=="NPC" & organism=="mouse" & sample!="VI035" & (experiment=="APH concentration" | grepl(".*\\((22|22/37|22/5|47/5|18/4|38/3)\\)", group)))
    # dplyr::filter(celltype=="NPC" & organism=="mouse" & sample!="VI035" & grepl("(Ctnna2|Nrxn1) promoter/enhancer|APH concentration$", experiment) | grepl(".*\\((22|22/37|22/5|47/5|18/4|38/3)\\)", group))
  #
  # Load information about sequencer for each run
  #
  batch_sequencer_kits_df = readr::read_tsv("data/batch_sequencer_kits.tsv")

  #
  # Find al the raw files
  #
  raw_df = samples_df %>%
    dplyr::group_by(run) %>%
    dplyr::do((function(z){
      zz<<-z
      # z = samples_df %>% dplyr::filter(run=="B400_015")
      run_path = paste0("/omics/groups/OE0574/internal/peggy/HTGTS/", z$run[1])
      metadata_path = Sys.glob(paste0(run_path, "/", z$run[1], "*metadata*"))
      if(any(file.exists(metadata_path))) {
        metadata_path = metadata_path[order(nchar(basename(metadata_path)))][1]
        metadata_df = z %>%
          dplyr::left_join(readr::read_tsv(metadata_path), by=c("sample"="Library")) %>%
          dplyr::mutate(metadata_exists=T) %>%
          dplyr::mutate(raw_path=paste0(run_path, "/preprocess/", sample, "_", run, "_R1.fq.gz"), raw_path_exists=file.exists(raw_path)) %>%
          dplyr::mutate(tlx_path=paste0(run_path, "/results/", sample, "_", run, "/", sample, "_", run, "_result.tlx"), tlx_path_exists=file.exists(tlx_path)) %>%
          dplyr::select_at(c(colnames(z), "metadata_exists", "raw_path", "raw_path_exists", "tlx_path", "tlx_path_exists"))
      } else {
        metadata_df = z %>% dplyr::bind_cols(metadata_exists=F, raw_path=NA_character_, raw_path_exists=F, tlx_path=NA_character_, tlx_path_exists=F)
      }

      metadata_df
    })(.)) %>%
    dplyr::ungroup()

  # library name	title	organism	tissue	cell line	cell type	genotype	treatment	molecule	single or paired-end	instrument model	description	processed data file 	processed data file 	raw file	raw file	raw file	raw file
  geo_samples_df = samples_df %>%
    dplyr::left_join(raw_df %>% dplyr::select(sample, raw_path), by="sample") %>%
    dplyr::left_join(batch_sequencer_kits_df, by="run") %>%
    dplyr::mutate(genotype=gsub("/", "-", ifelse(grepl("ESC-NPC;NXP[0-9]+;([0-9/-]+).*", description, ignore.case=T), gsub("ESC-NPC;NXP[0-9]+;([0-9/-]+).*", "\\1", description, ignore.case=T), ""))) %>%
    dplyr::group_by(organism, celltype, bait_name, treatment, genotype) %>%
    dplyr::mutate(replicate=paste0("rep", 1:dplyr::n())) %>%
    data.frame() %>%
    dplyr::mutate(
      `library name`=paste0("HTGTS_Sample", 1:dplyr::n()),
      sample_id=gsub("_+", "_", gsub(" ", "-", paste0(bait_name, "_", gsub("[ .]", "", gsub(" uM [0-9]+h$", "", treatment)), "_", `genotype`, "_", replicate))),
      `title`=gsub("; *", "; ", paste(description, replicate, sep="; ")),
      `organism`=organism,
      `cell line`=toupper(gsub("ESC-NPC;([^;]+);.*", "\\1", description)),
      `cell type`=gsub(";.*", "", description),
      `treatment`=treatment,
      `molecule`="genomic DNA",
      `single or paired-end`="paired-end",
      `instrument model`=sequencer_kit,
      description="",
      raw_extention=gsub("[^.]+(\\..*)$", "\\1", basename(raw_path)),
      raw_local_path1=raw_path,
      raw_local_path2=gsub("_R1", "_R2", raw_local_path1),
      raw_file1=paste0(sample_id, "_R1", raw_extention),
      raw_file2=paste0(sample_id, "_R2", raw_extention),
      tlx_local_path=path,
      tlx_file=paste0(sample_id, ".tlx")
    )


  #
  # Do a sanity check
  #
  if(T)
  {
    error = ""
    if(any(!complete.cases(geo_samples_df))) {
      current_error = paste0("NA values present in sample ids: ", paste(geo_samples_df$sample_id[!complete.cases(geo_samples_df)], collapse=", "))
      error = paste0(error, "\n", current_error)
    }

    if(any(duplicated(geo_samples_df$sample_id))) {
      current_error = paste0("Duplicate sample ids: ", paste(geo_samples_df$sample_id[duplicated(geo_samples_df$sample_id)], collapse=", "))
      error = paste0(error, "\n", current_error)
    }

    raw_local_files = c(geo_samples_df$raw_local_path1, geo_samples_df$raw_local_path2)
    if(any(duplicated(raw_local_files))) {
      current_error = paste0("Duplicate raw_local_path1 or raw_local_path2: ", paste(raw_local_files[duplicated(raw_local_files)], collapse=", "))
      error = paste0(error, "\n", unique(current_error))
    }

    file_sizes_df = dplyr::bind_rows(file.info(c(geo_samples_df$raw_local_path1, geo_samples_df$raw_local_path2, geo_samples_df$tlx_local_path))) %>%
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

    if(any(duplicated(geo_samples_df$raw_local_path2))) {
      current_error = paste0("Duplicate raw_local_path2: ", paste(geo_samples_df$raw_local_path2[duplicated(geo_samples_df$raw_local_path2)], collapse=", "))
      error = paste0(error, "\n", current_error)
    }

    if(any(duplicated(geo_samples_df$tlx_local_path))) {
      current_error = paste0("Duplicate tlx_local_path: ", paste(geo_samples_df$tlx_local_path[duplicated(geo_samples_df$tlx_local_path)], collapse=", "))
      error = paste0(error, "\n", current_error)
    }

    if(error!="") {
      stop(error)
    }
  }


  geo_samples_export_df = geo_samples_df %>%
    dplyr::select(`library name`, `title`, `organism`, `cell line`, `cell type`, `genotype`, `treatment`, `molecule`, `single or paired-end`, `instrument model`, `description`, tlx_file, raw_file1, raw_file2) %>%
    magrittr::set_names(c("library name", "title", "organism", "cell line", "cell type", "genotype", "treatment", "molecule", "single or paired-end", "instrument model", "description", "processed data file", "raw file", "raw file"))

  geo_paired_export_df = geo_samples_df %>%
    dplyr::select(`file name 1`=raw_file1, `file name 2`=raw_file2)

  geo_md5_df=dplyr::bind_rows(
    geo_samples_df %>% dplyr::select(local_path=raw_local_path1, file=raw_file1),
    geo_samples_df %>% dplyr::select(local_path=raw_local_path2, file=raw_file2),
    geo_samples_df %>% dplyr::select(local_path=tlx_local_path, file=tlx_file)) %>%
    dplyr::mutate(md5=tools::md5sum(local_path)) %>%
    dplyr::select(`file name`=file, `file checksum`=md5)


  excel_path = "reports/00-upload_ncbi/ncbi.xlsx"
  dir.create(dirname(excel_path), recursive=T, showWarnings=F)
  xlsx::write.xlsx(geo_samples_export_df, excel_path, sheetName="2-1. Metadata Template",  col.names=T, row.names=F, append=F)
  xlsx::write.xlsx(geo_paired_export_df, excel_path, sheetName="2-2. PAIRED-END EXPERIMENTS",  col.names=T, row.names=F, append=T)
  xlsx::write.xlsx(geo_md5_df, excel_path, sheetName="3. MD5 Checksums",  col.names=T, row.names=F, append=T)
}