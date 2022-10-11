library(rtracklayer)
library(readr)
library(dplyr)
library(cowplot)
library(stringr)
library(scico)
library(dplyr)
library(ggplot2)
library(tensorflow)
library(reticulate)
library(keras)
library(abind)
devtools::load_all('breaktools/')
source("tzNN_utils.R")
source("00-utils.R")

layer_activate = function(input, name, activation="selu", dropout=0, normalize=T) {
    output = input %>%
      keras::layer_batch_normalization(name=paste0(name, "_batchnorm")) %>%
      keras::layer_activation(activation=activation, name=paste0(name, "_", activation))
    if(normalize) return(output %>% keras::layer_spatial_dropout_2d(rate=dropout))
    ouptut
}

loss_dice = function(y_true, y_pred, smooth = 1.0) {
  y_true_f = keras::k_flatten(y_true)
  y_pred_f = keras::k_flatten(y_pred)
  intersection = keras::k_sum(y_true_f * y_pred_f)
  -(2 * intersection + smooth) / (keras::k_sum(y_true_f) + keras::k_sum(y_pred_f) + smooth)
}

pretrube_mirror = function(m) {
    array(m[,dim(m)[2]:1,], dim=dim(m))
}

metric_dice = keras::custom_metric("dice", function(y_true, y_pred, smooth = 1.0) -loss_dice(y_true, y_pred, smooth))
metric_intersection = keras::custom_metric("intersection", function(y_true, y_pred) keras::k_sum(keras::k_flatten(y_true) * keras::k_flatten(y_pred)))
metric_junction = keras::custom_metric("junction", function(y_true, y_pred) keras::k_sum(keras::k_flatten(y_true)) + keras::k_sum(keras::k_flatten(y_pred)))

get_deep = function(input_shape=c(20, 16, 1), dropout1=0.8, dropout2=0.2) {
    inputs = keras::layer_input(shape = input_shape) %>%
      keras::layer_batch_normalization(name="initial_batchnorm")


    layer01 = inputs %>%
      keras::layer_conv_2d(filters=32, kernel_size=c(3, 3), padding="valid", name="layer01_conv1") %>%
      layer_activate(name="layer01_act1")
    layer1pool = layer01 %>%
      keras::layer_average_pooling_2d(pool_size=c(5, 5), stride=c(1, 1), padding="valid", name="layer1pool_pool")
    layer10 = layer01 %>%
      keras::layer_conv_2d(filters=4, kernel_size=c(1, 1), padding="valid", name="layer10_conv") %>%
      layer_activate(name="layer10_act1", dropout=dropout2)

    layer11 = layer01 %>%
      keras::layer_conv_2d(filters=8, kernel_size=c(1, 1), padding="valid", name="layer11_conv1") %>%
      layer_activate(name="layer11_act1", dropout=dropout2) %>%
      keras::layer_depthwise_conv_2d(kernel_size=c(3, 3), padding="valid", name="layer11_dv_conv1") %>%
      keras::layer_conv_2d(filters=8, kernel_size=c(1, 1), padding="valid", name="layer11_conv2") %>%
      # keras::layer_conv_2d(filters=8, kernel_size=c(3, 3), padding="valid") %>%
      layer_activate(name="layer11_act2")

    layer12 = layer01 %>%
      keras::layer_conv_2d(filters=8, kernel_size=c(1, 1), padding="valid", name="layer12_conv1") %>%
      layer_activate(name="layer12_act1", dropout=dropout2) %>%
      keras::layer_depthwise_conv_2d(kernel_size=c(5, 5), padding="valid", name="layer12_dv_conv1") %>%
      keras::layer_conv_2d(filters=8, kernel_size=c(1, 1), padding="valid", name="layer12_conv2") %>%
      # keras::layer_conv_2d(filters=8, kernel_size=c(5, 5), padding="valid") %>%
      layer_activate(name="layer12_act2")

    layer13 = layer01 %>%
      keras::layer_conv_2d(filters=8, kernel_size=c(1, 1), padding="valid", name="layer13_conv1")  %>%
      layer_activate(name="layer13_act1", dropout=dropout2)%>%
      keras::layer_depthwise_conv_2d(kernel_size=c(7, 7), padding="valid", name="layer13_dv_conv1") %>%
      keras::layer_conv_2d(filters=8, kernel_size=c(1, 1), padding="valid", name="layer13_conv2") %>%
      # keras::layer_conv_2d(filters=8, kernel_size=c(7, 7), padding="valid") %>%
      layer_activate(name="layer13_act2")

    inputs_crop = inputs %>%
      keras::layer_conv_2d(filters=1, kernel_size=c(1, 3), padding="valid", name="inputs_conv1") %>%
      keras::layer_cropping_2d(cropping=c(4,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,0))
    layer10_crop = layer10 %>%
      keras::layer_cropping_2d(cropping=c(3,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,0))

    layer11_crop = layer11 %>%
      keras::layer_cropping_2d(cropping=c(2,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,1))
    layer12_crop = layer12 %>%
      keras::layer_cropping_2d(cropping=c(1,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,2))
    layer13_crop = layer13 %>%
      keras::layer_cropping_2d(cropping=c(0,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,3))
    layer1pool_crop = layer1pool %>%
      keras::layer_cropping_2d(cropping=c(1,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,2))


    concatenate1 = keras::layer_concatenate(list(layer11_crop, layer12_crop, layer13_crop), axis=3) %>%
      keras::layer_spatial_dropout_2d(rate = dropout1)
    classify = keras::layer_concatenate(list(inputs_crop, layer10_crop, layer1pool_crop, concatenate1), axis=3) %>%
      keras::layer_conv_2d(filters=1, kernel_size=c(1, 1), padding="valid", name="classify_conv1") %>%
      keras::layer_max_pooling_2d(pool_size=c(1, 14), stride=c(1, 1), padding="valid", name="classify_pool1") %>%
      layer_activate(name="classify_act1", activation="sigmoid") %>%
      keras::layer_zero_padding_2d(padding=c(4,0)) %>%
      keras::layer_conv_2d(filters=1, kernel_size=c(7, 1), padding="valid", name="classify_conv2") %>%
      layer_activate(name="classify_act2", activation="sigmoid") %>%
      keras::layer_zero_padding_2d(padding=c(3,0))

    model = keras::keras_model(inputs=inputs, outputs=classify)
    # model %>% keras::compile(optimizer=keras::optimizer_rmsprop(learning_rate=learning_rate), loss=keras::loss_binary_crossentropy)
    return(model)
}


tzNN_prepare_training_data = function()
{
    training_step = 5e5
    training_margin = 150e3
    training_width = 1e6
    training_binsize = 5e4

    #
    # Load data
    #
    repliseq1_df = readr::read_tsv("~/Workspace/Datasets/B400_RS_001/results/DMSO_repliseq50000.tsv") %>% dplyr::mutate(Run="RS1_mNPC8")
    repliseq2_df = readr::read_tsv("~/Workspace/Datasets/B400_RS_002/results/DMSO_repliseq50000.tsv") %>% dplyr::mutate(Run="RS2_mNPC8")
    repliseq3_df = readr::read_tsv("~/Workspace/Datasets/B400_RS_003/results/DMSO_repliseq50000.tsv") %>% dplyr::mutate(Run="RS3_mNPC8")
    repliseq_df = dplyr::bind_rows(
      readr::read_tsv("~/Workspace/Datasets/zhao_bmc_repliseq_2020/results/zhao_mESC_repliseq50000.tsv") %>% dplyr::mutate(Run="Zhao_mESC"),
      readr::read_tsv("~/Workspace/Datasets/zhao_bmc_repliseq_2020/results/zhao_mNPC_repliseq50000.tsv") %>% dplyr::mutate(Run="Zhao_mNPC")) %>%
      tidyr::crossing(reduction=1:2) %>% dplyr::mutate(Run=paste0(Run, 16/reduction)) %>%
      dplyr::mutate(repliseq_fraction=ceiling(repliseq_fraction/reduction)) %>%
      dplyr::group_by(Run, reduction, repliseq_chrom, repliseq_start, repliseq_end, repliseq_fraction) %>%
      dplyr::summarise(repliseq_value=sum(repliseq_value, na.rm=T), repliseq_value_abs=sum(repliseq_value_abs)/reduction[1]) %>%
      dplyr::group_by(Run, repliseq_chrom, repliseq_start, repliseq_end) %>%
      dplyr::mutate(repliseq_value=repliseq_value/max(repliseq_value)) %>%
      dplyr::ungroup()

    #
    # Combine all datasets and calculate average fraction for each locus
    #
    repliseq_compare_df = dplyr::bind_rows(repliseq_df, repliseq1_df, repliseq2_df, repliseq3_df)
    repliseqTime_compare_df = repliseq_compare_df %>%
      dplyr::group_by(Run, repliseq_chrom, repliseq_start, repliseq_end) %>%
      dplyr::summarise(repliseq_mean_fraction=weighted.mean(repliseq_fraction, tidyr::replace_na(repliseq_value, 0)), .groups="keep")

    #
    # Compare different Repli-Seq datasets
    #
    ggplot(repliseq_compare_df %>% dplyr::filter(repliseq_value_abs<=5000)) +
      ggridges::geom_density_ridges(aes(y=Run, x=repliseq_value_abs)) +
      labs(x="Number of reads per bin/fraction", y="")
    ggplot(repliseq_compare_df %>% dplyr::filter(repliseq_value<=1)) +
      ggridges::geom_density_ridges(aes(y=Run, x=repliseq_value)) +
      labs(x="Relative fraction abundance bin", y="")



    #
    # Load OK-seq and Amina's annotations
    #
    rfd_ranges = rtracklayer::import.bw("~/Workspace/Datasets/petryk_okseq_mESC/GSM3290342_OK-seq_smooth_results_w1000_s30_d30_z1_filter.bw") %>%
      as.data.frame() %>%
      dplyr::rename(rfd_chrom="seqnames", rfd_start="start", rfd_end="end", rfd_score="score") %>%
      dplyr::select(rfd_chrom, rfd_start, rfd_end, rfd_score) %>%
      df2ranges(rfd_chrom, rfd_start, rfd_end)
    iz_ranges = readr::read_tsv("~/Workspace/Datasets/petryk_okseq_mESC/GSM3290342_Ok_IZ.txt") %>%
      tidyr::extract(col=break_ID, into=c("iz_chrom", "iz_start", "iz_end"), regex="([^:]+):(\\d+)-(\\d+)") %>%
      dplyr::mutate(iz_start=as.numeric(iz_start), iz_end=as.numeric(iz_end)) %>%
      dplyr::select(iz_chrom, iz_start, iz_end) %>%
      df2ranges(iz_chrom, iz_start, iz_end)
    colfork_df = as.data.frame(rtracklayer::import.bed("data/zhao_colliding_forks.bed")) %>%
      dplyr::rename(colfork_chrom="seqnames", colfork_start="start", colfork_end="end", colfork_width="width") %>%
      dplyr::select(-strand)%>%
      dplyr::mutate(colfork_id=paste0(colfork_chrom, ":", colfork_start, "-", colfork_end))
    colfork_ranges = colfork_df %>%
      df2ranges(colfork_chrom, colfork_start, colfork_end)

    #
    # Find termination zones
    #
    tz_df = iz_ranges %>%
      innerJoinByOverlaps(colfork_ranges) %>%
      dplyr::group_by(colfork_id, colfork_chrom, colfork_start, colfork_end) %>%
      dplyr::filter(dplyr::n()==2) %>%
      dplyr::summarize(colfork_min_pos=min(iz_start)+25e3, colfork_max_pos=max(iz_end)-25e3, .groups='drop') %>%
      dplyr::ungroup() %>%
      df2ranges(colfork_chrom, colfork_min_pos, colfork_max_pos) %>%
      innerJoinByOverlaps(rfd_ranges) %>%
      dplyr::arrange(colfork_chrom, colfork_min_pos, colfork_max_pos, rfd_chrom, rfd_start) %>%
      dplyr::group_by(colfork_chrom, colfork_min_pos, colfork_max_pos) %>%
      dplyr::do((function(z) {
        zz <<- z
        # z = tz_df %>% dplyr::filter(grepl("99430908", colfork_id))

        rfd_model = loess(rfd_score ~ rfd_start, data=z, span=0.3)
        rfd_pred_df = data.frame(rfd_start=seq(min(z$rfd_start), max(z$rfd_start), 1e3)) %>%
          dplyr::mutate(rfd_chrom=z$colfork_chrom[1], rfd_end=rfd_start+1e3, rfd_score = predict(rfd_model, .))
        rfd_pred_tz_df = rfd_pred_df %>%
          dplyr::filter(abs(rfd_score) < 0.01)
        if(nrow(rfd_pred_tz_df)==0) return(data.frame())
        rfd_pred_tz_df = rfd_pred_tz_df %>%
          df2ranges(rfd_chrom, rfd_start, rfd_end) %>%
          GenomicRanges::reduce(min.gapwidth=1e4) %>%
          as.data.frame() %>%
          dplyr::select(rfd_reduced_chrom=seqnames, rfd_reduced_start=start, rfd_reduced_end=end) %>%
          df2ranges(rfd_reduced_chrom, rfd_reduced_start, rfd_reduced_end) %>%
          innerJoinByOverlaps(rfd_pred_df %>% df2ranges(rfd_chrom, rfd_start, rfd_end)) %>%
          dplyr::group_by(tz_chrom=rfd_chrom, rfd_reduced_start, rfd_reduced_end) %>%
          dplyr::summarise(tz_start=weighted.mean(rfd_start, w=max(abs(rfd_score))-abs(rfd_score)), tz_end=tz_start, tz_iz_left=z$colfork_min_pos[1], tz_iz_right=z$colfork_max_pos[1], .groups='drop') %>%
          dplyr::mutate(tz_start=tz_start-training_binsize/4, tz_end=tz_end+training_binsize/4)

        # tz_ggplot = data.frame(tz_start=GenomicRanges::start(rfd_pred_ranges), tz_end=GenomicRanges::end(rfd_pred_ranges)) %>% dplyr::mutate(tz_length=tz_end-tz_start)
        # ggplot(z, aes(x=rfd_start, y=rfd_score)) +
        #   geom_line() +
        #   geom_line(aes(x=rfd_start, y=rfd_score), data=rfd_pred_df) +
        #   geom_segment(aes(x=tz_start, xend=tz_end, y=0, yend=0), data=rfd_pred_tz_df, size=3, color="#FF0000") +
        #   geom_smooth(method="loess", span=.3)
        #   # geom_vline(xintercept=GenomicRanges::start(rfd_pred_ranges)[1])

        if(nrow(rfd_pred_tz_df)>1) {
            return(data.frame())
        }

        rfd_pred_tz_df %>% dplyr::select(dplyr::matches("tz_")) %>% dplyr::mutate(tz_alternative_n=1)
        # r = data.frame(tz_chrom=z$colfork_chrom[1], tz_start=GenomicRanges::start(rfd_pred_ranges), tz_end=GenomicRanges::end(rfd_pred_ranges), tz_alternative_n=nrow(rfd_pred_df))
      })(.)) %>%
      dplyr::ungroup() %>%
      dplyr::select(dplyr::starts_with("tz_"))

    tz_df = tz_df %>%
      df2ranges(tz_chrom, tz_start-100e3, tz_end+100e3) %>%
      innerJoinByOverlaps(repliseqTime_compare_df %>% dplyr::filter(Run=="Zhao_mESC16") %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end)) %>%
      dplyr::group_by(tz_chrom, tz_start, tz_end, tz_iz_left, tz_iz_right, tz_alternative_n) %>%
      dplyr::mutate(tz_peak_start=(repliseq_start+repliseq_end)[which.max(repliseq_mean_fraction)]/2) %>%
      dplyr::filter(repliseq_mean_fraction>=max(repliseq_mean_fraction)*0.99) %>%
      dplyr::do((function(tz) {
          tzz<<-tz
          #print(nrow(tz))
          #asdaD()
          tz_ranges = tz %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end)
          tz_reduced_ranges = tz_ranges %>%
            GenomicRanges::reduce() %>%
            as.data.frame() %>%
            # dplyr::mutate(tz_repliseq_peak_pos=(tz$repliseq_start+tz$repliseq_end)[which.max(tz$repliseq_mean_fraction)]/2, tz_extended_length=end-start) %>%
            dplyr::select(tz_extended_chrom=seqnames, tz_extended_start=start, tz_extended_end=end) %>%
            df2ranges(tz_extended_chrom, tz_extended_start, tz_extended_end)
          tz_reduced_ranges %>%
            innerJoinByOverlaps(tz_ranges)  %>%
            dplyr::filter(tz_extended_start<=tz_peak_start & tz_peak_start<=tz_extended_end) %>%
            dplyr::group_by(tz_extended_chrom, tz_extended_start, tz_extended_end) %>%
            dplyr::summarise(tz_peak_start=weighted.mean(repliseq_start/2+repliseq_end/2, repliseq_mean_fraction)-1, tz_peak_end=tz_peak_start+2, .groups="keep")

            # dplyr::select(tz_extended_start=start, tz_extended_end=end, tmp_corected_pos, tz_ragion_length)
      })(.)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(tz_okseq_start=tz_start, tz_okseq_end=tz_end, tz_start=tz_peak_start, tz_end=tz_peak_end)
      #dplyr::filter(tz_extended_start<=tz_peak_start & tz_peak_start<=tmp_corrected_end) %>%
      #dplyr::mutate(tz_start=tmp_corrected_start/2+tmp_corrected_end/2, tz_end=tmp_corrected_start/2+tmp_corrected_end/2) %>%
      #dplyr::select(dplyr::starts_with("tz_"))
    #table(tza_df$tz_extended_end-tza_df$tz_extended_start)/nrow(tza_df)

    readr::write_tsv(tz_df, file="data/tz.tsv")

    #
    # Split into training regions
    #
    training_regions_df = as.data.frame(GenomicRanges::reduce(colfork_ranges)) %>%
      dplyr::select(-strand) %>%
      dplyr::rename(colfork_reduced_chrom="seqnames", colfork_reduced_start="start", colfork_reduced_end="end", colfork_reduced_width="width") %>%
      dplyr::mutate(seq_start=colfork_reduced_start-training_margin, seq_end=colfork_reduced_end-training_width+training_margin) %>%
      dplyr::mutate(seq_start=round(seq_start/training_binsize)*training_binsize, seq_end=round(seq_end/training_binsize)*training_binsize) %>%
      dplyr::filter(seq_end>=seq_start) %>%
      dplyr::group_by(colfork_reduced_chrom, colfork_reduced_start, colfork_reduced_end) %>%
      dplyr::summarize(training_chrom=colfork_reduced_chrom[1], training_start=c(seq(seq_start, seq_end, by=training_step), seq_end-training_step), training_end=training_start+training_width) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(colfork_id=paste0(colfork_reduced_chrom, ":", colfork_reduced_start, "-", colfork_reduced_end), training_id=paste0(training_chrom, "_", format(training_start, scientific=F, justify="none", trim=T), "-", format(training_end, scientific=F, justify="none", trim=T)))

    training_regions_clean_df = training_regions_df %>%
      df2ranges(training_chrom, training_start, training_end) %>%
      leftJoinByOverlaps(tz_df %>% df2ranges(tz_chrom, tz_start, tz_end)) %>%
      dplyr::group_by(training_chrom, training_start, training_end) %>%
      dplyr::mutate(tz_margin=min(c(1e6, abs(training_start-tz_start), abs(training_end-tz_end)), na.rm=T)) %>%
      dplyr::filter(tz_margin>200e3) %>%
      dplyr::distinct_at(colnames(training_regions_df))

    #
    # Create training and label data.frames
    #
    data_df = training_regions_clean_df %>%
      df2ranges(training_chrom, training_start, training_end) %>%
      innerJoinByOverlaps(repliseq_compare_df %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end), minoverlap=training_binsize) %>%
      dplyr::group_by(Run, training_chrom, training_start, training_end) %>%
      dplyr::mutate(repliseq_rel_start=repliseq_start-training_start) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(Run, colfork_id, training_id, repliseq_chrom, repliseq_start, repliseq_end) %>%
      dplyr::mutate(repliseq_value=repliseq_value/max(repliseq_value)) %>%
      dplyr::ungroup()

    label_df = training_regions_clean_df %>%
      dplyr::group_by(colfork_id, training_id, training_chrom, training_start, training_end) %>%
      dplyr::summarize(training_tz_start=seq(training_start, training_end-1, training_binsize), training_tz_end=training_tz_start+training_binsize, .groups='drop') %>%
      dplyr::mutate(training_tz_rel_start=training_tz_start-training_start) %>%
      df2ranges(training_chrom, training_tz_start, training_tz_end) %>%
      leftJoinByOverlaps(tz_df %>% df2ranges(tz_chrom, tz_start, tz_end)) %>%
      dplyr::group_by(colfork_id, training_id, training_chrom, training_start, training_end, tz_start, tz_end, training_tz_rel_start, training_tz_start, training_tz_end) %>%
      dplyr::summarize(training_tz_present=as.integer(!is.na(tz_alternative_n)), .groups='drop') %>%
      dplyr::arrange(training_tz_rel_start) %>%
      dplyr::ungroup()
      # dplyr::inner_join(data_df %>% dplyr::distinct(repliseq_chrom, repliseq_rel_start, repliseq_) %>% dplyr::summarize(repliseq_value.mean=weighted.mean(repliseq_fraction, repliseq_value)), by=c("training_chrom"="repliseq_chrom", "training_tz_rel_start"="repliseq_rel_start")) %>%
      # dplyr::group_by(colfork_id, training_id, training_start, training_end, tz_start, tz_end) %>%
      # dplyr::mutate(training_tz_present=ifelse(training_width-training_tz_rel_start>=training_binsize*2, training_tz_present, 0)) %>%
      # dplyr::ungroup()


    if(debug==F)
    {
        #
        # Plot training samples distribution
        #
        debug_tz_df = training_regions_clean_df %>%
          df2ranges(training_chrom, training_start, training_end) %>%
          innerJoinByOverlaps(tz_df %>% df2ranges(tz_chrom, tz_start, tz_end)) %>%
          dplyr::group_by(training_chrom, training_start, training_end) %>%
          dplyr::mutate(tz_rel_start=tz_start-training_start, tz_rel_end=tz_end-training_start) %>%
          dplyr::ungroup()

        debug_rfd_df = training_regions_clean_df %>%
          df2ranges(training_chrom, training_start, training_end) %>%
          innerJoinByOverlaps(rfd_ranges) %>%
          dplyr::group_by(training_chrom, training_start, training_end) %>%
          dplyr::mutate(rfd_rel_start=rfd_start-training_start) %>%
          dplyr::ungroup()

        debug_repliseqTime_df = training_regions_clean_df %>%
          df2ranges(training_chrom, training_start, training_end) %>%
          innerJoinByOverlaps(repliseqTime_compare_df %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end)) %>%
          dplyr::group_by(training_chrom, training_start, training_end) %>%
          dplyr::mutate(repliseq_rel_start=repliseq_start-training_start, repliseq_rel_end=repliseq_end-training_start) %>%
          dplyr::ungroup()

        training_regions_ggplot = label_df %>%
          dplyr::group_by(training_id) %>%
          dplyr::summarize(training_tz_count=as.character(sum(training_tz_present, na.rm=T))) %>%
          dplyr::group_by(training_tz_count) %>%
          dplyr::summarize(training_samples_count=dplyr::n()) %>%
          dplyr::bind_rows(data.frame(training_tz_count="Total", training_samples_count=label_df %>% dplyr::distinct(training_id) %>% nrow()))
        ggplot(training_regions_ggplot) +
          geom_bar(aes(x=training_tz_count, y=training_samples_count), stat="identity") +
          labs(y="Training samples", x="No. of termination zones in sample") +
          scale_y_continuous(breaks = scales::pretty_breaks(n=20)) +
          ggpubr::theme_pubclean()

        # names(label_matrix)[rep(samples_training, 1)
        # reshape2::melt(predict(unet, array(x_train[1:(nrow(x_train)/2),,], dim=c(dim(x_train), 1)))[1:(nrow(x_train)/2),,1,1]) %>% dplyr::mutate(training_id=names(label_matrix)[rep(samples_training, 1)[Var1]], source="train", type="pred")

        #
        # Plot some examples
        #
        example_training_ids = training_regions_clean_df$training_id[c(1,500,1000,1500)]
        example_training_ids = with(training_regions_clean_df, training_id[grepl("chr18_628", training_id)])
        example_training_ids = with(training_regions_clean_df, training_id[grepl("chr1_183550000|chr1_88850000|chr1_134350000|chr1_117300000", training_id)])

        pdf("reports/tzNN_examples.pdf", width=8.27, height=8.27)
        for(example_training_ids in with(training_regions_clean_df, training_id[grepl("chr1_183550000|chr1_88850000|chr1_134350000|chr1_117300000", training_id)]))
        {
            debug_rfd_ggplot = debug_rfd_df %>%
              dplyr::filter(training_id %in% example_training_ids) %>%
              dplyr::arrange(training_id, rfd_rel_start) %>%
              dplyr::mutate(set="RFD", rfd_score=dplyr::case_when(rfd_score>0 & dplyr::lag(rfd_score)>0 | rfd_score<0 & dplyr::lag(rfd_score)<0~rfd_score, T~0))
            data_ggplot = dplyr::bind_rows(
              data_df %>% dplyr::mutate(set=Run, x=repliseq_rel_start, y=repliseq_fraction, value=repliseq_value) %>% dplyr::filter(Run=="Zhao_mESC16"),
              # as.matrix(as.data.frame(predict(unet, x), row.names=names(data_matrix$Zhao_mESC16), check.names=F)) %>%
              #     reshape2::melt() %>%
              #     dplyr::mutate(x=as.numeric(gsub("^V", "", as.character(Var2)))*training_binsize) %>%
              #     dplyr::mutate(training_id=as.character(Var1), colfork_id=training_id, set="Prediction") %>%
              #     dplyr::select(set, training_id, colfork_id, x, value) %>%
              #     tidyr::crossing(y=1:2),
              # df_all %>%
              #     dplyr::filter(validation_chrom=="chr1" & type=="pred" & source=="validation") %>%
              #     dplyr::mutate(x=Var2*training_binsize, colfork_id=training_id, set="Prediction") %>%
              #     dplyr::select(set, training_id, colfork_id, x, value) %>%
              #     tidyr::crossing(y=1:2),
              label_df %>%
                tidyr::crossing(data.frame(y=1)) %>%
                dplyr::mutate(set="Labels", x=training_tz_rel_start, value=training_tz_present)) %>%
              dplyr::filter(training_id %in% example_training_ids)

            p = ggplot(data_ggplot) +
                geom_tile(aes(x=as.numeric(x)+training_binsize/2, y=y, fill=value)) +
                geom_line(aes(x=repliseq_rel_start, y=repliseq_mean_fraction), color="#FFFFFF", data=debug_repliseqTime_df %>% dplyr::filter(training_id %in% example_training_ids & Run == "Zhao_mESC16") %>% dplyr::mutate(set=Run)) +
                # geom_area(aes(x=rfd_rel_start, y=rfd_score*4), fill="#FF0000", data=debug_rfd_ggplot %>% dplyr::filter(!is.na(rfd_score) & rfd_score>0)) +
                # geom_area(aes(x=rfd_rel_start, y=rfd_score*4), fill="#0000FF", data=debug_rfd_ggplot %>% dplyr::filter(!is.na(rfd_score) & rfd_score<=0)) +
                geom_area(aes(x=rfd_rel_start, y=rfd_score*4), data=debug_rfd_ggplot) +
                geom_segment(aes(x=tz_rel_start, xend=tz_rel_end, y=0, yend=0), data=debug_tz_df %>% dplyr::filter(training_id %in% example_training_ids) %>% dplyr::mutate(set="RFD"), color="#FF0000", size=2) +
                geom_segment(aes(y=yintercept, yend=yintercept, x=0, xend=0), data=data.frame(set="RFD", yintercept=0)) +
                facet_grid(set~training_id, scales="free", space="free", switch="y") +
                scico::scale_fill_scico(palette = "bilbao", guide="none") +
                scale_x_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="K"))  +
                theme_bw() +
                theme(
                  strip.text.y.left=element_text(angle=0),
                  strip.background = element_blank(),
                  axis.title=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  # panel.spacing=unit(0, "lines"),
                  panel.border=element_blank())
            print(p)
        }
        dev.off()
    }

    #
    # Final matrix
    #
    data_matrix = data_df %>%
      dplyr::mutate(repliseq_rel_start=format(repliseq_rel_start, scientific=F, justify="none", trim=F)) %>%
      split(.$Run) %>%
      lapply(function(y) {
        y %>%
          reshape2::dcast(training_id+repliseq_fraction ~ repliseq_rel_start, value.var="repliseq_value") %>%
          dplyr::arrange(dplyr::desc(repliseq_fraction)) %>%
          split(.$training_id) %>%
          lapply(FUN=function(z) z %>% `rownames<-`(NULL) %>% tibble::column_to_rownames("repliseq_fraction") %>% dplyr::select(-training_id) %>% as.matrix()) %>%
          `[`(training_regions_clean_df$training_id)
      })

    label_matrix = label_df %>%
      dplyr::mutate(training_tz_rel_start=format(training_tz_rel_start, scientific=F, justify="none", trim=T), training_tz_rel_start=factor(training_tz_rel_start, unique(training_tz_rel_start))) %>%
      reshape2::dcast(training_id ~ training_tz_rel_start, value.var="training_tz_present", drop=F) %>%
      split(.$training_id) %>%
      lapply(FUN=function(z) z %>% dplyr::select(-training_id) %>% as.matrix()) %>%
      `[`(training_regions_clean_df$training_id)

    #
    # Write final data to json files
    #
    dir.create("data/tz", recursive=T, showWarnings=F)
    save(data_matrix, file="data/tz/data_matrix.rda")
    save(label_matrix, file="data/tz/label_matrix.rda")

    if(debug==T) {
      #
      # PLoting final data for testing
      #
      i = x
      dplyr::bind_rows(
        reshape2::melt(data_matrix$Zhao_mESC[[i]]) %>% dplyr::mutate(set="Data"),
        reshape2::melt(label_matrix[[i]]) %>% dplyr::mutate(set="Labels")) %>%
        dplyr::mutate(Var1=factor(Var1), training_id=gsub("_", ":", names(data_matrix$Zhao_mESC8)[i])) %>%
        ggplot() +
          geom_tile(aes(x=Var2, y=Var1, fill=value)) +
          facet_grid(set~training_id, scales="free", space="free") +
          scico::scale_fill_scico(palette = "bilbao") +
          theme_bw() +
          theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    }


    #
    # Noise model
    #
    mtx_source <- do.call(cbind, data_matrix$Zhao_mESC8)
    mtx_target <- do.call(cbind, data_matrix$RS1_mNPC8)
    both_ldf = dplyr::bind_rows(
      reshape2::melt(unname(mtx_source)) %>% dplyr::mutate(dataset="source"),
      reshape2::melt(unname(mtx_target)) %>% dplyr::mutate(dataset="target")) %>%
      dplyr::rename(fraction="Var1", bin="Var2") %>%
      dplyr::group_by(dataset, bin) %>%
      dplyr::filter(all(!is.na(value))) %>%
      dplyr::mutate(fraction_highest=fraction[which.max(value)]) %>%
      dplyr::ungroup()

    noise_df = both_ldf %>%
      dplyr::group_by(dataset, fraction_highest, fraction) %>%
      dplyr::summarise(vsd=sd(value), value=mean(value)) %>%
      dplyr::group_by(fraction_highest, fraction) %>%
      dplyr::summarise(
        noise_mean=value[dataset=="target"]-value[dataset=="source"],
        noise_sd=pmax(0.01, vsd[dataset=="target"]-vsd[dataset=="source"]))

    if(debug==T) {
        i = 202
        data_noise = lapply(data_matrix$RS1_mNPC8, tzNN_noise_add, noise_df)
        dplyr::bind_rows(
          reshape2::melt(data_matrix$Zhao_mESC8[[i]]) %>% dplyr::mutate(set="Zhao"),
          reshape2::melt(data_noise[[i]]) %>% dplyr::mutate(set="Zhao + Noise"),
          reshape2::melt(data_matrix$RS1_mNPC[[i]]) %>% dplyr::mutate(set="Wei Run 01")) %>%
          dplyr::mutate(Var1=factor(Var1), training_id=gsub("_", ":", names(data_matrix)[i])) %>%
          ggplot() +
          geom_tile(aes(x=Var2, y=Var1, fill=value)) +
          facet_grid(set~training_id, scales="free", space="free") +
          scico::scale_fill_scico(palette = "bilbao") +
          theme_bw() +
          theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    }
}

tzNN_train = function() {
  load(file="data/tz/data_matrix.rda")
  load("data/tz/label_matrix.rda")
  x = aperm(array(unlist(data_matrix$Zhao_mESC16), dim=c(16, 20, length(data_matrix$Zhao_mESC16))), c(3,2,1))
  y = aperm(array(as.double(unlist(label_matrix)), dim=c(1, 20, length(data_matrix$Zhao_mESC16))), c(3,2,1))


  df_all = data.frame()
  for(chr in c("*", paste0("chr", 1:19)) {
      #sample_validation = grep("(chr12|chr19)_", names(label_matrix))
      sample_validation = grep(paste0("^", chr, "_"), names(label_matrix))
      samples_training = setdiff(1:nrow(x), sample_validation)
      print(paste0("Validation proportion (", chr, "): ", round(100*length(sample_validation)/(length(sample_validation)+length(samples_training)), 2), "%"))
      x_validate = x[sample_validation,,]
      y_validate = array(y[sample_validation,,], dim=c(length(sample_validation), 20, 1))
      x_train = x[samples_training,,]
      y_train = array(y[samples_training,,], dim=c(length(samples_training), 20, 1))
      x_train = abind::abind(x_train, pretrube_mirror(x_train), along=1)
      y_train = abind::abind(y_train, pretrube_mirror(y_train), along=1)

      set.seed(123)
      args = list(batch_size=32, learning_rate=1e-4, dropout1=0.9, dropout2=0.2)
      unet = get_deep(input_shape=c(dim(x_train)[2:3], 1),  dropout1=args$dropout1, dropout2=args$dropout2)
      unet = unet %>% keras::compile(
            optimizer=keras::optimizer_adam(learning_rate=args$learning_rate),
            loss=loss_dice,
            #loss=keras::loss_binary_crossentropy(),
            #loss=keras::loss_poisson(),
            metrics=list(metric_dice, keras::metric_binary_crossentropy, keras::metric_poisson, metric_intersection, metric_junction)
          )
      #print(unet)
      model_math = paste0("logs/", format(Sys.time(), "%m-%d_%H%M"), "-", stringr::str_glue("INCEPTION_3x5_doA{do1}_doB{do2}_bs{bs}_lr{lr}_DICE_{chr}", chr=chr, do1=args$dropout1, do2=args$dropout2, bs=args$batch_size, lr=args$learning_rate))
      history = keras::fit(unet, x_train, y_train, validation_data=list(x_val=x_validate, y_val=y_validate),
        batch_size=args$batch_size, epochs=10000, verbose=1, shuffle=T,
        callbacks=list(keras::callback_tensorboard(model_math, histogram_freq=5, embeddings_freq=5), keras::callback_progbar_logger(count_mode="steps"), keras::callback_early_stopping(monitor="val_loss", min_delta=0.000001, patience=1000, restore_best_weights=T)))
      print(unet %>% evaluate(x_validate, y_validate))


      keras::save_model_tf(unet, paste0(model_math, "/model"))
      # unet = keras::load_model_tf(paste0(model_math, "/model/saved_model.pb"))
      pdf(paste0(model_math, "/evaluation.pdf"), width=8.27, height=11.69)
      plot(unet)

      predictions_val = predict(unet, x_validate)
      rockr_val_pred = ROCR::prediction(as.numeric(predictions_val), as.numeric(y_validate))
      rockr_val_acc = ROCR::performance(rockr_val_pred, measure = "prec", x.measure="rec")
      rockr_val_perf = ROCR::performance(rockr_val_pred, measure = "tpr", x.measure = "fpr")
      plot(rockr_val_acc, colorize=T, main="ROC curve (validation data)")
      plot(rockr_val_perf, colorize=T, main="Precision/recall curve (validation data)"); abline(a=0, b=1)


      predictions_train = predict(unet, x_train)
      rockr_train_pred = ROCR::prediction(as.numeric(predictions_train), as.numeric(y_train))
      rockr_train_acc = ROCR::performance(rockr_train_pred, measure = "prec", x.measure="rec")
      rockr_train_perf = ROCR::performance(rockr_train_pred, measure = "tpr", x.measure = "fpr")
      plot(rockr_train_acc, colorize=T, main="ROC curve (training data)")
      plot(rockr_train_perf, colorize=T, main="Precision/recall curve (training data)"); abline(a=0, b=1)


      df = dplyr::bind_rows(
        as.matrix(as.data.frame(predict(unet, x), row.names=names(data_matrix$Zhao_mESC16), check.names=F)) %>%
          reshape2::melt() %>%
          dplyr::mutate(Var2=as.numeric(gsub("^V", "", as.character(Var2)))) %>%
          dplyr::mutate(training_id=as.character(Var1), type="pred"),
        as.matrix(as.data.frame(y, row.names=names(data_matrix$Zhao_mESC16), check.names=F)) %>%
          reshape2::melt() %>%
          dplyr::mutate(Var2=as.numeric(gsub("^V", "", as.character(Var2))))  %>%
          dplyr::mutate(training_id=as.character(Var1), type="true")
      ) %>%
        dplyr::mutate(source=dplyr::case_when(training_id %in% names(data_matrix$Zhao_mESC16)[samples_training] ~ "training", T~"validation")) %>%
        dplyr::mutate(value=pmax(0, value)) %>%
        dplyr::group_by(source, training_id) %>%
        dplyr::filter(any(value>=0.1)) %>%
        dplyr::group_by(source) %>%
        dplyr::filter(training_id %in% unique(training_id)[sample(length(unique(training_id)), pmin(length(unique(training_id)), 100))]) %>%
        dplyr::ungroup()
      df_all = dplyr::bind_rows(df_all, df %>% dplyr::mutate(validation_chrom=chr))
      # df = df_all %>% dplyr::filter(validation_chrom=="chr1")
      df2 = df %>%
        dplyr::distinct(Var2, training_id, type, source, value) %>%
        tidyr::crossing(data.frame(shift=c(-0.4001, -0.4, 0.4, 0.4001), value_mod=c(0, 1, 1, 0))) %>%
        dplyr::mutate(Var2=Var2+shift, value=value*value_mod, diff=0) %>%
        dplyr::arrange(Var2, training_id, type, source, value) %>%
        dplyr::bind_rows(df)
      ggplot(df2) +
        ggridges::geom_ridgeline(aes(x=Var2, y=training_id, height=value, color=type), alpha=0.0, scale=0.5) +
        facet_wrap(~source, scales="free_y") +
        scale_color_manual(values=c("pred"="#333333", "true"="#FF0000")) +
        theme_paper()
      dev.off()
  }

  # df_all %>%
  #   readr::write_tsv("reports/tzNN_crossvalidation_results.tsv")
  # df_all = readr::read_tsv("reports/tzNN_crossvalidation_results.tsv")

  pdf("reports/tzNN_performance.pdf", width=8.27, height=8.27)
  rockr_all_data = df_all %>%
    dplyr::arrange(source, validation_chrom, training_id, Var2) %>%
    dplyr::filter(source=="validation") %>%
    reshape2::dcast(source+validation_chrom+training_id+Var2~type, value.var="value") %>%
    dplyr::group_by(source, validation_chrom, training_id) %>%
    dplyr::mutate(swap=dplyr::case_when(
      true==1 & !is.na(dplyr::lag(pred)) & dplyr::lag(pred)>pred & dplyr::lag(pred)>dplyr::lead(pred) ~ "lag",
      true==1 & !is.na(dplyr::lead(pred)) & dplyr::lead(pred)>pred & dplyr::lag(pred)<=dplyr::lead(pred)  ~ "lead",
      T ~ "orig"
    )) %>% plyr::mutate(pred=dplyr::case_when(
      swap=="lag" ~ dplyr::lag(pred),
      swap=="lead"  ~ dplyr::lead(pred),
      true==0 & dplyr::lag(swap)=="lead" ~ 0,
      true==0 & dplyr::lead(swap)=="lag" ~ 0,
      # true==0 & (!is.na(dplyr::lag(true)) & dplyr::lag(true)>true | !is.na(dplyr::lead(true)) & dplyr::lead(true)>true) ~ pmin(dplyr::lag(true), dplyr::lead(true), na.rm=T),
      T ~ pred
    )) %>%
    dplyr::ungroup()

  rockr_all_pred = ROCR::prediction(rockr_all_data$pred, rockr_all_data$true)
  # rockr_all_pred = ROCR::prediction(rockr_all_data %>% dplyr::group_split(validation_chrom) %>% lapply(function(l) l$pred), rockr_all_data %>% dplyr::group_split(validation_chrom) %>% lapply(function(l) l$true))
  rockr_all_acc = ROCR::performance(rockr_all_pred, measure = "prec", x.measure="rec")
  rockr_all_perf = ROCR::performance(rockr_all_pred, measure = "tpr", x.measure = "fpr")
  plot(rockr_all_acc, colorize=T, main="Precision/recall curve (validation data)")
  plot(rockr_all_perf, colorize=T, main="ROC curve (validation data)"); abline(a=0, b=1)
  dev.off()
}