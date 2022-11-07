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


# 1. Connect small gaps
# 2. Adjust using lower line of repliseq
# 3. Better model for correct TZ selection
# 4 Allow some left/right margin when traaining
# 5. Return denser steps for training

save_model = function(model, path, include_model=T, include_weights=T){
    dir.create(dirname(path), recursive=T, showWarnings=F)
    # Save model and weights. JSON can be visualized with https://Netron.app
    if(include_model) {
        write(keras::model_to_json(model), paste0(path, "_model.json"))
    }
    if(include_weights) {
        keras::save_model_weights_hdf5(model, paste0(path, "_weights.h5"))
    }
}

load_model = function(path) {
    model_path = paste0(path, "_model.json")
    weights_path = paste0(path, "_weights.h5")
    if(!file.exists(model_path)) {
        stop(paste0("Model file '", weights_path, "' not found"))
    }

    model = keras::model_from_json(readLines(model_path)) %>%
        keras::compile(optimizer=keras::optimizer_adam(), loss=loss_dice, metrics=list(metric_dice, keras::metric_binary_crossentropy, keras::metric_poisson))

    if(file.exists(weights_path)) {
        model = model %>% keras::load_model_weights_hdf5(weights_path)
    } else {
        writeLines(paste0("Weights file '", weights_path, "' not found. Skipping..."))
    }

    model
}

create_model = function(input_shape, dropout1=0.8, dropout2=0.2) {
    inputs = keras::layer_input(shape = input_shape) %>%
      keras::layer_batch_normalization(name="initial_batchnorm")

    layer01 = inputs %>%
      keras::layer_conv_2d(filters=32, kernel_size=c(3, 3), padding="valid", name="layer01_conv1") %>%
      layer_activate(name="layer01_act1")
    # layer1pool = layer01 %>%
    #   keras::layer_average_pooling_2d(pool_size=c(5, 5), stride=c(1, 1), padding="valid", name="layer1pool_pool")


    layer0pool = inputs %>%
      keras::layer_average_pooling_2d(pool_size=c(5, 16), stride=c(1, 1), padding="valid", name="layer0pool_pool")
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
    layer14 = layer01 %>%
      keras::layer_conv_2d(filters=8, kernel_size=c(1, 1), padding="valid", name="layer14_conv1")  %>%
      layer_activate(name="layer14_act1", dropout=dropout2)%>%
      keras::layer_depthwise_conv_2d(kernel_size=c(9, 7), padding="valid", name="layer14_dv_conv1") %>%
      keras::layer_conv_2d(filters=8, kernel_size=c(1, 1), padding="valid", name="layer14_conv2") %>%
      # keras::layer_conv_2d(filters=8, kernel_size=c(7, 7), padding="valid") %>%
      layer_activate(name="layer14_act2")


    layer21 = layer11 %>%
      keras::layer_conv_2d(filters=4, kernel_size=c(1, 1), padding="valid", name="layer21_conv1")  %>%
      layer_activate(name="layer21_act1", dropout=dropout2) %>%
      keras::layer_depthwise_conv_2d(kernel_size=c(5, 5), padding="valid", name="layer21_dv_conv1") %>%
      keras::layer_conv_2d(filters=4, kernel_size=c(1, 1), padding="valid", name="layer21_conv2") %>%
      # keras::layer_conv_2d(filters=8, kernel_size=c(7, 7), padding="valid") %>%
      layer_activate(name="layer21_act2")
    layer22 = layer11 %>%
      keras::layer_conv_2d(filters=4, kernel_size=c(1, 1), padding="valid", name="layer22_conv1")  %>%
      layer_activate(name="layer22_act1", dropout=dropout2)%>%
      keras::layer_depthwise_conv_2d(kernel_size=c(7, 5), padding="valid", name="layer22_dv_conv1") %>%
      keras::layer_conv_2d(filters=4, kernel_size=c(1, 1), padding="valid", name="layer22_conv2") %>%
      # keras::layer_conv_2d(filters=8, kernel_size=c(7, 7), padding="valid") %>%
      layer_activate(name="layer22_act2")
    layer23 = layer11 %>%
      keras::layer_conv_2d(filters=4, kernel_size=c(1, 1), padding="valid", name="layer23_conv1")  %>%
      layer_activate(name="layer23_act1", dropout=dropout2)%>%
      keras::layer_depthwise_conv_2d(kernel_size=c(9, 5), padding="valid", name="layer23_dv_conv1") %>%
      keras::layer_conv_2d(filters=4, kernel_size=c(1, 1), padding="valid", name="layer23_conv2") %>%
      # keras::layer_conv_2d(filters=8, kernel_size=c(7, 7), padding="valid") %>%
      layer_activate(name="layer23_act2")


    inputs_crop = inputs %>%
      keras::layer_conv_2d(filters=1, kernel_size=c(1, 3), padding="valid", name="inputs_conv1") %>%
      keras::layer_cropping_2d(cropping=c(4,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,0))
    layer10_crop = layer10 %>%
      keras::layer_cropping_2d(cropping=c(3,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,0))
    layer0pool_crop = layer0pool %>%
      keras::layer_cropping_2d(cropping=c(2,0)) %>%
      keras::layer_zero_padding_2d(padding=list(c(0,0), c(0,13)))


    layer11_crop = layer11 %>%
      keras::layer_cropping_2d(cropping=c(2,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,1))
    layer12_crop = layer12 %>%
      keras::layer_cropping_2d(cropping=c(1,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,2))
    layer13_crop = layer13 %>%
      keras::layer_cropping_2d(cropping=c(0,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,3))
    layer14_crop = layer14 %>%
      keras::layer_cropping_2d(cropping=c(0,0)) %>%
      keras::layer_zero_padding_2d(padding=c(1,3))


    layer21_crop = layer21 %>%
      keras::layer_cropping_2d(cropping=c(0,0)) %>%
      keras::layer_zero_padding_2d(padding=c(0,3))
    layer22_crop = layer22 %>%
      keras::layer_cropping_2d(cropping=c(0,0)) %>%
      keras::layer_zero_padding_2d(padding=c(1,3))
    layer23_crop = layer23 %>%
      keras::layer_cropping_2d(cropping=c(0,0)) %>%
      keras::layer_zero_padding_2d(padding=c(2,3))


    concatenate1 = keras::layer_concatenate(list(layer11_crop, layer12_crop, layer21_crop, layer22_crop, layer23_crop), axis=3) %>%
      keras::layer_spatial_dropout_2d(rate = dropout1)
    classify = keras::layer_concatenate(list(inputs_crop, layer10_crop, layer0pool_crop, concatenate1), axis=3) %>%
      keras::layer_conv_2d(filters=4, kernel_size=c(1, 1), padding="valid", name="classify_conv1") %>%
      layer_activate(name="classify_act1") %>%
      keras::layer_conv_2d(filters=2, kernel_size=c(1, 1), padding="valid", name="classify_conv2") %>%
      keras::layer_max_pooling_2d(pool_size=c(1, 14), stride=c(1, 1), padding="valid", name="classify_pool1") %>%
      layer_activate(name="classify_act2", activation="sigmoid") %>%
      keras::layer_zero_padding_2d(padding=c(4,0)) %>%
      keras::layer_conv_2d(filters=1, kernel_size=c(9, 1), padding="valid", name="classify_conv3") %>%
      layer_activate(name="classify_act3", activation="sigmoid") %>%
      keras::layer_zero_padding_2d(padding=c(4,0))

    model = keras::keras_model(inputs=inputs, outputs=classify)
    # model %>% keras::compile(optimizer=keras::optimizer_rmsprop(learning_rate=learning_rate), loss=keras::loss_binary_crossentropy)
    return(model)
}

repliseq_zoomout = function(repliseq_df, binsize=50e3) {
  repliseq_zoom_df = repliseq_df %>%
    dplyr::arrange(repliseq_chrom, repliseq_start, repliseq_fraction) %>%
    dplyr::mutate(repliseq_group=floor(repliseq_start/(2*binsize))) %>%
    dplyr::group_by(repliseq_chrom, repliseq_fraction) %>%
    dplyr::mutate(
      repliseq_end=dplyr::lead(repliseq_end),
      repliseq_value=ifelse(is.na(dplyr::lead(repliseq_value)), repliseq_value, repliseq_value/2+dplyr::lead(repliseq_value)/2)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(repliseq_chrom, repliseq_group, repliseq_fraction, .keep_all=T) %>%
    dplyr::mutate(repliseq_zoomout=2) %>%
    dplyr::select(-repliseq_group)

  dplyr::bind_rows(repliseq_df %>% dplyr::mutate(repliseq_zoomout=1), repliseq_zoom_df)
}

prepare_data = function(repliseq_df, binsize=50e3, binsmargin=6, binspersample=30, fractions=16) {
  r_df = repliseq_zoomout(repliseq_df, binsize=binsize) %>%
    dplyr::group_by(repliseq_zoomout, repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::mutate(repliseq_value=repliseq_value/max(repliseq_value, na.rm=T)) %>%
    dplyr::ungroup()

  prediction_regions_df = r_df %>%
    dplyr::distinct(repliseq_zoomout, repliseq_chrom, repliseq_start) %>%
    dplyr::arrange(repliseq_zoomout, repliseq_chrom, repliseq_start) %>%
    dplyr::group_by(repliseq_zoomout, repliseq_chrom) %>%
    dplyr::summarise(repliseq_prediction_start=seq(min(repliseq_start, na.rm=T), max(repliseq_start, na.rm=T), by=(binspersample-binsmargin)*binsize*repliseq_zoomout[1]), repliseq_prediction_end=repliseq_prediction_start+binsize*binspersample*repliseq_zoomout[1], .groups="keep")

  r_pred_df = prediction_regions_df %>%
    # dplyr::filter(repliseq_zoomout==1) %>%
    dplyr::group_by(repliseq_zoomout, repliseq_chrom, repliseq_prediction_start, repliseq_prediction_end) %>%
    dplyr::summarise(
      repliseq_start=min(repliseq_prediction_start, na.rm=T)+((1:binspersample)-1)*binsize*repliseq_zoomout[1],
      repliseq_end=repliseq_start+binsize*repliseq_zoomout[1],
      repliseq_rel_start=repliseq_start-repliseq_prediction_start,
      repliseq_rel_end=repliseq_end-repliseq_prediction_start,
      repliseq_rel_step=repliseq_rel_start/(binsize*repliseq_zoomout),
      .groups="keep") %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(r_df, by=c("repliseq_zoomout", "repliseq_chrom", "repliseq_start", "repliseq_end")) %>%
    dplyr::arrange(repliseq_zoomout, repliseq_chrom, repliseq_rel_step, repliseq_fraction) %>%
    dplyr::mutate(repliseq_rel_start=format(repliseq_rel_start, scientific=F, justify="none", trim=T)) %>%
    dplyr::mutate(repliseq_rel_start=factor(repliseq_rel_start, unique(repliseq_rel_start))) %>%
    dplyr::mutate(repliseq_rel_step=factor(repliseq_rel_step, unique(repliseq_rel_step))) %>%
    reshape2::dcast(repliseq_zoomout+repliseq_chrom+repliseq_prediction_start+repliseq_prediction_end+repliseq_fraction ~ repliseq_rel_step, value.var="repliseq_value") %>%
    dplyr::arrange(repliseq_zoomout, repliseq_chrom, repliseq_prediction_start, repliseq_prediction_end, dplyr::desc(repliseq_fraction)) %>%
    dplyr::mutate(prediction_id=paste0(repliseq_zoomout, "|", repliseq_chrom, ":", format(repliseq_prediction_start, scientific=F, justify="none", trim=T), "-", format(repliseq_prediction_end, scientific=F, justify="none", trim=T))) %>%
    dplyr::select(-repliseq_zoomout, -repliseq_chrom, -repliseq_prediction_start, -repliseq_prediction_end)

  r_pred_split_matrix = r_pred_df %>%
    dplyr::select(-prediction_id) %>%
    base::split(f=r_pred_df$prediction_id) %>%
    lapply(FUN=function(z) {
      z %>%
        `rownames<-`(NULL) %>%
        tibble::column_to_rownames("repliseq_fraction") %>%
        as.matrix()
    })

  r_pred_matrix = base::aperm(array(unlist(r_pred_split_matrix), dim=c(fractions, binspersample, length(r_pred_split_matrix))), c(3,2,1))
  dimnames(r_pred_matrix) = list(names(r_pred_split_matrix), 1:ncol(r_pred_matrix), as.character(fractions:1))

  r_pred_matrix
}

predict_model = function(model, repliseq_df, binsize=50e3, binspersample=30, binsmargin=6, fractions=16, threshold=0.5) {
  repliseq_matrix = prepare_data(repliseq_df, binsize=binsize, binsmargin=binsmargin, binspersample=binspersample, fractions=fractions)
  prediction_regions_df = as.data.frame(stringr::str_split(rownames(repliseq_matrix), pattern="[|:-]", simplify=T)) %>%
    dplyr::rename(prediction_zoomout="V1", prediction_chrom="V2", prediction_start="V3", prediction_end="V4") %>%
    dplyr::mutate(prediction_zoomout=as.numeric(prediction_zoomout), prediction_start=as.numeric(prediction_start), prediction_end=as.numeric(prediction_end))

  repliseqTime_df = repliseq_df %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::summarise(repliseq_mean_fraction=weighted.mean(repliseq_fraction, tidyr::replace_na(repliseq_value, 0)), .groups="keep") %>%
    dplyr::mutate(repliseq_pos=round(repliseq_start/2+repliseq_end/2))

  chromsizes_df = repliseq_df %>%
    dplyr::group_by(chromsize_chrom=repliseq_chrom) %>%
    dplyr::summarise(chromsize_length=max(repliseq_start))

  predictions_df = as.matrix(predict(model, repliseq_matrix)[,,,1], check.names=F) %>%
    `dimnames<-`(dimnames(repliseq_matrix)[1:2]) %>%
    reshape2::melt(value.name="tz_pred") %>%
    tidyr::separate(Var1, sep="[|:-]", into=c("tz_zoomout", "tz_chrom", "prediction_start", "prediction_end"), remove=F) %>%
    # dplyr::filter(grepl("(1|2)\\|chr6:60000000", Var1) | grepl("(1|2)\\|chr6:48000000", Var1)) %>%
    dplyr::mutate(
      tz_zoomout=as.numeric(tz_zoomout),
      tz_step=as.numeric(gsub("^V", "", as.character(Var2))),
      prediction_start=as.numeric(prediction_start),
      prediction_end=as.numeric(prediction_end),
      # tz_chrom=gsub(":.*", "", as.character(Var1)),
      # prediction_start=as.numeric(gsub(".* ", "", as.character(Var1))),
      tz_start=prediction_start+(tz_step-1)*binsize*tz_zoomout,
      tz_end=tz_start+binsize*tz_zoomout
    ) %>%
    dplyr::arrange(tz_zoomout, tz_chrom, tz_start, dplyr::desc(tz_pred)) %>%
    dplyr::select(tz_zoomout, tz_chrom, tz_start, tz_end, tz_step, tz_pred, tz_prediction_start=prediction_start) %>%
    dplyr::distinct(tz_zoomout, tz_chrom, tz_start, tz_end, .keep_all=T)

  predictions_significant_df = predictions_df %>%
    dplyr::filter(tz_pred>=threshold) %>%
    df2ranges(tz_chrom, tz_start, tz_end) %>%
    GenomicRanges::reduce() %>%
    as.data.frame() %>%
    dplyr::select(tz_chrom="seqnames", tz_start="start", tz_end="end") %>%
    df2ranges(tz_chrom, tz_start, tz_end) %>%
    innerJoinByOverlaps(predictions_df %>% dplyr::select(raw_zoomout=tz_zoomout, raw_chrom=tz_chrom, raw_start=tz_start, raw_end=tz_end, raw_pred=tz_pred) %>% df2ranges(raw_chrom, raw_start, raw_end)) %>%
    dplyr::arrange(raw_zoomout, dplyr::desc(raw_pred)) %>%
    dplyr::distinct(tz_chrom, tz_start, tz_end, .keep_all=T) %>%
    dplyr::mutate(tz_start=raw_start, tz_end=raw_end, tz_pred=raw_pred) %>%
    dplyr::select(dplyr::starts_with("tz_")) %>%
    dplyr::mutate(tz_pos=round(tz_start/2+tz_end/2))


  forks_df = dplyr::bind_cols(
    predictions_significant_df,
    predictions_significant_df[predictions_significant_df %>% df2ranges(tz_chrom, tz_pos, tz_pos) %>% GenomicRanges::follow(),] %>% dplyr::select(tz_left=tz_pos),
    predictions_significant_df[predictions_significant_df %>% df2ranges(tz_chrom, tz_pos, tz_pos) %>% GenomicRanges::precede(),] %>% dplyr::select(tz_right=tz_pos)
  ) %>%
    # dplyr::filter(tz_chrom=="chr13", tz_start>=119850000) %>%
    # dplyr::filter(tz_chrom=="chr4" & tz_start>=39.2e6 & tz_end<=40.6e6) %>%
    inner_join(chromsizes_df, by=c("tz_chrom"="chromsize_chrom")) %>%
    df2ranges(tz_chrom, ifelse(!is.na(tz_left), tz_left, 1), ifelse(!is.na(tz_right), tz_right, chromsize_length)) %>%
    innerJoinByOverlaps(repliseqTime_df %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end)) %>%
    dplyr::group_by(tz_chrom, tz_pos, tz_left, tz_right, tz_pred) %>%
    dplyr::mutate(is_left=repliseq_pos<tz_pos & repliseq_pos>tz_left, is_right=repliseq_pos>tz_pos & repliseq_pos<tz_right) %>%
    dplyr::summarise(
      tz_iz_left=ifelse(!is.na(tz_left[1]), repliseq_pos[is_left][which.min(repliseq_mean_fraction[is_left])], 1),
      tz_iz_right=ifelse(!is.na(tz_right[1]), repliseq_pos[is_right][which.min(repliseq_mean_fraction[is_right])], chromsize_length[1]),
      tz_iz_left=ifelse(!is.na(tz_iz_left), tz_iz_left, tz_left[1]/2+tz_pos[1]/2),
      tz_iz_right=ifelse(!is.na(tz_iz_right), tz_iz_right, tz_right[1]/2+tz_pos[1]/2)
    ) %>%
    # dplyr::mutate(tz_iz_left=round(tz_start/2+tz_left/2), tz_iz_right=round(tz_end/2+tz_right/2)) %>%
    dplyr::ungroup() %>%
    reshape2::melt(measure.vars=c("tz_iz_left", "tz_iz_right"), variable.name="fork_direction", value.name="fork_iz") %>%
    dplyr::mutate(
      fork_chrom=tz_chrom,
      fork_direction=dplyr::case_when(fork_direction=="tz_iz_left"~"right", T~"left"),
      fork_tz=tz_pos,
      fork_start=pmin(fork_iz, fork_tz),
      fork_end=pmax(fork_tz, fork_iz),
      fork_pred=tz_pred
    ) %>%
    dplyr::select(dplyr::starts_with("fork_")) %>%
    dplyr::arrange(fork_chrom, fork_tz, fork_start)

  list(all=predictions_df, significant=predictions_significant_df, forks=forks_df, regions=prediction_regions_df)
}

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

tzNN_prepare_training_data = function()
{
  debug = F
  training_margin=300e3
  training_width=1.5e6
  training_binsize = 5e4
  # training_step = 5e5
  # Coverage once
  training_step = 9e5

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
    chromsizes_df = repliseq_df %>%
      dplyr::group_by(chromsize_chrom=repliseq_chrom) %>%
      dplyr::summarise(chromsize_length=max(repliseq_start))

    #
    # Combine all datasets and calculate average fraction for each locus
    #
    repliseq_compare_df = dplyr::bind_rows(repliseq_df, repliseq1_df, repliseq2_df, repliseq3_df)
    repliseqTime_compare_df = repliseq_compare_df %>%
      dplyr::group_by(Run, repliseq_chrom, repliseq_start, repliseq_end) %>%
      dplyr::summarise(repliseq_mean_fraction=weighted.mean(repliseq_fraction, tidyr::replace_na(repliseq_value, 0)), .groups="keep") %>%
      dplyr::ungroup()

    repliseqTime_compare_df = repliseq_compare_df %>%
      dplyr::filter(Run=="Zhao_mESC16") %>%
      dplyr::arrange(Run, repliseq_chrom, repliseq_start, repliseq_end, repliseq_value) %>%
      dplyr::group_by(Run, repliseq_chrom, repliseq_start, repliseq_end) %>%
      dplyr::summarise(
        repliseq_mean_fraction=weighted.mean(
          repliseq_fraction,
          tidyr::replace_na(ifelse(repliseq_fraction>repliseq_fraction[which.max(repliseq_value)]+1, 0, repliseq_value), 0)), .groups="keep") %>%
      dplyr::ungroup()
    repliseqTime_compare_df %>%
      dplyr::ungroup() %>%
      dplyr::filter(Run=="Zhao_mESC16") %>%
      dplyr::select(repliseq_chrom, repliseq_start, repliseq_end, repliseq_mean_fraction) %>%
      readr::write_tsv("data/tz/tz_annotation_repliseq2.bedgraph", col_names=F)

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

    colfork_df = readr::read_tsv("data/zhao_colliding_forks.bed", col_names=c("colfork_chrom", "colfork_start", "colfork_end", "colfork_name")) %>%
      dplyr::mutate(colfork_id=paste0(colfork_chrom, ":", colfork_start, "-", colfork_end))

  #
  #
  # HERE
  #
  #
  find_tz_between_iz = function(colfork_ranges, iz_ranges, rfd_ranges, debug=F, loess_span=0.2) {
    #
    # Find termination zones
    #
    ret = iz_ranges %>%
      innerJoinByOverlaps(colfork_ranges) %>%
      dplyr::group_by(colfork_id, colfork_chrom, colfork_start, colfork_end) %>%
      dplyr::filter(dplyr::n()>=2) %>%
      dplyr::slice(c(1, dplyr::n())) %>%
      dplyr::summarize(colfork_min_pos=min(iz_start)+5e3, colfork_max_pos=max(iz_end)-5e3, .groups='drop') %>%
      dplyr::ungroup() %>%
      df2ranges(colfork_chrom, colfork_min_pos, colfork_max_pos) %>%
      innerJoinByOverlaps(rfd_ranges) %>%
      dplyr::filter(rfd_start>colfork_min_pos & rfd_end<colfork_max_pos) %>%
      dplyr::arrange(colfork_chrom, colfork_min_pos, colfork_max_pos, rfd_chrom, rfd_start) %>%
      dplyr::group_by(colfork_chrom, colfork_min_pos, colfork_max_pos) %>%
      dplyr::do((function(z) {
        zz <<- z
        # z = tz_df %>% dplyr::filter(grepl("99430908", colfork_id))

        z = zz
        rfd_model = loess(rfd_score ~ rfd_start, data=z, span=loess_span)
        rfd_pred_df = data.frame(rfd_start=seq(min(z$rfd_start), max(z$rfd_start), 1e3)) %>%
          dplyr::mutate(rfd_chrom=z$colfork_chrom[1], rfd_end=rfd_start+1e3, rfd_score = predict(rfd_model, .))
        rfd_pred_significant_df = rfd_pred_df %>%
          dplyr::filter(abs(rfd_score) < 0.01)
        if(nrow(rfd_pred_significant_df)==0) return(data.frame())
        rfd_pred_tz_df = rfd_pred_significant_df %>%
          df2ranges(rfd_chrom, rfd_start, rfd_end) %>%
          GenomicRanges::reduce(min.gapwidth=1e4) %>%
          as.data.frame() %>%
          dplyr::select(rfd_reduced_chrom=seqnames, rfd_reduced_start=start, rfd_reduced_end=end) %>%
          dplyr::mutate(rfd_reduced_length=200e3) %>%
          dplyr::group_by(rfd_reduced_chrom, rfd_reduced_start, rfd_reduced_end) %>%
          dplyr::summarise(
            colfork_min_pos=z$colfork_min_pos[1],
            colfork_max_pos=z$colfork_max_pos[1],
            rfd_region_offset=-1:1,
            rfd_region_start=c(pmax(z$colfork_min_pos[1], rfd_reduced_start-rfd_reduced_length), rfd_reduced_start, rfd_reduced_end),
            rfd_region_end=c(rfd_reduced_start, rfd_reduced_end, pmin(z$colfork_max_pos[1], rfd_reduced_end+rfd_reduced_length)), .groups='drop') %>%
          dplyr::ungroup() %>%
          data.frame() %>%
          df2ranges(rfd_reduced_chrom, rfd_region_start, rfd_region_end) %>%
          innerJoinByOverlaps(rfd_pred_df %>% df2ranges(rfd_chrom, rfd_start, rfd_end)) %>%
          dplyr::group_by(tz_chrom=rfd_chrom, rfd_reduced_start, rfd_reduced_end) %>%
          dplyr::mutate(
            rfd_left_score1=max(rfd_score[rfd_region_offset==-1], na.rm=T),
            rfd_right_score1=-min(rfd_score[rfd_region_offset==1], na.rm=T),
            rfd_left_score2=sum(rfd_score[rfd_region_offset==-1]>0, na.rm=T),
            rfd_right_score2=sum(rfd_score[rfd_region_offset==1]<0, na.rm=T)) %>%
          dplyr::filter(rfd_region_offset==0) %>%
          # dplyr::summarise(tz_start=min(rfd_start, na.rm=T), tz_end=max(rfd_end, na.rm=T), tz_iz_left=z$colfork_min_pos[1], tz_iz_right=z$colfork_max_pos[1], rfd_trend_score=rfd_right_score[1]-rfd_middle_score[1], .groups='drop') %>%
          dplyr::summarise(
            tz_start=weighted.mean(rfd_start, w=max(abs(rfd_score))-abs(rfd_score)),
            tz_end=tz_start, tz_iz_left=z$colfork_min_pos[1],
            tz_iz_right=z$colfork_max_pos[1],
            rfd_trend_score1=rfd_right_score1[1]+rfd_left_score1[1],
            rfd_trend_score2=rfd_right_score2[1]+rfd_left_score2[1],
            .groups='drop') %>%
          # dplyr::summarise(tz_start=rfd_start[which.min(abs(rfd_score))], tz_end=tz_start, tz_iz_left=z$colfork_min_pos[1], tz_iz_right=z$colfork_max_pos[1], rfd_trend_score=rfd_right_score[1]-rfd_middle_score[1], .groups='drop') %>%
          dplyr::arrange(dplyr::desc(rfd_trend_score1), dplyr::desc(rfd_trend_score2)) %>%
          dplyr::slice(1) %>%
          dplyr::mutate(tz_start=tz_start-training_binsize/4, tz_end=tz_end+training_binsize/4) %>%
          dplyr::select(dplyr::starts_with("tz_"))

        if(debug) {
          print(ggplot(z, aes(x=rfd_start, y=rfd_score)) +
            geom_line() +
            geom_area(aes()) +
            geom_line(aes(x=rfd_start, y=rfd_score), data=rfd_pred_df) +
            geom_segment(aes(x=tz_start, xend=tz_end, y=0, yend=0), data=rfd_pred_tz_df, size=3, color="#FF0000") +
            geom_smooth(method="loess", span=loess_span) +
            coord_cartesian(ylim=c(-max(abs(z$rfd_score), na.rm=T), max(abs(z$rfd_score), na.rm=T))))
        }

        if(nrow(rfd_pred_tz_df)>1) {
            return(data.frame())
        }

        rfd_pred_tz_df %>% dplyr::select(dplyr::matches("tz_")) %>% dplyr::mutate(tz_alternative_n=nrow(rfd_pred_tz_df))
        # r = data.frame(tz_chrom=z$colfork_chrom[1], tz_start=GenomicRanges::start(rfd_pred_ranges), tz_end=GenomicRanges::end(rfd_pred_ranges), tz_alternative_n=nrow(rfd_pred_df))
      })(.)) %>%
      dplyr::ungroup() %>%
      dplyr::select(dplyr::starts_with("tz_"))
  }

    # data.frame(colfork_id=c("chr2:31,211,044-31,650,066")) %>%
    #   dplyr::mutate(colfork_id=gsub(",", "", colfork_id)) %>%
    #   tidyr::separate(colfork_id, into=c("colfork_chrom", "colfork_start", "colfork_end"), sep="[:-]", remove=F) %>%
    #   dplyr::mutate(colfork_start=as.numeric(colfork_start), colfork_end=as.numeric(colfork_end)) %>%
    #   df2ranges(colfork_chrom, colfork_start, colfork_end) %>%
    #   find_tz_between_iz(iz_ranges, rfd_ranges, loess_span=0.2, debug=T)

    tz_step1_df = colfork_df %>%
      df2ranges(colfork_chrom, colfork_start, colfork_end) %>%
      find_tz_between_iz(iz_ranges, rfd_ranges, debug=F)

    #
    # Adjust final TZ position to match Repliseq data (+/- 300kb)
    #
    max_adjustment = 1e6
    tz_df = tz_step1_df %>%
      dplyr::mutate(tz_region_left=round(pmax(tz_iz_left, tz_start-max_adjustment)), tz_region_right=round(pmin(tz_iz_right, tz_end+max_adjustment))) %>%
      dplyr::mutate(tz_start=ceiling(tz_start/training_binsize)*training_binsize, tz_end=ceiling(tz_end/training_binsize)*training_binsize) %>%
      df2ranges(tz_chrom, tz_region_left, tz_region_right) %>%
      innerJoinByOverlaps(repliseqTime_compare_df %>% dplyr::filter(Run=="Zhao_mESC16") %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end)) %>%
      dplyr::group_by(tz_chrom, tz_start, tz_end, tz_iz_left, tz_iz_right, tz_alternative_n) %>%
      # dplyr::mutate(tz_peak_start=(repliseq_start+repliseq_end)[which.max(repliseq_mean_fraction)]/2) %>%
      # dplyr::filter(repliseq_mean_fraction>=max(repliseq_mean_fraction)*0.99) %>%
      dplyr::do((function(tz) {
          tzz<<-tz

        # tz = tza_df %>% dplyr::filter(tz_chrom=="chr1" & tz_region_right==145402000)
        writeLines(with(tz[1,], paste0(tz_chrom, ":", tz_iz_left, "-", tz_iz_right, " => ", nrow(tz))))
          #print(nrow(tz))
          # asdaD()

        if(nrow(tz)>10) {
          tz_model = loess(repliseq_mean_fraction ~ repliseq_start, data=tz, span=0.6)
          tz = tz %>% dplyr::mutate(repliseq_mean_fraction=predict(tz_model, data=.))
        }
        tz_prepared = tz %>%
          dplyr::mutate(
            tuning_direction=dplyr::case_when(tz_start>=repliseq_start & tz_end<=repliseq_end~0, tz_start-repliseq_end>=0~-1, repliseq_start-tz_end>=0~1, T~0),
            repliseq_mean_fraction_diff=dplyr::case_when(tuning_direction==0~0, tuning_direction==1~repliseq_mean_fraction-dplyr::lag(repliseq_mean_fraction), T~repliseq_mean_fraction-dplyr::lead(repliseq_mean_fraction)),
            repliseq_mean_fraction_difftype=dplyr::case_when(tuning_direction==0~"increase",repliseq_mean_fraction_diff>0~"increase", T~"decrease")
          ) %>%
          dplyr::select(tz_chrom, tz_start, tz_end, repliseq_start, repliseq_end, tuning_direction, repliseq_mean_fraction, repliseq_mean_fraction_diff, repliseq_mean_fraction_difftype)

        tuning_shift = tz_prepared %>%
          dplyr::group_by(tz_chrom, tz_start, tz_end, tuning_direction) %>%
          dplyr::do((function(z) {
            zz<<-z
            # if(z$tuning_direction[1]==0) asdasD()
            if(z$tuning_direction[1]==-1) z = z[nrow(z):1,]
            z.rle = base::rle(z$repliseq_mean_fraction_difftype)
            data.frame(tuning_distance=z.rle$lengths, repliseq_mean_fraction_difftype=z.rle$values, repliseq_mean_fraction=sapply(split(z$repliseq_mean_fraction, rep(1:length(z.rle$lengths),  times=z.rle$lengths)), max))
          })(.)) %>%
          dplyr::filter(tuning_direction==0 | dplyr::n()>=2 & repliseq_mean_fraction_difftype[1]=="increase" & repliseq_mean_fraction_difftype[2]=="decrease") %>%
          dplyr::slice(1) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(dplyr::desc(repliseq_mean_fraction)) %>%
          dplyr::slice(1)

        if(nrow(tuning_shift)==0) {
          stop("wrong")
          return(tz %>% dplyr::distinct(tz_peak_start=tz_start, tz_peak_end=tz_end) %>% dplyr::mutate(tz_is_tuned=F))
        }

        tz_tuned = tuning_shift %>%
          dplyr::mutate(tz_peak_start=tz_start+tuning_direction*tuning_distance*training_binsize, tz_peak_end=tz_end+tuning_direction*tuning_distance*training_binsize, tz_is_tuned=T) %>%
          dplyr::select(tz_peak_chrom=tz_chrom, tz_peak_start, tz_peak_end, tz_is_tuned, repliseq_mean_fraction)
        tz_final = tz %>%
          dplyr::filter(
            tuning_shift$tuning_direction<0 & repliseq_end>=(tz_start-tuning_shift$tuning_distance*training_binsize) & repliseq_start<=tz_end |
            tuning_shift$tuning_direction>0 & repliseq_start<=(tz_end+tuning_shift$tuning_distance*training_binsize) & repliseq_end>=tz_start |
            tuning_shift$tuning_direction==0 & repliseq_end>=tz_start & repliseq_start<=tz_end) %>%
          dplyr::filter(tz_tuned$repliseq_mean_fraction-repliseq_mean_fraction<=0.01 & tz_tuned$repliseq_mean_fraction-repliseq_mean_fraction>=0) %>%
          df2ranges(tz_chrom, repliseq_start, repliseq_end) %>%
          GenomicRanges::reduce() %>%
          as.data.frame() %>%
          dplyr::select(tz_chrom=seqnames, tz_start=start, tz_end=end) %>%
          df2ranges(tz_chrom, tz_start, tz_end) %>%
          innerJoinByOverlaps(tz_tuned %>% df2ranges(tz_peak_chrom, tz_peak_start, tz_peak_end)) %>%
          as.data.frame() %>%
          dplyr::mutate(tz_peak_start=tz_start, tz_peak_end=tz_end)

        if(nrow(tz_final)==0) {
          stop("wrong")
        }

        # ggplot(tz) +
        #   geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=repliseq_start, xmax=repliseq_end, fill=as.factor(tuning_direction)), data=tz_prepared) +
        #   # geom_line(aes(x=repliseq_start+training_binsize/2, y=6*(repliseq_mean_fraction_difftype=="increase"), color="difftype"), data=tz_prepared) +
        #   # geom_line(aes(x=repliseq_start+training_binsize/2, y=repliseq_mean_fraction_diff, color="diff"), data=tz_prepared) +
        #   geom_line(aes(x=repliseq_start+training_binsize/2, y=repliseq_mean_fraction), color="#CC0033") +
        #   geom_point(aes(x=repliseq_start+training_binsize/2, y=repliseq_mean_fraction), color="#CC0033") +
        #   geom_rect(aes(xmin=tz_peak_start, xmax=tz_peak_end, ymin=11.9, ymax=12), data=tz_tuned) +
        #   geom_rect(aes(xmin=tz_peak_start, xmax=tz_peak_end, ymin=12.1, ymax=12.2), data=tz_final) +
        #   geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=tz_start, xmax=tz_end, color="Original"), alpha=0.4, data=tz %>% dplyr::distinct(tz_start, tz_end)) +
        #   geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=tz_peak_start, xmax=tz_peak_end, color="Tuned"), alpha=0.4, fill="#FFFFFF00", data=tz_final) +
        #   geom_hline(yintercept=0) +
        #   theme_gray()

        # print(1)
        #   geom_line(aes(x=repliseq_start, y=repliseq_mean_fraction_diff*32, color="diff")) +
        #   geom_line(aes(x=repliseq_start, y=repliseq_mean_fraction_increase*16, color="increase")) +
        #   geom_line(aes(x=repliseq_start, y=shift_direction*16, color="shift"))
        # diff(tz$repliseq_mean_fraction)
        #
        #   tz_ranges = tz %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end)
        #   tz_reduced_ranges = tz_ranges %>%
        #     GenomicRanges::reduce() %>%
        #     as.data.frame() %>%
        #     # dplyr::mutate(tz_repliseq_peak_pos=(tz$repliseq_start+tz$repliseq_end)[which.max(tz$repliseq_mean_fraction)]/2, tz_extended_length=end-start) %>%
        #     dplyr::select(tz_extended_chrom=seqnames, tz_extended_start=start, tz_extended_end=end) %>%
        #     df2ranges(tz_extended_chrom, tz_extended_start, tz_extended_end)
        #   tz_reduced_ranges %>%
        #     innerJoinByOverlaps(tz_ranges)  %>%
        #     dplyr::filter(tz_extended_start<=tz_peak_start & tz_peak_start<=tz_extended_end) %>%
        #     dplyr::group_by(tz_extended_chrom, tz_extended_start, tz_extended_end) %>%
        #     dplyr::summarise(tz_peak_start=weighted.mean(repliseq_start/2+repliseq_end/2, repliseq_mean_fraction)-1, tz_peak_end=tz_peak_start+2, .groups="keep")

            # dplyr::select(tz_extended_start=start, tz_extended_end=end, tmp_corected_pos, tz_ragion_length)
        tz_final
      })(.)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(tz_okseq_start=tz_start, tz_okseq_end=tz_end, tz_start=tz_peak_start, tz_end=tz_peak_end)
      #dplyr::filter(tz_extended_start<=tz_peak_start & tz_peak_start<=tmp_corrected_end) %>%
      #dplyr::mutate(tz_start=tmp_corrected_start/2+tmp_corrected_end/2, tz_end=tmp_corrected_start/2+tmp_corrected_end/2) %>%
      #dplyr::select(dplyr::starts_with("tz_"))
    #table(tza_df$tz_extended_end-tza_df$tz_extended_start)/nrow(tza_df)

    readr::write_tsv(tz_df, file="data/tz/tz_annotation.tsv")

    tz_step1_df %>%
      dplyr::mutate(tz_start=round(tz_start), tz_end=round(tz_end)) %>%
      dplyr::select(tz_chrom, tz_start, tz_end) %>%
      readr::write_tsv(file="data/tz/tz_annotation_step1.bed", col_names=F)
    tz_step1_df %>%
      dplyr::mutate(tz_start=round(tz_start), tz_end=round(tz_end)) %>%
      dplyr::mutate(thinStart=as.integer(tz_iz_left), thinEnd=as.integer(tz_iz_right), name="TZ", strand="*", score=0, color="0,0,0") %>%
      dplyr::select(tz_chrom, thinStart, thinEnd, name, score,strand, tz_start, tz_end, color) %>%
      readr::write_tsv(file="data/tz/tz_annotation_step1ext.bed", col_names=F)
    tz_df %>%
      dplyr::mutate(tz_start=as.integer(tz_start), tz_end=as.integer(tz_end)) %>%
      dplyr::select(tz_chrom, tz_start, tz_end) %>%
      readr::write_tsv(file="data/tz/tz_annotation_step2.bed", col_names = F)
    tz_df %>%
      dplyr::mutate(tz_start=as.integer(tz_start), tz_end=as.integer(tz_end)) %>%
      dplyr::mutate(thinStart=as.integer(pmax(tz_iz_left, tz_start-max_adjustment)), thinEnd=as.integer(pmin(tz_iz_right, tz_end+max_adjustment)), name="TZ", strand="*", score=0, color="0,0,0") %>%
      dplyr::select(tz_chrom, thinStart, thinEnd, name, score,strand, tz_start, tz_end, color) %>%
      readr::write_tsv(file="data/tz/tz_annotation_step2ext.bed", col_names = F)

    #
    # Split into training regions
    #
    training_long_regions_df = as.data.frame(GenomicRanges::reduce(colfork_df %>% df2ranges(colfork_chrom, colfork_start, colfork_end), min.gapwidth=500e3)) %>%
      dplyr::select(-strand) %>%
      dplyr::rename(colfork_reduced_chrom="seqnames", colfork_reduced_start="start", colfork_reduced_end="end", colfork_reduced_width="width")
    training_starts_df = training_long_regions_df %>%
      dplyr::inner_join(chromsizes_df, by=c("colfork_reduced_chrom"="chromsize_chrom")) %>%
      # dplyr::filter(colfork_reduced_chrom=="chr3" & colfork_reduced_start>=96e6 & colfork_reduced_start<=100e6) %>%
      dplyr::mutate(seq_start=colfork_reduced_start-training_margin, seq_end=pmin(colfork_reduced_end, chromsize_length-training_margin)-training_width+training_margin) %>%
      dplyr::mutate(seq_start=ceiling(seq_start/training_binsize)*training_binsize, seq_end=ceiling(seq_end/training_binsize)*training_binsize)
    training_regions_df = training_starts_df %>%
      dplyr::filter(seq_end>=seq_start) %>%
      dplyr::group_by(colfork_reduced_chrom, colfork_reduced_start, colfork_reduced_end) %>%
      dplyr::summarize(training_chrom=colfork_reduced_chrom[1], training_start=unique(c(seq(seq_start, seq_end, by=training_step), seq_end)), training_end=training_start+training_width) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(colfork_id=paste0(colfork_reduced_chrom, ":", colfork_reduced_start, "-", colfork_reduced_end), training_id=paste0(training_chrom, "_", format(training_start, scientific=F, justify="none", trim=T), "-", format(training_end, scientific=F, justify="none", trim=T)))


    #
    # Export intermediate regions
    #
    training_starts_df %>%
      dplyr::select(colfork_reduced_chrom, seq_start, seq_end) %>%
      readr::write_tsv(file="data/tz/training_starts.bed", col_names = F)
    training_long_regions_df %>%
      dplyr::select(colfork_reduced_chrom, colfork_reduced_start, colfork_reduced_end) %>%
      readr::write_tsv(file="data/tz/training_long_regions.bed", col_names = F)
    training_regions_df %>%
      dplyr::group_by(training_chrom, training_start, training_end) %>%
      dplyr::summarise(name="training", strand="*", score=0, region_start=training_start, region_end=training_end, thickStart=training_start+training_margin, thickEnd=training_end-training_margin) %>%
      dplyr::ungroup() %>%
      dplyr::select(training_chrom, region_start, region_end, name, score, strand, thickStart, thickEnd) %>%
      readr::write_tsv(file="data/tz/training_regions.bed", col_names = F)



    #
    # Create training and label data.frames
    #
    data_df = training_regions_df %>%
      df2ranges(training_chrom, training_start, training_end) %>%
      innerJoinByOverlaps(repliseq_compare_df %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end), minoverlap=training_binsize) %>%
      dplyr::group_by(Run, training_chrom, training_start, training_end) %>%
      dplyr::mutate(repliseq_rel_start=repliseq_start-training_start) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(Run, colfork_id, training_id, repliseq_chrom, repliseq_start, repliseq_end) %>%
      dplyr::mutate(repliseq_value=repliseq_value/max(repliseq_value)) %>%
      dplyr::ungroup()

    #
    # TODO: remove and fix DISTINCT TZ!!!!
    #
    label_df = training_regions_df %>%
      dplyr::group_by(colfork_id, training_id, training_chrom, training_start, training_end) %>%
      dplyr::summarize(training_tz_start=seq(training_start, training_end-1, training_binsize), training_tz_end=training_tz_start+training_binsize, .groups='drop') %>%
      dplyr::mutate(training_tz_rel_start=training_tz_start-training_start) %>%
      df2ranges(training_chrom, training_tz_start, training_tz_end) %>%
      leftJoinByOverlaps(tz_df %>% df2ranges(tz_chrom, tz_start, tz_end)) %>%
      dplyr::distinct(training_id, training_tz_rel_start, .keep_all=T) %>%
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
        debug_tz_df = training_regions_df %>%
          df2ranges(training_chrom, training_start, training_end) %>%
          innerJoinByOverlaps(tz_df %>% df2ranges(tz_chrom, tz_start, tz_end)) %>%
          dplyr::group_by(training_chrom, training_start, training_end) %>%
          dplyr::mutate(tz_rel_start=tz_start-training_start, tz_rel_end=tz_end-training_start) %>%
          dplyr::ungroup()

        debug_rfd_df = training_regions_df %>%
          df2ranges(training_chrom, training_start, training_end) %>%
          innerJoinByOverlaps(rfd_ranges) %>%
          dplyr::group_by(training_chrom, training_start, training_end) %>%
          dplyr::mutate(rfd_rel_start=rfd_start-training_start) %>%
          dplyr::ungroup()

        debug_repliseqTime_df = training_regions_df %>%
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
        example_training_ids = training_regions_df$training_id[c(1,500,1000,1500)]
        example_training_ids = with(training_regions_df, training_id[grepl("chr18_628", training_id)])
        example_training_ids = with(training_regions_df, training_id[grepl("chr4:31900000|chr4_36000000|chr4_33700000|chr4_34600000|chr4_35500000", training_id)])


        pdf("reports/tzNN_examples.pdf", width=8.27, height=8.27)
        for(example_training_ids in with(training_regions_df, training_id[grepl("chr1_183550000|chr1_88850000|chr1_134350000|chr1_117300000", training_id)]))
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
              #     dplyr::filter(validation_k==1 & type=="pred" & source=="validation") %>%
              #     dplyr::mutate(x=Var2*training_binsize, colfork_id=training_id, set="Prediction") %>%
              #     dplyr::select(set, training_id, colfork_id, x, value) %>%
              #     tidyr::crossing(y=1:2),
              label_df %>%
                tidyr::crossing(data.frame(y=1)) %>%
                dplyr::mutate(set="Labels", x=training_tz_rel_start, value=training_tz_present)) %>%
              dplyr::filter(training_id %in% example_training_ids)

            p = ggplot(data_ggplot) +
                geom_tile(aes(x=as.numeric(x)+training_binsize/2, y=y, fill=value)) +
                geom_tile(aes(x=as.numeric(x)+training_binsize/2, y=y), fill="#CCCCCC", data=data.frame(x=c(0:3, 16:19)*50000, y=1, value=1, set="Labels", training_id=example_training_ids)) +
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
          `[`(training_regions_df$training_id)
      })

    label_matrix = label_df %>%
      dplyr::mutate(training_tz_present=ifelse(training_margin<=training_tz_rel_start & training_tz_rel_start<training_width-training_margin, training_tz_present, 0)) %>%
      dplyr::mutate(training_tz_rel_start=format(training_tz_rel_start, scientific=F, justify="none", trim=T), training_tz_rel_start=factor(training_tz_rel_start, unique(training_tz_rel_start))) %>%
      reshape2::dcast(training_id ~ training_tz_rel_start, value.var="training_tz_present", drop=F) %>%
      split(.$training_id) %>%
      lapply(FUN=function(z) z %>% dplyr::select(-training_id) %>% as.matrix()) %>%
      `[`(training_regions_df$training_id)

  table(sapply(label_matrix, function(z) colnames(z)[which.max(z)]))[order(as.numeric(names(table(sapply(label_matrix, function(z) colnames(z)[which.max(z)])))))]
  table(sapply(data_matrix$Zhao_mESC16, function(z) any(is.na(z))))

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
      i = i+1
      dplyr::bind_rows(
        reshape2::melt(data_matrix$Zhao_mESC16[[i]]) %>% dplyr::mutate(set="Data"),
        reshape2::melt(label_matrix[[i]]) %>% dplyr::mutate(set="Labels", value=ifelse(training_margin<=Var2 & Var2<training_width-training_margin, value, 0.5))) %>%
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
  #reticulate::virtualenv_create("r-tensorflow", python="/usr/bin/python3.8")
  #reticulate::py_install("tensorflow-gpu")
  # tensorflow::install_tensorflow(envname="r-tensorflow")
  #keras::install_keras(envname="r-tensorflow", tensorflow="gpu")

  # gpu_options <- tf$GPUOptions(per_process_gpu_memory_fraction = 0.3)
  # config <- tf$ConfigProto(gpu_options = gpu_options)
  # k_set_session(tf$Session(config = config))

  model_dir = "reports/09-replication_fork_nn/models_new6/"
  dir.create(model_dir, recursive=T, showWarnings=F)

  load(file="data/tz/data_matrix.rda")
  load("data/tz/label_matrix.rda")

  x_width = ncol(data_matrix$Zhao_mESC16[[1]])
  x_height = nrow(data_matrix$Zhao_mESC16[[1]])

 #data_matrix$Zhao_mESC16[sapply(data_matrix$Zhao_mESC16, function(z) any(is.na(z)))]

  x = aperm(array(unlist(data_matrix$Zhao_mESC16), dim=c(x_height, x_width, length(data_matrix$Zhao_mESC16))), c(3,2,1))
  y = aperm(array(as.double(unlist(label_matrix)), dim=c(1, x_width, length(data_matrix$Zhao_mESC16))), c(3,2,1))

  set.seed(123)
  crossvalidation = split(sample(1:length(label_matrix), length(label_matrix)), cut(1:length(label_matrix), 10))
  df_all = data.frame()
  for(chr in 0:length(crossvalidation)) {
      # sample_validation = grep(paste0("^", chr, "_"), names(label_matrix))
      if(chr==0) {
        sample_validation = numeric(0)
      } else {
        sample_validation = crossvalidation[[chr]]
      }
      samples_training = setdiff(1:nrow(x), sample_validation)
      print(paste0("Validation proportion (", chr, "): ", round(100*length(sample_validation)/(length(sample_validation)+length(samples_training)), 2), "%"))
      x_validate = x[sample_validation,,]
      y_validate = array(y[sample_validation,,], dim=c(length(sample_validation), x_width, 1))
      x_train = x[samples_training,,]
      y_train = array(y[samples_training,,], dim=c(length(samples_training), x_width, 1))
      x_train = abind::abind(x_train, pretrube_mirror(x_train), along=1)
      y_train = abind::abind(y_train, pretrube_mirror(y_train), along=1)
      x_weights = as.matrix(sapply(1:nrow(x_train), function(i) { ii<<-i
        m = max(apply(x_train[i,,], 1, function(w) weighted.mean(1:length(w), w))*y_train[i,,1])
          dplyr::case_when(m<4~0.5, m<10~0.7, T~1)
        }))

      set.seed(123)
      args = list(batch_size=32, learning_rate=1e-3, dropout1=0.9, dropout2=0.2)
      model = create_model(input_shape=c(NA_integer_, dim(x_train)[3], 1),  dropout1=args$dropout1, dropout2=args$dropout2) %>%
        keras::compile(
            optimizer=keras::optimizer_adam(learning_rate=args$learning_rate),
            loss=loss_dice,
            #loss=keras::loss_binary_crossentropy(),
            #loss=keras::loss_poisson(),
            metrics=list(metric_dice, keras::metric_binary_crossentropy, keras::metric_poisson)
          )

      model_math = paste0("logs/", format(Sys.time(), "%m-%d_%H%M"), "-", stringr::str_glue("WEIGHTED_3x5_doA{do1}_doB{do2}_bs{bs}_lr{lr}_DICE_cv{chr}", chr=chr, do1=args$dropout1, do2=args$dropout2, bs=args$batch_size, lr=args$learning_rate))
      keras::fit(model, x_train, y_train, validation_data=list(x_val=x_validate, y_val=y_validate),
        batch_size=args$batch_size, sample_weight=x_weights,  epochs=10000, verbose=1, shuffle=T,
        callbacks=list(keras::callback_tensorboard(model_math, histogram_freq=5, embeddings_freq=5), keras::callback_progbar_logger(count_mode="steps"), keras::callback_early_stopping(monitor="val_loss", min_delta=0.000001, patience=1000, restore_best_weights=T)))

      save_model(model, file.path(gsub("/$", "", model_dir), paste0("cv", chr)))

      df = dplyr::bind_rows(
        as.matrix(as.data.frame(predict(model, x), row.names=names(data_matrix$Zhao_mESC16), check.names=F)) %>%
          reshape2::melt() %>%
          dplyr::mutate(Var2=as.numeric(gsub("^V", "", as.character(Var2)))) %>%
          dplyr::mutate(training_id=as.character(Var1), type="pred"),
        as.matrix(as.data.frame(y, row.names=names(data_matrix$Zhao_mESC16), check.names=F)) %>%
          reshape2::melt() %>%
          dplyr::mutate(Var2=as.numeric(gsub("^V", "", as.character(Var2))))  %>%
          dplyr::mutate(training_id=as.character(Var1), type="true")
      ) %>%
        dplyr::mutate(source=dplyr::case_when(training_id %in% names(data_matrix$Zhao_mESC16)[samples_training] ~ "training", T~"validation")) %>%
        dplyr::mutate(value=pmax(0, value))
      df_all = dplyr::bind_rows(df_all, df %>% dplyr::mutate(validation_k=chr))

      pdf(paste0(model_dir, "/evaluation.pdf"), width=8.27, height=11.69)
      plot(model)
      if(dim(x_validate)[1]>0) {
          print(model %>% evaluate(x_validate, y_validate))
          predictions_val = predict(model, x_validate)
          rockr_val_pred = ROCR::prediction(as.numeric(predictions_val), as.numeric(y_validate))
          rockr_val_acc = ROCR::performance(rockr_val_pred, measure = "prec", x.measure="rec")
          rockr_val_perf = ROCR::performance(rockr_val_pred, measure = "tpr", x.measure = "fpr")
          plot(rockr_val_acc, colorize=T, main="ROC curve (validation data)")
          plot(rockr_val_perf, colorize=T, main="Precision/recall curve (validation data)")
          abline(a=0, b=1)
      }

      predictions_train = predict(model, x_train)
      rockr_train_pred = ROCR::prediction(as.numeric(predictions_train), as.numeric(y_train))
      rockr_train_acc = ROCR::performance(rockr_train_pred, measure = "prec", x.measure="rec")
      rockr_train_perf = ROCR::performance(rockr_train_pred, measure = "tpr", x.measure = "fpr")
      plot(rockr_train_acc, colorize=T, main="ROC curve (training data)")
      plot(rockr_train_perf, colorize=T, main="Precision/recall curve (training data)")
      abline(a=0, b=1)

      df_random = df %>%
        dplyr::group_by(source, training_id) %>%
        dplyr::filter(any(value>=0.1))%>%
        dplyr::ungroup() %>%
        dplyr::group_by(source) %>%
        dplyr::filter(training_id %in% unique(training_id)[sample(length(unique(training_id)), pmin(length(unique(training_id)), 50))])  %>%
        dplyr::ungroup() %>%
        dplyr::distinct(Var2, training_id, type, source, value)
      df2 = df_random %>%
        tidyr::crossing(data.frame(shift=c(-0.4001, -0.4, 0.4, 0.4001), value_mod=c(0, 1, 1, 0))) %>%
        dplyr::mutate(Var2=Var2+shift, value=value*value_mod, diff=0) %>%
        dplyr::arrange(Var2, training_id, type, source, value) %>%
        dplyr::bind_rows(df_random)
      ggplot(df2) +
        ggridges::geom_ridgeline(aes(x=Var2, y=training_id, height=value, color=type), alpha=0.0, scale=0.5) +
        facet_wrap(~source, scales="free_y") +
        scale_color_manual(values=c("pred"="#333333", "true"="#FF0000")) +
        theme_paper()
      dev.off()
  }

  df_all %>%
    readr::write_tsv(file.path(gsub("/$", "", model_dir), "crossvalidation_data_2022-10-30.tsv"))
}

if (!interactive()) {
  tzNN_train()
}

tzNN_evaluate = function() {
  df_all = readr::read_tsv("reports/09-replication_fork_nn/crossvalidation_data_2022-10-30.tsv")

  pdf("reports/09-replication_fork_nn/crossvalidation_performance_2022-10-30.pdf", width=8.27, height=8.27)
  rockr_all_data = df_all %>%
    dplyr::arrange(source, validation_k, training_id, Var2) %>%
    dplyr::filter(source=="validation" & validation_k!=0) %>%
    reshape2::dcast(source+validation_k+training_id+Var2~type, value.var="value") %>%
    dplyr::group_by(source, validation_k, training_id) %>%
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
  # rockr_all_pred = ROCR::prediction(rockr_all_data %>% dplyr::group_split(validation_k) %>% lapply(function(l) l$pred), rockr_all_data %>% dplyr::group_split(validation_k) %>% lapply(function(l) l$true))
  rockr_all_acc = ROCR::performance(rockr_all_pred, measure = "prec", x.measure="rec")
  rockr_all_perf = ROCR::performance(rockr_all_pred, measure = "tpr", x.measure = "fpr")
  plot(rockr_all_acc, colorize=T, main="Precision/recall curve (validation data)")
  plot(rockr_all_perf, colorize=T, main="ROC curve (validation data)")
  abline(a=0, b=1)
  dev.off()
}

tzNN_predict = function() {
  reticulate::virtualenv_python("r-tensorflow")
  model = load_model("reports/09-replication_fork_nn/models_new6/cv0")

  celltype = "NPC"
  repliseq_df = readr::read_tsv(paste0("~/Workspace/Datasets/zhao_bmc_repliseq_2020/results/zhao_m", celltype, "_repliseq50000.tsv"))
  p = predict_model(model, repliseq_df, binsize=50e3, threshold=0.5)

  #
  # Export predictions
  #
  p$forks %>%
    dplyr::mutate(
      replication_chrom=fork_chrom,
      replication_strand=dplyr::case_when(fork_direction=="right"~"+", T~"-"),
      replication_direction=dplyr::case_when(fork_direction=="right"~"telomeric", T~"centromeric"),
      replication_length=fork_end-fork_start,
      replication_start=dplyr::case_when(fork_direction=="right"~fork_start, T~fork_end),
      replication_end=dplyr::case_when(fork_direction=="right"~fork_end, T~fork_start)) %>%
    dplyr::select(dplyr::matches("replication_")) %>%
    readr::write_tsv("data/replication_reduced_subsets_predicted.tsv")

  #
  # Export predictions (BED)
  #
  path_suffix = "2022-11-01"
  fork_bed_path = paste0("reports/09-replication_fork_nn/forks-", celltype, "-", path_suffix, ".bed")
  writeLines('track itemRgb=On visibility=2 colorByStrand="255,0,0 0,0,255"', con=fork_bed_path)
  p$forks %>%
    dplyr::mutate(fork_strand=dplyr::case_when(fork_direction=="right"~"+", fork_direction=="left"~"-"), thickStat=fork_start, thickEnd=fork_end, fork_color=dplyr::case_when(fork_direction=="right"~"#FF0000", fork_direction=="left"~"#0000FF")) %>%
    dplyr::select(fork_chrom, fork_start, fork_end, fork_direction, fork_pred, fork_strand, thickStat, thickEnd, fork_color) %>%
    readr::write_tsv(fork_bed_path, append=T, col_names=F)
  p$significant %>%
    dplyr::mutate(score=-log10(tz_pred), name=as.character(tz_start), strand="*") %>%
    dplyr::select(tz_chrom, tz_start, tz_end, name, score, strand) %>%
    readr::write_tsv(paste0("reports/09-replication_fork_nn/predicted_tz_", celltype, "-", path_suffix, ".bed"), col_names=F)
  p$regions %>%
    dplyr::select(prediction_chrom, prediction_start, prediction_end) %>%
    readr::write_tsv(paste0("reports/09-replication_fork_nn/predicted_regions_tz_", celltype, "-", path_suffix, ".bed"), col_names=F)
  p$all %>%
    dplyr::group_by(tz_zoomout) %>%
    dplyr::do((function(z){
      z %>%
        dplyr::select(tz_chrom, tz_start, tz_end, tz_pred) %>%
        readr::write_tsv(paste0("reports/09-replication_fork_nn/predicted_tz_", celltype, "-", path_suffix, "-zoomout", z$tz_zoomout[1], ".bedgraph"), col_names=F)
    })(.))




  #
  # Export annotations
  #
  tz_df = readr::read_tsv("data/tz/tz_annotation.tsv")
  fork_annotated_path = "reports/09-replication_fork_nn/forks_annotated.bed"
  writeLines('track name="Replication forks" itemRgb=On visibility=2 colorByStrand="255,0,0 0,0,255"', con=fork_annotated_path)
  tz_df %>%
    dplyr::select(tz_train_chrom=tz_chrom, tz_train_end=tz_start, tz_train_start=tz_end, tz_iz_left, tz_iz_right) %>%
    df2ranges(tz_train_chrom, tz_train_start, tz_train_start) %>%
    innerJoinByOverlaps(predictions_significant_df %>% df2ranges(tz_chrom, tz_start, tz_end)) %>%
    reshape2::melt(measure.vars=c("tz_iz_left", "tz_iz_right"), variable.name="fork_iz_location", value.name="fork_iz") %>%
    dplyr::mutate(
      fork_tz=tz_end/2+tz_start/2,
      fork_chrom=tz_chrom, fork_name=dplyr::case_when(fork_iz_location=="tz_iz_left"~"right fork", fork_iz_location=="tz_iz_right"~"left fork"),
      fork_start=pmin(fork_iz, fork_tz), fork_end=pmax(fork_tz, fork_iz), fork_score=-log10(tz_pred),
      fork_strand=dplyr::case_when(fork_iz_location=="tz_iz_left"~"+", fork_iz_location=="tz_iz_right"~"-"),
      thickStat=fork_start, thickEnd=fork_end, fork_color=dplyr::case_when(fork_strand=="+"~"#FF0000", fork_strand=="-"~"#0000FF")
    ) %>%
    dplyr::select(fork_chrom, fork_start, fork_end, fork_name, fork_score, fork_strand, thickStat, thickEnd, fork_color) %>%
    readr::write_tsv(fork_annotated_path, append=T, col_names=F)

}