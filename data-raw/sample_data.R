devtools::load_all()


modes <- c("HILIC_neg", "HILIC_pos", "RP_neg", "RP_pos")

objects <- list()
for (mode in modes) {
  dada <- read_from_excel(file = system.file("extdata", paste0(mode, "_sample.xlsx"), package = "amp"),
                          name = mode, id_prefix = "Sample_")
  pd <- dada$pheno_data
  pd <- tidyr::separate(pd, col = "Class", into = c("Group", "Time"), sep = 1)
  pd$Group[pd$Group == "Q"] <- "QC"
  pd$Time[pd$Time == "C"] <- "QC"

  object <- construct_MetaboSet(exprs= dada$exprs, pheno_data = pd, feature_data = dada$feature_data,
                                group_col = "Group", time_col = "Time")
  objects[[mode]] <- object[[1]]
}

merged_sample <- merge_metabosets(objects)

hilic_neg_sample <- objects[[1]]
usethis::use_data(hilic_neg_sample, overwrite = TRUE)
hilic_pos_sample <- objects[[2]]
usethis::use_data(hilic_pos_sample, overwrite = TRUE)
rp_neg_sample <- objects[[3]]
usethis::use_data(rp_neg_sample, overwrite = TRUE)
rp_pos_sample <- objects[[4]]
usethis::use_data(rp_pos_sample, overwrite = TRUE)

usethis::use_data(merged_sample)
