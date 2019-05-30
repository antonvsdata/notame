
library(ggplot2)

devtools::document()

hilic_neg <- read_from_excel(file = "~/amp/inst/extdata/HILIC neg final.xlsx",
                       sheet = 1,
                       corner_row = 5, corner_column = "U", name = "HILIC_neg")

hilic_pos <- read_from_excel(file = "~/amp/inst/extdata/HILIC pos final.xlsx",
                             sheet = 1,
                             corner_row = 4, corner_column = "U",
                             name = "HILIC_pos")

rp_neg <- read_from_excel(file = "~/amp/inst/extdata/RP neg final.xlsx",
                             sheet = 1,
                             corner_row = 4, corner_column = "U",
                             name = "RP_neg")

rp_pos <- read_from_excel(file = "~/amp/inst/extdata/RP pos final.xlsx",
                             sheet = 1,
                             corner_row = 4, corner_column = "U",
                             name = "RP_pos")

all_modes <- list(hilic_neg, hilic_pos, rp_neg, rp_pos)

lcms <- list()

for (mode in all_modes) {
  lcms <- c(lcms, construct_MetaboSet(assay_data = mode$assay_data,
                              pheno_data = mode$pheno_data,
                              feature_data = mode$feature_data,
                              group_col = "Class"))
}

merged <- merge_metabosets(lcms)
identical(rbind(fData(lcms[[1]]), fData(lcms[[2]]), fData(lcms[[3]]), fData(lcms[[4]])), fData(merged))
identical(rbind(exprs(lcms[[1]]), exprs(lcms[[2]]), exprs(lcms[[3]]), exprs(lcms[[4]])), exprs(merged))


hilics <- merge_metabosets(lcms[[1]], lcms[[2]])


lol <- lcms[[1]]

dim(lol)

group_col(lol) <- "Class"


normalized <- inverse_normalize(lol)

xd <- correct_drift(lol)

lol <- assess_quality(lol)
xd <- assess_quality(xd)

quality(xd)[2:5] - quality(lol)[2:5]


ins <- inspect_dc(lol, xd, condition = "RSD_r < 0 & D_ratio_r < 0")

exprs(ins)

xd <- assess_quality(xd)
quality(xd)

quality(ins)


plot_dc(lol[1:10], xd[1:10], file = "lolxd.pdf")

pl <- plot_pca(lol, color = "Injection_order")
class(pl)


lol <- drop_qcs(lol)


save_subject_line_plots(drop_qcs(example_set), "xd.pdf")

save_group_lineplots(drop_qcs(example_set), file = "oho.pdf")

save_group_boxplots(drop_qcs(example_set), file = "nais.pdf")
