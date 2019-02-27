
library(ggplot2)

devtools::document()
doParallel::registerDoParallel(cores = 3)

res <- read_from_excel(file = "~/amp/inst/extdata/HILIC neg final.xlsx",
                       sheet = 1,
                       corner_row = 5, corner_column = "V",
                       split_by = c("Mode"))


lcms <- construct_MetaboSet(assay_data = res$assay_data,
                            pheno_data = res$pheno_data,
                            feature_data = res$feature_data)

lol <- lcms[[1]]

dim(lol)

group_col(lol) <- "Class"




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

plot_pca(lol, color = "Injection_order")
