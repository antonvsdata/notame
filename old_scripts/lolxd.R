library(ggplot2)


devtools::document()
doParallel::registerDoParallel(cores = 3)
res <- read_from_excel(file = "~/amp/inst/extdata/split_data.xlsx",
                       sheet = 2,
                       corner_row = 4, corner_column = "G",
                       split_by = c("Column"))


lcms <- construct_MetaboSet(assay_data = res$assay_data,
                            pheno_data = res$pheno_data,
                            feature_data = res$feature_data)

lol <- lcms[[1]]

group_col(lol) <- "Group"


xd <- correct_drift(lol)

lol <- assess_quality(lol)
xd <- assess_quality(xd)

quality(xd)[2:5] - quality(lol)[2:5]


ins <- inspect_dc(lol, xd, condition = "RSD_r < 0 & D_ratio_r < 0")

exprs(ins)

xd <- assess_quality(xd)
quality(xd)

quality(ins)


plot_dc(lol, xd, file = "lolxd.pdf")
