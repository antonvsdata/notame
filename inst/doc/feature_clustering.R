## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, include=FALSE------------------------------------------
devtools::load_all()

## ---- out.width = "400px", echo=FALSE, fig.align='center'----------------
knitr::include_graphics("algo_1.jpg")

## ---- fig.show = "hold", out.width = "40%", fig.align = "default", echo = FALSE----

knitr::include_graphics("clust_1_1.jpg")

knitr::include_graphics("clust_1_2.jpg")


## ---- fig.show = "hold", out.width = "40%", fig.align = "default", echo=FALSE----

knitr::include_graphics("clust_1_3.jpg")

knitr::include_graphics("clust_1_5.jpg")


## ---- fig.show = "hold", out.width = "40%", fig.align = "default", echo=FALSE----

knitr::include_graphics("clust_1_6.jpg")

knitr::include_graphics("clust_3_1.jpg")


## ---- out.width = "400px", echo=FALSE, fig.align='center'----------------
knitr::include_graphics("algo_2.jpg")

## ------------------------------------------------------------------------
clustered <- cluster_features(example_set, rt_window = 2, corr_thresh = 0.4, d_thresh = 0.6)


## ------------------------------------------------------------------------
colnames(results(clustered))
head(results(clustered)$Cluster_ID)

## ------------------------------------------------------------------------
compressed <- compress_clusters(clustered)

## ------------------------------------------------------------------------
data <- combined_data(example_set)
features <- fData(example_set)

## ------------------------------------------------------------------------
conn <- find_connections(data = data, features = features,
                 corr_thresh = 0.4, rt_window = 2,
                 name_col = "Feature_ID", mz_col = "Mass", rt_col = "RetentionTime")
head(conn)

## ------------------------------------------------------------------------
clusters <- find_clusters(connections = conn, d_thresh = 0.6)

## ------------------------------------------------------------------------
features_clustered <- assign_cluster_id(data, clusters, features, name_col = "Feature_ID")

## ------------------------------------------------------------------------
pulled <- pull_clusters(data, features_clustered, name_col = "Feature_ID")
cluster_data <- pulled$cdata
cluster_features <- pulled$cfeatures

