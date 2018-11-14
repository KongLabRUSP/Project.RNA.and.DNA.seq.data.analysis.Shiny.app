getwd()
source("source/app_v6/data_prep.R")
ct <- fread("data/featurescounts_uvb-skin_dedup_renyi_2-9-2018.csv", skip = 1)
info <- fread("data/design_table.csv", header = TRUE)
ct <- as.data.frame(ct)
info <- as.data.frame(info)
sample_label_selected_expl <- c("02w_CON_1", "02w_CON_0", "02w_UAA_1",  "02w_UAA_0",  "02w_SFN_1", "02w_SFN_0")
covariate_selected_expl <- "trt"
dt <- match.ct.info(ct, info, sample_label_selected_expl , "Geneid")

# selected design table
selected_info <- info[match(sample_label_selected_expl, info$`Sample Label`), c("Sample Label", covariate_selected_expl)]
colnames(selected_info) <- c("sample", "covariate")
# log transform
dt_log2 <- dt
dt_log2[, 2:ncol(dt_log2)] <- log2(dt_log2[, 2:ncol(dt_log2)] +1)
dt_log2

# need to do a cluster before heatmap
tmp <- dt_log2[, -1]
sample_dist <- as.matrix(dist(t(as.matrix(tmp))))
p3 <- plot_ly(
  x = colnames(sample_dist), y = rownames(sample_dist),
  z = sample_dist, type = "heatmap")
p3

# pca
matrix_pca <- t(dt_log2[, -1])
matrix_pca <- matrix_pca[, apply(matrix_pca, 2, var) != 0]
res_pca <- prcomp(matrix_pca,
                  center = TRUE,
                  scale. = TRUE)
df_out <- as.data.frame(res_pca$x)

percentage <- round(res_pca$sdev^2 / sum(res_pca$sdev^2) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste0( as.character(percentage), "%"), ")")

df_out$sample <- colnames(dt_log2)[-1]
df_out <- merge(df_out, selected_info, by = "sample")

pca <- ggplot(data = df_out,
                 aes(x = PC1,
                     y = PC2)) +
  geom_point(aes(fill = df_out[, ncol(df_out)]),
             shape = 21,
             size = 3,
             alpha = 0.5) +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  geom_text(data = df_out,
            aes(x = PC1 + 20,
                y = PC2,
                label = df_out$sample),
            size = 2,
            hjust = 0.5) +
  labs(fill = "legend")

p4 <- ggplotly(pca)
p4


# test cige
# do DEGexp


