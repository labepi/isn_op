bandt_pompe_path <- "~/ufrn/bandt_pompe/"

source(paste(bandt_pompe_path, "bandt_pompe.R", sep = ""))
source(paste(bandt_pompe_path, "measures.R", sep = ""))
source(paste(bandt_pompe_path, "helpers.R", sep = ""))

# for plotting
library(ggplot2)
library(ggrepel)

# for pvclust
library(pvclust)
library(dendextend)

#randomForest
library(randomForest)

# parameters
tau <- 1
D_l <- c(2, 3, 4, 5, 6, 7)

# directory of samples
dir_samples <- "../tcp-isn/data/samples/"

# getting the list of files
files <- list.files(dir_samples)

# the names of files, removing the final part
names <- gsub("-60k.csv", "", files)

# features totais
feats <- matrix(nrow = 0, ncol = (12 + 12 + 12 + 1))
colnames(feats) <- c(
    "H2", "C2", "F2",
    "H2diff", "C2diff", "F2diff",
    "H3", "C3", "F3",
    "H3diff", "C3diff", "F3diff",
    "H4", "C4", "F4",
    "H4diff", "C4diff", "F4diff",
    "H5", "C5", "F5",
    "H5diff", "C5diff", "F5diff",
    "H6", "C6", "F6",
    "H6diff", "C6diff", "F6diff",
    "H7", "C7", "F7",
    "H7diff", "C7diff", "F7diff",
    "label"
)

# looping nas amostras
for (i in 1:length(files)) {
    filename <- paste(dir_samples, files[i], sep = "")

    x <- read.csv(filename)

    # features for each serie
    feats_i <- c()

    # computing features for all D's
    for (D in D_l)
    {
        # bandt-pomping the series
        bp <- bandt_pompe_distribution(x[, 1], D = D)
        HP <- shannon_entropy(bp$probabilities, normalized = TRUE)
        SC <- complexity(bp$probabilities, HP)
        FI <- fisher_information(bp$probabilities)
        feats_i <- c(feats_i, HP, SC, FI)

        # bandt-pomping the diff of series
        bp <- bandt_pompe_distribution(diff(x[, 1]), D = D)
        HP <- shannon_entropy(bp$probabilities, normalized = TRUE)
        SC <- complexity(bp$probabilities, HP)
        FI <- fisher_information(bp$probabilities)
        feats_i <- c(feats_i, HP, SC, FI)
    }

    # adding the label
    feats_i <- c(feats_i, names[i])

    # adding the features to the matrix
    feats <- rbind(feats, feats_i)
}


# ajuste das features
feats_df <- as.data.frame(cbind(apply(feats[, 1:36], 2, as.numeric)))
feats_df$label <- as.factor(feats[, 37])
rownames(feats_df) = names

# ja sabemos que estas sao as melhores features
feats_list = c("F2", "F3", "F4", "F5", "F6", "F7", "F2diff", "F3diff", "F4diff", "F5diff", "F6diff", "F7diff")

# 1.
# pvclust

# testando todas as features
# feats_only = feats_df[, -c(37)]
#
# set.seed(123)
#
# # res_pv <- pvclust(t(feats_only), method.dist = "cor", method.hclust = "average", nboot = 10)
# res_pv <- pvclust(t(feats_only), method.dist = "euclid", method.hclust = "average", nboot = 1000)
# plot(res_pv, labels=names)
# pvrect(res_pv)


# testando best features
feats_best = feats_df[, feats_list]

set.seed(123)

res_pv_best <- pvclust(t(feats_best), method.dist = "euclid", method.hclust = "average", nboot = 1000)

# organizando com o dendextend

png("images/dendro/dendrogram_all_pvclust.png", width=1500, height=600)

par(mar = c(15,2,2,2))

res_pv_best %>% as.dendrogram %>% 
    set("branches_k_color", k = 6) %>% 
    set("branches_lwd", 4) %>%
    set("labels_cex", 2) %>%
    set("labels_col", k=6) %>%
    # set("nodes_pch", 21) %>% # set all nodes to be pch 21
    set("leaves_pch", 19) %>%
    set("leaves_cex", 1.2) %>%
    set("leaves_col", 1) %>%
    # set("hang", 0.0) %>%
    pvclust_show_signif(res_pv_best, alpha=0.05) %>%
    # pvclust_show_signif_gradient(result, alpha=0.01) %>%
    plot()

res_pv_best %>% as.dendrogram %>% 
    rect.dendrogram(k=6, lty = 5, lwd = 0, col=rgb(0.1, 0.2, 0.4, 0.1), upper_rect=0, lower_rect=-2.5)

dev.off()


# subclusters
dend_list <- get_subdendrograms(as.dendrogram(res_pv_best), 6)
dend_linux = dend_list[[4]]


png("images/dendro/dendrogram_linux_pvclust.png", width=1500, height=600)

par(mar = c(20,2,2,2))

dend_linux %>%  
    set("branches_k_color", k = 2) %>% 
    set("branches_lwd", 4) %>%
    set("labels_cex", 2) %>%
    set("labels_col", k=2) %>%
    # set("nodes_pch", 21) %>% # set all nodes to be pch 21
    set("leaves_pch", 19) %>%
    set("leaves_cex", 1.2) %>%
    set("leaves_col", 1) %>%
    # set("hang", 0.0) %>%
    # pvclust_show_signif(res_pv_best, alpha=0.05) %>%
    # pvclust_show_signif_gradient(result, alpha=0.01) %>%
    plot()

dend_linux %>% rect.dendrogram(k=2, lty = 5, lwd = 0, col=rgb(0.1, 0.2, 0.4, 0.1), upper_rect=0, lower_rect=-0.25)

dev.off()

#
#
#
# png("plot_strategy_4.png", width=1000, height=800)
#
# # with pvclust
# plot(res_pv_best, labels=names, print.pv=c("au"), print.num=F, hang=.05, cex=1.05, cex.pv=1.05)
# pvrect(res_pv_best, alpha=0.95, pv="au", type="geq")
#
#
# # with dendextend
# dend <- as.dendrogram(res_pv_best)
# res_pv_best %>% as.dendrogram %>% 
#    plot(main = "Cluster dendrogram with AU/BP values (%)\n reproduced plot with dendrogram")
# res_pv_best %>% text
# res_pv_best %>% pvrect
#
# dend %>% pvclust_show_signif(res_pv_best) %>% 
#    plot(main = "Cluster dendrogram \n bp values are highlighted by signif")
#
# result = res_pv_best
#
# dend %>% pvclust_show_signif(result, show_type = "lwd") %>% 
#    plot(main = "Cluster dendrogram with AU/BP values (%)\n bp values are highlighted by signif")
# result %>% text
# result %>% pvrect(alpha=0.95)
#
#
# dend %>% pvclust_show_signif_gradient(result) %>% 
#    plot(main = "Cluster dendrogram with AU/BP values (%)\n bp values are colored by signif")
#
# # organizando com o dendextend
#
# png("dendrogram_all_pvclust.png", width=1500, height=400)
#
# par(mar = c(15,2,2,2))
#
# res_pv_best %>% as.dendrogram %>% 
#     set("branches_k_color", k = 6) %>% 
#     set("branches_lwd", 4) %>%
#     set("labels_cex", 1.6) %>%
#     set("labels_col", k=6) %>%
#     # set("nodes_pch", 21) %>% # set all nodes to be pch 21
#     set("leaves_pch", 19) %>%
#     set("leaves_cex", 1.2) %>%
#     set("leaves_col", 1) %>%
#     # set("hang", 0.0) %>%
#     pvclust_show_signif(result, alpha=0.05) %>%
#     # pvclust_show_signif_gradient(result, alpha=0.01) %>%
#     plot()
#
# res_pv_best %>% as.dendrogram %>% 
#     rect.dendrogram(k=6, lty = 5, lwd = 0, col=rgb(0.1, 0.2, 0.4, 0.1), upper_rect=0, lower_rect=-2.5)
#
# #print.pv=c("au")) #horiz = TRUE, hang=1.5)
# # res_pv_best %>% text(print.num=F, alpha=0.99)
#
# # res_pv_best %>% pvrect(alpha=0.99, border=1)
#
# # res_pv_best %>% pvrect2(lower_rect=-0.7, border=2, alpha = 0.95, lty = 2, lwd=3, max.only=T)
#
# # res_pv_best %>% pvrect2(lower_rect=-0.7, border=1, alpha = 0.95, lty = 2, lwd=3, max.only=T)
#
# # rect.dendrogram(dend, 6, border = 1:6)
#
#
# dev.off()




# obtendo os clusters
clusters = pvpick(res_pv_best)$clusters

# organizando as suas classes
class_l = matrix(nrow = 0, ncol = 2)
for (i in 1:length(clusters))
{
    class_l = rbind(class_l, do.call(rbind, lapply(clusters[[i]], c, i)))
}

print("names e suas classes")
print(class_l)

# to save the names and classes in a csv file
write.table(class_l, "classes_list.csv", sep = ",",
    row.names = FALSE, col.names = FALSE)


# TODO:
# - testar clusters com 10k, .., 60k (da os mesmos)
# - ver os scatter plots das diferentes classes (class_l)
# - fazer novo split de datasets de series de 10k
# - adicionar ao final da serie a sua classe (ver class_l)
# - fazer uma classificacao com varios algoritmos de classif.
# - ver como imprimir pvclust legal (ggplot2?)



# 2.
# hclust

res_dd = dist(feats_best, method="euclid")
res_hc = hclust(res_dd, method = "average")


# plot(as.phylo(res_hc), type = "fan")
# hc_colors = c("red", "blue", "green", "black", "orange", "")

# sabemos que tem 6 clusters 
clus6 = cutree(res_hc, 6)

png("plot_strategy_5.png", width=1000, height = 1000)

# plot(as.phylo(res_hc), type = "fan", tip.color = hc_colors[clus6])
plot(as.phylo(res_hc), type = "fan", tip.color = clus6, label.offset=0.1, edge.color=1, edge.lty=1, edge.width=1, font=2)

dev.off()












#############################################################################
#############################################################################
#############################################################################
####################      FIM AQUI      #####################################
#############################################################################
#############################################################################
#############################################################################



# analise de plots (antigo)

# par(mfrow = c(2, 4))
#
# plot.ccep(D = 3, H = feats_df$H3, SC = feats_df$C3)
# plot.ccep(D = 3, H = feats_df$H3diff, SC = feats_df$C3diff)
# plot.ccep(D = 4, H = feats_df$H4, SC = feats_df$C4)
# plot.ccep(D = 4, H = feats_df$H4diff, SC = feats_df$C4diff)
# plot.ccep(D = 5, H = feats_df$H5, SC = feats_df$C5)
# plot.ccep(D = 5, H = feats_df$H5diff, SC = feats_df$C5diff)
# plot.ccep(D = 6, H = feats_df$H6, SC = feats_df$C6)
# plot.ccep(D = 6, H = feats_df$H6diff, SC = feats_df$C6diff)

gplot_plane <- function(feats_df, D = 6, x = "H6diff", y = "C6diff", xlim = c(0, 1.05), ylim = c(0, 0.5), title = "") {
    # gerando o plot
    p <- gplot.ccep(D = D, xlim = xlim, ylim = ylim, main = title)

    # changing the labs
    p <- p + xlab(expression(Normalized ~ Permutation ~ Entropy ~ (H[S]))) +
        ylab(expression(Statistical ~ Complexity ~ (C[JS])))

    # adjusting layout and centering title
    p <- p + theme_bw(base_size = 24) + theme(plot.title = element_text(hjust = 0.5))

    # adding points
    p <- p + geom_point(aes(x = {{x}}, y = {{y}}, colour = label), size = 4, shape = 19, data = feats_df) + scale_colour_discrete(name = "OS")

    # adding text
    p <- p + geom_label_repel(aes(x = {{x}}, y = {{y}}, label = label, color = label),
        segment.alpha = 0.6,
        size = 7,
        label.size = NA,
        point.padding = 0.2,
        box.padding = 0.2,
        show.legend = FALSE,
        max.overlaps = 20,
        data = feats_df
    )

    return(p)
}

#+ scale_shape_discrete(name="OS")

# , values=c('blue', 'red', 'green', 'darkgreen')) +
#        scale_shape_manual(name='Types', values=c(19, 15, 17, 1))




# clustering with kmeans
########################

# resk <- kmeans(feats_df[, 1:16], centers = 5)
# cbind(names, resk$cluster)

resk <- kmeans(feats_df[, c("H3diff", "C3diff", "H4diff", "C4diff", "H5diff", "C5diff", "H6diff", "C6diff")], centers = 5)
cbind(names, resk$cluster)

resk <- kmeans(feats_df[, c("H6diff", "C6diff")], centers = 6, iter.max = 50, nstart = 25)
cbind(names, resk$cluster)

resk <- kmeans(feats_df[, c("H2diff", "C2diff")], centers = 6, iter.max = 50, nstart = 25)
cbind(names, resk$cluster)

resk <- kmeans(feats_df[, c("H7diff", "C7diff")], centers = 6, iter.max = 50, nstart = 25)
cbind(names, resk$cluster)


resk <- kmeans(feats_df[, feats_list], centers = 6, iter.max = 50, nstart = 25)
cbind(names, resk$cluster)

###############


# so teste
# cbind(feats_df_small[,c(15,16)], feats_df_large[,c(15,16)], names)


# Hclust
########


plot(hclust(dist(feats_df[, c("H6diff", "C6diff")]), method = "single"), labels = names)

plot(hclust(dist(feats_df[, c("H6diff", "C6diff")]), method = "complete"), labels = names)

plot(hclust(dist(feats_df[, c("H2diff", "C2diff")]), method = "complete"), labels = names)

plot(hclust(dist(feats_df[, c("H7diff", "C7diff")]), method = "complete"), labels = names)

plot(hclust(dist(feats_df[, c("H3diff", "C3diff", "H4diff", "C4diff", "H5diff", "C5diff", "H6diff", "C6diff")]), method = "complete"), labels = names)


###############

library(plotly)

data <- read.csv(paste(dir_samples,files[3], sep=""))

data_diff <- diff(data[, 1])
x <- 2:length(data_diff)

plot_ly(x=data_diff[x], y=data_diff[x-1], type='scatter', mode='markers')


# create dataset of entropy/complexity
##############
dataset <- matrix(nrow = 0, ncol = 25)
colnames(dataset) <- c(
    "H2", "C2",
    "H2diff", "C2diff",
    "H3", "C3",
    "H3diff", "C3diff",
    "H4", "C4",
    "H4diff", "C4diff",
    "H5", "C5",
    "H5diff", "C5diff",
    "H6", "C6",
    "H6diff", "C6diff",
    "H7", "C7",
    "H7diff", "C7diff",
    "label"
)

file <- read.csv("../tcp-isn/data/dataset.csv")
for (i in 1:nrow(file)) {
    print(i)
    row <- t(as.numeric(file[i,]))

    rows <- head(row[1, ], -1)

    hc2 <- complexity_entropy(rows, 2, tau)
    hc2diff <- complexity_entropy(diff(rows), 2, tau)
    hc3 <- complexity_entropy(rows, 3, tau)
    hc3diff <- complexity_entropy(diff(rows), 3, tau)
    hc4 <- complexity_entropy(rows, 4, tau)
    hc4diff <- complexity_entropy(diff(rows), 4, tau)
    hc5 <- complexity_entropy(rows, 5, tau)
    hc5diff <- complexity_entropy(diff(rows), 5, tau)
    hc6 <- complexity_entropy(rows, 6, tau)
    hc6diff <- complexity_entropy(diff(rows), 6, tau)
    hc7 <- complexity_entropy(rows, 7, tau)
    hc7diff <- complexity_entropy(diff(rows), 7, tau)

    dataset_i <- c(
        hc2[1], hc2[2],
        hc2diff[1], hc2diff[2],
        hc3[1], hc3[2],
        hc3diff[1], hc3diff[2],
        hc4[1], hc4[2],
        hc4diff[1], hc4diff[2],
        hc5[1], hc5[2],
        hc5diff[1], hc5diff[2],
        hc6[1], hc6[2],
        hc6diff[1], hc6diff[2],
        hc7[1], hc7[2],
        hc7diff[1], hc7diff[2],
        file$label[i]
    )
    dataset <- rbind(dataset, dataset_i)
}

dataset_df <- as.data.frame(cbind(apply(dataset[, 1:24], 2, as.numeric)))
dataset_df$label <- as.factor(dataset[, 25])
