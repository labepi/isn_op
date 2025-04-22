bandt_pompe_path <- "~/ufrn/bandt_pompe/"

source(paste(bandt_pompe_path, "bandt_pompe.R", sep = ""))
source(paste(bandt_pompe_path, "measures.R", sep = ""))
source(paste(bandt_pompe_path, "helpers.R", sep = ""))

# for plotting
library(ggplot2)
library(ggrepel)

# for nice dendrogram
library(ape)

# for pvclust
library(pvclust)
library(dendextend)

# for feature selection
library(caret)

# randomForest
library(randomForest)


# parameters
tau <- 1
D_l <- c(2, 3, 4, 5, 6, 7)

# directory of samples
dir_samples <- "../tcp-isn/data/switches/"

# getting the list of files
files <- list.files(dir_samples)

# the names of files, removing the final part
names <- gsub("-60k.csv", "", files)

# getting names short
names_short <- unlist(lapply(strsplit(names, "-"), "[", 1))

# features totais
feats <- matrix(nrow = 0, ncol = (12 + 12 + 12 + 2))
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
    "label", "class"
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

    # adding the label and class
    feats_i <- c(feats_i, names[i], names_short[i])

    # adding the features to the matrix
    feats <- rbind(feats, feats_i)
}

# ajuste das features
feats_df <- as.data.frame(cbind(apply(feats[, 1:36], 2, as.numeric)))
feats_df$label <- as.factor(feats[, 37])
feats_df$class <- as.factor(feats[, 38])
rownames(feats_df) = names

# ESTRATEGIA:
# 1. Estima-se a quantidade de clusters com o pvclust
# 2. Utiliza-se esta quantidade para verificar as features mais importantes
# 3. Faz uma classificaçao com essa quantidade de classes e features


# 1. estimando a quantidade de features

feats_only = feats_df[, -c(37,38)]

set.seed(123)

res_pv <- pvclust(t(feats_only), method.dist = "euclid", method.hclust = "average", nboot = 1000)
# res_pv <- pvclust(t(feats_only), method.dist = "cor", method.hclust = "average", nboot = 10)

png("images/dendro/dendrogram_switches_pvclust.png", width=1500, height=500)

par(mar = c(10,2,2,2))

res_pv_dend = as.dendrogram(res_pv)

# atualizando/simplificando a lista de labels
lab_names = labels(res_pv_dend)
lab_short <- unlist(lapply(strsplit(lab_names, "-"), "[", 1))
labels(res_pv_dend) = lab_short

res_pv_dend  %>% 
    set("branches_k_color", k = 3) %>% 
    set("branches_lwd", 4) %>%
    set("labels_cex", 2) %>%
    set("labels_col", k=3) %>%
    # set("nodes_pch", 21) %>% # set all nodes to be pch 21
    set("leaves_pch", 19) %>%
    set("leaves_cex", 1.2) %>%
    set("leaves_col", 1) %>%
    # set("hang", 0.0) %>%
    pvclust_show_signif(res_pv, alpha=0.05) %>%
    # pvclust_show_signif_gradient(result, alpha=0.01) %>%
    plot()

res_pv_dend %>% 
    rect.dendrogram(k=3, lty = 5, lwd = 0, col=rgb(0.1, 0.2, 0.4, 0.1), upper_rect=0, lower_rect=-2)

#res_pv %>% text()

dev.off()



# png("plot_strategy_1.png", width=1000, height=800)
#
# plot(res_pv, labels=names_short)
# pvrect(res_pv)
#
# dev.off()


# sabemos que sao 3 classes, entao, fazendo este corte
class_pv = cutree(res_pv, k=3)


clusters = pvpick(res_pv)$clusters

# organizando as suas classes
class_l = matrix(nrow = 0, ncol = 2)
for (i in 1:length(clusters))
{
    class_l = rbind(class_l, do.call(rbind, lapply(clusters[[i]], c, i)))
}

print("names e suas classes")
print(class_l)



# 2. 

# criando um novo dataset, apenas com as features e as 3 classes descobertas
feats_imp = feats_only
feats_imp$class_pv = class_pv

# nao deu muito certo aqui
# control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# results <- rfe(class_pv ~ ., data=feats_imp, sizes=c(1:36), rfeControl=control)

# using randomForest

set.seed(123)

rf1 = randomForest(class_pv ~ ., data=feats_imp, importance=TRUE, ntree=1000)
imp = importance(rf1)

# ordering by IncNodePurity
imp_order = imp[order(imp[,2], decreasing=T),c(2,1)]
imp_order

# TODO: formatar uma tabela aqui
#F3        3.84262870 11.192388
#F4        3.78024889 11.301750
#F2        3.73211059 11.154605
#F5        3.47780622 10.451854
#F7        3.47233059 10.466777
#F6        3.42555357 10.764874
#F2diff    1.36398844  7.981063
#F7diff    0.94910366  7.075326
#F5diff    0.79794566  6.725263
#F3diff    0.66551643  6.199401
#F6diff    0.65947811  5.961829
#F4diff    0.65140719  7.069170
#C3diff    0.42942925  4.096629
#C5diff    0.41344464  4.273094
#C4diff    0.37292135  3.805934
#C7diff    0.24487603  4.755358
#C7        0.16106920  3.197786
#C2diff    0.14941777  6.926048
#
#C6diff    0.14510587  3.820731
#H7        0.14107573  4.615293
#C6        0.13852316  3.955866
#H6diff    0.12907386  5.279817
#H3        0.12578292  4.745957
#H3diff    0.12219956  5.046966
#H4diff    0.10811311  4.923966
#H7diff    0.09907142  4.505349
#H5diff    0.09831714  5.419994
#H4        0.08859440  4.165168
#H2diff    0.08512328  6.094627
#H5        0.05966909  3.155293
#C3        0.05184454  3.051156
#C2        0.03366424  2.986155
#C4        0.02935758  2.460076
#H2        0.02847134  3.997762
#C5        0.02785580  3.219042
#H6        0.02516636  3.562481

# gerando a imagem da impureza nos nós

bpd.df = data.frame(x=row.names(imp_order), y=imp_order[,1])
row.names(bpd.df) = 1:36

p = ggplot(aes(x=reorder(x, y, decreasing=T), y=y), data=bpd.df) +
        geom_bar(stat='identity', width=0.7, 
                 fill='darkgray', 
                colour='black') + 
        ylab('Inc. na Impureza') + xlab('Atributos') + 
        geom_hline(yintercept = 0.5, colour="red") +
        theme_bw(base_size=28) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              plot.title=element_text(hjust=0.5))
        
# ggtitle(title)

ggsave("images/strategy/feature_sel_nodeimpurity.png", p)

#, width=800, height=200, units="px")




png("plot_strategy_2.png", width=1000, height=800)

vimp = varImp(rf1)
varImpPlot(rf1)

dev.off()


# melhor feature?
boxplot(F4 ~ class_pv, data = feats_imp)
# pior feature?
boxplot(C2 ~ class_pv, data = feats_imp)






# 3. testando uma classificacao

# feats_list = c("F2", "F3", "F4", "F5", "F6", "F7")
feats_list = c("F2", "F3", "F4", "F5", "F6", "F7", "F2diff", "F3diff", "F4diff", "F5diff", "F6diff", "F7diff")

# feats_best = feats_imp[, c(feats_list, "class_pv")]
feats_best = feats_imp[, c(feats_list, "class_pv")]

# 3

# treinando o model
set.seed(123)
rfmodel = randomForest(class_pv ~ ., data=feats_best, ntree=1000)

# fazendo a classificacao (apenas pra testar)
cbind(round(predict(rfmodel, feats_df[, feats_list])), class_pv)



#############################################################################
#############################################################################
#############################################################################
####################      FIM AQUI      #####################################
#############################################################################
#############################################################################
#############################################################################



# feature importance
####################

# TODO: fazer feature importance aqui

feats_only = feats_df[, -c(37,38)]

# testing correlation
corm = cor(feats_only)
corlist = findCorrelation(corm, cutoff=0.5)


# testing fests importance
feats_class = feats_df[,-c(37)]

# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3)

# train the model
TG = expand.grid(k=1:3,size=seq(5,20,by=5))
# model <- train(x=feats_only, y=factor(names_short), method="lvq", preProcess="scale", trControl=control, tuneGrid=TG)

model <- train(class ~ ., data=feats_class, method="lvq", preProcess="scale", trControl=control, tuneGrid=TG)

# estimate variable importance
importance <- varImp(model, scale=FALSE)


# testing RFE

control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(class ~ ., data=feats_class, sizes=c(1:36), rfeControl=control)
# summarize the results


# clustering with kmeans
########################

# resk <- kmeans(feats_df[, 1:16], centers = 5)
# cbind(names, resk$cluster)

resk <- kmeans(feats_df[, c("H3diff", "C3diff", "H4diff", "C4diff", "H5diff", "C5diff", "H6diff", "C6diff")], centers = 4)
cbind(names, resk$cluster)

resk <- kmeans(feats_df[, c("H6diff", "C6diff")], centers = 4, iter.max = 50, nstart = 25)
cbind(names, resk$cluster)

resk <- kmeans(feats_df[, c("H2diff", "C2diff")], centers = 4, iter.max = 50, nstart = 25)
cbind(names, resk$cluster)

resk <- kmeans(feats_df[, c("H7diff", "C7diff")], centers = 4, iter.max = 50, nstart = 25)
cbind(names, resk$cluster)


resk <- kmeans(feats_df[, -c(37,38)], centers = 3, iter.max = 50, nstart = 25)
cbind(names, resk$cluster)



resk <- kmeans(feats_df[, feats_list2], centers = 3, iter.max = 50, nstart = 25)
cbind(names, resk$cluster)


###############


# Hclust
########


plot(hclust(dist(feats_df[, c("H6diff", "C6diff", "FIdiff")]), method = "single"), labels = names_short)

plot(hclust(dist(feats_df[, c("H6diff", "C6diff")]), method = "complete"), labels = names_short)

plot(hclust(dist(feats_df[, c("H2diff", "C2diff")]), method = "complete"), labels = names_short)

plot(hclust(dist(feats_df[, c("H7diff", "C7diff")]), method = "complete"), labels = names_short)

plot(hclust(dist(feats_df[, c("H3diff", "C3diff", "H4diff", "C4diff", "H5diff", "C5diff", "H6diff", "C6diff")]), method = "complete"), labels = names)


# testando outros modelos de dendrogram

dd <- dist(scale(feats_df[, c("H6diff", "C6diff")]), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

colors <- c("red", "blue", "green", "black")
clus4 <- cutree(hc, 4)

plot(as.phylo(hc),
    type = "fan", tip.color = colors[clus4],
    label.offset = 1, cex = 0.7
)


# pvclust
#########

set.seed(1234)

# result <- pvclust(t(feats_df[,-c(37, 38)]), method.dist = "euclid", method.hclust = "average", nboot = 10)
# result <- pvclust(t(feats_df[,-c(37, 38)]), method.dist = "euclid", method.hclust = "centroid", nboot = 10)
result <- pvclust(t(feats_df[,-c(37, 38)]), method.dist = "cor", method.hclust = "average", nboot = 10)

png("plot_pvclust_allfeats.png", width=1800, height="800")

plot(result, labels=names_short)
pvrect(result)

dev.off()



############################
############################

# funcao para imprimir o ccep
# analise de plots (antigo)
gplot_plane <- function(feats_df, D = 6, x = "H6diff", y = "C6diff", xlim = c(0, 1.05), ylim = c(0, 0.5), title = "") {
    # gerando o plot
    p <- gplot.ccep(D = D, xlim = xlim, ylim = ylim, main = title)

    # changing the labs
    p <- p + xlab(expression(Normalized ~ Permutation ~ Entropy ~ (H[S]))) +
        ylab(expression(Statistical ~ Complexity ~ (C[JS])))

    # adjusting layout and centering title
    p <- p + theme_bw(base_size = 24) + theme(plot.title = element_text(hjust = 0.5))

    # adding points
    p <- p + geom_point(aes_string(x = x, y = y, colour = "label"), size = 4, shape = 19, data = feats_df) + scale_colour_discrete(name = "OS")

    # adding text
    p <- p + geom_label_repel(aes_string(x = x, y = y, label = "label", color = "label"),
        segment.alpha = 0.6,
        size = 7,
        label.size = NA,
        point.padding = 0.2,
        box.padding = 0.2,
        show.legend = FALSE,
        max.overlaps = 50,
        data = feats_df
    )

    return(p)
}

#############################
#############################




###############

library(plotly)

data <- read.csv(paste(dir_samples, files[3], sep = ""))

data_diff <- diff(data[, 1])
x <- 2:length(data_diff)

plot_ly(x = data_diff[x], y = data_diff[x - 1], type = "scatter", mode = "markers")


# create dataset of entropy/complexity
##############
tau <- 1
D_l <- c(2, 3, 4, 5, 6, 7)
dataset <- matrix(nrow = 0, ncol = (12 + 12 + 12 + 2))
colnames(dataset) <- c(
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
    "label", "class"
)


file <- read.csv("../tcp-isn/data/switches_dataset.csv")
classes <- unlist(lapply(strsplit(file$label, "-"), "[", 1))
for (i in 1:nrow(file)) {
    print(i)
    row <- t(as.numeric(file[i, ]))

    rows <- head(row[1, ], -1)

    dataset_i <- c()

    # computing features for all D's
    for (D in D_l)
    {
        # bandt-pomping the series
        bp <- bandt_pompe_distribution(rows, D = D)
        HP <- shannon_entropy(bp$probabilities, normalized = TRUE)
        SC <- complexity(bp$probabilities, HP)
        FI <- fisher_information(bp$probabilities)

        dataset_i <- c(dataset_i, HP, SC, FI)

        # bandt-pomping the diff of series
        bp <- bandt_pompe_distribution(diff(rows), D = D)
        HP <- shannon_entropy(bp$probabilities, normalized = TRUE)
        SC <- complexity(bp$probabilities, HP)
        FI <- fisher_information(bp$probabilities)

        dataset_i <- c(dataset_i, HP, SC, FI)
    }

    # adding the label and class
    dataset_i <- c(dataset_i, file$label[i], classes[i])

    # adding the features to the matrix
    dataset <- rbind(dataset, dataset_i)
}

# ajuste das features
dataset_df <- as.data.frame(cbind(apply(dataset[, 1:36], 2, as.numeric)))
dataset_df$label <- as.factor(dataset[, 37])
dataset_df$class <- as.factor(dataset[, 38])
