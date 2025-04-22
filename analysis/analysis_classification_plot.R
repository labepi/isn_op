library(ggplot2)
library(rgl)
library(scatterplot3d)
library(ggsci)
library(scales)

suppressMessages(library(ggforce))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(factoextra))
suppressMessages(library(ggraph))

# loading data
x = read.csv("classification_results.csv")

df = as.data.frame(x)

df$num = as.factor(df$num)

df$acc = round(100*df$acc, digits=2)
df$f1 = round(100*df$f1, digits=2)

# accuracy

p = ggplot(df, aes(x=reorder(num, acc, decreasing=T) ,y=acc,group=Classificador,fill=Classificador)) +
        geom_bar(stat='identity',position=position_dodge(), width=0.5) +
        # geom_errorbar(aes(ymin=acc-sd, ymax=acc+sd), width=.4,
                 # position=position_dodge(.9)) +
        geom_text(aes(label=acc), position=position_dodge(width=0.55), hjust=0, vjust=-0.2, angle=60, size=7) +
        theme_bw(base_size=24) + scale_shape_manual(values=15:24) +
        xlab('Tamanho das Séries') +
        ylab('Acurácia (%)') +
        scale_color_d3(palette="category20") +
        scale_fill_d3(palette="category20") +
        coord_cartesian(ylim=c(0.4, 1.00)) + 
        theme(legend.position="bottom") + 
        guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
        coord_cartesian(ylim=c(50, 110))

ggsave("images/strategy/classification_devices.png", p, width=20, height=6)


# f1-score

p = ggplot(df, aes(x=reorder(num, f1, decreasing=T) ,y=f1,group=Classificador,fill=Classificador)) +
        geom_bar(stat='identity',position=position_dodge(), width=0.5) +
        # geom_errorbar(aes(ymin=acc-sd, ymax=acc+sd), width=.4,
                 # position=position_dodge(.9)) +
        geom_text(aes(label=f1), position=position_dodge(width=0.9), hjust=0, vjust=-0.2, angle=60, size=7) +
        theme_bw(base_size=24) + scale_shape_manual(values=15:24) +
        xlab('Tamanho das Séries') +
        ylab('F1-Score') +
        scale_color_d3(palette="category20") +
        scale_fill_d3(palette="category20") +
        coord_cartesian(ylim=c(0.4, 1.00)) + 
        theme(legend.position="bottom") + 
        guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
        coord_cartesian(ylim=c(50, 110))

        # scale_y_continuous(limits = c(0.8, 1.0), oob = rescale_none)

ggsave("images/strategy/f1score_devices.png", p, width=20, height=6)



#device,num,acc,f1


