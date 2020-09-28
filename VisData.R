library(tidyverse)

pairwise <- read_csv("Data/pairwise.csv")
blups <- read_csv("Data/blups.csv")

pairwise %>%
      ggplot() +
      geom_histogram(aes(x=btilde - b1), binwidth=0.5) + 
      theme_bw() + xlab("Difference") +
  theme(panel.grid.minor=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'lines'),
        panel.grid.major = element_line(colour = "gray80"),
        axis.title.x=element_text(size=16, vjust=0),
        axis.title.y=element_text(size=16, vjust=1, angle=90),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        strip.text.y=element_text(size=12),
        strip.text.x=element_text(size=12))
ggsave("HistDiff.pdf")

pairwise %>%
      filter(n==50) %>%
      ggplot() +
      geom_point(aes(x=btilde , y=b1), alpha=0.5) + 
      theme_bw() + xlab("Pairwise") + ylab("Model") +
  theme(panel.grid.minor=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'lines'),
        panel.grid.major = element_line(colour = "gray80"),
        axis.title.x=element_text(size=16, vjust=0),
        axis.title.y=element_text(size=16, vjust=1, angle=90),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        strip.text.y=element_text(size=12),
        strip.text.x=element_text(size=12))
ggsave("PairwiseDiff.pdf")



###### SAMPLING 

stacked <- blups %>%
  group_by(cluster) %>%
  summarise("Cluster_SRS"=mean(Cluster_SRS),
            "Cluster_SRSP"=mean(Cluster_SRSP),
            "Cluster_Size"=mean(Cluster_Size),
            "Cluster_Latent"=mean(Cluster_Latent),
            "V2_SRS" = mean(V2_SRS),
            "V2_SRSP" = mean(V2_SRSP),
            "V2_Size" = mean(V2_Size),
            "V2_Latent" = mean(V2_Latent),
            "NoC_SRS" = mean(NoC_SRS),
            "NoC_SRSP" = mean(NoC_SRSP),
            "NoC_Size" = mean(NoC_Size),
            "NoC_Latent" = mean(NoC_Latent))

full_output <- stacked %>%
  merge(cbind(BLUPS, "cluster"=1:N), by="cluster") ### BLUPS are from blupscluster doc
#  merge(cbind(pairwise_BLUPS, "cluster"=1:N), by="cluster")

full_output %>%
  ggplot(aes(x=NoC_SRS, y=BLUPS)) +
  geom_point() + theme_bw() + xlab("No cluster weights SRS") +
  theme(panel.grid.minor=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'lines'),
        panel.grid.major = element_line(colour = "gray80"),
        axis.title.x=element_text(size=16, vjust=0),
        axis.title.y=element_text(size=16, vjust=1, angle=90),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        strip.text.y=element_text(size=12),
        strip.text.x=element_text(size=12))
ggsave("SRS.pdf")


full_output %>%
  ggplot(aes(x=NoC_Size, y=BLUPS)) +
  geom_point() + theme_bw()+ xlab("No cluster weights Sample wrt Size") +
  theme(panel.grid.minor=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'lines'),
        panel.grid.major = element_line(colour = "gray80"),
        axis.title.x=element_text(size=16, vjust=0),
        axis.title.y=element_text(size=16, vjust=1, angle=90),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        strip.text.y=element_text(size=12),
        strip.text.x=element_text(size=12))

ggsave("SRSP.pdf")

full_output %>%
  ggplot(aes(x=NoC_Latent, y=BLUPS)) +
  geom_point() + theme_bw()   + xlab("No cluster weights Sample wrt outcome") +
  theme(panel.grid.minor=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'lines'),
        panel.grid.major = element_line(colour = "gray80"),
        axis.title.x=element_text(size=16, vjust=0),
        axis.title.y=element_text(size=16, vjust=1, angle=90),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        strip.text.y=element_text(size=12),
        strip.text.x=element_text(size=12))
ggsave("latent.pdf")
##### Version 2 VS BLUPS
full_output %>%
  ggplot(aes(x=V2_SRS, y=BLUPS)) +
  geom_point() + theme_bw()+ xlab("Version 2 SRS") +
  theme(panel.grid.minor=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'lines'),
        panel.grid.major = element_line(colour = "gray80"),
        axis.title.x=element_text(size=16, vjust=0),
        axis.title.y=element_text(size=16, vjust=1, angle=90),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        strip.text.y=element_text(size=12),
        strip.text.x=element_text(size=12))
ggsave("v2SRS.pdf")


full_output %>%
  ggplot(aes(x=V2_Size, y=BLUPS)) +
  geom_point() + theme_bw()+ xlab("Version 2 Sample wrt Size") +
  theme(panel.grid.minor=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'lines'),
        panel.grid.major = element_line(colour = "gray80"),
        axis.title.x=element_text(size=16, vjust=0),
        axis.title.y=element_text(size=16, vjust=1, angle=90),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        strip.text.y=element_text(size=12),
        strip.text.x=element_text(size=12))
ggsave("v2SRSP.pdf")

full_output %>%
  ggplot(aes(x=V2_Latent, y=BLUPS)) +
  geom_point() + theme_bw() + xlab("Version 2 Sample wrt outcome") +
  theme(panel.grid.minor=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'lines'),
        panel.grid.major = element_line(colour = "gray80"),
        axis.title.x=element_text(size=16, vjust=0),
        axis.title.y=element_text(size=16, vjust=1, angle=90),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.key=element_blank(),
        strip.text.y=element_text(size=12),
        strip.text.x=element_text(size=12))
ggsave("v2latent.pdf")
