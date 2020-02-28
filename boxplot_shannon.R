library(tidyverse)
library(readxl)

#########################################################
## ID group in fonction of 1st component of PCA on HAP ##
#########################################################
read_excel("PCA HAP T1.xlsx") %>%
  select(Subjid, Group2) %>%
  transmute(Individual = as.character(Subjid), Group=Group2) -> group_df


#######################
## Load Shannon data ##
#######################
Lee_table <- read_excel("AF17Lee.xlsx", col_names = TRUE, range = "A3:D407")

shannon_bact <- read.table("estimates_shannon_bacteria.txt", header=T) %>%
  left_join(Lee_table, by="Sample") %>%
  mutate(Individual = as.character(Individual)) %>%
  left_join(group_df, by="Individual") %>%
  filter(Group != "NA")


shannon_bact %>% head()

shannon_fungi <- read.table("estimates_shannon_fungi.txt", header=T) %>%
  left_join(Lee_table, by="Sample") %>%
  mutate(Individual = as.character(Individual)) %>%
  left_join(group_df, by="Individual") %>%
  filter(Group != "NA")

shannon_fungi <- shannon_fungi[!duplicated(shannon_fungi$Sample),] 

shannon_fungi %>% head()

#############
## Boxplot ##
#############

## Bacteria Cheek
shannon_bact %>%
  filter(Site == "Cheek") %>%
  ggplot(aes(x = factor(Group), y = shannon_est)) +
    geom_boxplot(col="grey30", outlier.shape = NA, width=.5, fill="grey90") +
    geom_point(aes(col = City), position = position_jitter(w = 0.1, h = 0)) +
    scale_color_manual(values = c("#7B58A4", "#8DC63F")) + 
    theme_bw() + ylab("Shannon diversity index") + xlab("Pollution group") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size=12), legend.position = "bottom") + 
    geom_smooth(method = "lm", se=TRUE, aes(group=1)) +
    ggtitle("Repartition of 16S shannon diversity index in pollution group")

## Fungi Cheek
shannon_fungi %>%
  filter(Site == "Cheek") %>%
  ggplot(aes(x = factor(Group), y = shannon_est)) +
    geom_boxplot(col="grey30", outlier.shape = NA, width=.5, fill="grey90") +
    geom_point(aes(col = City), position = position_jitter(w = 0.1, h = 0)) +
    scale_color_manual(values = c("#7B58A4", "#8DC63F")) + 
    theme_bw() + ylab("Shannon diversity index") + xlab("Pollution group") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size=12), legend.position = "bottom") + 
    geom_smooth(method = "lm", se=TRUE, aes(group=1)) +
    ggtitle("Repartition of ITS shannon diversity index in pollution group")

## Bacteria Scalp
shannon_bact %>%
  filter(Site == "Scalp") %>%
  ggplot(aes(x = factor(Group), y = shannon_est)) +
  geom_boxplot(col="grey30", outlier.shape = NA, width=.5, fill="grey90") +
    geom_point(aes(col = City), position = position_jitter(w = 0.1, h = 0)) +
    scale_color_manual(values = c("#7B58A4", "#8DC63F")) + 
    theme_bw() + ylab("Shannon diversity index") + xlab("Pollution group") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size=12), legend.position = "bottom") + 
    geom_smooth(method = "lm", se=TRUE, aes(group=1)) +
    ggtitle("Repartition of 16S shannon diversity index in pollution group")

## Fungi Scalp
shannon_fungi %>%
  filter(Site == "Scalp") %>%
  ggplot(aes(x = factor(Group), y = shannon_est)) +
  geom_boxplot(col="grey30", outlier.shape = NA, width=.5, fill="grey90") +
  geom_point(aes(col = City), position = position_jitter(w = 0.1, h = 0)) +
  scale_color_manual(values = c("#7B58A4", "#8DC63F")) + 
  theme_bw() + ylab("Shannon diversity index") + xlab("Pollution group") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size=12), legend.position = "bottom") + 
  geom_smooth(method = "lm", se=TRUE, aes(group=1)) +
  ggtitle("Repartition of ITS shannon diversity index in pollution group")



#################
## Correlation ##
#################
HAP_score <- read_excel("PCA HAP T1.xlsx") %>%
  mutate(Individual = as.character(Subjid)) %>%
  select(Individual, City, PC1)

## Bacteria Scalp
HAP_score %>%
  left_join(select(shannon_bact, shannon_est,Individual, Site), by="Individual") %>%
  filter(shannon_est != "NA") %>%
  filter(Site == "Scalp") %>%
  ggpubr::ggscatter(x = "PC1", y = "shannon_est", color = "City",
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE # Add confidence interval
  ) + ggpubr::stat_cor(method = "pearson", label.x = 5, label.y =6) + 
  scale_color_manual(values = c("#7B58A4", "#8DC63F")) + 
  ylab("Shannon diversity index") + xlab("Pollution score") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size=12), legend.position = "bottom") +
  ggtitle("Correlation between ITS shannon diversity and HAP score (Bacteri scalp")

## Bacteria Cheek
HAP_score %>%
  left_join(select(shannon_bact, shannon_est, Individual, Site), by="Individual") %>%
  filter(shannon_est != "NA") %>%
  filter(Site == "Cheek") %>%
  ggpubr::ggscatter(x = "PC1", y = "shannon_est", color = "City",
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE # Add confidence interval
  ) + ggpubr::stat_cor(method = "pearson", label.x = 5, label.y =6) + 
  scale_color_manual(values = c("#7B58A4", "#8DC63F")) + 
  ylab("Shannon diversity index") + xlab("Pollution score") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size=12), legend.position = "bottom") +
  ggtitle("Correlation between ITS shannon diversity and HAP score (Bacteria cheek")

## fungi Cheek
HAP_score %>%
  left_join(select(shannon_fungi, shannon_est, Individual, Site), by="Individual") %>%
  filter(shannon_est != "NA") %>%
  filter(Site == "Cheek") %>%
  ggpubr::ggscatter(x = "PC1", y = "shannon_est", color = "City",
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE # Add confidence interval
  ) + ggpubr::stat_cor(method = "pearson", label.x = 5, label.y =6) + 
  scale_color_manual(values = c("#7B58A4", "#8DC63F")) + 
  ylab("Shannon diversity index") + xlab("Pollution score") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size=12), legend.position = "bottom") +
  ggtitle("Correlation between ITS shannon diversity and HAP score (Fungi cheek")

## fungi Cheek
HAP_score %>%
  left_join(select(shannon_fungi, shannon_est, Individual, Site), by="Individual") %>%
  filter(shannon_est != "NA") %>%
  filter(Site == "Scalp") %>%
  ggpubr::ggscatter(x = "PC1", y = "shannon_est", color = "City",
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE # Add confidence interval
  ) + ggpubr::stat_cor(method = "pearson", label.x = 5, label.y =6) + 
  scale_color_manual(values = c("#7B58A4", "#8DC63F")) + 
  ylab("Shannon diversity index") + xlab("Pollution score") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size=12), legend.position = "bottom") +
  ggtitle("Correlation between ITS shannon diversity and HAP score (Fungi cheek")



