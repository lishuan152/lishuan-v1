rm(list = ls())
setwd()
otu <- read.delim(' .txt', row.names = 1, stringsAsFactors = FALSE)
head(otu)
otu <- data.frame(t(otu))
group <- read.delim(' .txt', stringsAsFactors = FALSE)
head(group)
install.packages("vegan")
library(vegan)
seed <- sample(.Machine$integer.max, size=1)
seed
set.seed(seed)
adonis_group <- adonis2(otu~group, group, distance = 'bray', permutations = 999)
adonis_group
bray_dis <- vegdist(otu, method = 'bray')
pcoa <- cmdscale(bray_dis, k = 3, eig = TRUE)
pcoa_exp <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])
pcoa1 <- paste0('PCoA1: ', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste0('PCoA2: ', round(100*pcoa_exp[2], 2), '%')
site <- data.frame(pcoa$point)[1:2]
site$sample <- rownames(site)
bray_dis
beta_disp <- betadisper(bray_dis, group = group)
print(beta_disp)
distances_to_centroid <- beta_disp$distances
head(distances_to_centroid)
anova(beta_disp)
TukeyHSD(beta_disp)
site <- merge(site, group, by = 'sample')
names(site)[2:3] <- c('pcoa1', 'pcoa2')
install.packages("ggplot2")
install.packages("ggrepel")
library(ggplot2)
library(ggrepel)
p <- ggplot(data = site) + 
  geom_point(aes(x = pcoa1, y = pcoa2, color = group), size = 2) + 
  geom_text_repel(aes(x = pcoa1, y = pcoa2, label = sample, color = group), size = 4,
  box.padding = unit(0.4, 'lines'), max.overlaps = 40,show.legend = FALSE) + 
  scale_color_manual(limits = c('A', 'B','C'), values = c("#EB7C68","#F8BA4F","#47BBC6")) + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color = 'black'), legend.key = element_blank()) + 
        labs(x = pcoa1, y = pcoa2, color =group)
p
#The following code calculates the beta distances between the three samples within each group based on the distance matrix calculated by the aforementioned code.
dist_matrix <- matrix(0, nrow = 9, ncol = 9)
rownames(dist_matrix) <- colnames(dist_matrix) <- c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
dist_matrix[2,1] <- 0.12936
dist_matrix[3,1] <- 0.11852; dist_matrix[3,2] <- 0.16081
dist_matrix[4,1] <- 0.90208; dist_matrix[4,2] <- 0.89221; dist_matrix[4,3] <- 0.90067
dist_matrix[5,1] <- 0.90013; dist_matrix[5,2] <- 0.88636; dist_matrix[5,3] <- 0.89048; dist_matrix[5,4] <- 0.16959
dist_matrix[6,1] <- 0.89937; dist_matrix[6,2] <- 0.89568; dist_matrix[6,3] <- 0.90403; dist_matrix[6,4] <- 0.17523; dist_matrix[6,5] <- 0.18239
dist_matrix[7,1] <- 0.97994; dist_matrix[7,2] <- 0.97788; dist_matrix[7,3] <- 0.98699; dist_matrix[7,4] <- 0.91987; dist_matrix[7,5] <- 0.92138; dist_matrix[7,6] <- 0.9242
dist_matrix[8,1] <- 0.98178; dist_matrix[8,2] <- 0.97875; dist_matrix[8,3] <- 0.9884; dist_matrix[8,4] <- 0.93049; dist_matrix[8,5] <- 0.93255; dist_matrix[8,6] <- 0.93461; dist_matrix[8,7] <- 0.16461
dist_matrix[9,1] <- 0.97799; dist_matrix[9,2] <- 0.97517; dist_matrix[9,3] <- 0.98493; dist_matrix[9,4] <- 0.9177; dist_matrix[9,5] <- 0.92052; dist_matrix[9,6] <- 0.92323; dist_matrix[9,7] <- 0.18868; dist_matrix[9,8] <- 0.1889
dist_matrix <- dist_matrix + t(dist_matrix)
diag(dist_matrix) <- 0
bray_dis <- as.dist(dist_matrix)
group <- factor(rep(c("A", "B", "C"), each = 3))
names(group) <- c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
beta_disp <- betadisper(bray_dis, group = group)
print(beta_disp)
set.seed(123)
perm_result <- permutest(beta_disp, permutations = 9999)
print(perm_result)
anova_result <- anova(beta_disp)
print(anova_result)
tukey_result <- TukeyHSD(beta_disp)
print(tukey_result)
distances <- data.frame(
  Sample = names(group),
  Group = group,
  Distance_to_centroid = beta_disp$distances
)
dist_matrix <- matrix(0, nrow = 9, ncol = 9)
rownames(dist_matrix) <- colnames(dist_matrix) <- c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
dist_matrix[2,1] <- 0.73872
dist_matrix[3,1] <- 0.74242; dist_matrix[3,2] <- 0.47845
dist_matrix[4,1] <- 0.98754; dist_matrix[4,2] <- 0.98013; dist_matrix[4,3] <- 0.98182
dist_matrix[5,1] <- 0.98418; dist_matrix[5,2] <- 0.97407; dist_matrix[5,3] <- 0.97912; dist_matrix[5,4] <- 0.80774
dist_matrix[6,1] <- 0.97778; dist_matrix[6,2] <- 0.96465; dist_matrix[6,3] <- 0.96599; dist_matrix[6,4] <- 0.7734; dist_matrix[6,5] <- 0.65017
dist_matrix[7,1] <- 0.97946; dist_matrix[7,2] <- 0.95286; dist_matrix[7,3] <- 0.95623; dist_matrix[7,4] <- 0.97879; dist_matrix[7,5] <- 0.97138; dist_matrix[7,6] <- 0.95522
dist_matrix[8,1] <- 0.98552; dist_matrix[8,2] <- 0.96599; dist_matrix[8,3] <- 0.97071; dist_matrix[8,4] <- 0.98418; dist_matrix[8,5] <- 0.98384; dist_matrix[8,6] <- 0.97576; dist_matrix[8,7] <- 0.75286
dist_matrix[9,1] <- 0.98081; dist_matrix[9,2] <- 0.95421; dist_matrix[9,3] <- 0.95724; dist_matrix[9,4] <- 0.97946; dist_matrix[9,5] <- 0.97677; dist_matrix[9,6] <- 0.95993; dist_matrix[9,7] <- 0.54411; dist_matrix[9,8] <- 0.75859
dist_matrix <- dist_matrix + t(dist_matrix)
diag(dist_matrix) <- 0
bray_dis <- as.dist(dist_matrix)
group <- factor(rep(c("A", "B", "C"), each = 3))
names(group) <- c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
beta_disp <- betadisper(bray_dis, group = group)
print(beta_disp)
set.seed(123)
perm_result <- permutest(beta_disp, permutations = 9999)
print(perm_result)
anova_result <- anova(beta_disp)
print(anova_result)
tukey_result <- TukeyHSD(beta_disp)
print(tukey_result)
distances <- data.frame(
  Sample = names(group),
  Group = group,
  Distance_to_centroid = beta_disp$distances
)