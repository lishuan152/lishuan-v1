rm(list = ls())
getwd()
setwd("")
otu_table<-read.csv('C.csv',header=T,row.names=1)
n
rel_abundance <- apply(otu_table, 2, function(x) x/sum(x))
rel_abundance
write.csv(rel_abundance, 'non-GC_rel_abundance.csv')
mean_rel_abundance <- rowMeans(rel_abundance)
low_rel_abundance_otu <- rownames(otu_table)[mean_rel_abundance < 0.001]
otu_table_filtered <- otu_table[!(rownames(otu_table) %in% low_rel_abundance_otu), ]
freq <- apply(otu_table_filtered, 1, function(x) sum(x > 0)/length(x))
keep <- freq >= 1/3
otu_table_filt <- otu_table_filtered[keep, ]
write.csv(otu_table_filt, 'non-GC_otu_table_filt.csv')
otu<-otu_table_filt
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
n
BiocManager::install("impute")
BiocManager::install("preprocessCore")
n
library(impute)
library(GO.db)
library(WGCNA)
library(psych)
library(reshape2)
library(igraph)
cor = corAndPvalue(t(otu),y=NULL,use = "pairwise.complete.obs", 
                   alternative='two.sided',method='spearman')
r = cor$cor
p = cor$p
p = p.adjust(p, method = 'BH')
r[p > 0.001 | abs(r) < 0.60] = 0
write.csv(data.frame(r, check.names = FALSE), 'corr.matrix.csv')
g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) 
g
g = delete_vertices(g, names(degree(g)[degree(g) == 0]))
E(g)$corr = E(g)$weight
E(g)$weight = abs(E(g)$weight)
tax = read.csv('C tax.csv', row.names=1, header=T)
tax = tax[as.character(V(g)$name), ] 
V(g)$Kingdom = tax$Kingdom 
V(g)$Phylum = tax$Phylum 
V(g)$Class = tax$Class
V(g)$Order = tax$Order
V(g)$Family = tax$Family
V(g)$Genus = tax$Genus
V(g)$Species = tax$Species
node_list = data.frame(
  label = names(V(g)),
  kingdom = V(g)$Kingdom,
  phylum = V(g)$Phylum,
  class = V(g)$Class,
  order = V(g)$Order,
  family = V(g)$Family,
  genus=V(g)$Genus,
  species = V(g)$Species)
head(node_list)
write.csv(node_list, 'network.node_list.csv')
edge = data.frame(as_edgelist(g))
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)
head(edge_list)
write.csv(edge_list, 'network.edge_list.csv')
write_graph(g, 'network.graphml', format = 'graphml')
nodes_num = length(V(g))
edges_num = length(E(g))
positive.cor_num = sum(E(g)$corr>0)
negative.cor_num = sum(E(g)$corr<0)
average_degree = mean(degree(g))
average_path_length = average.path.length(g, directed = FALSE)
network_diameter = diameter(g, directed = FALSE)
network_density = edge_density(g)
clustering_coefficient = transitivity(g)
network_parameter = data.frame(nodes_num,
                               edges_num,
                               positive.cor_num,
                               negative.cor_num,
                               average_degree,
                               average_path_length,
                               network_diameter,
                               network_density,
                               clustering_coefficient
)
network_parameter
write.csv(network_parameter, 'network_parameter.csv')
otu1 = otu
otu1[otu1>0] = 1
write.csv(otu1, 'adjacent_matrix.csv')