
library(tidyverse)
library(dplyr)
library(ggiraph)
library(clusterProfiler)
library(enrichplot)
library(igraph)



# Leer network 

BBC <- read_graph("BasalBreastCancer.graphml", format = "graphml")


# Definir organism para humano
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)

  
#------ LOUVAIN CLUSTERING -------
lc.BBC <- cluster_louvain(BBC) # Encontrar las comunidades
communities(lc.BBC) # Para ver todas las listas de comunidades = 9 comunidades

memb.BBC <- membership(lc.BBC) # Para ver a que comunidad pertence cada gen


# Universo de todos los genes de la red para usar como background
universe_genes <- V(BBC)$name


head(V(BBC)$name)


#---- COMUNIDAD 1 ------ Sí 

genes_c1 <- V(BBC)$name [memb.BBC == 1] # Genes solo de comunidad 1


# Encrich comunidad 1

go_enrich.1 <- enrichGO(gene = genes_c1,
                      universe = universe_genes,
                      OrgDb = organism,
                      keyType = "SYMBOL",
                      readable = TRUE,
                      ont = "MF",
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)


plot_com.1 <- dotplot(go_enrich.1) # Dotplot para visualizar


# Lo mismo con todas las comunidades:

 
#---- COMUNIDAD 2 ------ Sí

genes_c2 <- V(BBC)$name [memb.BBC == 2]

go_enrich.2 <- enrichGO(gene = genes_c2,
                        universe = universe_genes,
                        OrgDb = organism,
                        keyType = "SYMBOL",
                        readable = TRUE,
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

plot_com.2 <- dotplot(go_enrich.2)



#---- COMUNIDAD 3 ------ No

genes_c3 <- V(BBC)$name [memb.BBC == 3]

go_enrich.3 <- enrichGO(gene = genes_c3,
                        universe = universe_genes,
                        OrgDb = organism,
                        keyType = "SYMBOL",
                        readable = TRUE,
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

plot_com.3 <- dotplot(go_enrich.3)


#---- COMUNIDAD 4 ------ No

genes_c4 <- V(BBC)$name [memb.BBC == 4]

go_enrich.4 <- enrichGO(gene = genes_c4,
                        universe = universe_genes,
                        OrgDb = organism,
                        keyType = "SYMBOL",
                        readable = TRUE,
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

plot_com.4 <- dotplot(go_enrich.4)


#---- COMUNIDAD 5 ------ No

genes_c5 <- V(BBC)$name [memb.BBC == 5]

go_enrich.5 <- enrichGO(gene = genes_c5,
                        universe = universe_genes,
                        OrgDb = organism,
                        keyType = "SYMBOL",
                        readable = TRUE,
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

plot_com.5 <- dotplot(go_enrich.5)


#---- COMUNIDAD 6 ------ Sí

genes_c6 <- V(BBC)$name [memb.BBC == 6]

go_enrich.6 <- enrichGO(gene = genes_c6,
                        universe = universe_genes,
                        OrgDb = organism,
                        keyType = "SYMBOL",
                        readable = TRUE,
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

plot_com.6 <- dotplot(go_enrich.6)


#---- COMUNIDAD 7 ------ Sí = 8

genes_c7 <- V(BBC)$name [memb.BBC == 7]

go_enrich.7 <- enrichGO(gene = genes_c7,
                        universe = universe_genes,
                        OrgDb = organism,
                        keyType = "SYMBOL",
                        readable = TRUE,
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

plot_com.7 <- dotplot(go_enrich.7)


#---- COMUNIDAD 8 ------ NO

genes_c8 <- V(BBC)$name [memb.BBC == 8]

go_enrich.8 <- enrichGO(gene = genes_c8,
                        universe = universe_genes,
                        OrgDb = organism,
                        keyType = "SYMBOL",
                        readable = TRUE,
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

plot_com.8 <- dotplot(go_enrich.8)



#---- COMUNIDAD 9 ------ NO

genes_c9 <- V(BBC)$name [memb.BBC == 9]

go_enrich.9 <- enrichGO(gene = genes_c9,
                        universe = universe_genes,
                        OrgDb = organism,
                        keyType = "SYMBOL",
                        readable = TRUE,
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

plot_com.9 <- dotplot(go_enrich.9)













