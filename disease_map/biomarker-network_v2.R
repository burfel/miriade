# old file name : network_graph_attempt_v2 (2021-06-18)

library(OmnipathR)
library(dplyr)
library(igraph)

op <- OmnipathR::import_omnipath_interactions()
#op_all <- OmnipathR::import_all_interactions()

#######################################################
# Prepare the edges (using the names still gets also  #
# the uniprots. It's just the uniprots are in cols 1  #
# and 2, and the names in cols 3 and 4)               #
#######################################################

#######################################################
# Alternative 1 - With gene names                     #
#######################################################

# exchanged U-PAR by PLAUR (gene)
genes <- c("PEBP1", "DDAH1", "SPON1", "SLITRK2", "CLSTN3", "SEZ6L", "PLAUR", "DDC", "CRH", "MMP1", "ENO2")
# All Genes
genes_ad <- c("PEBP1", "DDAH1", "SPON1")
# for AD
genes_ftd <- c("SLITRK2", "CLSTN3", "SEZ6L", "PLAUR", "SPON1")
# for FTD
genes_dlb <- c("DDC", "CRH", "MMP1", "ENO2")
# for DLB

opout <- op[op$source_genesymbol %in% genes | op$target_genesymbol %in% genes,]

ig <- graph_from_edgelist(as.matrix(opout[,3:4]))

#######################################################
# Alternative 2 - With uniprots instead of names      #
#######################################################

genes <- c("P30086", "O94760", "Q9HCB6", "Q9H156", "Q9BQT9", "Q9BYH1", "Q03405", "P20711", "P06850", "P03956", "P09104")
# All Genes
genes_ad <- c("P30086", "O94760", "Q9HCB6")
# for AD
genes_ftd <- c("Q9H156", "Q9BQT9", "Q9BYH1", "Q03405", "Q9HCB6")
# for FTD
genes_dlb <- c("P20711", "P06850", "P03956", "P09104")
# for DLB

# 25NOV2022
genes_by_type <- list("AD" = genes_ad, "FTD"=genes_ftd, "DLB"=genes_dlb)
# 25NOV2022
type_list <- list("AD","FTD","DLB")

opout <- op[op$source %in% genes | op$target %in% genes,]

ig <- graph_from_edgelist(as.matrix(opout[,1:2]))

## NEW 25NOV2022 - no longer working
#draw_igraph(ig, genes_by_type, type_list)

# 6DEC2022
genes_by_types <- generate_all_vertex_type_combinations(genes_by_type)
color_seq <- generate_color_sequence(length(genes_by_types) + 1)

#######################################################
# Legend for colors - color definitions               #
#######################################################

# Define the column names.
colnames = c("name", "color")
color_map<-c("AD", "red", "FTD", "green", "DLB", "cyan", "AD+FTD", "yellow")
names_and_colors <- matrix(color_map, nrow = 4, byrow = TRUE, 
                           dimnames = list(NULL, colnames))
legend_mat <- paste(c(names_and_colors[,1]), c(names_and_colors[,2]))

#######################################################
# Creating graph and adjusting node sizes and colors  #
#######################################################

#V(ig)$colour <- "lightgrey"
V(ig)$color <- "lightgrey"
V(ig)$size <- 16

V(ig)[V(ig)$name %in% genes_ad]$color <- "red"
V(ig)[V(ig)$name %in% genes_ad]$size <- 20
#V(ig)$size <- 30
V(ig)[V(ig)$name %in% genes_ftd]$color <- "green"
V(ig)[V(ig)$name %in% genes_ftd]$size <- 20

V(ig)[V(ig)$name %in% genes_dlb]$color <- "cyan"
V(ig)[V(ig)$name %in% genes_dlb]$size <- 20

V(ig)[V(ig)$name %in% genes_ad & V(ig)$name %in% genes_ftd]$color <- "yellow"

#######################################################
# Creating graph layout and plotting it               #
#######################################################

#plot(ig)

# this ensures the starting random position is the same
# for the layouts that use a random starting position
set.seed(1492) 

#minC <- rep(-Inf, vcount(ig))
#maxC <- rep(Inf, vcount(ig))
#minC[1] <- maxC[1] <- 0
#l <- layout.fruchterman.reingold(ig, niter=5000, 
 #                                minx=minC, maxx=maxC,
  #                               miny=minC, maxy=maxC,
   #                              weights=rep(0.01,vcount(ig)))

l <- layout_with_dh(ig, coords = NULL, maxiter = 10, 
               fineiter = max(10, log2(vcount(ig))), 
               cool.fact = 0.75, weight.node.dist = 19,
               weight.border = 0, 
               #weight.edge.lengths = edge_density(ig)/17,
               weight.edge.lengths = edge_density(ig)/22,
               weight.edge.crossings = 100 - sqrt(edge_density(ig)),
               weight.node.edge.dist = 0.2* (1 - edge_density(ig)))

plot.igraph(ig, layout=l, 
            edge.width=1, edge.arrow.size = 0.05, #edge.length = 10, 
            vertex.label.cex=0.6, 
            vertex.label.family="Helvetica",
            vertex.label.font=1,
            vertex.shape="circle",)
            #edge.arrow.mode=2,)
            #edge.curved=TRUE,)


legend('topleft', legend = legend_mat)

write_graph(ig, "bm-network.graphml", format = c("graphml"))

##############################################################################
##############################################################################
# Getting 2nd degree interactions within the existing vertices               #
##############################################################################
##############################################################################

#
# pseudocode:
# opnew <- edges in op where both vertices are in V(ig)
#
vertices <- V(ig)$name
#######################################################
# Alternative 1 - With gene names                     #
#######################################################
op_seconddeg <- op[op$source_genesymbol %in% vertices & op$target_genesymbol %in% vertices,]
ig_seconddeg <- graph_from_edgelist(as.matrix(op_seconddeg[,3:4]))

#######################################################
# Alternative 2 - With uniprots instead of names      #
#######################################################
op_seconddeg <- op[op$source %in% vertices & op$target %in% vertices,]
ig_seconddeg <- graph_from_edgelist(as.matrix(op_seconddeg[,1:2]))

#######################################################
# Creating graph and adjusting node sizes and colors  #
#######################################################

V(ig_seconddeg)$color <- "lightgrey"
V(ig_seconddeg)$size <- 4

V(ig_seconddeg)[V(ig_seconddeg)$name %in% genes_ad]$color <- "red"
V(ig_seconddeg)[V(ig_seconddeg)$name %in% genes_ad]$size <- 5
#V(ig_seconddeg)$size <- 30
V(ig_seconddeg)[V(ig_seconddeg)$name %in% genes_ftd]$color <- "green"
V(ig_seconddeg)[V(ig_seconddeg)$name %in% genes_ftd]$size <- 5

V(ig_seconddeg)[V(ig_seconddeg)$name %in% genes_dlb]$color <- "blue"
V(ig_seconddeg)[V(ig_seconddeg)$name %in% genes_dlb]$size <- 5

V(ig_seconddeg)[V(ig_seconddeg)$name %in% genes_ad & V(ig_seconddeg)$name %in% genes_ftd]$color <- "yellow"

#######################################################
# Creating graph layout and plotting it               #
#######################################################

# this ensures the starting random position is the same
# for the layouts that use a random starting position
set.seed(1492) 

l2 <- layout_with_dh(ig_seconddeg, coords = NULL, maxiter = 10, 
                    fineiter = max(10, log2(vcount(ig_seconddeg))), 
                    cool.fact = 0.75, weight.node.dist = 19,
                    weight.border = 0, 
                    #weight.edge.lengths = edge_density(ig_seconddeg)/17,
                    weight.edge.lengths = edge_density(ig_seconddeg)/22,
                    weight.edge.crossings = 100 - sqrt(edge_density(ig_seconddeg)),
                    weight.node.edge.dist = 0.2* (1 - edge_density(ig_seconddeg)))

plot.igraph(ig_seconddeg, layout=l2, 
            edge.width=1, edge.arrow.size = 0.02, #edge.length = 10, 
            vertex.label.cex=0.2, 
            vertex.label.family="Helvetica",
            vertex.label.font=1,
            vertex.shape="circle",)

write_graph(ig_seconddeg, "bm-network-second-degree-light.graphml", format = c("graphml"))

##############################################################################
##############################################################################
# Getting 2nd degree interactions in general                                 #
##############################################################################
##############################################################################

#
# pseudocode:
# opnew <- edges in op where one or both vertices are in V(ig)
#
#######################################################
# Alternative 1 - With gene names                     #
#######################################################
op_seconddeg_hardcore <- op[op$source_genesymbol %in% vertices | op$target_genesymbol %in% vertices,]
ig_seconddeg_hardcore <- graph_from_edgelist(as.matrix(op_seconddeg_hardcore[,3:4]))

#######################################################
# Alternative 2 - With uniprots instead of names      #
#######################################################
opseconddeg_hardcore <- op[op$source %in% vertices | op$target %in% vertices,]
ig_seconddeg_hardcore <- graph_from_edgelist(as.matrix(op_seconddeg_hardcore[,1:2]))


#######
# FTD
#SLITRK2
#CLSTN3
#SEZ6L
#U-PAR
#(SPON1)


# DLB
#DDC, CRH, MMP-1, ENO2

# all dementias
#CHIT1, AQP4


#%%---------
#  Biomarker

#ABL1
#THOP1
#ITGB2
#CYR61

#SP29
#KYAT1
#MMP-10
#ITGAM


#NUDT5
#SOD2
#BLM hydrolase
#GLO1
#GH
#STX8