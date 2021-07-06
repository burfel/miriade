# visualisation of interactons of biomarkers, stringDB (does not include omnipath)
# Author: Felicia Burtscher
# Date:   28-06-2021

library(dplyr)
library(igraph)
library(httr)

#######################################################
# Getting the json/tsv frin STRING-db with a simple   #
# get request. The species=9606 tells it to only      #
# return human interactions. It uses the names        #
# separated by %0d, and required score of 40          #
#######################################################

# Using network does by default only interactions between the list you provide
# using the add_nodes parameter allow you to extend the interaction neighborhood
STRINGdbJSON <- GET('http://string-db.org/api/json/network?identifiers=PEBP1%0dDDAH1%0dSPON1%0dSLITRK2%0dCLSTN3%0dSEZ6L%0dPLAUR%0dDDC%0dCRH%0dMMP1%0dENO2&required_score=40&species=9606')
# using interaction_partners does by default interactions between the list you
# provide and EVERYTHING THEY HAVE
STRINGdbJSON <- GET('http://string-db.org/api/json/interaction_partners?identifiers=PEBP1%0dDDAH1%0dSPON1%0dSLITRK2%0dCLSTN3%0dSEZ6L%0dPLAUR%0dDDC%0dCRH%0dMMP1%0dENO2&required_score=40&species=9606')
content(STRINGdbJSON,'text')

#STRINGdbTSV <- GET('https://string-db.org/api/tsv/interaction_partners?identifiers=PEBP1%0dDDAH1%0dSPON1%0dSLITRK2%0dCLSTN3%0dSEZ6L%0dPLAUR%0dDDC%0dCRH%0dMMP1%0dENO2&required_score=950&species=9606')
#
# Adjustable request: Adjust the score_parameter, paste it into the request url and send the request
# can do a similar thing with the limit parameter
#
score_parameter <- 950
limit_parameter <- 100
request_string <- gsub(" ", "", 
                       paste('https://string-db.org/api/tsv/interaction_partners?identifiers=PEBP1%0dDDAH1%0dSPON1%0dSLITRK2%0dCLSTN3%0dSEZ6L%0dPLAUR%0dDDC%0dCRH%0dMMP1%0dENO2&required_score=',
                             score_parameter,'&limit=',limit_parameter,'&species=9606'))
STRINGdbTSV <- GET(request_string)
#content(STRINGdbTSV,'text')
STRINGdbTSVParsed <- content(STRINGdbTSV,'parsed')

#######################################################
# We shall extract the number of edges and nodes      #
# per required_score                                  #
#######################################################

# manually first
# Depending on the required_score, here
# are the number of edges that are returned for all sorts
# of rquired_score parameter:
# 950: 81
# 40,100,150: 10699
# 400: 1619
# 200: 6489
# 180: 8179
#
scoreToEdges <- matrix(c(950, 81, 40, 10699, 100, 10699, 150, 10699,  400, 1619, 200, 6489, 180, 8179,
                         160, 9812,155,10248,0,10699,250,4068,300,2846), 
                       byrow = TRUE, ncol = 2,
                       dimnames = list(NULL, c("required_score", "num_of_edges")))
plot(scoreToEdges, xaxp=c(0, 1000, 20))

# manually again
required_score <- c(950, 150, 400, 200, 180, 
                    160, 155,250,300)
num_of_edges <- c(81, 10699, 1619, 6489, 8179, 
                  9812, 10248,4068,2846)
lo <- loess(num_of_edges~required_score)
plot(required_score,num_of_edges)
lines(predict(lo), col='red', lwd=2)

# It seems to me obvious that it starts to get interesting at 150.
# We will create a vector for required_scores, and try to extract
# the num of edges and nodes for it..
required_score <- c(150,155,160,170,180,190,200,225,250,275,300,350,400,
                    450,550,650,750,850,950,1000)

num_of_edges <- character(length(required_score))
num_of_nodes <- character(length(required_score))
for (i in 1:length(required_score)) {
  # We use paste to concatenate the url with the score, and gsub to remove the whitespace added by paste
  request_string <- gsub(" ", "", 
                         paste('https://string-db.org/api/tsv/interaction_partners?identifiers=PEBP1%0dDDAH1%0dSPON1%0dSLITRK2%0dCLSTN3%0dSEZ6L%0dPLAUR%0dDDC%0dCRH%0dMMP1%0dENO2&required_score=',
                               required_score[i],'&species=9606'))
  tempTSV <- content(GET(request_string),'parsed')
  num_of_edges[i] <- nrow(tempTSV)
  # Attempting to extract the unique values.
  #unique(data.frame(output = Reduce(c, (apply(tempTSV[,3:4], 1, c)))))
  # So, googling and some trial and error got us the following:
  # Reduce(c, (apply(tempTSV[,3:4], 1, c))) will take
  # the 3rd and 4th columns of tempTSV, and will reduce them into
  # a long vector. unique will extract the unique values.
  # And finally, length gets the actual number
  num_of_nodes[i] <- length(unique(Reduce(c, (apply(tempTSV[,3:4], 1, c)))))
}

#######################################################
# Plot the number of edges per required_score         #
#######################################################

lo <- loess(num_of_edges~required_score)
plot(required_score,num_of_edges, xaxp=c(150, 1000, 17))
lines(predict(lo), col='red', lwd=2)

#######################################################
# Plot the number of nodes per required_score         #
#######################################################

lo <- loess(num_of_nodes~required_score)
plot(required_score,num_of_nodes, xaxp=c(150, 1000, 17))
lines(predict(lo), col='red', lwd=2)

# BIG ISSUE: It seems that STRING-db returns edges without
# information about directedness. igraph assumes it is
# directed, and graphs it with outgoing arrows
# from our set of proteins

#######################################################
# With gene names                                     #
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

sdbGraph <- graph_from_edgelist(as.matrix(STRINGdbTSVParsed[,3:4]))
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

V(sdbGraph)$color <- "lightgrey"
V(sdbGraph)$size <- 16

V(sdbGraph)[V(sdbGraph)$name %in% genes_ad]$color <- "red"
V(sdbGraph)[V(sdbGraph)$name %in% genes_ad]$size <- 20

V(sdbGraph)[V(sdbGraph)$name %in% genes_ftd]$color <- "green"
V(sdbGraph)[V(sdbGraph)$name %in% genes_ftd]$size <- 20

V(sdbGraph)[V(sdbGraph)$name %in% genes_dlb]$color <- "cyan"
V(sdbGraph)[V(sdbGraph)$name %in% genes_dlb]$size <- 20

V(sdbGraph)[V(sdbGraph)$name %in% genes_ad & V(sdbGraph)$name %in% genes_ftd]$color <- "yellow"

#######################################################
# Creating graph layout and plotting it               #
#######################################################

# this ensures the starting random position is the same
# for the layouts that use a random starting position
set.seed(1492) 

l <- layout_with_dh(sdbGraph, coords = NULL, maxiter = 10, 
                    fineiter = max(10, log2(vcount(sdbGraph))), 
                    cool.fact = 0.75, weight.node.dist = 19,
                    weight.border = 0,
                    weight.edge.lengths = edge_density(sdbGraph)/22,
                    weight.edge.crossings = 100 - sqrt(edge_density(sdbGraph)),
                    weight.node.edge.dist = 0.2* (1 - edge_density(sdbGraph)))

plot.igraph(sdbGraph, layout=l, 
            edge.width=1, edge.arrow.size = 0.05, 
            vertex.label.cex=0.6, 
            vertex.label.family="Helvetica",
            vertex.label.font=1,
            vertex.shape="circle",)

legend('topleft', legend = legend_mat)

write_graph(sdbGraph, "bm-network-stringdb.graphml", format = c("graphml"))

