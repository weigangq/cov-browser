library(igraph)
library(qgraph)
library(plotrix)


snp <- read.table("ResearchProjects/covid19/edges_and_mutation.tsv", header = T) 
edges <- as.matrix(snp[,1:2]) #extract edge list
g <- graph_from_edgelist(edges) #build the graph
degree_out <- as.matrix(degree(g, mode = "out"))

#design color gradient of nodes
color_grad <- color.scale(out_degree,0,c(0.1,1),1,color.spec = "hsv") #red
color_grad <- color.scale(out_degree,c(0.2,0),c(0.4,1),1,color.spec = "hsv") #yellow to red

#assign edge color
edge_col <- c(rep('grey',120))
edge_col[which(snp$syn_nonsyn==0)] <- 'blue'

#plot the graph, color nodes by outgoing degree, color edges by aa change
pdf("ResearchProjects/covid19/hap_network_covid.png")
tkplot(g,layout=layout_with_fr,
       vertex.size=15,
       vertex.label.cex=.8,
       vertex.color=color_grad,
       edge.color = edge_col,
       edge.arrow.size=.7,
       edge.arrow.width=.4)

#plot using qgraph
qgraph(edges, color=color_grad, vsize=5, esize=2, edge.color= edge_col)

#make a color palette as the legend 
plot(0:10,type="n",axes = F)
gradient.rect(0,0,3,3,col = color.scale(1:10,c(0.2,0),c(0.4,1),1,color.spec = "hsv"),gradient="x")

library(reshape2)
library(ggplot2)

###Regression
snp_sub <- snp[,c("from","to","syn_nonsyn")]
mutation <- dcast(snp_sub, from~syn_nonsyn)

vertex_list <- get.vertex.attribute(g)$name
nodes <- data.frame(node=vertex_list, degreeOUT= as.vector(degree_out))

#mutation destribution in different loci
locus <- as.data.frame(table(snp$locus))
d <- dcast(snp, locus~syn_nonsyn)
locus <- cbind(locus,d[,2:4])
colnames(locus) <- c("locus", "mut_num","nonsyn","syn","unknown")
mut_loc <- read.table("ResearchProjects/covid19/mutation_locus.txt",header = T)
loc_len <- as.vector(c(200,200,200,700,1300,80,21000,1000,200,4000))
loc_len <- data.frame(locus=locus$locus, length=loc_len)
mut_loc <- merge(loc_len, mut_loc, by="locus")
mut_loc$normal <- mut_loc$number/mut_loc$length
ggplot(data = mut_loc, aes(x=locus, y=normal, fill=mutation))+
  geom_bar(stat = "identity")


