install.packages("plotrix")
###plot CoronaVirus Network- color by degree
###Regression analysis

install.packages("RColorBrewer")

library(RColorBrewer)
library(igraph)



#######Define the Graph######

edgeTable<- read.table("C:/Users/brian/Dropbox/LabShared/edges.tsv",header=T)
dim(edgeTable)



###there exist 85 edges for 86 nodes. is this correct? Yes, it's correct for since there are no loops. 


edgeT<-unique(edgeTable[,1:2])  ###take out repeated rows using unique() function

edgeT<-as.matrix(edgeT)
dim(edgeT)
g<-graph_from_edgelist(edgeT)

#####Degree Of Out Edges for eac vertex######


degree(g, mode="out")                          #gives the degree of each vertice for those that are going out. 
uniquedegrees<-unique(degree(g, mode="out"))   #Use the mode="out" to get the out going nodes.  
#######Use heat.colors to get color assignments, 1 for each unique degree####
dg<-sort(uniquedegrees)                       

ht<-heat.colors(30,rev=T)    #ht a vector consisting of 30 colors.  
#ht
uniqueColors<-ht[2.5*dg+1]  #uniqueColors-one color for each degree. 8 colors total. 
uniqueColors[8]<-ht[25]      # a little adjustment to colors
#length(uniqueColors)

df<-data.frame(sort(uniquedegrees),uniqueColors)    #dataFrame matching degree to color


##Vcolors- assigns a color to each vertice.
Vcolors<-c()   ##colors by degree

for (i in 1:length(V(g))){
  c<-which(df[,1]==as.vector(degree(g,mode="out"))[i])
  Vcolors[i]<-as.character(df[c,2])
}
Vcolors


#####plot with degree by color####
##igraph plot-- it's bad
plot.igraph(g, layout=layout_with_lgl(g),vertex.size=6,xlim=c(-.3,.3),vertex.label.cex=.8,vertex.color=Vcolors,edge.arrow.size=.4,edge.arrow.width=.4)
  legend('topleft',legend=sort(uniquedegrees),pt.cex=2,col='black',pch=21, pt.bg=uniqueColors,title="deg")

##tkplot- it's better but no legend 
  
tkplot(g,vertex.size=15,vertex.label.cex=.8,vertex.color=Vcolors,edge.arrow.size=.7,edge.arrow.width=.4)



  
###########Regression Analysis#########


######Prep The Data########
edgeSplit<-split(edgeTable,edgeTable[,1])
st<-edgeSplit[[6]]


edgeSplit[[6]]
edgeSplitOut<-split(edgeTable,edgeTable[,1])
length(edgeSplitOut)

edgeSplitIn<-split(edgeTable,edgeTable[,2])
edgeSplitIn

length(edgeSplitIn)  #corresponds to the number of vertice take out bat vertex

####3 functions for  4 lapply functions i'll use later####

NonSyn<-function(data){          ###counts nonSynonymous mutations
  
  data[is.na(data)]<-0
  degO<-nrow(data)
  nonSynDeg<-nrow(data)-sum(data$syn_nonsyn)
  return(nonSynDeg)
}

syn<-function(data){          ###counts nonSynonymous mutations
  data[is.na(data)]<-0
  degO<-nrow(data)
  SynDeg<-sum(data$syn_nonsyn)
  return(SynDeg)
}

degO<-function(data){                 ######counts number of out edges
  data[is.na(data)]<-0
  data<-unique(data[,c(1,2)])
  return(nrow(data))
}



NumMutation<-function(data){         #####counts number of mutations for each node
  data[is.na(data)]<-0
  numM<-nrow(data)
  return(numM)
}




#####3 lapply functions

numMut<-unlist(lapply(edgeSplitOut,NumMutation))
degreeOut<-unlist(lapply(edgeSplitOut,degO))
NonS<-unlist(lapply(edgeSplitOut,NonSyn))
Syn<-unlist(lapply(edgeSplitOut,syn))

numMutIn<-unlist(lapply(edgeSplitIn,NumMutation))
#degreeIn<-unlist(lapply(edgeSplitIn,degI))
NonSIn<-unlist(lapply(edgeSplitIn,NonSyn))
SynIn<-unlist(lapply(edgeSplitIn, syn))
####df2-  dataframe for regression models


length(edgeSplitOut)
edgeSplit[[1]]
df2<-data.frame(degreeOut,NonS,numMut,Syn)
df2
df2$diffMut<-df2$numMut-df2$NonS
df2$nodeNames<-names(edgeSplit)
df2$mutRatio<-df2$NonS/df2$numMut
nrow(df2)


df3<-data.frame(NonSIn,numMutIn,SynIn)
df3
df3$diffMut<-df2$numMut-df2$NonS
df3$nodeNames<-names(edgeSplitIn)
df3$mutRatio<-df2$NonS/df2$numMut
df4<-merge(df3,df2,by="nodeNames")
dim(df4)
df4<-df4[,-9]
head(df4)
head(df4)
head(df2)

cordf<-cor(df4[,2:8])
cordf[cordf==1]<-NA
cordf
hist(cordf,ylab="Correlation Values",main="histogram of Correlation Values")
?cor

length(df4[,1])

df2<-df2[-1,]     #take the ST1 -the bat out of regression data set. 
plot(df2$degreeOut,df2$mutRatio)

df2
library(ggplot2)

x <- df2$degreeOut
y <- df2$NonS
lm3 <- lm(y ~x)

xc <- seq(1, 18, by=0.1)
dfc <- data.frame(x=xc)
dfc$lmPred<-predict(lm3,dfc)


###linear regression  ---- NonSynonymousMutationCount~degreeOut
lm2<-lm(NonS~degreeOut,data=df2)

summary(lm2)




###ggplot of linear model

library(ggplot2)

pl<-ggplot(df2[,1:2], aes(x=degreeOut, y=NonS)) +
  geom_point(shape=1) + geom_jitter()   # Use hollow circles
  
pl+geom_smooth(method="lm")+ 
  ylab("#nonsynonymousMutations") + xlab("Degrees Out of a Node")# Add a loess smoothed fit curve with confidence region
head(df2)





###regression experiments - (messy)
glm1<-glm(cbind(NonS,diffMut)~degreeOut,family ="binomial",data=df2)
summary(glm1)
help.search("Lemeshow")
df2

poiss<-glm(NonS~degreeOut,family="poisson",data=df2)

head(df2)
plot(df2$degreeOut,df2$NonS, ylab="#non-synonymous Mutations")

plot(fitted(lm1),residuals(lm1))

qqnorm(residuals(lm1),ylab="Residuals",main="QQnorm For Residuals")
cook<-cooks.distance(lm1)
library(faraway)
halfnorm(cook,2,df2$nodeNames)
lmCooks<-lm(NonS~degreeOut,data=df2, subset=(cook<0.5))
summary(lmCooks)

qqline(residuals(lm1))

