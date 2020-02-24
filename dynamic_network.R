library(vegan)
library(network)
library(networkDynamic)
require(ndtv)

# On charge les données 
data(varespec)
data(varechem)

spe <- varespec # matrice d'espèce
chem <- varechem # matrice de variables environnementales


# Fenêtre glissante -------------------------------------------------------

# ici on calcule des réseau le long d'une fenêtre glissante, et on stock la pile de réseaux dans une liste : dynnet
dynnet <- list() 
window <- 10 
newWtdorder <- order(chem$pH)
EnvM <- min(chem$pH)
EnvMs <- vector()

for(i in 1:(nrow(chem)- window +1)){
  indexes <- i:(i + window -1) 
  cat(mean(chem$pH[newWtdorder][indexes]), "\n")
  EnvMs[i] <- mean(chem$pH[newWtdorder][indexes])
  subspe <- spe[,indexes]
  adjmat <- cor(subspe) # on calcule une matrice d'adjacence sur la base des corrélations. Il y a de bien meilleures façon de faire mais c'est pour l'exemple
  adjmat[abs(adjmat) >= 0.6] <- 1 #on filtre
  adjmat[abs(adjmat) < 0.6] <- 0
  
  
  colnames(adjmat) <- rownames(adjmat) <- colnames(subspe)
  
  g <- network(adjmat, directed = FALSE)

  set.vertex.attribute(x = g, attrname = "PID", value = network.vertex.names(g))
  dynnet[[i]] <- g
}

# on crée un réseau dynamique
adjmat.dyn <- networkDynamic(network.list=dynnet,onsets=EnvMs*10, termini=append(EnvMs[-1],3.15)*10)
get.change.times(adjmat.dyn)

#on fait une vido qui tabasse
render.animation(adjmat.dyn)
saveVideo(ani.replay(),video.name="networks_movie_step10.mp4",
          ani.width=1680,ani.height=1050)


timeline(adjmat.dyn, xlab="WTD", plot.vertex.spells = FALSE, e.label=FALSE, edge.col="blue", ylab="Edges", main="Edges spells phase plot over WTD gradient",edge.lwd=0.1)
filmstrip(adjmat.dyn,displaylabels=FALSE, frames = 12, mfrow=c(4,3))
