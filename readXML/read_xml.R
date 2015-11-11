install.packages("XML")
# setwd("~/Research15/intellij/bnkit/micro_data")
library(XML)
library(stringr)
library(bnlearn)
library(reshape)
library(ggplot2)

bnFile <- xmlParse("ExampleNetwork.new")
root = xmlRoot(bnFile)
nodes <- xmlSize(root[["def"]])
nodeVars <- list()

##Extract node information##
##Record network structure##
network <- ""
# n3NetSrc <- c()
# n3NetTarget <- c()
for (i in seq(1, nodes, 1)) {
  
  cur <- root[[1]][[i]]
  curName <- xmlAttrs(cur)["name"]
  if (xmlAttrs(cur)["type"] == "Boolean") {
    params = "True;False"
  } else {
    params = paste(xmlAttrs(cur)["params"], "", sep = "")
  }
  node <- root[[i+1]]
  if (xmlSize(node) == 1) { #root node - no parents
    par = "NP"
    cpt <- strsplit(xmlValue(node[[1]]), "\n")[[1]]
    cpt <- cpt[2:length(cpt)] 
    network <- paste(network, "[", xmlAttrs(cur)["name"], "]", sep="")
  } else if (xmlSize(node) == 2) { #Single parent
    par <- paste(xmlAttrs(node[[1]])["var"], "", sep = "")
    if (xmlAttrs(node)["type"] == "GDT") { #Can't handle GDTs
      cpt <- "GDT"
    } else {
      cpt <- strsplit(xmlValue(node[[2]]), "\n")[[1]]
      cpt <- cpt[2:length(cpt)] 
    }
    network <- paste(network, "[", xmlAttrs(cur)["name"], "|", par, "]", sep = "")
  } else {
    #multiple parents
    par = ""
    network <- paste(network, "[", xmlAttrs(cur)["name"], "|", sep = "")
    for (j in seq(1, xmlSize(node)-1, 1)) {
      print(j)
      par <- paste(par, xmlAttrs(node[[j]])["var"], ";", sep = "")
      if (j < xmlSize(node)-1) {
        network <- paste(network, xmlAttrs(node[[j]])["var"], ":", sep = "")
      } else if (j == xmlSize(node)-1) {
        network <- paste(network, xmlAttrs(node[[j]])["var"], "]", sep = "")
      }
    }
    
    if (xmlAttrs(node)["type"] == "GDT") { #Can't handle GDTs
      cpt <- "GDT"
    } else {
      cpt <- strsplit(xmlValue(node[[j+1]]), "\n")[[1]]
      cpt <- cpt[2:length(cpt)] 
    }
  }
  nodeVars[[curName]] <- list("parent" = par, 
                              "name" = paste(xmlAttrs(cur)["name"], "", sep = ""), 
                              "params" = params, 
                              "cpt" = cpt)
}

##Print and create network##
print(network)
testNet <- model2network(network)

##Create CPT matrices##
cptMats <- list()
for (name in names(nodeVars)) {

  cpt <- nodeVars[[name]]$cpt
#   if (nodeVars[[name]]$cpt == "GDT") { #Can't display GDTs
#     next
#   }
  if ("GDT" %in% nodeVars[[name]]$cpt) { #Can't display GDTs
    ##Store structure to pass to bn.fit
    parParams <- strsplit(nodeVars[[nodeVars[[name]]$parent]]$params, ";")[[1]]
    parName <- nodeVars[[name]]$parent
    nodeParams <- c("Mean", "Variance")
    x <- runif(2)
    x <- x/sum(x)
    start <- c(rep(x, length(parParams)))
    dim(start) <- c(length(nodeParams),length(parParams))
    listNames <- list()
    listNames[[name]] <- nodeParams
    listNames[[parName]] <- parParams
    dimnames(start) <- listNames
    cptMats[[name]] <- start
    
    #Store dataframe to attempt independent plotting
    
    next
  }
  if (length(cpt) == 1) { #root node - CPT = prior
    cptVals <- str_replace_all(cpt, ";", "")
    cptVals <- as.numeric(strsplit(cptVals, ",")[[1]])
    params <- strsplit(nodeVars[[name]]$params, ";")[[1]]
    cptMat <- matrix(cptVals, ncol = length(cptVals), dimnames = list(NULL, params))
    cptMats[[name]] <- cptMat
  } 
  else if (! grepl(";", nodeVars[[name]]$parent)) { #single Parent
    parParams <- strsplit(nodeVars[[nodeVars[[name]]$parent]]$params, ";")[[1]]
    parName <- nodeVars[[name]]$parent
    nodeParams <- strsplit(nodeVars[[name]]$params, ";")[[1]]
    cptVals <- str_replace_all(cpt, "[0-9]+: ", "")
    cptVals <- str_replace_all(cptVals, "; (.*)", "")
    start <- c()
    for (i in seq(1, length(cptVals), 1)) {
      start <- c(start, as.numeric(strsplit(cptVals[i], ",")[[1]]))
    }
    dim(start) <- c(length(nodeParams),length(parParams))
    listNames <- list()
    listNames[[name]] <- nodeParams
    listNames[[parName]] <- parParams
    dimnames(start) <- listNames
    cptMats[[name]] <- start
  } else { #multi parent
    nodeParams <- strsplit(nodeVars[[name]]$params, ";")[[1]]
    parents <- strsplit(nodeVars[[name]]$parent, ";")[[1]]
    parParams <- list()
    for (k in seq(1, length(parents), 1)) {
      params <- strsplit(nodeVars[[parents[k]]]$params, ";")[[1]]
      parParams[[parents[k]]] <- params
    }
    states <- expand.grid(parParams)
    
    cptVals <- str_replace_all(cpt, "[0-9]+: ", "")

    start <- c()
    for (y in seq(1, dim(states)[1], 1)) {
      stateString <- "("
      for (x in seq(1, length(states), 1)) {
        if (x < length(states)) {
          stateString <- paste(stateString, states[y,x], ", ", sep = "")
        } else {
          stateString <- paste(stateString, states[y,x], ")", sep = "")
        }
      }
      for (z in seq(1, length(cptVals), 1)) {
        if (grepl(stateString, cptVals[z])) {
          curCPT <- str_replace_all(cptVals[z], "; (.*)", "")
          start <- c(start, as.numeric(strsplit(curCPT, ",")[[1]]))
          if (sum(as.numeric(strsplit(curCPT, ",")[[1]])) > 2) {
            print(stateString)
            print(as.numeric(strsplit(curCPT, ",")[[1]]))
          }
          break
        } else if (z == length(cptVals)) { #Tested all options
          x <- runif(length(nodeParams))
          x <- x/sum(x) #Random Distribution
          start <- c(start, x) #not included in CPT so create empty result
        }
      }
    }
    dimList <- list()
    dimList[[name]] <- nodeParams
    dimList <- append(dimList, parParams)
    dims <- c()
    for (val in dimList) {
      dims <- c(dims, length(val))
    }
    dim(start) <- dims
    dimnames(start) <- dimList

    cptMats[[name]] <- start
  }
}

##Load network with CPTs##
netFit = custom.fit(testNet, dist = cptMats)

bn.fit.barchart(netFit$Sequence)[2] #can access individual subplots
bn.fit.barchart(netFit$RepeatSeq)
pdf("testBarchartVL.pdf", width = 10, height = 10)
bn.fit.barchart(netFit$VarianceL)
# par(cex=3, cex.axis=3, cex.lab=3, cex.main=5, cex.sub=3)
dev.off()


##Data frame testing##
matSeq <- cptMats$Sequence
dfSeq <- melt(matSeq, id = "row.names")
ggplot(dfSeq, aes(x = Sequence, y = value, group = RepeatSeq, fill = RepeatSeq)) +
  geom_bar(stat = "identity") + coord_flip() + theme_bw()

ggplot(dfSeq, aes(x = Sequence, y = value, group = RepeatSeq, fill = RepeatSeq)) +
  geom_bar(stat = "identity") + coord_flip() + 
  facet_wrap(~RepeatSeq) + theme_bw()

ggplot(dfSeq, aes(x = Sequence, y = value, group = RepeatSeq, colour = RepeatSeq)) +
  geom_point(size = 5) + 
  facet_wrap(~RepeatSeq) + theme_bw()

ggplot(dfSeq, aes(x = Sequence, y = value, group = RepeatSeq, colour = RepeatSeq)) +
  geom_point(size = 5) + 
  theme_bw()

matVL <- cptMats$VarianceL
dfVL <- melt(matVL)

ggplot(dfVL, aes(x = factor(VarianceL), y = value, group = RepeatSeq, fill = RepeatSeq)) +
  geom_bar(stat = "identity") + coord_flip() + 
  facet_wrap( ~ Genomic.Location + Sequence) + theme_bw()

ggplot(dfVL, aes(x = factor(VarianceL), y = value, group = RepeatSeq, fill = Sequence)) +
  geom_bar(stat = "identity") + coord_flip() + 
  facet_wrap( ~ Genomic.Location + RepeatSeq) + theme_bw()

ggplot(dfVL, aes(x = factor(VarianceL), y = value)) +
  geom_bar(stat = "identity", fill = "deepskyblue") + coord_flip() + 
  facet_wrap( ~ Genomic.Location + RepeatSeq + Sequence) + theme_bw()



