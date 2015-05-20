#Preliminary TRY trait analysis
require(parallel)
require(geiger)
require(OUwie)
require(caper)
require(phytools)
source("~/Code/old_data/lib/numeric_TRY_data_cleaning.R")

drip.node.labels <- function(tree){
    #Assertions and argument checking
    if(is.null(tree$node.label)) stop("Phylogeny must have internally-labelled nodes")
    if(!inherits(tree, "phylo")) stop("Phylogeny must be an ape::phylo object")
    
    #Order node labels according to age
    named.nodes <- which(tree$node.label != "")
    names(named.nodes) <- tree$node.label[named.nodes]
    ages.nodes <- sapply(named.nodes, function(x) nodeheight(tree, x))
    named.nodes <- named.nodes[order(ages.nodes)]

    #Loop through nodes and label down the tree (slow, I know...)
    node.labels <- character(tree$Nnode)
    for(i in seq_along(named.nodes)){
        next.edges <- numeric(0)
        curr.edges <- tree$edge[tree$edge[,1]==named.nodes[i]+length(tree$tip.label),2]
        curr.edges <- curr.edges[curr.edges > length(tree$tip.label)]
        while(length(curr.edges) > 1){
            for(j in seq_along(curr.edges)){
                node.labels[curr.edges[j]] <- names(named.nodes[i])
                t <- tree$edge[tree$edge[,1]==curr.edges[j],2]
                next.edges <- append(next.edges, t[t > length(tree$tip.label)])
            }
            curr.edges <- next.edges
            next.edges <- numeric(0)
        }
    }

    return(node.labels[-seq_along(tree$tip.label)])
}

#Trim TRY (will get more later once taxon merge is complete)
data <- read.table("~/Desktop/clean-ish-try.txt")
species.data <- rownames(data) <- data$species
data$species <- NULL
data <- numericize.df(data, 0.25)
data <- data[,apply(data, 2, function(x) sum(!is.na(x))/nrow(data)) > 0.05]

#Get numerical traits (crudely) and transform
n.data <- data
n.data[,c(2,4,5,6,7,8,9,10,11,13,15,19,20,21,22,25,27,28)] <- log10(n.data[,c(2,4,5,6,7,8,9,10,11,13,15,19,20,21,22,25,27,28)])
n.data <- n.data[,-c(1,14)]
rownames(n.data) <- species.data

#Load phylogeny and subset
tree <- read.tree("~/Dropbox/Common/Vascular_Plants_rooted.dated.tre")
tree$tip.label <- tolower(gsub("_", " ", tree$tip.label))
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(n.data)))
n.data <- n.data[rownames(n.data) %in% tree$tip.label,]

#Paint clades with useful names
tree$node.label <- ifelse(tree$node.label %in% c("Fabidae", "Fabales", "Magnoliidae", "Rosaceae", "Monocotyledoneae", "Superrosidae", "Superasteridae", "angio", "Angiospermae", "Poaceae"), tree$node.label, "")
tree$node.label <- drip.node.labels(tree)
tree$nodel.label[is.na(tree$node.label)] <- ""
n.data <- n.data[match(tree$tip.label, rownames(n.data)),]

#Calculate Brownian motion
regimes <- rep(NA, nrow(tree$edge))
for(i in seq(2, nrow(tree$edge)))
    if(any(tree$edge[,2] == tree$edge[i,1]))
        regimes[i] <- which(tree$edge[,2] == tree$edge[i,1])

regimes <- tree$node.label[regimes[tree$edge[,2] <= length(tree$tip.label)]]
regimes[is.na(regimes)] <- "background"

wrap.data <- wrap.data <- function(trait.index){
    t <- data.frame(matrix(c(rownames(n.data), regimes, n.data[,2], rep(NA, nrow(n.data))), ncol=4), stringsAsFactors=FALSE)
    t[,3] <- as.numeric(t[,3])
    t <- t[!is.na(t[,3]),]
    t[,4] <- rep(NA, nrow(t))
    t.tree <- drop.tip(tree, setdiff(tree$tip.label, t[,1]))
    saveRDS(OUwie(t.tree, t, "BMS"), paste0(trait.index, "_bms", ".RDS"))
}

models <- mcmapply(wrap.data, trait.index=rep(seq_len(ncol(n.data))), mc.cores=24)
