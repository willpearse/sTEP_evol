#Preliminary TRY trait analysis
require(parallel)
require(geiger)
require(mvMORPH)

#Trim TRY (will get more later once taxon merge is complete)
data <- read.table("clean-ish-try.txt", header=TRUE)
data[data==-9999] <- NA
t <- sapply(data, function(x) sum(!is.na(x)))
data <- data[,t > 1000]
data$species <- as.character(data$species)

#Get numerical traits (crudely) and transform
n.data <- data
for(i in seq(2,ncol(n.data)))
    try(n.data[,i] <- as.numeric(as.character(n.data[,i])))
n.data <- n.data[, sapply(n.data, function(x) sum(!is.na(x))) > 1000]
n.data[,c(2:3,5:15,18:21,24,26,27)] <- log10(n.data[,c(2:3,5:15,18:21,24,26,27)])
n.data <- n.data[,-2]

#Load phylogeny and subset
tree <- read.tree("Vascular_Plants_rooted.dated.tre")
tree$tip.label <- tolower(gsub("_", " ", tree$tip.label))
tree$node.label <- NULL
tree <- drop.tip(tree, setdiff(tree$tip.label, n.data$species))
n.data <- n.data[n.data$species %in% tree$tip.label,]

#Calcualte phylogenetic signal!
wrap.data <- function(trait.index, model){
    trait <- setNames(n.data[,trait.index], n.data$species)
    trait <- trait[!is.na(trait)]
    tree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
    saveRDS(fitContinuous(tree, trait, model=model), paste0(trait.index, "_", model, ".RDS"))
}

#Multivariate
t.data <- data[,c(3,54)] #Picked by coverage...
rownames(t.data) <- data$species
for(i in seq_len(ncol(t.data)))
    t.data[,i] <- as.numeric(as.character(t.data[,i]))
t.data <- na.omit(t.data)
t.tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(t.data)))
model <- mvBM(t.tree, t.data)
