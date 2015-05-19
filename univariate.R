#Preliminary TRY trait analysis
require(parallel)
require(geiger)
require(mvMORPH)

#Trim TRY (will get more later once taxon merge is complete)
data <- read.delim("~/Code/try_clean/try_trait.txt")
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
tree <- read.tree("~/Dropbox/Common/Vascular_Plants_rooted.dated.tre")
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

models <- mcmapply(wrap.data, trait.index=rep(seq(2,ncol(n.data)),each=8), model=rep(c("BM","OU","EB","trend","lambda","kappa","delta","white"),ncol(n.data)-1), mc.cores=24)

#Load the RDS files
setwd("rds")

files <- list.files()
files <- files[grepl(".RDS", files)]
aics <- data.frame(aic=numeric((ncol(n.data)-1)*8), model=character((ncol(n.data)-1)*8), trait=numeric((ncol(n.data)-1)*8), stringsAsFactors=FALSE)
for(i in seq_along(files)){
    t <- readRDS(files[i])
    name <- strsplit(files[i], "_|\\.")[[1]]
    aics[i,] <- c(t$opt$aic, name[2], as.numeric(name[1]))
}

aics$aic <- as.numeric(aics$aic)
aics <- aics[aics$aic < 20000,]
aics <- aics[aics$model != "",]
min.aic <- tapply(aics$aic, aics$trait, min)
aics$diff <- aics$aic - min.aic[aics$trait]
with(aics, boxplot(diff ~ model))

png("aic_boxplot.png", width=1000, height=1000)
t.aics <- aics[!aics$model %in% c("OU", "trend"),]
ranks <- unlist(with(t.aics, tapply(setNames(aic,model), trait, rank)))
ranks.model <- unname(sapply(names(ranks), function(x) gsub("[0-9]*\\.", "", x)))
ranks.model <- ranks.model[!ranks.model %in% c("OU", "trend")]
boxplot(ranks ~ factor(ranks.model), cex.axis=2)
dev.off()

lam.files <- files[grepl("lambda", files)]
lam <- setNames(rep(NA, length(lam.files)), lam.files)
for(i in seq_along(lam.files)){
    t <- readRDS(lam.files[i])
    lam[i] <- t$opt$lambda
}
names(lam) <- sapply(names(lam), function(x) strsplit(x, "_")[[1]][1])

png("lam_hist.png", width=1000, height=1000)
hist(lam, breaks=20, col="grey", xlab=expression(lambda), cex.lab=3, ylab="", main="")
dev.off()
