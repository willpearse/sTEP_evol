library(ape)

#' Trims young tips
#'
#' @param phy phylo objects
#' @param min_age_exclusive Minimum species age in million years.
#' @return a list that includes a trimmed phylo objects and the names of the trimmed tips
#' @example
#' \dontrun{
#'
#' # Load TRY data and Zanne's tree
#' dat <- read.delim("~/Desktop/trydat.txt", sep = " ")
#' tre <- ape::read.tree("~/Desktop/Vascular_Plants_rooted.dated.tre")
#'
#' # Make species names compatible between datasets
#' tre$tip.label  <- gsub("_", " ", tolower(tre$tip.label))
#'
#' # Subset datasets
#' tre <- ape::drop.tip(tre, setdiff(tre$tip.label, dat$species))
#' dat <- dat[dat$species %in% tre$tip.label, ]
#'
#' # Trim tips older than 50k years (0.05MY)
#' tre_new <- trim_young_tips(tre, min_age_exclusive = 0.05)
#' # If you want to use the trimmed phylogeny, it lives in
#' tre_new$phy
#' # But if you just want the names of the excluded tips, just do
#' tre_new$tips_exluded
#' }
trim_young_tips = function(phy, min_age_exclusive) {
    if(!ape::is.ultrametric(phy))
        stop("Tree must be ultrametric")

    trash_names <- NULL ## Yeah, you know what this fuckhead is getting ready for...
    repeat {
        tip   <- phy$edge[ , 2] <= ape::Ntip(phy)
        young <- phy$edge.length < min_age_exclusive

        if(any(tip & young)){
            young_edges <- phy$edge[tip & young, ]
            dupes       <- duplicated(young_edges[ , 1])
            trash       <- young_edges[dupes, 2]
            trash_names <- c(trash_names, phy$tip.label[trash]) ## grow a vector inside a loop!
            phy         <- ape::drop.tip(phy, trash)
        } else {
            break
        }
    }
    r <- list("phy" = phy, "tips_exluded" = trash_names)
    return(r)
}
