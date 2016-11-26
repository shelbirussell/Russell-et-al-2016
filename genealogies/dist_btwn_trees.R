library(ape)
library(phangorn)

symtree <- read.tree("RAxML_bipartitions.SvsymGenomeOG_mauve_raxmltree-renamed")
symrooted <- root(symtree, c("Svelesiana", "Spervernic", "Selarraich"))
symcollapse <- pruneTree(symrooted, 50)

mitotree <- read.tree("RAxML_bipartitions.SvmitoGenomeOG_mauve_raxmltree-renamed")
mitorooted <- root(mitotree, c("Svelesiana", "Spervernic", "Selarraich"))
mitocollapse <- pruneTree(mitorooted, 50)


#The topological distance is defined as twice the number of internal branches defining different bipartitions of the tips (Penny and Hendy 1985). Rzhetsky and Nei (1992) proposed a modification of the original formula to take multifurcations into account.
treedist(symcollapse, mitocollapse)
RF.dist(symcollapse, mitocollapse, rooted=TRUE)
dist.topo(symcollapse, mitocollapse)

