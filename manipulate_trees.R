
tree\3 <- read.tree("\2\3\4", tree.names="\3")
rooted\3 <- root(tree\3, c("Svelesiana", "Spervernic"))
collapse\3 <- pruneTree(rooted\3, 50)
trees[[\1]] <- collapse\3
rmv\3 <- drop.tip(collapse\3, root.edge=1, c("Svelesiana", "Spervernic", "Selarraich"))
rmOGtrees[[\1]] <- rmv\3\n

trees <- vector("list")
class(trees) <- "multiPhylo"
rmbrnchtrees <- vector("list")
class(rmbrnchtrees) <- "multiPhylo"

treescaffold1_1_100000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_1_100000_OG_mauve_raxmltree", tree.names="scaffold1_1_100000")
rootedscaffold1_1_100000 <- root(treescaffold1_1_100000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[1]] <- rootedscaffold1_1_100000
rmbrnchscaffold1_1_100000 <- rootedscaffold1_1_100000
rmbrnchscaffold1_1_100000$edge.length <- NULL
rmbrnchtrees[[1]] <- rmbrnchscaffold1_1_100000

treescaffold1_100001_200000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_100001_200000_OG_mauve_raxmltree", tree.names="scaffold1_100001_200000")
rootedscaffold1_100001_200000 <- root(treescaffold1_100001_200000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[2]] <- rootedscaffold1_100001_200000
rmbrnchscaffold1_100001_200000 <- rootedscaffold1_100001_200000
rmbrnchscaffold1_100001_200000$edge.length <- NULL
rmbrnchtrees[[2]] <- rmbrnchscaffold1_100001_200000

treescaffold1_200001_300000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_200001_300000_OG_mauve_raxmltree", tree.names="scaffold1_200001_300000")
rootedscaffold1_200001_300000 <- root(treescaffold1_200001_300000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[3]] <- rootedscaffold1_200001_300000
rmbrnchscaffold1_200001_300000 <- rootedscaffold1_200001_300000
rmbrnchscaffold1_200001_300000$edge.length <- NULL
rmbrnchtrees[[3]] <- rmbrnchscaffold1_200001_300000

treescaffold1_300001_400000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_300001_400000_OG_mauve_raxmltree", tree.names="scaffold1_300001_400000")
rootedscaffold1_300001_400000 <- root(treescaffold1_300001_400000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[4]] <- rootedscaffold1_300001_400000
rmbrnchscaffold1_300001_400000 <- rootedscaffold1_300001_400000
rmbrnchscaffold1_300001_400000$edge.length <- NULL
rmbrnchtrees[[4]] <- rmbrnchscaffold1_300001_400000

treescaffold1_400001_500000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_400001_500000_OG_mauve_raxmltree", tree.names="scaffold1_400001_500000")
rootedscaffold1_400001_500000 <- root(treescaffold1_400001_500000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[5]] <- rootedscaffold1_400001_500000
rmbrnchscaffold1_400001_500000 <- rootedscaffold1_400001_500000
rmbrnchscaffold1_400001_500000$edge.length <- NULL
rmbrnchtrees[[5]] <- rmbrnchscaffold1_400001_500000

treescaffold1_500001_600000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_500001_600000_OG_mauve_raxmltree", tree.names="scaffold1_500001_600000")
rootedscaffold1_500001_600000 <- root(treescaffold1_500001_600000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[6]] <- rootedscaffold1_500001_600000
rmbrnchscaffold1_500001_600000 <- rootedscaffold1_500001_600000
rmbrnchscaffold1_500001_600000$edge.length <- NULL
rmbrnchtrees[[6]] <- rmbrnchscaffold1_500001_600000

treescaffold1_600001_700000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_600001_700000_OG_mauve_raxmltree", tree.names="scaffold1_600001_700000")
rootedscaffold1_600001_700000 <- root(treescaffold1_600001_700000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[7]] <- rootedscaffold1_600001_700000
rmbrnchscaffold1_600001_700000 <- rootedscaffold1_600001_700000
rmbrnchscaffold1_600001_700000$edge.length <- NULL
rmbrnchtrees[[7]] <- rmbrnchscaffold1_600001_700000

treescaffold1_700001_800000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_700001_800000_OG_mauve_raxmltree", tree.names="scaffold1_700001_800000")
rootedscaffold1_700001_800000 <- root(treescaffold1_700001_800000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[8]] <- rootedscaffold1_700001_800000
rmbrnchscaffold1_700001_800000 <- rootedscaffold1_700001_800000
rmbrnchscaffold1_700001_800000$edge.length <- NULL
rmbrnchtrees[[8]] <- rmbrnchscaffold1_700001_800000

treescaffold1_800001_900000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_800001_900000_OG_mauve_raxmltree", tree.names="scaffold1_800001_900000")
rootedscaffold1_800001_900000 <- root(treescaffold1_800001_900000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[9]] <- rootedscaffold1_800001_900000
rmbrnchscaffold1_800001_900000 <- rootedscaffold1_800001_900000
rmbrnchscaffold1_800001_900000$edge.length <- NULL
rmbrnchtrees[[9]] <- rmbrnchscaffold1_800001_900000

treescaffold1_900001_1000000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_900001_1000000_OG_mauve_raxmltree", tree.names="scaffold1_900001_1000000")
rootedscaffold1_900001_1000000 <- root(treescaffold1_900001_1000000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[10]] <- rootedscaffold1_900001_1000000
rmbrnchscaffold1_900001_1000000 <- rootedscaffold1_900001_1000000
rmbrnchscaffold1_900001_1000000$edge.length <- NULL
rmbrnchtrees[[10]] <- rmbrnchscaffold1_900001_1000000

treescaffold1_1000001_1100000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_1000001_1100000_OG_mauve_raxmltree", tree.names="scaffold1_1000001_1100000")
rootedscaffold1_1000001_1100000 <- root(treescaffold1_1000001_1100000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[11]] <- rootedscaffold1_1000001_1100000
rmbrnchscaffold1_1000001_1100000 <- rootedscaffold1_1000001_1100000
rmbrnchscaffold1_1000001_1100000$edge.length <- NULL
rmbrnchtrees[[11]] <- rmbrnchscaffold1_1000001_1100000

treescaffold1_1100001_1200000 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_1100001_1200000_OG_mauve_raxmltree", tree.names="scaffold1_1100001_1200000")
rootedscaffold1_1100001_1200000 <- root(treescaffold1_1100001_1200000, "Spervernic", resolve.root=TRUE)
trees[[12]] <- rootedscaffold1_1100001_1200000
rmbrnchscaffold1_1100001_1200000 <- rootedscaffold1_1100001_1200000
rmbrnchscaffold1_1100001_1200000$edge.length <- NULL
rmbrnchtrees[[12]] <- rmbrnchscaffold1_1100001_1200000

treescaffold1_1200001_1213831 <- read.tree("RAxML_bipartitions.Svsym_scaffold1_1200001_1213831_OG_mauve_raxmltree", tree.names="scaffold1_1200001_1213831")
rootedscaffold1_1200001_1213831 <- root(treescaffold1_1200001_1213831, c("Svelesiana", "Selarraichensis"), resolve.root=TRUE)
trees[[13]] <- rootedscaffold1_1200001_1213831
rmbrnchscaffold1_1200001_1213831 <- rootedscaffold1_1200001_1213831
rmbrnchscaffold1_1200001_1213831$edge.length <- NULL
rmbrnchtrees[[13]] <- rmbrnchscaffold1_1200001_1213831

treescaffold2_1_100000 <- read.tree("RAxML_bipartitions.Svsym_scaffold2_1_100000_OG_mauve_raxmltree", tree.names="scaffold2_1_100000")
rootedscaffold2_1_100000 <- root(treescaffold2_1_100000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[14]] <- rootedscaffold2_1_100000
rmbrnchscaffold2_1_100000 <- rootedscaffold2_1_100000
rmbrnchscaffold2_1_100000$edge.length <- NULL
rmbrnchtrees[[14]] <- rmbrnchscaffold2_1_100000

treescaffold2_100001_200000 <- read.tree("RAxML_bipartitions.Svsym_scaffold2_100001_200000_OG_mauve_raxmltree", tree.names="scaffold2_100001_200000")
rootedscaffold2_100001_200000 <- root(treescaffold2_100001_200000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[15]] <- rootedscaffold2_100001_200000
rmbrnchscaffold2_100001_200000 <- rootedscaffold2_100001_200000
rmbrnchscaffold2_100001_200000$edge.length <- NULL
rmbrnchtrees[[15]] <- rmbrnchscaffold2_100001_200000

treescaffold2_200001_300000 <- read.tree("RAxML_bipartitions.Svsym_scaffold2_200001_300000_OG_mauve_raxmltree", tree.names="scaffold2_200001_300000")
rootedscaffold2_200001_300000 <- root(treescaffold2_200001_300000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[16]] <- rootedscaffold2_200001_300000
rmbrnchscaffold2_200001_300000 <- rootedscaffold2_200001_300000
rmbrnchscaffold2_200001_300000$edge.length <- NULL
rmbrnchtrees[[16]] <- rmbrnchscaffold2_200001_300000

treescaffold2_300001_400000 <- read.tree("RAxML_bipartitions.Svsym_scaffold2_300001_400000_OG_mauve_raxmltree", tree.names="scaffold2_300001_400000")
rootedscaffold2_300001_400000 <- root(treescaffold2_300001_400000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[17]] <- rootedscaffold2_300001_400000
rmbrnchscaffold2_300001_400000 <- rootedscaffold2_300001_400000
rmbrnchscaffold2_300001_400000$edge.length <- NULL
rmbrnchtrees[[17]] <- rmbrnchscaffold2_300001_400000

treescaffold2_400001_500000 <- read.tree("RAxML_bipartitions.Svsym_scaffold2_400001_500000_OG_mauve_raxmltree", tree.names="scaffold2_400001_500000")
rootedscaffold2_400001_500000 <- root(treescaffold2_400001_500000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[18]] <- rootedscaffold2_400001_500000
rmbrnchscaffold2_400001_500000 <- rootedscaffold2_400001_500000
rmbrnchscaffold2_400001_500000$edge.length <- NULL
rmbrnchtrees[[18]] <- rmbrnchscaffold2_400001_500000

treescaffold2_500001_600000 <- read.tree("RAxML_bipartitions.Svsym_scaffold2_500001_600000_OG_mauve_raxmltree", tree.names="scaffold2_500001_600000")
rootedscaffold2_500001_600000 <- root(treescaffold2_500001_600000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[19]] <- rootedscaffold2_500001_600000
rmbrnchscaffold2_500001_600000 <- rootedscaffold2_500001_600000
rmbrnchscaffold2_500001_600000$edge.length <- NULL
rmbrnchtrees[[19]] <- rmbrnchscaffold2_500001_600000

treescaffold2_600001_700000 <- read.tree("RAxML_bipartitions.Svsym_scaffold2_600001_700000_OG_mauve_raxmltree", tree.names="scaffold2_600001_700000")
rootedscaffold2_600001_700000 <- root(treescaffold2_600001_700000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[20]] <- rootedscaffold2_600001_700000
rmbrnchscaffold2_600001_700000 <- rootedscaffold2_600001_700000
rmbrnchscaffold2_600001_700000$edge.length <- NULL
rmbrnchtrees[[20]] <- rmbrnchscaffold2_600001_700000

treescaffold2_700001_800000 <- read.tree("RAxML_bipartitions.Svsym_scaffold2_700001_800000_OG_mauve_raxmltree", tree.names="scaffold2_700001_800000")
rootedscaffold2_700001_800000 <- root(treescaffold2_700001_800000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[21]] <- rootedscaffold2_700001_800000
rmbrnchscaffold2_700001_800000 <- rootedscaffold2_700001_800000
rmbrnchscaffold2_700001_800000$edge.length <- NULL
rmbrnchtrees[[21]] <- rmbrnchscaffold2_700001_800000

treescaffold2_800001_892555 <- read.tree("RAxML_bipartitions.Svsym_scaffold2_800001_892555_OG_mauve_raxmltree", tree.names="scaffold2_800001_892555")
rootedscaffold2_800001_892555 <- root(treescaffold2_800001_892555, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[22]] <- rootedscaffold2_800001_892555
rmbrnchscaffold2_800001_892555 <- rootedscaffold2_800001_892555
rmbrnchscaffold2_800001_892555$edge.length <- NULL
rmbrnchtrees[[22]] <- rmbrnchscaffold2_800001_892555

treescaffold3_1_100000 <- read.tree("RAxML_bipartitions.Svsym_scaffold3_1_100000_OG_mauve_raxmltree", tree.names="scaffold3_1_100000")
rootedscaffold3_1_100000 <- root(treescaffold3_1_100000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[23]] <- rootedscaffold3_1_100000
rmbrnchscaffold3_1_100000 <- rootedscaffold3_1_100000
rmbrnchscaffold3_1_100000$edge.length <- NULL
rmbrnchtrees[[23]] <- rmbrnchscaffold3_1_100000

treescaffold3_100001_200000 <- read.tree("RAxML_bipartitions.Svsym_scaffold3_100001_200000_OG_mauve_raxmltree", tree.names="scaffold3_100001_200000")
rootedscaffold3_100001_200000 <- root(treescaffold3_100001_200000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[24]] <- rootedscaffold3_100001_200000
rmbrnchscaffold3_100001_200000 <- rootedscaffold3_100001_200000
rmbrnchscaffold3_100001_200000$edge.length <- NULL
rmbrnchtrees[[24]] <- rmbrnchscaffold3_100001_200000

treescaffold3_200001_300000 <- read.tree("RAxML_bipartitions.Svsym_scaffold3_200001_300000_OG_mauve_raxmltree", tree.names="scaffold3_200001_300000")
rootedscaffold3_200001_300000 <- root(treescaffold3_200001_300000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[25]] <- rootedscaffold3_200001_300000
rmbrnchscaffold3_200001_300000 <- rootedscaffold3_200001_300000
rmbrnchscaffold3_200001_300000$edge.length <- NULL
rmbrnchtrees[[25]] <- rmbrnchscaffold3_200001_300000

treescaffold3_300001_400000 <- read.tree("RAxML_bipartitions.Svsym_scaffold3_300001_400000_OG_mauve_raxmltree", tree.names="scaffold3_300001_400000")
rootedscaffold3_300001_400000 <- root(treescaffold3_300001_400000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[26]] <- rootedscaffold3_300001_400000
rmbrnchscaffold3_300001_400000 <- rootedscaffold3_300001_400000
rmbrnchscaffold3_300001_400000$edge.length <- NULL
rmbrnchtrees[[26]] <- rmbrnchscaffold3_300001_400000

treescaffold3_400001_500000 <- read.tree("RAxML_bipartitions.Svsym_scaffold3_400001_500000_OG_mauve_raxmltree", tree.names="scaffold3_400001_500000")
rootedscaffold3_400001_500000 <- root(treescaffold3_400001_500000, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[27]] <- rootedscaffold3_400001_500000
rmbrnchscaffold3_400001_500000 <- rootedscaffold3_400001_500000
rmbrnchscaffold3_400001_500000$edge.length <- NULL
rmbrnchtrees[[27]] <- rmbrnchscaffold3_400001_500000

treescaffold3_500001_537613 <- read.tree("RAxML_bipartitions.Svsym_scaffold3_500001_537613_OG_mauve_raxmltree", tree.names="scaffold3_500001_537613")
rootedscaffold3_500001_537613 <- root(treescaffold3_500001_537613, c("Svelesiana", "Spervernic"), resolve.root=TRUE)
trees[[28]] <- rootedscaffold3_500001_537613
rmbrnchscaffold3_500001_537613 <- rootedscaffold3_500001_537613
rmbrnchscaffold3_500001_537613$edge.length <- NULL
rmbrnchtrees[[28]] <- rmbrnchscaffold3_500001_537613

treescaffold4_1_28016 <- read.tree("RAxML_bipartitions.Svsym_scaffold4_1_28016_OG_mauve_raxmltree", tree.names="scaffold4_1_28016")
rootedscaffold4_1_28016 <- root(treescaffold4_1_28016, "Svelesiana", resolve.root=TRUE)
trees[[29]] <- rootedscaffold4_1_28016
rmbrnchscaffold4_1_28016 <- rootedscaffold4_1_28016
rmbrnchscaffold4_1_28016$edge.length <- NULL
rmbrnchtrees[[29]] <- rmbrnchscaffold4_1_28016
