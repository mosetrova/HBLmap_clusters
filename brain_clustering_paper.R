library(ggplot2)
library(ape)
library(ComplexHeatmap)
library(tidyr)
library(reshape2)
library(dplyr)
library(Rtsne)
library(gridExtra)
library(ggrepel)
library(cowplot)

args = c("MS-75reg", "AB_RNA-98reg")
labels = c("mRNA, 35 regions (our data)", "mRNA, 98 regions (Allen Institute)")
tsne.plots = list()

trans = read.delim("ab_to_ours.txt", check.names = FALSE)
assoc = melt(as.matrix(trans))
assoc = assoc[!is.na(assoc$value),]

assoc = data.frame(ms = assoc$Var2, ab = assoc$Var1, link = assoc$value)

anat.colors = read.delim("anatomical_colors.txt", comment.char = "", header = FALSE, stringsAsFactors = FALSE)
named.colors = anat.colors$V2
names(named.colors) = anat.colors$V1

clust.colors = read.table("cluster_colors.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE, comment.char = "")$V1

regions = read.delim("new_75_regions.txt", header = TRUE, stringsAsFactors = FALSE)

data = as.data.frame(read.table(paste0(args[1], "-data.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE))
peaks = as.data.frame(read.table(paste0(args[1], "-peaks.txt"), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE))

peaks = peaks[,data$V8 != ""]
data = data[data$V8 != "",]

norm = peaks

for (b in unique(data$V6)){
  avg.b = c()
  for (i in unique(data$V8)){
    m = norm[,(data$V8==i)&(data$V6==b)]
    s = sum((data$V8==i)&(data$V6==b))
    if(s>1){
      m = rowMeans(m)
    }
    if(s>0){
      avg.b = cbind(avg.b,m)
    }
  }
  coeff = rowMeans(avg.b)
  norm[,data$V6==b] = apply(norm[,data$V6==b],2,function (x) x-coeff)
}

res.man = manova(as.matrix(t(norm)) ~ data$V8)
pvals.r = as.matrix(as.data.frame(lapply(summary.aov(res.man),function (x) x[["Pr(>F)"]][[1]])))
length(pvals.r)
BH.r = p.adjust(pvals.r,method="BH")
sum(BH.r<0.00001)
sum(pvals.r<0.00001)

norm = norm[BH.r<0.00001,]

ms.peaks = peaks
ms.norm = norm
ms.data = data

reg = c()

for (r in unique(data$V8)) {
  if(sum(data$V8 == r) > 1) {
    reg = cbind(reg, rowMeans(norm[,data$V8 == r]))
  }
  else {
    reg = cbind(reg, norm[,data$V8 == r])
  }
}

colnames(reg) = unique(data$V8)

ms.reg = reg
ms.reginfo = data.frame(region = colnames(ms.reg), area = sapply(colnames(ms.reg), function(x) regions$region.area[regions$region == x][1]))

data = as.data.frame(read.table(paste0(args[2], "-data.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE))
peaks = as.data.frame(read.table(paste0(args[2], "-peaks.txt"), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE))

peaks = peaks[,data$V8 != ""]
data = data[data$V8 != "",]

norm = peaks

for (b in unique(data$V6)){
  avg.b = c()
  for (i in unique(data$V8)){
    m = norm[,(data$V8==i)&(data$V6==b)]
    s = sum((data$V8==i)&(data$V6==b))
    if(s>1){
      m = rowMeans(m)
    }
    if(s>0){
      avg.b = cbind(avg.b,m)
    }
  }
  coeff = rowMeans(avg.b)
  norm[,data$V6==b] = apply(norm[,data$V6==b],2,function (x) x-coeff)
}

res.man = manova(as.matrix(t(norm)) ~ data$V8)
pvals.r = as.matrix(as.data.frame(lapply(summary.aov(res.man),function (x) x[["Pr(>F)"]][[1]])))
length(pvals.r)
BH.r = p.adjust(pvals.r,method="BH")
sum(BH.r<0.00001)
sum(pvals.r<0.00001)

norm = norm[BH.r<0.00001,]

reg = c()

for (r in unique(data$V8)) {
  if(sum(data$V8 == r) > 1) {
    reg = cbind(reg, rowMeans(norm[,data$V8 == r]))
  }
  else {
    reg = cbind(reg, norm[,data$V8 == r])
  }
}

colnames(reg) = unique(data$V8)

ab.reg = reg
ab.reginfo = data.frame(region = colnames(ab.reg), area = sapply(colnames(ab.reg), function(x) data$V9[data$V8 == x][1]))

### hierarchical clustering

knum = 6

ms.hc = hclust(dist(t(ms.reg)))
ab.hc = hclust(dist(t(ab.reg)))

ms.clusters = data.frame(region = names(cutree(ms.hc, 4)), cluster = cutree(ms.hc, 4))
ms.clusters = left_join(ms.clusters, ms.reginfo, by = "region")

ab.clusters = data.frame(region = names(cutree(ab.hc, 6)), cluster = cutree(ab.hc, 6))
ab.clusters = left_join(ab.clusters, ab.reginfo, by = "region")

h1 = Heatmap(ms.clusters[,c("cluster")],
             rect_gp = gpar(col = "white", lwd = 2),
             col = clust.colors[3:length(clust.colors)],
             name = "cluster",
             cluster_rows = as.dendrogram(ms.hc),
             row_dend_width = unit(10, "cm"))
h2 = Heatmap(ms.clusters[,c("area")],
             rect_gp = gpar(col = "white", lwd = 2),
             col = named.colors,
             name = "area")

ht_list = h1 + h2
draw(ht_list, ht_gap = unit(0, "mm"))

h1 = Heatmap(ab.clusters[,c("cluster")],
             rect_gp = gpar(col = "white", lwd = 2),
             col = clust.colors[3:length(clust.colors)],
             name = "cluster",
             cluster_rows = as.dendrogram(ab.hc),
             row_dend_width = unit(10, "cm"))
h2 = Heatmap(ab.clusters[,c("area")],
             rect_gp = gpar(col = "white", lwd = 2), 
             col = named.colors,
             name = "area")

ht_list = h1 + h2
draw(ht_list, ht_gap = unit(0, "mm"))

### draw trees and tanglegrams

ms.tree = as.phylo(ms.hc)
ab.tree = as.phylo(ab.hc)

plot(ms.tree)
plot(ab.tree)

assoc = c()

for (m in colnames(ms.reg)[colnames(ms.reg) %in% colnames(trans)]) {
  for (a in colnames(ab.reg)[colnames(ab.reg) %in% rownames(trans)]) {
    if (!is.na(trans[a, m])) {
      assoc = rbind(assoc, c(m, a))
    }
  }
}

cophyloplot(ms.tree, ab.tree, assoc = assoc,
            #lwd = assoc[,3], 
            length.line = 0,
            #gap = 0,
            space = 150,
            rotate = FALSE,
            type = "phylogram",
            show.tip.label = FALSE,
            use.edge.length = TRUE)

tsne = Rtsne(t(norm), dims = 2, perplexity=10, verbose=TRUE, max_iter = 500)
tsne.df = as.data.frame(tsne$Y)
colnames(tsne.df) = c("x", "y")
tsne.df$brain = data$V6
tsne.df$region = data$V8
tsne.df$area = data$V9
tsne.df$cluster = sapply(tsne.df$region, function(x) ab.clusters$cluster[ab.clusters$region == x])

t1 = ggplot(tsne.df, aes(x, y, fill = area)) +
  geom_point(size = 3, shape = 21, show.legend = FALSE) +
  scale_fill_manual(name = "brain region", values = named.colors) +
  #ggtitle("PCA (brain norm)") +
  theme(aspect.ratio = 1,
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank()) +
  xlab("") +
  ylab("") 

t2 = ggplot(tsne.df, aes(x, y, fill = factor(cluster))) +
  geom_point(size = 3, shape = 21, show.legend = FALSE) +
  scale_fill_manual(name = "brain region", values = clust.colors) +
  #ggtitle("PCA (brain norm)") +
  theme(aspect.ratio = 1,
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank()) +
  xlab("") +
  ylab("")

grid.arrange(t1, t2, nrow = 2)

tsne = Rtsne(t(ms.norm), dims = 2, perplexity=10, verbose=TRUE, max_iter = 500)
tsne.df = as.data.frame(tsne$Y)
colnames(tsne.df) = c("x", "y")
tsne.df$brain = ms.data$V6
tsne.df$region = ms.data$V8
tsne.df$area = sapply(ms.data$V8, function(x) ms.reginfo$area[ms.reginfo$region == x])
tsne.df$cluster = sapply(tsne.df$region, function(x) ms.clusters$cluster[ms.clusters$region == x])

ms1 = ggplot(tsne.df, aes(x, y, fill = area)) +
  geom_point(size = 3, shape = 21, show.legend = FALSE) +
  scale_fill_manual(name = "brain region", values = named.colors) +
  #ggtitle("PCA (brain norm)") +
  theme(aspect.ratio = 1,
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank()) +
  xlab("") +
  ylab("") 

ms2 = ggplot(tsne.df, aes(x, y, fill = factor(cluster))) +
  geom_point(size = 3, shape = 21, show.legend = FALSE) +
  scale_fill_manual(name = "brain region", values = clust.colors) +
  #ggtitle("PCA (brain norm)") +
  theme(aspect.ratio = 1,
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank()) +
  xlab("") +
  ylab("")

grid.arrange(ms1, ms2, nrow = 2)


### test clustering stability for lipids

pMatrix.min <- function(A, B) { 
  #finds the permutation P of A such that ||PA - B|| is minimum in Frobenius norm 
  # Uses the linear-sum assignment problem (LSAP) solver in the "clue" package 
  
  # Returns P%*%A and the permutation vector `pvec' such that 
  # A[pvec, ] is the permutation of A closest to B 
  n <- nrow(A) 
  D <- matrix(NA, n, n) 
  for (i in 1:n) { 
    for (j in 1:n) { 
      D[j, i] <- (sum((B[j, ] - A[i, ])^2)) 
    }
  } 
  vec <- c(solve_LSAP(D)) 
  list(A=A[vec,], pvec=vec) 
} 

require(clue)  # need this package to solve the LSAP 

knum = 4

ms.clusters = data.frame(region = names(cutree(ms.hc, knum)), cluster = cutree(ms.hc, knum))
ms.clusters = arrange(ms.clusters, cluster)

ind.clusters = c()
for (ind in unique(ms.data$V6)) {
  temp.norm = ms.norm[, ms.data$V6 != ind]
  temp.data = ms.data[ms.data$V6 != ind,]
  temp.reg = c()
  for (r in unique(temp.data$V8)) {
    if(sum(temp.data$V8 == r) > 1) {
      temp.reg = cbind(temp.reg, rowMeans(temp.norm[,temp.data$V8 == r]))
    }
    else {
      temp.reg = cbind(temp.reg, temp.norm[,temp.data$V8 == r])
    }
  }
  colnames(temp.reg) = unique(temp.data$V8)
  temp.hc = hclust(dist(t(temp.reg)))
  ms.clusters = left_join(ms.clusters, data.frame(region = names(cutree(temp.hc, knum)), cluster = cutree(temp.hc, knum)), by = "region", suffix = c("", paste0("_no_", ind)))
  ind.clusters = cbind(ind.clusters, pMatrix.min(table(ms.clusters[,c("cluster", paste0("cluster_no_", ind))]), diag(1, knum))$pvec[ms.clusters[,paste0("cluster_no_", ind)]])
}

ms.clusters[,paste0("cluster_no_", unique(ms.data$V6))] = ind.clusters

ind.log = as.data.frame(apply(ind.clusters, 2, function(x) x == ms.clusters$cluster))
rownames(ind.log) = ms.clusters$region
ind.log$cluster = ms.clusters$cluster
ind.log$prob = (table(ind.log$cluster)/nrow(ms.clusters))[ind.log$cluster]
ind.log$pval = apply(ind.log, 1, function(x) sum(dbinom(sum(x[1:4]):4, 4, x[colnames(ind.log) == "prob"])))
ind.log$region = rownames(ind.log)
ind.log$cluster.name = paste0("cluster", ind.log$cluster)

ind.log$region = factor(ind.log$region, levels = rev(labels(as.dendrogram(ms.hc))))
ind.log$cluster.name = factor(ind.log$cluster.name, levels = rev(unique(ind.log$cluster.name[order(ind.log$region)])))

p2 = ggplot(ind.log, aes(region, pval)) +
  geom_point() +
  geom_hline(yintercept = 0.05, lty = 3, col = "black") +
  coord_flip() +
  xlab("") +
  ylab("binomial p-value") +
  theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
  #scale_y_log10() +
  facet_grid(cluster.name~., scales = "free_y", space = "free")

ind.log %>%
  group_by(cluster.name) %>%
  summarize(sign.prop = sum(pval < 0.05)/n())

p1 = melt(ms.clusters) %>%
  ggplot(aes(region, variable, fill = factor(value))) +
  geom_tile() +
  theme(panel.background = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip() +
  scale_fill_manual(name = "cluster", values = clust.colors[3:length(clust.colors)]) +
  scale_x_discrete(limits = rev(labels(as.dendrogram(ms.hc)))) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab("")

grid.arrange(p1, p2, nrow = 1)

ind.log$sign = ifelse(ind.log$pval < 0.05, "*", "")

clusters.melt = melt(ms.clusters)
clusters.melt = left_join(clusters.melt, ind.log[,c("region", "sign")])
clusters.melt$value = paste0("cluster", clusters.melt$value)
clusters.melt = rbind(clusters.melt, data.frame(region = regions$region,
                                                variable = 'anatomy',
                                                value = regions$region.area,
                                                #area = regions$region.area,
                                                sign = ""))
names(clust.colors) = paste0("cluster", 1:length(clust.colors))

ggplot(clusters.melt, aes(region, variable, fill = factor(value))) +
  scale_fill_manual(name = "cluster", values = c(named.colors, clust.colors)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_tile() +
  theme(panel.background = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        aspect.ratio = 75/16,
        axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_text(aes(x = region, y = "cluster", label = sign), col = "black") +
  coord_flip() +
  scale_x_discrete(limits = rev(labels(as.dendrogram(ms.hc))), position = "top") +
  guides(fill = "none") +
  xlab("") +
  ylab("")

### let's plot regions

back.reg = read.delim("background_coords_df.tsv")
reg.to.plot = read.delim("region_coords_df.tsv", stringsAsFactors = FALSE)

regions$num = as.numeric(gsub("^(\\d*)\\s.*", "\\1", regions$region))

reg.to.plot = left_join(reg.to.plot, regions[,c("num", "region", "region.area")])

back.reg$x = as.numeric(as.character(back.reg$x))
back.reg$y = as.numeric(as.character(back.reg$y))

reg.to.plot$x = as.numeric(as.character(reg.to.plot$x))
reg.to.plot$y = as.numeric(as.character(reg.to.plot$y))

b1 = ggplot(data = back.reg) + 
  geom_polygon(aes(x = x, y = y, group = group), fill = "lightgrey", color = "black", size = 0.5) + 
  geom_polygon(data = reg.to.plot[!is.na(reg.to.plot$region.area),], aes(x = x, y = y, fill = region.area, group = group), color = "black", size = 0.5) + 
  #geom_text(aes(x = x_poi, y = y_poi, label = col)) + # background numbers
  #geom_text(data = reg, aes(x = x_poi, y = y_poi, label = num)) + # region numbers
  #guides(fill=FALSE) +
  scale_fill_manual(name = "", values = named.colors, na.value = "transparent") +
  theme(aspect.ratio = (max(back.reg$y) - min(back.reg$y))/(max(back.reg$x) - min(back.reg$x)), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = NA), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

reg.to.plot$cluster = paste0("cluster", sapply(reg.to.plot$region, function(x) ms.clusters$cluster[ms.clusters$region == x]))

b2 = ggplot(data = back.reg) + 
  geom_polygon(aes(x = x, y = y, group = group), fill = "lightgrey", color = "black", size = 0.5) + 
  geom_polygon(data = reg.to.plot[!is.na(reg.to.plot$region.area),], aes(x = x, y = y, fill = cluster, group = group), color = "black", size = 0.5) + 
  #geom_text(aes(x = x_poi, y = y_poi, label = col)) + # background numbers
  #geom_text(data = reg, aes(x = x_poi, y = y_poi, label = num)) + # region numbers
  #guides(fill=FALSE) +
  scale_fill_manual(name = "", values = clust.colors, na.value = "transparent") +
  theme(aspect.ratio = (max(back.reg$y) - min(back.reg$y))/(max(back.reg$x) - min(back.reg$x)), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = NA), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

plot_grid(b1, b2, nrow = 2, align = "hv")

### cluster analysis for 35 regions

args = c("MS-75reg", "ourRNA-35reg")
labels = c("mRNA, 35 regions (our data)", "mRNA, 98 regions (Allen Institute)")
tsne.plots = list()

region.35 = unique(read.table(paste0(args[2], "-data.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V8)

data = as.data.frame(read.table(paste0(args[1], "-data.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE))
peaks = as.data.frame(read.table(paste0(args[1], "-peaks.txt"), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE))

peaks = peaks[,data$V8 != ""]
data = data[data$V8 != "",]

peaks = peaks[,data$V8 %in% region.35]
data = data[data$V8 %in% region.35,]

norm = peaks

for (b in unique(data$V6)){
  avg.b = c()
  for (i in unique(data$V8)){
    m = norm[,(data$V8==i)&(data$V6==b)]
    s = sum((data$V8==i)&(data$V6==b))
    if(s>1){
      m = rowMeans(m)
    }
    if(s>0){
      avg.b = cbind(avg.b,m)
    }
  }
  coeff = rowMeans(avg.b)
  norm[,data$V6==b] = apply(norm[,data$V6==b],2,function (x) x-coeff)
}

res.man = manova(as.matrix(t(norm)) ~ data$V8)
pvals.r = as.matrix(as.data.frame(lapply(summary.aov(res.man),function (x) x[["Pr(>F)"]][[1]])))
length(pvals.r)
BH.r = p.adjust(pvals.r,method="BH")
sum(BH.r<0.00001)
sum(pvals.r<0.00001)

norm = norm[BH.r<0.00001,]

ms.peaks = peaks
ms.norm = norm
ms.data = data

reg = c()

for (r in unique(data$V8)) {
  if(sum(data$V8 == r) > 1) {
    reg = cbind(reg, rowMeans(norm[,data$V8 == r]))
  }
  else {
    reg = cbind(reg, norm[,data$V8 == r])
  }
}

colnames(reg) = unique(data$V8)

ms.reg = reg
ms.reginfo = data.frame(region = colnames(ms.reg), area = sapply(colnames(ms.reg), function(x) regions$region.area[regions$region == x][1]))

data = as.data.frame(read.table(paste0(args[2], "-data.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE))
peaks = as.data.frame(read.table(paste0(args[2], "-peaks.txt"), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE))

peaks = peaks[,data$V8 != ""]
data = data[data$V8 != "",]

norm = peaks

for (b in unique(data$V6)){
  avg.b = c()
  for (i in unique(data$V8)){
    m = norm[,(data$V8==i)&(data$V6==b)]
    s = sum((data$V8==i)&(data$V6==b))
    if(s>1){
      m = rowMeans(m)
    }
    if(s>0){
      avg.b = cbind(avg.b,m)
    }
  }
  coeff = rowMeans(avg.b)
  norm[,data$V6==b] = apply(norm[,data$V6==b],2,function (x) x-coeff)
}

res.man = manova(as.matrix(t(norm)) ~ data$V8)
pvals.r = as.matrix(as.data.frame(lapply(summary.aov(res.man),function (x) x[["Pr(>F)"]][[1]])))
length(pvals.r)
BH.r = p.adjust(pvals.r,method="BH")
sum(BH.r<0.00001)
sum(pvals.r<0.00001)

norm = norm[BH.r<0.00001,]

reg = c()

for (r in unique(data$V8)) {
  if(sum(data$V8 == r) > 1) {
    reg = cbind(reg, rowMeans(norm[,data$V8 == r]))
  }
  else {
    reg = cbind(reg, norm[,data$V8 == r])
  }
}

colnames(reg) = unique(data$V8)

our.reg = reg
our.reginfo = data.frame(region = colnames(our.reg), area = sapply(colnames(our.reg), function(x) regions$region.area[regions$region == x][1]))

### hierarchical clustering

assoc = cbind(region.35, region.35)
knum = 4

ms.hc = hclust(dist(t(ms.reg)))
our.hc = hclust(dist(t(our.reg)))

ms.clusters = data.frame(region = names(cutree(ms.hc, 4)), cluster = cutree(ms.hc, 4))
ms.clusters = left_join(ms.clusters, ms.reginfo, by = "region")

our.clusters = data.frame(region = names(cutree(our.hc, 5)), cluster = cutree(our.hc, 5))
our.clusters = left_join(our.clusters, ms.reginfo, by = "region")

h1 = Heatmap(ms.clusters[,c("cluster")],
             #row_order = 1:nrow(ms.clusters),
             rect_gp = gpar(col = "white", lwd = 2),
             name = "cluster",
             cluster_rows = as.dendrogram(ms.hc),
             #row_split = knum,
             row_dend_width = unit(10, "cm"))
h2 = Heatmap(ms.clusters[,c("area")],
             #row_order = 1:nrow(ms.clusters),
             rect_gp = gpar(col = "white", lwd = 2),
             col = named.colors,
             name = "area")

ht_list = h1 + h2
draw(ht_list, ht_gap = unit(0, "mm"))

h1 = Heatmap(our.clusters[,c("cluster")],
             #row_order = 1:nrow(our.clusters),
             rect_gp = gpar(col = "white", lwd = 2),
             name = "cluster",
             cluster_rows = as.dendrogram(our.hc),
             #row_split = knum,
             row_dend_width = unit(10, "cm"))
h2 = Heatmap(our.clusters[,c("area")],
             #row_order = 1:nrow(our.clusters),
             rect_gp = gpar(col = "white", lwd = 2), 
             col = named.colors,
             name = "area")

ht_list = h1 + h2
draw(ht_list, ht_gap = unit(0, "mm"))

ms.tree = as.phylo(ms.hc)
our.tree = as.phylo(our.hc)

plot(ms.tree)
plot(our.tree)

cophyloplot(ms.tree, our.tree, assoc = assoc,
            #lwd = assoc[,3], 
            length.line = 0,
            #gap = 0,
            space = 150,
            rotate = FALSE,
            type = "phylogram",
            show.tip.label = FALSE,
            use.edge.length = TRUE)

tsne = Rtsne(t(norm), dims = 2, perplexity=10, verbose=TRUE, max_iter = 500)
tsne.df = as.data.frame(tsne$Y)
colnames(tsne.df) = c("x", "y")
tsne.df$brain = data$V6
tsne.df$region = data$V8
tsne.df$area = our.reginfo$area
tsne.df$cluster = paste0("cluster", sapply(tsne.df$region, function(x) our.clusters$cluster[our.clusters$region == x]))

t3 = ggplot(tsne.df, aes(x, y, fill = area)) +
  geom_point(size = 3, shape = 21, show.legend = FALSE) +
  scale_fill_manual(name = "brain region", values = named.colors) +
  #ggtitle("PCA (brain norm)") +
  theme(aspect.ratio = 1,
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank()) +
  xlab("") +
  ylab("")

t4 = ggplot(tsne.df, aes(x, y, fill = factor(cluster))) +
  geom_point(size = 3, shape = 21, show.legend = FALSE) +
  scale_fill_manual(name = "brain region", values = clust.colors) +
  #ggtitle("PCA (brain norm)") +
  theme(aspect.ratio = 1,
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank()) +
  xlab("") +
  ylab("")

grid.arrange(t3, t4, nrow = 2)

grid.arrange(t1, t3, t2, t4, nrow = 2)

### compare lipidome and transcriptome clusters to functional networks

networks = read.delim("networks.txt")
net.colors = read.delim("networks_colors.txt")

net.distances = c()
for (network in colnames(networks)) {
  temp.regions = rownames(networks)[networks[,network] == 1]
  temp.ms.dist = as.matrix(dist(t(ms.reg)))
  diag(temp.ms.dist) = NA
  
  temp.our.dist = as.matrix(dist(t(our.reg)))
  diag(temp.our.dist) = NA
  
  temp.regions = temp.regions[temp.regions %in% colnames(temp.ms.dist)]
  
  temp.distances = c(lip.in = median(temp.ms.dist[temp.regions, temp.regions], na.rm = TRUE),
                     lip.out = median(temp.ms.dist[temp.regions, !(colnames(temp.ms.dist) %in% temp.regions)], na.rm = TRUE),
                     expr.in = median(temp.our.dist[temp.regions, temp.regions], na.rm = TRUE),
                     expr.out = median(temp.our.dist[temp.regions, !(colnames(temp.ms.dist) %in% temp.regions)], na.rm = TRUE))
  
  net.distances = rbind(net.distances, temp.distances)
}

net.distances = as.data.frame(net.distances)
net.distances$network = gsub("\\.", " ", gsub(".network", "", colnames(networks)))
net.distances$lip.ratio = net.distances$lip.out/net.distances$lip.in
net.distances$expr.ratio = net.distances$expr.out/net.distances$expr.in

ggplot(net.distances, aes(expr.ratio, lip.ratio, label = network)) +
  geom_point(pch = 21, size = 6, fill = "lavender") +
  geom_label_repel(fill = "lightgrey") +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(limits = c(0.5, 2.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 2.5), expand = c(0, 0)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab("Transcriptome") +
  ylab("Lipidome") +
  ggtitle("Ratio between median distances\noutside and inside of functional networks")

wilcox.test(net.distances$lip.ratio, net.distances$expr.ratio, alternative = "greater")

net.plots = list()
for (network in colnames(networks)) {
  temp.reg.to.plot = reg.to.plot
  temp.networks = data.frame(region = rownames(networks) , network = networks[,network])
  temp.reg.to.plot = left_join(temp.reg.to.plot, temp.networks)
  p = ggplot(data = back.reg) + 
    geom_polygon(aes(x = x, y = y, group = group), fill = NA, color = "black", size = 0.2) + 
    geom_polygon(data = temp.reg.to.plot[!is.na(temp.reg.to.plot$region.area),], aes(x = x, y = y, fill = factor(network), group = group), color = "black", size = 0.2) + 
    #geom_text(aes(x = x_poi, y = y_poi, label = col)) + # background numbers
    #geom_text(data = reg, aes(x = x_poi, y = y_poi, label = num)) + # region numbers
    guides(fill=FALSE) +
    scale_fill_manual(name = "", values = c(NA, as.character(net.colors$color[net.colors$network == gsub("\\.", " " , network)])), na.value = "transparent") +
    theme(aspect.ratio = (max(back.reg$y) - min(back.reg$y))/(max(back.reg$x) - min(back.reg$x)), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, color = NA), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    ggtitle(paste0(gsub("\\.", " " , network), " (", sum(temp.networks$network), " reg.)"))
  net.plots[[network]] = p
}

plot_grid(plotlist = net.plots, nrow = 3, align = "hv")

### analyze correlation to spatial distances and SC/FC

data = as.data.frame(read.table(paste0(args[1], "-data.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE))
peaks = as.data.frame(read.table(paste0(args[1], "-peaks.txt"), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE))

peaks = peaks[,data$V8 != ""]
data = data[data$V8 != "",]

peaks = peaks[,data$V8 %in% region.35]
data = data[data$V8 %in% region.35,]

norm = peaks

for (b in unique(data$V6)){
  avg.b = c()
  for (i in unique(data$V8)){
    m = norm[,(data$V8==i)&(data$V6==b)]
    s = sum((data$V8==i)&(data$V6==b))
    if(s>1){
      m = rowMeans(m)
    }
    if(s>0){
      avg.b = cbind(avg.b,m)
    }
  }
  coeff = rowMeans(avg.b)
  norm[,data$V6==b] = apply(norm[,data$V6==b],2,function (x) x-coeff)
}

ms.peaks = peaks
ms.norm = norm
ms.data = data

reg = c()

for (r in unique(data$V8)) {
  if(sum(data$V8 == r) > 1) {
    reg = cbind(reg, rowMeans(norm[,data$V8 == r]))
  }
  else {
    reg = cbind(reg, norm[,data$V8 == r])
  }
}

colnames(reg) = unique(data$V8)

ms.reg = reg
ms.reginfo = data.frame(region = colnames(ms.reg), area = sapply(colnames(ms.reg), function(x) regions$region.area[regions$region == x][1]))

data = as.data.frame(read.table(paste0(args[2], "-data.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE))
peaks = as.data.frame(read.table(paste0(args[2], "-peaks.txt"), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE))

peaks = peaks[,data$V8 != ""]
data = data[data$V8 != "",]

norm = peaks

for (b in unique(data$V6)){
  avg.b = c()
  for (i in unique(data$V8)){
    m = norm[,(data$V8==i)&(data$V6==b)]
    s = sum((data$V8==i)&(data$V6==b))
    if(s>1){
      m = rowMeans(m)
    }
    if(s>0){
      avg.b = cbind(avg.b,m)
    }
  }
  coeff = rowMeans(avg.b)
  norm[,data$V6==b] = apply(norm[,data$V6==b],2,function (x) x-coeff)
}

reg = c()

for (r in unique(data$V8)) {
  if(sum(data$V8 == r) > 1) {
    reg = cbind(reg, rowMeans(norm[,data$V8 == r]))
  }
  else {
    reg = cbind(reg, norm[,data$V8 == r])
  }
}

colnames(reg) = unique(data$V8)

our.reg = reg
our.reginfo = data.frame(region = colnames(our.reg), area = sapply(colnames(our.reg), function(x) regions$region.area[regions$region == x][1]))


ctx = read.table("56_regions_cortex.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
neoctx = read.table("49_regions_neocortex.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

ctx = colnames(our.reg)[colnames(our.reg) %in% ctx$V1]
neoctx = colnames(our.reg)[colnames(our.reg) %in% neoctx$V1]

mri.coords = read.table("75_regions_mri_coord.txt", sep = "\t", header = TRUE)

rownames(mri.coords) = mri.coords$name

lip.dist = dist(t(ms.reg[,neoctx]))
expr.dist = dist(t(our.reg[,neoctx]))
eucl.dist = dist(mri.coords[neoctx,3:5])

cor(lip.dist, eucl.dist)
cor(expr.dist, eucl.dist)
cor(lip.dist, expr.dist)

lip.dist.mtx = as.matrix(lip.dist)
diag(lip.dist.mtx) = NA
expr.dist.mtx = as.matrix(expr.dist)
diag(expr.dist.mtx) = NA
eucl.dist.mtx = as.matrix(eucl.dist)
diag(eucl.dist.mtx) = NA

dist.cor = data.frame(
  #lip.eucl = diag(cor(as.matrix(lip.dist), as.matrix(eucl.dist))),
  #expr.eucl = diag(cor(as.matrix(expr.dist), as.matrix(eucl.dist))),
  #lip.expr = diag(cor(as.matrix(lip.dist), as.matrix(expr.dist)))
  lip.eucl = diag(apply(lip.dist.mtx, 1, function(x) apply(eucl.dist.mtx, 1, function(y) cor(x, y, use = "complete.obs")))),
  expr.eucl = diag(apply(expr.dist.mtx, 1, function(x) apply(eucl.dist.mtx, 1, function(y) cor(x, y, use = "complete.obs")))),
  lip.expr = diag(apply(lip.dist.mtx, 1, function(x) apply(expr.dist.mtx, 1, function(y) cor(x, y, use = "complete.obs"))))
)

str.conn = as.matrix(read.table("str_conn_59reg_dist.txt", sep = "\t", check.names = FALSE))
func.conn = as.matrix(read.table("func_conn_59reg_dist.txt", sep = "\t", check.names = FALSE))

diag(str.conn) = NA
diag(func.conn) = NA

str.dist = max(str.conn[neoctx,neoctx], na.rm = TRUE) - str.conn[neoctx,neoctx]
func.dist = max(func.conn[neoctx,neoctx], na.rm = TRUE) - func.conn[neoctx,neoctx]

str.cor = data.frame(
  lip.str = diag(apply(lip.dist.mtx, 1, function(x) apply(str.dist, 1, function(y) cor(x, y, use = "complete.obs")))),
  expr.str = diag(apply(expr.dist.mtx, 1, function(x) apply(str.dist, 1, function(y) cor(x, y, use = "complete.obs"))))
)

func.cor = data.frame(
  lip.func = diag(apply(lip.dist.mtx, 1, function(x) apply(func.dist, 1, function(y) cor(x, y, use = "complete.obs")))),
  expr.func = diag(apply(expr.dist.mtx, 1, function(x) apply(func.dist, 1, function(y) cor(x, y, use = "complete.obs"))))
)

str.cor$group = str.cor$lip.str > str.cor$expr.str
func.cor$group = func.cor$lip.func > func.cor$expr.func

p1 = ggplot(dist.cor, aes(x = expr.eucl, y = lip.eucl, label = num, fill = group)) +
  geom_vline(xintercept = 0, col = "darkgrey") +
  geom_hline(yintercept = 0, col = "darkgrey") +
  geom_abline(slope = 1, intercept = 0, lty = 6) +
  geom_point(pch = 21, col = "black", size = 4) +
  scale_fill_manual(values = c("lightgreen", "lightgrey")) +
  #geom_text(size = 3) +
  scale_x_continuous(limits = c(-0.8, 0.85)) +
  scale_y_continuous(limits = c(-0.8, 0.85)) +
  guides(fill = "none") +
  theme(aspect.ratio = 1, 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black")) +
  xlab("Correlataion with transcriptome") +
  ylab("Correlataion with lipidome") +
  ggtitle("Correlation to spatial distances")

p2 = ggplot(str.cor, aes(x = expr.str, y = lip.str, fill = group)) +
  geom_vline(xintercept = 0, col = "darkgrey") +
  geom_hline(yintercept = 0, col = "darkgrey") +
  geom_abline(slope = 1, intercept = 0, lty = 6) +
  geom_point(pch = 21, col = "black", size = 4) +
  scale_fill_manual(values = c("lightgreen", "lightgrey")) +
  #geom_text(size = 3) +
  scale_x_continuous(limits = c(-0.8, 0.85)) +
  scale_y_continuous(limits = c(-0.8, 0.85)) +
  guides(fill = "none") +
  theme(aspect.ratio = 1, 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black")) +
  xlab("Correlataion with transcriptome") +
  ylab("Correlataion with lipidome") +
  ggtitle("Correlation to structural connectivity")

p3 = ggplot(func.cor, aes(x = expr.func, y = lip.func, fill = group)) +
  geom_vline(xintercept = 0, col = "darkgrey") +
  geom_hline(yintercept = 0, col = "darkgrey") +
  geom_abline(slope = 1, intercept = 0, lty = 6) +
  geom_point(pch = 21, col = "black", size = 4) +
  scale_fill_manual(values = c("lightgreen", "lightgrey")) +
  #geom_text(size = 3) +
  scale_x_continuous(limits = c(-0.8, 0.85)) +
  scale_y_continuous(limits = c(-0.8, 0.85)) +
  guides(fill = "none") +
  theme(aspect.ratio = 1, 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black")) +
  xlab("Correlataion with transcriptome") +
  ylab("Correlataion with lipidome") +
  ggtitle("Correlation to functional connectivity")

grid.arrange(p1, p2, p3, nrow = 1)

a1 = melt(dist.cor[,1:2]) %>%
  ggplot(aes(variable, value)) +
  geom_boxplot(width = 0.3, fill = "lightgrey") +
  geom_jitter(width = 0.12) +
  theme_bw() +
  scale_x_discrete(labels = c("Lipidome", "Transcriptome")) +
  xlab("") +
  ylab("Correlation with spatial distances") +
  theme(aspect.ratio = 1) +
  ggtitle(paste0("p = ", round(wilcox.test(dist.cor$lip.eucl, dist.cor$expr.eucl, paired = TRUE, alternative = "less")$p.value, 3)))

a2 = melt(str.cor[,1:2]) %>%
  ggplot(aes(variable, value)) +
  geom_boxplot(width = 0.3, fill = "lightgrey") +
  geom_jitter(width = 0.12) +
  theme_bw() +
  scale_x_discrete(labels = c("Lipidome", "Transcriptome")) +
  xlab("") +
  ylab("Correlation with structural connectivity") +
  theme(aspect.ratio = 1) +
  ggtitle(paste0("p = ", round(wilcox.test(str.cor$lip.str, str.cor$expr.str, paired = TRUE, alternative = "greater")$p.value, 3)))

a3 = melt(func.cor[,1:2]) %>%
  ggplot(aes(variable, value)) +
  geom_boxplot(width = 0.3, fill = "lightgrey") +
  geom_jitter(width = 0.12) +
  theme_bw() +
  scale_x_discrete(labels = c("Lipidome", "Transcriptome")) +
  xlab("") +
  ylab("Correlation with functional connectivity") +
  theme(aspect.ratio = 1) +
  ggtitle(paste0("p = ", round(wilcox.test(func.cor$lip.func, func.cor$expr.func, paired = TRUE, alternative = "greater")$p.value, 3)))


plot_grid(a1, a2, a3, nrow = 1)