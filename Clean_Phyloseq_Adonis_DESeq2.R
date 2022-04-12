library("ggplot2")
library("phyloseq")
library("DESeq2")
library("vegan")

##Loading OTU and taxa tables
my_data.frame = read.csv(file = "Shwachman_otutable.csv", header = T, sep = ",")
rownames(my_data.frame) <- paste0("OTU", 1:nrow(my_data.frame))
my_otu_table = otu_table(my_data.frame, taxa_are_rows = TRUE)

my_tax_data.frame = read.csv(file = "Shwachman_taxtable.csv", header=F, sep=",")
rownames(my_tax_data.frame) <- paste0("OTU", 1:nrow(my_tax_data.frame))
my_tax_tab = tax_table(my_tax_data.frame, errorIfNULL = TRUE)
rownames(my_tax_tab) <- paste0("OTU", 1:nrow(my_tax_tab))
colnames(my_tax_tab) <- c("Genus", "Species")

##Create phyloseq class experiment-level object
physeq = phyloseq(my_otu_table, my_tax_tab)

## Loading Sample data and labeling metadata
sampledata = read.csv("meta_Shwachman.csv", header = F, sep = ",")
sampledata <- as.data.frame(sampledata)
rownames(sampledata) <- sampledata[,1]
sampledata[,1] <- NULL
colnames(sampledata) <- c("ETF", "PtID", "PK", "BSvsNS")
sampledata = sample_data(sampledata)
#```

#########PHYLOSEQ############


##Merge phyloseq and sampledata to create new phylo-class experiment-level object with taxa, samples, and variables
ps = merge_phyloseq(sampledata, physeq)
ps

##We are now ready to use phyloseq.

theme_set(theme_bw())

##**Visualize alpha-diversity**:
##  Use `plot_richness`

##```{r richness, warning=FALSE}
apl = plot_richness(ps, 
              x="ETF",
              measures=c("Observed", "Shannon", "InvSimpson"),
              color="PtID")

newSTorder = c("E", "T", "F")


apl$data$ETF <-as.character(apl$data$ETF)
apl$data$ETF <- factor(apl$data$ETF, levels = newSTorder)
print (apl)
##```

##**Ordinate**:

#Ordinate ETF
##  ```{r ordinate}
ord.nmds.bray <- ordinate(ps,
                          method="NMDS",
                          distance="bray")
bdl = plot_ordination(physeq = ps,
                ordination = ord.nmds.bray,
                color="ETF",
                title="Bray NMDS") + 
  geom_point(size = 6) 

newSTorder = c("E", "T", "F")

bdl$data$ETF <-as.character(bdl$data$ETF)
bdl$data$ETF <- factor(bdl$data$ETF, levels = newSTorder)

#Add ellipse with t distribution
ebdl = bdl + stat_ellipse(type = "t")
print (ebdl)

##Ordinate PK
ord.nmds.bray <- ordinate(ps,
                          method="NMDS",
                          distance="bray")
PKplot = plot_ordination(physeq = ps,
                      ordination = ord.nmds.bray,
                      color="PK",
                      title="Bray NMDS PK/PD") + 
  geom_point(size = 6) 

#Add ellipse with t distribution
EPKplot = PKplot + stat_ellipse(type = "t")
print (EPKplot)

##Ordinate BSvsNS
ord.nmds.bray <- ordinate(ps,
                          method="NMDS",
                          distance="bray")
ASplot = plot_ordination(physeq = ps,
                ordination = ord.nmds.bray,
                color="BSvsNS",
                title="Bray NMDS Antibiotic Spectrum") + 
  geom_point(size = 6) 

#Add ellipse with t distribution
EASplot = ASplot + stat_ellipse(type = "t")
print (EASplot)

##**Bar plot**:

##  ```{r bar-plot}
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
rapl = plot_bar(ps.top20, x="PtID", 
         fill="ETF", title = "Relative Abundance Top 20 Species") + facet_wrap(~Species, scales="free_x") + theme(axis.text.x = element_text(size = 5))
##```
newSTorder = c("E", "T", "F")

rapl$data$ETF <-as.character(rapl$data$ETF)
rapl$data$ETF <- factor(rapl$data$ETF, levels = newSTorder)
print (rapl)


#######PERMANOVA################

##Start of Permanova
##My phyloseq object called physeq2; tutorial called erie

##PK

#Calculate bray curtis distance matrix
erie_bray <- phyloseq::distance(ps, method = "bray")
##can also do using "jaccard"

#Make a data frame from the sample data
sampledf <- data.frame(sample_data(ps))

#Adonis test 
adonis(erie_bray ~ PK, data=sampledf, strata = sampledf$PtID)

#Homogeneity of dispersion test
beta <- betadisper(erie_bray, sampledf$PK)
permutest(beta)

##Antibiotic Spectrum
#Adonis test 
adonis(erie_bray ~ BSvsNS, data=sampledf, strata = sampledf$PtID)

#Homogeneity of dispersion test
beta <- betadisper(erie_bray, sampledf$BSvsNS)
permutest(beta)

######DESEQ2###############

##Merge phyloseq and sampledata to create new phylo-class experiment-level object with taxa, samples, and variables
physeq2 = merge_phyloseq(sampledata, physeq)

##PK

# Changing phyloseq object into DEseq object
patho.des = phyloseq_to_deseq2(physeq2, ~ PK)

# calculate geometric means prior to estimate size factors - this is because we have zeros
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(patho.des), 1, gm_mean)

#Size factor and dispersion estimates
patho.des = estimateSizeFactors(patho.des, geoMeans = geoMeans)
patho.des <- estimateDispersions(patho.des)
patho.des.var <- getVarianceStabilizedData(patho.des)
patho.des.deseq = DESeq(patho.des, fitType="local")


# Filtering to find the top most differentially abundant OTUs, ordering by p value, and removing taxa with NA value
#contrast - "variable", "numerator", "denominator", so positive log fold change would be higher in numerator 
#and negative log fold change would be higher in the denominator
res = results(patho.des.deseq, cooksCutoff = FALSE, contrast = c("PK", "Yes", "No"))
res = res[order(res$padj, na.last=NA), ]
res 

#export res table into excel
write.csv(x=res, file="PK_res.csv")

#set significance
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq2)[rownames(sigtab), ], "matrix"))

head(sigtab)

dim(sigtab)

#export significant log2 fold change into excel
write.csv(x=sigtab, file="PK_sigtab.csv")

#### Visualizing differentially abundant OTUs 

#Manhattan Plot, p<0.1
plotMA(res, ylim=c(-30,30))

### Genus (Genus) order ###
x = tapply(sigtab$log2FoldChange, sigtab$PK, function(x) max(x))
x = sort(x, TRUE)
sigtab$PK = factor(as.character(sigtab$PK), levels=names(x))

# Colored by Species
patho.sp.genus.l2fc_plot <- ggplot(sigtab, aes(x=log2FoldChange, y=Genus, color=Species)) + geom_vline(xintercept = 0.0, color = "gray", size = 1) + geom_point(size=3) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
patho.sp.genus.l2fc_plot
ggsave("patho.sp.genus.l2fc.pdf", patho.sp.genus.l2fc_plot)


##Antibiotic Spectrum

# Changing phyloseq object into DEseq object
patho.des = phyloseq_to_deseq2(physeq2, ~ BSvsNS)

# calculate geometric means prior to estimate size factors - this is because we have zeros
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(patho.des), 1, gm_mean)

#Size factor and dispersion estimates
patho.des = estimateSizeFactors(patho.des, geoMeans = geoMeans)
patho.des <- estimateDispersions(patho.des)
patho.des.var <- getVarianceStabilizedData(patho.des)
patho.des.deseq = DESeq(patho.des, fitType="local")


# Filtering to find the top most differentially abundant OTUs, ordering by p value, and removing taxa with NA value
#contrast - "variable", "numerator", "denominator", so positive log fold change would be higher in numerator 
#and negative log fold change would be higher in the denominator
res = results(patho.des.deseq, cooksCutoff = FALSE, contrast = c("BSvsNS", "BS", "NS"))
res = res[order(res$padj, na.last=NA), ]
res 

#export res table into excel
write.csv(x=res, file="BSvsNS_res.csv")

#set significance
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq2)[rownames(sigtab), ], "matrix"))

head(sigtab)

dim(sigtab)

#export significant log2 fold change into excel
write.csv(x=sigtab, file="BSvsNS_sigtab.csv")

#### Visualizing differentially abundant OTUs 

#Manhattan Plot, p<0.05
plotMA(res, ylim=c(-30,30))

### Genus (Genus) order ###
x = tapply(sigtab$log2FoldChange, sigtab$BSvsNS, function(x) max(x))
x = sort(x, TRUE)
sigtab$BSvsNS = factor(as.character(sigtab$BSvsNS), levels=names(x))

# Colored by Species
patho.sp.genus.l2fc_plot <- ggplot(sigtab, aes(x=log2FoldChange, y=Genus, color=Species)) + geom_vline(xintercept = 0.0, color = "gray", size = 1) + geom_point(size=3) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
patho.sp.genus.l2fc_plot
ggsave("patho.sp.genus.l2fc.pdf", patho.sp.genus.l2fc_plot)
