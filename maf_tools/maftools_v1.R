library(maftools)
library(glue)
library(ggplot2)
library(readr)
library(glue)




file_prefixes <- c("Tumour_Ade", "Tumour_AMN", "Metastasis","TumorOnly")


for (file_prefix in file_prefixes) {

    laml.maf =glue("{file_prefix}_noGermline_NoCommonvariant_VAFGT0.01_LT0.4.tsv")  # _SelectedGenes #"TumorOnly_noGermline_NoCommonvariant_VAFGT0.01_LT0.4_SelectedGenes.tsv" "TumorOnly_noGermline_NoCommonvariant.tsv"
    laml.clin=glue("{file_prefix}_clinical_details.tsv")
    prefix=glue("{file_prefix}_Complete_Genes_Group_category2_") #Complete_ #Selected
    clin_df <- read_delim(laml.clin , delim = "\t")

    #MAF object
    laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
    #Plotting MAF summary.
    pdf(file = glue("{prefix}Variant_summary.pdf"))
    plotmafSummary(maf = laml, rmOutlier = FALSE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()

    #Color coding 
    fabcolors = RColorBrewer::brewer.pal(n = 1,name = 'Spectral')
    names(fabcolors) = c(file_prefix)
    fabcolors = list(Group = fabcolors)

    #Drawing oncoplots
    pdf(file = glue("{prefix}Genes_Oncoplot.pdf"))
    oncoplot(maf = laml, top = 20,clinicalFeatures = "Group",draw_titv = TRUE,   sortByAnnotation = TRUE)
    dev.off()

    #Drawing oncoplots  with pathway
    pdf(file = glue("{prefix}Genes_Oncoplot_Pathway.pdf"), width = 10, height = 6)
    oncoplot(maf = laml, top = 20,clinicalFeatures = "Group",  sortByAnnotation = TRUE,
    annotationColor = fabcolors,pathways = "auto",fontSize = 0.4,gene_mar = 8)
    dev.off()

    }


## Signature analysis

library("BSgenome.Hsapiens.UCSC.hg38")
library('NMF')
library('pheatmap')
library(maftools)
library(glue)
library(ggplot2)
library(readr)
library(glue)


file_prefixes <- c("Tumour_Ade", "Tumour_AMN", "Metastasis","TumorOnly")

for (file_prefix in file_prefixes) {

file_prefix="TumorOnly"
laml.maf =glue("{file_prefix}_noGermline_NoCommonvariant_VAFGT0.01_LT0.4.tsv")  # _SelectedGenes #"TumorOnly_noGermline_NoCommonvariant_VAFGT0.01_LT0.4_SelectedGenes.tsv" "TumorOnly_noGermline_NoCommonvariant.tsv"
laml.clin=glue("{file_prefix}_clinical_details.tsv")
prefix=glue("{file_prefix}_Complete_Genes_Group_category2_") #Complete_ #Selected
clin_df <- read_delim(laml.clin , delim = "\t")

#MAF object
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

laml.tnm = trinucleotideMatrix(maf = laml,add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6,pConstant = 1e-9)

pdf(file = glue("{prefix}Cophenetic_plot.pdf"), width = 10, height = 6)
plotCophenetic(res = laml.sign)
dev.off()

laml.sig = extractSignatures(mat = laml.tnm, n = 3,pConstant = 1e-9)

#Compate against original 30 signatures 
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")

pdf(file = glue("{prefix}cosinesimilarity_against_original_30_signatures_plot.pdf"), width = 10, height = 6)
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
dev.off()

laml.og30_df=as.data.frame(unlist(laml.og30.cosm$best_match))
write.csv(laml.og30_df,glue("{prefix}cosinesimilarity_against_original_30_signatures.csv",row.name=None))

#Compate against updated version3 60 signatures 
laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")

pdf(file = glue("{prefix}cosinesimilarity_against_updatedveersion3_60_signatures_plot.pdf"), width = 10, height = 6)
pheatmap::pheatmap(mat = laml.v3.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
dev.off()

laml.og60_df=as.data.frame(unlist(laml.v3.cosm$best_match))
write.csv(laml.og60_df,glue("{prefix}cosinesimilarity_against_updatedveersion3_60_signatures.csv",row.name=None))

#plot signatures

pdf(file = glue("{prefix}signatures_plot.pdf"), width = 10, height = 6)
maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "SBS")
dev.off()


########## Somatic Interactions
pdf(file = glue("{prefix}somaticInteration_plot.pdf"), width = 10, height = 6)
par(mar = c(8, 4.1, 4.1, 4.1))

tryCatch({somaticInteractions_df=somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1),fontSize=0.6)
}, error = function(e) {
  # Error handling code
  cat("An error occurred:", conditionMessage(e), "\n")})

dev.off()
write.csv(somaticInteractions_df,glue("{prefix}somaticInteration_plot.csv",row.name=None))


#Detecting cancer driver genes based on positional clustering
laml.sig = oncodrive(maf = laml, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
head(laml.sig)

pdf(file = glue("{prefix}DriverGenes_plot.pdf"), width = 10, height = 6)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
dev.off()
write.csv(laml.sig,glue("{prefix}DriverGenes_plot.csv",row.name=None))


##Adding and summarizing pfam domains
#laml.pfam = pfamDomains(maf = laml, AACol = 'HGVSp_Short', top = 10)


#Oncogenic Signaling Pathways
#pws = pathways(maf = laml, plotType = 'treemap')

}
