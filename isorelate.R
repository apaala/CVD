library(isoRelate)
#library(snpStats)
args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args)
pf=args[1]
mf=args[2]
outdir=args[3]
if(length(args)>3){
mafcut=args[4]
}else{
mafcut=0.001
}
### setting default values, can make these params user inputs later
imm=0.1
smm=0.1
chr=NULL
cores=10
minsnp=20
minbp=50000
err=0.001
#for plot settings (non exhaustive)
alpha=0.1
annotation_genes=NULL

####Read plink files
print("Reading plink files")
ped=read.table(pf, header=F, as.is=T)
map=read.table(mf, header=F, as.is=T)
pdmp=list(ped,map)  

###For test data set only!!!setting it to 2 instead of 0 as it wont pass filter otherwise. It can take 1 for haploid, 2 for diploid
pdmp[[1]][["V5"]]=2

#### Get Genotypes
print("Getting genotypes")
my_genotypes<-getGenotypes(ped.map=pdmp,reference.ped.map = NULL,maf=mafcut,isolate.max.missing = imm,snp.max.missing = smm,chromosomes = chr,input.map.distance = "cM",reference.map.distance = "cM")

#### Get Params
print("Generating parameters")
my_parameters <- getIBDparameters(ped.genotypes = my_genotypes, number.cores = cores)
write.table(my_parameters,paste0(outdir, "/ibdparameters.txt"),sep = "\t", row.names=FALSE, col.names = F, quote=FALSE)

####Get IBD segments
print("Getting IBD segments")
my_ibd <- getIBDsegments(ped.genotypes = my_genotypes,parameters = my_parameters, number.cores = cores, minimum.snps = minsnp,minimum.length.bp = minbp, error = err)

###Write result
write.table(my_ibd,paste0(outdir, "/ibdSegments.txt"),sep = "\t", row.names=FALSE, col.names = F, quote=FALSE)

####Get IBD summary
print("Generating summary")
summary=capture.output(getIBDsummary(ped.genotypes = my_genotypes, ibd.segments = my_ibd))
write.table(summary,paste0(outdir, "/ibdsummary.txt"),sep = "\t", row.names=FALSE, col.names = F, quote=FALSE)
write.table(getIBDsummary(ped.genotypes = my_genotypes, ibd.segments = my_ibd),paste0(outdir, "/ibdsummary.txt"),sep = "\t", row.names=FALSE, col.names = F, quote=FALSE)

####Can include if needed, was useful for my testing
#ped_1<-read.pedfile(p)
#geno=(ped_1$genotypes)
#ped_sum=col.summary(geno)
#Write the summary!
#ggplot(ped_sum, aes(x=ped_sum$MAF)) + geom_density()+labs(x="Normalized Counts", y="Density",fill="")+geom_hline(yintercept = 0)

####get posterior
print("Get IBDposterior")
my_posterior <- getIBDposterior(ped.genotypes = my_genotypes,parameters = my_parameters, number.cores = 1, error = 0.001)
write.table(my_posterior,paste0(outdir, "/ibdposterior.txt"),sep = "\t", row.names=FALSE, col.names = F, quote=FALSE)


####Plot IDB segments
print("Generating plot")
plot1=plotIBDsegments(ped.genotypes = my_genotypes, ibd.segments = my_ibd, interval = NULL,annotation.genes = NULL,
                annotation.genes.color = NULL, highlight.genes = NULL, highlight.genes.labels = FALSE, highlight.genes.color = NULL,
                highlight.genes.alpha = 0.1, segment.height = 0.6, number.per.page = NULL, fid.label = FALSE, 
                iid.label = FALSE, ylabel.size = 9, add.rug = FALSE, plot.title = "Distribution of IBD segments in PNG", 
                add.legend = TRUE, segment.color = NULL)

###write plot to outdir!
pdf(file=paste0(outdir,"IBD_Plot1.pdf"), width=11,height = 10)

plotIBDsegments(ped.genotypes = my_genotypes, ibd.segments = my_ibd, interval = NULL,annotation.genes = NULL,
                annotation.genes.color = NULL, highlight.genes = NULL, highlight.genes.labels = FALSE, highlight.genes.color = NULL,
                highlight.genes.alpha = 0.1, segment.height = 0.6, number.per.page = NULL, fid.label = FALSE,
                iid.label = FALSE, ylabel.size = 9, add.rug = FALSE, plot.title = "Distribution of IBD segments in PNG",
                add.legend = TRUE, segment.color = NULL)
dev.off()

### Binary matrix
my_matrix <- getIBDmatrix(ped.genotypes = my_genotypes, 
                          ibd.segments = my_ibd)

my_iR <- getIBDiR(ped.genotypes = my_genotypes, 
                  ibd.matrix = my_matrix, 
                  groups = NULL)

pdf(file=paste0(outdir,"IBD_Plot2.pdf"), width=11,height = 10)

plotIBDiR(ibd.iR = my_iR, 
          interval = NULL, 
          annotation.genes = NULL,
          annotation.genes.color = NULL,
          highlight.genes = NULL,
          highlight.genes.labels = FALSE,
          highlight.genes.color = NULL,
          highlight.genes.alpha = 0.1,
          point.size = 1,
          point.color = NULL,
          add.rug = FALSE, 
          plot.title = "Significance of IBD sharing", 
          add.legend = FALSE,
          facet.label = TRUE, 
          facet.scales = "fixed")

dev.off()

####Clusters

#my_p_clusters <- getIBDpclusters(ped.genotypes = my_genotypes, 
#                                 ibd.segments = my_ibd, 
#                                 prop=0.9, 
#                                 hi.clust = FALSE)




####Plot with interval, needs to be user input based.
#plot2=plotIBDsegments(ped.genotypes = my_genotypes, ibd.segments = my_ibd, interval = interval_list, annotation.genes = annotation_genes, annotation.genes.color = NULL, highlight.genes = highlight_genes, highlight.genes.labels = FALSE, highlight.genes.color = NULL, highlight.genes.alpha = 0.1, segment.height = 0.8, number.per.page = NULL, fid.label = FALSE, iid.label = FALSE, ylabel.size = 9, add.rug = TRUE, plot.title = "Distribution of IBD segments in PNG", add.legend = TRUE, segment.color = c("purple","green"))

####


