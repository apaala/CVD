library(moimix)
library(SeqArray)
args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args)
vcf=args[1]
outdir=args[2]
#varid=args[3]

#####Requires gzipped vcf and tabix

fname=paste0(outdir,"/",basename(vcf),".gds")

seqVCF2GDS(vcf,fname, storage.option="LZ4_RA")

isolates=seqOpen(fname)

seqSummary(isolates)

####save sample IDs
sample.id <- seqGetData(isolates, "sample.id")

####get genomic coords for all variants
coords<- getCoordinates(isolates)

####Filter variant ids not
#seqSetFilter(isolates, variant.id = coords$variant.id[coords$chromosome != varid])

# matrix of read counts for ref and alt alleles
#counts <- alleleCounts(isolates)
set.seed(2002)

####Estimate BAF matrix
#isolate_baf <- bafMatrix(isolates)
#for(sampleid in sample.id){
#plot(isolate_baf, sampleid)

# fit a model on sample, multiple clone infection-find a way to define k
#m1 <- binommix(counts, sample.id = sampleid, k = 2)
#summary(m1)

####Estimates moi with binommix
#param.estimates <- getTheta(m1)
#param.estimates
#pdf(paste0(outdir,"/",basename(vcf),sampleid
#plot(isolate_baf, sampleid)
#abline(h = param.estimates$mu.hat)

#####Estimate moi with Fws
fws_all <- getFws(isolates)
hist(fws_all)
pedname=fname=paste0(outdir,"/moimix_out_",basename(vcf))
extractPED(isolates, moi.estimates=fws_all, outfile=pedname)
# see if our sample that we estimated is multiclonal according to fws

for(sampleid in sample.id){
fws_all[sampleid] < 0.95
}


