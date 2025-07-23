#cn.mops.R
install.packages(cn.mops)
BiocManager::install("cn.mops") #either this or line 2
library(cn.mops)
BAMFiles <- "/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/Bam/OICRM4CA-07-01-P.bam"
bamDataRanges <- getReadCountsFromBAM(BAMFiles, parallel=8) #parallel = number of cores
# Run with specified region and window length
res <- cn.mops(bamDataRanges)
plot(res,which=1)
resCNMOPS <- cn.mops(XRanges) # run cn.mops on the GRanges object
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS) # calculate integer copy number
(resCNMOPS)
cnvs(resCNMOPS)[1:5]
cnvr(resCNMOPS)[1,1:5]
segplot(resCNMOPS,sampleIdx=13) #chromosome plots
plot(resCNMOPS,which=1) #cnv region plots


