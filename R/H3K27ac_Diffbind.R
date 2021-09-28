#################################################################################################################################################
# Metadata Preparaton
bamsToCount_Rep1 <- dir("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/bowtie_out", full.names = TRUE, pattern = "*.25M.bam")
bamsToCount_Rep2 <- dir("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_Rep2_analysis/bowtie_out", full.names = TRUE, pattern = "*.dup.filtered.srt.q10.bam$")
H3K27ac.bam = c(bamsToCount_Rep1,bamsToCount_Rep2)
H3K27ac_metadata = data.frame(SampleID=c(gsub(".25M.bam|.Rep2.dup.filtered.srt.q10.bam","",basename(H3K27ac.bam))))
H3K27ac_metadata$Tissue = c(rep("DCs",28))
H3K27ac_metadata$Factor = c(rep("H3K27ac",28))
H3K27ac_metadata$Condition = factor(H3K27ac_metadata$SampleID)
H3K27ac_metadata$Treatment = factor(gsub(".25M.bam|.Rep2.dup.filtered.srt.q10.bam","",basename(H3K27ac.bam)))
H3K27ac_metadata$Replicate = factor(c(rep(1,14),rep(2,14)))

# Replcae "NCoR1.2hr.CpG" to "Emp.0hr" and vice versa
H3K27ac.bam[c(15,24)] = H3K27ac.bam[c(24,15)] #c("NCoR1.2hr.CpG","Emp.0hr")
H3K27ac_metadata$bamReads = H3K27ac.bam
H3K27ac_metadata$ControlID = c(rep("Input",28))
H3K27ac_metadata$bamControl = c(rep("/home/imgsb/Gyan/NCOR1/Analysis_result_of_NCoR1_ChIPSeq/NCOR1_mm10_alignment/Uniquley_mapped_reads/sorted_NCoR1_Ch_Input_R1.bam",28))
H3K27ac_metadata$SampleID = paste0(H3K27ac_metadata$SampleID,"_",H3K27ac_metadata$Replicate)

peak_Rep1 <- dir("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks", full.names = TRUE, pattern = "*_summits_2kb.bed")[-8]
peak_Rep2 <- dir("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_Rep2_analysis/macs_peak/", full.names = TRUE, pattern = "*_summits_2kb.bed")
H3K27ac.peak = c(peak_Rep1,peak_Rep2)
H3K27ac.peak[c(15,24)] = H3K27ac.peak[c(24,15)]
H3K27ac_metadata$Peaks  = H3K27ac.peak
H3K27ac_metadata$PeakCaller = c(rep("bed",28))
################################################################################################################################################
H3K27ac.diffbind = dba(sampleSheet = H3K27ac_metadata)
plot(H3K27ac.diffbind)




H3K27ac.diffbind <- dba.count(H3K27ac.diffbind, summits=500,bParallel=TRUE)
H3K27ac.diffbind <- dba.normalize(H3K27ac.diffbind, bRetrieve=TRUE)


H3K27ac.diffbind = dba.contrast(H3K27ac.diffbind, categories=c(DBA_CONDITION))
H3K27ac.diffbind = dba.analyze(H3K27ac.diffbind,bParallel=TRUE)
H3K27ac.diffbind.DB <- dba.report(H3K27ac.diffbind)
DiffBind::dba.plotHeatmap(H3K27ac.diffbind, contrast=1, correlations=FALSE)


dba.plotMA(H3K27ac.diffbind,  contrast=c("Condition","NCoR1_2hr_pIC",""),
            bNormalized=FALSE, sub="Non-Normalized")
dba.plotMA(H3K27ac.diffbind,bXY=FALSE,bSmooth=FALSE,bNormalized=TRUE)
