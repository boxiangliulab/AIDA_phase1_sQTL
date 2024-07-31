source('smr/example/plot/plot_SMR.mod.r')
source('smr/example/plot/image.scale.r')

# Read the data file in R:
SMRData = ReadSMRData("plot/myplot.ENSG00000143889.txt")
# Plot the SMR results in a genomic region centred around a probe:
pdf('plot.SMRLocusPlot.pdf', width=8,height=8)
SMRLocusPlot(data=SMRData, smr_thresh=8.4e-6, heidi_thresh=0.05, plotWindow=1000, max_anno_probe=16)
dev.off()
# smr_thresh: genome-wide significance level for the SMR test.
# heidi_thresh: threshold for the HEIDI test. The default value is 0.05.
# cis_wind: size of a window centred around the probe to select cis-eQTLs for plot. The default value is 2000Kb.
# max_anno_probe: maximum number of probe names to be displayed on the figure. The default value is 16.

pdf('plot.SMREffectPlot.pdf', width=7,height=7)
SMREffectPlot(data=SMRData, trait_name="lymphocyteCount") 
dev.off()
