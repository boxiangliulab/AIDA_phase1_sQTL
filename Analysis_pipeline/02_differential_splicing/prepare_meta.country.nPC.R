library(stringr)
library(dplyr)
library(data.table)

dat <- fread('/ebs1/zhangyuntian/project/aida/aida_perind_numers.counts.gz')
meta <- read.delim('~/project/proj_aida_sqtl/scripts/batch_v1_0/30_DS_pseudobulk/AIDA_fist_phase_metadata.ind.txt')

#ancestry_PC <- read.csv('/ebs1/users/tianchi/project/proj_aida_sqtl/scripts/batch_v1_0/21_spliZ_post/linear_model/2987_withou_mix_blood_pca_ethnicity_result.csv')
#ID_convert <- read.table('/ebs1/users/tianchi/project/proj_aida_sqtl/scripts/batch_v1_0/01_STARsolo/conversion_DCPID_sampleID.csv')
#ancestry_PC <- merge(ancestry_PC[, c(1,3,4,5,6,7)], ID_convert, by = 'V1')
#colnames(ancestry_PC) <- c('GenotypeID', 'gPC1', 'gPC2', 'gPC3', 'gPC4', 'gPC5', 'DCP_ID')

#meta.all <- merge(meta, ancestry_PC, by='DCP_ID')



celltypes <- read.table('/ebs1/users/tianchi/project/proj_aida_sqtl/scripts/batch_v1_0/summary_cell_type/cell_type_use.txt')$V1

for (ct in celltypes){
  #print(ct)
  idx <- which(ct == gsub('.*\\.(.*)\\..*', '\\1', colnames(dat)))
  #tb.tmp <- dat[, ..idx]
  #rownames(tb.tmp) <- dat$V1
  #write.table(tb.tmp, paste0('/ebs2/aida-sqtl/FreezeV1_0/31_DS_pseudobulk_leafcutter/aida_perind_numers.', ct, '.counts'), quote=F, sep=' ')
  #system(paste0('gzip', ' /ebs2/aida-sqtl/FreezeV1_0/31_DS_pseudobulk_leafcutter/aida_perind_numers.', ct, '.counts'))
  
  pc <- as.data.frame(t(read.delim(paste0('/ebs1/zhangyuntian/project/aida/PC_file/new3_21_PC_freeeze/', ct, '_PC.txt'), row.names=1)))
  colnames(pc) <- c('pPC1', 'pPC2', 'pPC3', 'pPC4', 'pPC5', 'pPC6', 'pPC7', 'pPC8', 'gPC1', 'gPC2', 'gPC3', 'gPC4', 'gPC5', 'sex', 'age')
  pc$DCP_ID <- rownames(pc)

  dat.tmp <- data.frame(V1 = colnames(dat)[idx])
  dat.tmp[c('DCP_ID','celltype')] <- str_split_fixed(dat.tmp$V1, '\\.', 3)[,1:2]
  merge <- merge(dat.tmp, pc, by='DCP_ID')
  merge <- merge(merge, meta[, c('DCP_ID', 'ethnicity')], by='DCP_ID')

  #df.tmp.merge = merge %>% mutate(
  #  ancestry = case_when(
  #    ethnicity == 'Chinese' ~ "East",
  #    ethnicity == 'Indian' ~ "South",
  #    ethnicity == 'Japanese' ~ "East",
  #    ethnicity == 'Korean' ~ "East",
  #    ethnicity == 'Malay' ~ "Southeast"
  #  )
  #)
  #country_list <- sort(unique(merge$ethnicity))
  country_list <- c('Chinese', 'Indian', 'Malay')
  for (i in 1:length(country_list)){
    for (j in 1:length(country_list)){
      if (j > i){
        countryA <- country_list[i]
	countryB <- country_list[j]
        df.tmp.merge <- subset(merge, ethnicity == countryA | ethnicity == countryB)
	#df.tmp.merge$country <- ifelse(df.tmp.merge$ethnicity == countryA, 1, 0)

	pcs <- c('pPC1', 'pPC2', 'pPC3', 'pPC4', 'pPC5', 'pPC6', 'pPC7', 'pPC8')
  	for (n in 1:8){
    	  nPC <- pcs[1:n]
	  write.table(df.tmp.merge[, c('V1', 'ethnicity', 'sex', 'age', nPC)], paste0('/ebs2/aida-sqtl/FreezeV1_0/31_DS_pseudobulk_leafcutter_include_splicing_PC/groups_file.country.', countryA, '_vs_', countryB, '.', n, 'PC.', ct, '.txt'), quote=F, row.names=F, col.names=F, sep='\t')
	}
      }
    }
  }
}
