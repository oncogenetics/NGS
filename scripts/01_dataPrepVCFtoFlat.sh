# 08/09/2019
# Workspace -------------------------------------------------------------------
module load bcftools/1.5

cd /data/rds/DGE/DUDGE/OGENETIC/Data/VCF_repository/Filtered_Reannotated_VEP

folderVCF=/data/rds/DGE/DUDGE/OGENETIC/Data/VCF_repository/Filtered_Reannotated_VEP
folderOut=/data/rds/DGE/DUDGE/OGENETIC/Data/VCF_repository/Filtered_Reannotated_VEP/ShinyAppInput

# vcf to flat -----------------------------------------------------------------
for i in `ls *.vcf.gz`; do \
fileVCF=${folderVCF}/${i}
fileOut=${folderOut}/`basename ${fileVCF}`
bcftools query -f'%CHROM %POS %ID %REF %ALT[ %GT]\n' ${fileVCF} > ${fileOut}.GT
echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tCLNSIG" > ${fileOut}.CLNSIG
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/CLNSIG\n' ${fileVCF} >> ${fileOut}.CLNSIG
echo -e "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|GENE_PHENO|SIFT|PolyPhen|LoF|LoF_filter|LoF_flags|LoF_info|ExAC_AF|ExAC_AF_AFR|ExAC_AF_AMR|ExAC_AF_Adj|ExAC_AF_CONSANGUINEOUS|ExAC_AF_EAS|ExAC_AF_FEMALE|ExAC_AF_FIN|ExAC_AF_MALE|ExAC_AF_NFE|ExAC_AF_OTH|ExAC_AF_POPMAX|ExAC_AF_SAS|CADD_PHRED|CADD_RAW|REVEL|gnomAD_exomes|gnomAD_exomes_AF_nfe" > ${fileOut}.CSQ
bcftools query -f '%INFO/CSQ\n' ${fileVCF} >> ${fileOut}.CSQ
bcftools query -l ${fileVCF} > ${fileOut}.SAMPLE
done


#zip temp_file_transfer.zip ShinyAppInput/*

# make ped map ----------------------------------------------------------------
#module load plink/1.90.beta2k
#plink --vcf ${fileVCF} --no-sex --no-pheno --double-id --vcf-half-call m --keep-allele-order --recode 12 --out ${fileOut}


