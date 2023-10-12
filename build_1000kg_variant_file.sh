ml bcftools/1.18

hapmap3_dbsnp="/gpfs/space/GI/eQTLGen/public_data/dbSNP/GCF_000001405.filtered.hm3.chrSorted.40.vcf.gz"

for chr in {1..22}; do
  filename=ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

  bcftools query ${filename} --include 'INFO/EUR_AF[0]>0.05&INFO/EUR_AF[0]<0.95' \
  -f '%CHROM\t%POS\t%POS\t%ID\n' > processed/variants.chr${chr}.filtered.txt
done

ml purge

cat processed/variants.chr*.filtered.txt > processed/variants.chrALL.M_5_50.txt

ml any/python/3.9.9
#pip3 install CrossMap

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

CrossMap.py bed processed/hg19ToHg38.over.chain.gz processed/variants.chrALL.M_5_50.txt processed/variants.chrALL.hg38.M_5_50.txt

ml bcftools/1.18
bcftools query ${hapmap3_dbsnp} -R processed/variants.chrALL.hg38.M_5_50.txt \
  -f '%CHROM\t%POS\t%POS\t%ID\n' > processed/variants.chrALL.hg38.w_hm3.M_5_50.bed