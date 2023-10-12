
rm variants.chrALL.hg38.w_hm3.M_5_50.bed
touch variants.chrALL.hg38.w_hm3.M_5_50.bed

for chr in {1..22}; do
  echo $chr
  zcat /gpfs/space/GI/eQTLGen/public_data/ldsc/eur_w_ld_chr/${chr}.l2.ldscore.gz \
  | awk 'BEGIN{FS="\t"; OFS=FS} $5 > 0.05 && $5 < 0.95 { print $1,$3,$3,$2 }' >> variants.chrALL.hg38.w_hm3.M_5_50.bed
done