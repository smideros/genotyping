#file: ~\workdir\genotyping\genotypes.sh
#date: 5november2015

# Installed vcftools:

#vcftools --gzvcf settu_28a_unique.vcf --freq --out settu_28a_unique_analysis
#vcftools --gzvcf settu_28a_unique.vcf --freq2 --out settu_28a_unique_analysis2
#vcftools --gzvcf /Users/smideros/Box\ Sync/Pathogenesis/Seto/data\ \&\ analysis\ 2015/settu_28a_unique.vcf --counts --out settu_28a_unique_analysis
#vcftools --gzvcf /Users/smideros/Box\ Sync/Pathogenesis/Seto/data\ \&\ analysis\ 2015/settu_28a_unique.vcf --depth --out settu_28a_unique_analysis
#vcftools --gzvcf /Users/smideros/Box\ Sync/Pathogenesis/Seto/data\ \&\ analysis\ 2015/settu_28a_unique.vcf --site-depth --out settu_28a_unique_analysis
#vcftools --gzvcf /Users/smideros/Box\ Sync/Pathogenesis/Seto/data\ \&\ analysis\ 2015/settu_28a_unique.vcf --site-mean-depth --out settu_28a_unique_analysis
#vcftools --gzvcf /Users/smideros/Box\ Sync/Pathogenesis/Seto/data\ \&\ analysis\ 2015/settu_28a_unique.vcf --geno-depth --out settu_28a_unique_analysis
#vcftools --gzvcf /Users/smideros/Box\ Sync/Pathogenesis/Seto/data\ \&\ analysis\ 2015/settu_28a_unique.vcf --het --out settu_28a_unique_analysis

# Installed vcflib


vcfhetcount /Users/smideros/Box\ Sync/Pathogenesis/Seto/data\ \&\ analysis\ 2015/settu_28a_unique.vcf > settu_28a_unique.vcf.hetcount
vcfhethomratio /Users/smideros/Box\ Sync/Pathogenesis/Seto/data\ \&\ analysis\ 2015/settu_28a_unique.vcf > settu_28a_unique.vcf.hetratio
vcf2tsv -g /Users/smideros/Box\ Sync/Pathogenesis/Seto/data\ \&\ analysis\ 2015/settu_28a_unique.vcf > settu_28a_unique.vcf.tsv
