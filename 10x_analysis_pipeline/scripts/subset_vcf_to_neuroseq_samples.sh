# This script prepares a VCF file for use with demuxlet.
# Note that the VCF produced is likely not useful for other purposes.

INDIR=/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Filtered
OUTDIR=/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/genotypes

#### Subsetting VCFs ####

# only need a subset of all hipsci genotypes
# note - check that I only use open access genotypes here - OUTDIR can only contain open access data
SAMPLES=/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/metadata/neuroseq_cell_line_ids.txt

# select only relatively common exome SNPs
REGIONS=/hps/nobackup/hipsci/scratch/singlecell_endodiff/resources/references/gnomad.exomes.r2.0.2.sites.allele-freq-only.common.biallelic.vcf.gz

# loop over each autosomal chromosome, subset to regions and samples, and ensure the output is sorted
#### TEST - don't subset samples - we give a sample list to demuxlet anyway

for chr in {1..22}; do
    INFILE=$INDIR/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr.$chr.norm.renamed.recode.vcf.gz
    OUTFILE=$OUTDIR/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr.$chr.norm.renamed.recode.gnomad.exomes.common.biallelic.neuroseq_donors.vcf.gz
    OUTFILE2=$OUTDIR/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr.$chr.norm.renamed.recode.gnomad.exomes.common.biallelic.neuroseq_donors.sorted.vcf.gz
#    bcftools view -R $REGIONS -S $SAMPLES -Oz -o $OUTFILE $INFILE
    bcftools view -R $REGIONS -Oz -o $OUTFILE $INFILE
    bcftools sort -Oz -o $OUTFILE2 $OUTFILE
    rm $OUTFILE
done


### Concatenating VCFs ####

# remove any previous file with this name
rm $OUTDIR/sorted_list_of_chr_vcfs.txt
for chr in 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9; do
echo $OUTDIR/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr.$chr.norm.renamed.recode.gnomad.exomes.common.biallelic.neuroseq_donors.sorted.vcf.gz >> $OUTDIR/sorted_list_of_chr_vcfs.txt ;
done

OUTFILE3=$OUTDIR/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.allchr.norm.renamed.recode.gnomad.exomes.common.biallelic.neuroseq_donors.vcf.gz
bcftools concat -o $OUTFILE3 -f $OUTDIR/sorted_list_of_chr_vcfs.txt -O z
tabix $OUTFILE3


#### Final adjustment of VCF for use with demuxlet ####
#remove unnecessary contigs from header that crash demuxlet

OUTFILE4=$OUTDIR/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.allchr.norm.renamed.recode.gnomad.exomes.common.biallelic.neuroseq_donors.vcf
zcat $OUTFILE3 | grep -v -e '##contig=<ID=X' -e '##contig=<ID=Y' -e '##contig=<ID=MT' > $OUTFILE4
