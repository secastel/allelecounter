# allelecounter
Counts the number of reads which map to either the reference or alternate allele at each heterozygous SNP.
Runs on Python2.7.x and has the following dependencies: [pyvcf](https://github.com/jamescasbon/PyVCF), [Samtools](http://www.htslib.org) (1.2+), awk.

#Usage
Requires a BAM, VCF, and Reference, produces read counts for each allele at each heterozygous SNP. Output is in the same format as the [GATK ASEReadCounterTool](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php) for cross compatibility. By default this tool will discard duplicate reads, and only count reads whose mates overlap once. For more flexibility please consider using the GATK tool.

##Arguments
* **--vcf** - VCF file containing genotype for the sample. Must be gzipped and Tabix indexed. To improve runtime this file should be pre-processed to only include bi-allelic heterozygous SNPs.
* **--sample** - Name of sample to use in VCF file.
* **--bam** - BAM file containing reads. Duplicates should be marked, and file must be indexed with samtools index.
* **--ref** - Reference sequence in FASTA format.
* **--min_cov** - Minimum coverage for a SNP before it is included in output. Minimum is 2 becuase samtools will not generate pileups for sites with only 1 read.
* **--min_baseq** - Minimum base quality at the SNP required for reads to be counted.
* **--min_mapq** - Mimimum mapping qualityfor reads to be counted. Note that due to limitations with Samtools mpileup the maximum value is 93. Any value higher than 93 will be set to 93.
* **--o** - Output file name.
