# allelecounter
Counts the number of reads which map to either the reference or alternate allele at each heterozygous SNP.
Runs on Python2.7.x and has the following dependencies: [pyvcf](https://github.com/jamescasbon/PyVCF), [Samtools](http://www.htslib.org) (1.2+), awk.

# Citation
Castel, S. E., Levy-Moonshine, A., Mohammadi, P., Banks, E. & Lappalainen, T. Tools and best practices for data processing in allelic expression analysis. Genome Biol. 16, 195 (2015).

# Usage
Requires a BAM, VCF, and Reference, produces read counts for each allele at each heterozygous SNP. Output is in the same format as the [GATK ASEReadCounterTool](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php) for cross compatibility. If you would like to exclude duplicate reads they must first be marked using a separate tool, as this script does not mark duplicate reads itself. I suggest using [Picard](https://broadinstitute.github.io/picard/) MarkDuplicates. If duplicate reads have been marked, this tool will discard them. Reads whose mates overlap will only be counted once. For more flexibility please consider using the GATK tool.

## Arguments
* **--vcf** - VCF file containing genotype for the sample. Must be gzipped and Tabix indexed. To improve runtime this file should be pre-processed to only include bi-allelic heterozygous SNPs.
* **--sample** - Name of sample to use in VCF file.
* **--bam** - BAM file containing reads. Duplicates should be marked, and file must be indexed with Samtools index.
* **--ref** - Reference sequence in FASTA format.
* **--min_cov** - Minimum coverage for a SNP before it is included in output. Minimum is 2 becuase samtools will not generate pileups for sites with only 1 read.
* **--min_baseq** - Minimum base quality at the SNP required for reads to be counted.
* **--min_mapq** - Mimimum mapping qualityfor reads to be counted. Note that due to limitations with Samtools mpileup the maximum value is 93. Any value higher than 93 will be set to 93.
* **--max_depth** _(100000)_ - Maximum depth to report at each heterozygous variant when calling samtools mpileup. This value must be greater than the most covered site in your data. Any sites that are covered by more reads than this depth could have false positive signals of allelic imbalance because samtools does not randomly sample reads.
* **--o** - Output file name.
