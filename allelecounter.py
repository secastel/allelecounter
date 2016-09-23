import tempfile;
import argparse;
import subprocess;
import calendar;
import time;
import random;
import os;
import pysam;
import gzip;

def main():
	#Arguments passed 
	parser = argparse.ArgumentParser()
	parser.add_argument("--vcf", required=True, help="Genotype VCF for sample (gzipped with index).")
	parser.add_argument("--sample", required=True, help="Sample name")
	parser.add_argument("--bam", required=True, help="BAM file containing RNA-seq reads (duplicates marked with index)")
	parser.add_argument("--ref", required=True, help="Refrence genome file (fasta with index)")
	parser.add_argument("--min_cov", required=True, type=int, help="Minimum total coverage to report read counts")
	parser.add_argument("--min_baseq", required=True, type=int, help="Minimum base quality to count reads")
	parser.add_argument("--min_mapq", required=True, type=int, help="Minimum map quality to count reads")
	parser.add_argument("--o", required=True, help="Output file")
	args = parser.parse_args()
	
	# deal with STAR alignments having mapq of 255
	# the highest mapq that can be encoded by samtools mpileup is 93
	# so if min_mapq > 93 set it to 93
	if args.min_mapq > 93: args.min_mapq = 93;
	
	start_timestamp = calendar.timegm(time.gmtime());
	
	version = "0.5";
	print("");
	print("##################################################")
	print("          Welcome to allelecounter v%s"%(version));
	print("  Author: Stephane Castel (scastel@nygenome.org)")
	print("##################################################");
	print("");
	
	print("NOTE: chromosome names must be consistent between VCF and BAM for proper function.\n");
	
	#1 perform the pileup using samtools
	print("Generating pileup...");
	
	# instantiate temp file to dump to
	pileup_out = tempfile.NamedTemporaryFile(delete=False);
	
	#A unzip the VCF and convert to BED use samtools to produce a read pileup
	if args.vcf.endswith(".gz"):
		subprocess.check_output("gunzip -c "+args.vcf+" | awk 'BEGIN { OFS=\"\t\"; FS=\"\t\"; } { if (index($0, \"#\") == 0) { print($1,$2-1,$2,$3,$6,$4,$5,$7,$8,$9); } }' | samtools mpileup -I -B -q 0 -Q 0 -s -l - -f "+args.ref+" "+args.bam+" > "+pileup_out.name, shell=True);
	else:
		print("FATAL ERROR: input VCF must be gzipped and indexed.")
		quit();
	
	# 2 process the mpileup result
	# A load the VCF
	#vcf_reader = vcf.Reader(filename=args.vcf);
	vcf_reader = pysam.Tabixfile(args.vcf,"r");
	vcf_map = sample_column_map(args.vcf);
	
	# B go through the pileup result line by line
	out_stream = open(args.o, "w");
	
	print("Processing pileup...");
	out_stream.write("contig	position	variantID	refAllele	altAllele	refCount	altCount	totalCount	lowMAPQDepth	lowBaseQDepth	rawDepth	otherBases\n");
	
	stream_in = open(pileup_out.name, "r");
	
	for line in stream_in:
		cols = line.replace("\n","").split("\t");
		#chr	pos	REF	count	reads
		chr = cols[0];
		pos = int(cols[1]);
		ref = "";
		alt = "";
		rsid = ".";
		# first retrieve the VCF record for this SNP
		records = vcf_reader.fetch(chr,pos-1,pos);
		
		for record in records:
			vcf_cols = record.rstrip().split("\t");
			if int(vcf_cols[1]) == pos:
				alleles = [vcf_cols[3]] + vcf_cols[4].split(",");
				gt_index = vcf_cols[8].split(":").index("GT");
				gt = vcf_cols[vcf_map[args.sample]].split(":")[gt_index].replace("|","").replace("/","");
				# only want heterozygous sites
				if len(set(gt)) > 1:
					allele_indices = map(int, list(gt));
					# set ref to minimum allele index
					ref = alleles[min(allele_indices)];
					# and alt to max
					alt = alleles[max(allele_indices)];
					
					rsid = vcf_cols[2];
					
		# if the site is biallelic and we have ref / alt genotype data go ahead and do the read counts
		if ref != "" and alt != "":
			reads = list(cols[4]);
			reads = [x for x in reads if (x != "$" and x != "^" and x != "~" and x != "\"" and x != "!")]
			block = 0;
			out_reads = [];
			for read in reads:
				if block == -1:
					if ord(read) >= 48 and ord(read) <= 57:
						block_str += read;
					else:
						block = int(block_str) - 1;
				elif read == "+" or read == "-":
					block = -1;
					block_str = "";
				elif block > 0:
					block -= 1;
				elif block == 0:
					out_reads.append(read);
			reads = out_reads;
			baseqs = list(cols[5]);
			mapqs = list(cols[6]);
			
			low_baseq = 0;
			low_mapq = 0;
			ref_count = 0;
			alt_count = 0;
			other_count = 0;
			raw_depth = len([x for x in reads if (x != ">" and x != "<")]);
			
			for read,baseq,mapq in zip(reads,baseqs,mapqs):
				if read != "<" and read != ">":
					basequal = ord(baseq)-33;
					mapqual = ord(mapq)-33;
					if basequal >= args.min_baseq and mapqual >= args.min_mapq:
						# count this base
						if read == "." or read == ",":
							ref_count += 1;
						elif read == alt.upper() or read == alt.lower():
							alt_count += 1;
						else:
							other_count += 1;
					if basequal < args.min_baseq:
						low_baseq += 1;
					if mapqual < args.min_mapq:
						low_mapq += 1;
			totalCount = ref_count+alt_count;
			
			if totalCount >= args.min_cov:
				out_stream.write("\t".join([cols[0],cols[1],rsid,ref,alt,str(ref_count),str(alt_count),str(totalCount),str(low_mapq),str(low_baseq),str(raw_depth),str(other_count),"\n"]));
	
	out_stream.close();
	stream_in.close();
	
	os.remove(pileup_out.name);
	
	stop_timestamp = calendar.timegm(time.gmtime());
	print("Total time to complete: %d seconds"%(stop_timestamp-start_timestamp));

def sample_column_map(path, start_col=9, line_key="#CHR"):
	stream_in = gzip.open(path, "r");
	
	out_map = {};
	for line in stream_in:
		if line_key in line:
			line = line.rstrip().split("\t");
			for i in range(start_col,len(line)):
				out_map[line[i]] = i;
		
			break;
	
	stream_in.close();
	
	return(out_map);

if __name__ == "__main__":
	main();
