configfile: "configs/config.yaml"

"""
A snakemake script for downsampling bams and seeing if we can still find specific variants.

"""

# For each downsample rate specified in the config
for down_sample in config["down_samples"]:
	# For each replciate specified in the config
	for replicate in config["down_sample_replicates"]:

		rule:
			input:
				bam="input/{sample}.bam",
			output:
				"output/bams/{sample}_" + down_sample + "_" + str(replicate) + ".bam"
			params:
				prob = down_sample,
				walltime = "1:00:00",
				ncpus = "1"
			shell:
				#Downsample the bam at a particular rate
				"picard DownsampleSam I={input} O={output} STRATEGY=Chained P={params.prob} RANDOM_SEED=null"
# Index bams
rule index_bams:
	input:
		"output/bams/{sample}.bam"
	output:
		"output/bams/{sample}.bam.bai"
	params:
		walltime = "1:00:00",
		ncpus = "1"
	shell:
		"samtools index {input}"

# Use sambamba to calculate depth
rule calculate_depth:
	input:
		bam ="output/bams/{sample}.bam",
		index="output/bams/{sample}.bam.bai"
	output:
		"output/depth/{sample}.txt"
	params:
		bed_file = config["bed_file"],
		walltime = "1:00:00",
		ncpus = "1"
	shell:
		"sambamba depth region -L {params.bed_file} {input.bam} > {output}"

# Call variants using mutect2
rule call_variants:
	input:
		bam = "output/bams/{sample}_{down_sample}_{replicate}.bam",
		index = "output/bams/{sample}_{down_sample}_{replicate}.bam.bai"
	output:
		vcf = "output/vcfs/{sample}_{down_sample}_{replicate}.vcf.gz",
		vcf_index = "output/vcfs/{sample}_{down_sample}_{replicate}.vcf.gz.tbi",
		bam = "output/mutect_bam/{sample}_{down_sample}_{replicate}_mutect.bam",
	params:
		java_options = config["gatk_java_options"],
		gatk_ref_genome = config["gatk_ref_genome"],
		gatk_interval_bed = config["gatk_interval_bed"],
		gatk_interval_padding = config["gatk_interval_padding"],
		gatk_max_population_af = config["gatk_max_population_af"],
		gatk_germline_resource = config["gatk_germline_resource"],
		gatk_af_of_alleles_not_in_resource = config["gatk_af_of_alleles_not_in_resource"],
		sample_id = lambda wildcards: wildcards.sample.split("_")[-1],
		walltime = "1:00:00",
		ncpus = "1"
	shell:
		"gatk Mutect2 "
		"--reference {params.gatk_ref_genome} "
		"--input {input.bam} "
		"--tumor {params.sample_id} "
		"--genotype-germline-sites true "
		"--genotyping-mode DISCOVERY "
		"--intervals {params.gatk_interval_bed} "
		"--interval-padding {params.gatk_interval_padding} "
		"--max-population-af {params.gatk_max_population_af} "
		"--output-mode EMIT_ALL_SITES "
		"--germline-resource {params.gatk_germline_resource} "
		"--af-of-alleles-not-in-resource {params.gatk_af_of_alleles_not_in_resource} "
		"--output {output.vcf} "
		"--bamout {output.bam} "
		"--verbosity ERROR "
		"--QUIET true "

# Split multiallelics
rule vt_decompose:
	input:
		vcf = "output/vcfs/{sample}_{down_sample}_{replicate}.vcf.gz",
		index = "output/vcfs/{sample}_{down_sample}_{replicate}.vcf.gz.tbi"
	output:
		vcf = "output/vcfs_decompose/{sample}_{down_sample}_{replicate}.vcf.gz",
		index = "output/vcfs_decompose/{sample}_{down_sample}_{replicate}.vcf.gz.tbi",
	params:
		walltime = "1:00:00",
		ncpus = "1"
	shell:
		"vt decompose {input.vcf} -o {output.vcf} && tabix {output.vcf}"

# decompose blocks
rule vt_decompose_block:
	input:
		vcf = "output/vcfs_decompose/{sample}_{down_sample}_{replicate}.vcf.gz",
		index = "output/vcfs_decompose/{sample}_{down_sample}_{replicate}.vcf.gz.tbi",
	output:
		vcf = "output/vcfs_decompose_block/{sample}_{down_sample}_{replicate}.vcf.gz",
		index = "output/vcfs_decompose_block/{sample}_{down_sample}_{replicate}.vcf.gz.tbi"
	params:
		walltime = "1:00:00",
		ncpus = "1"
	shell:
		"vt decompose_blocksub -a {input.vcf} -o {output.vcf} && tabix {output.vcf}"

# Convert to vcfs
rule vcf_to_csv:
	input:
		vcf = "output/vcfs_decompose_block/{sample}_{down_sample}_{replicate}.vcf.gz",
		index = "output/vcfs_decompose_block/{sample}_{down_sample}_{replicate}.vcf.gz.tbi"
	output:
		"output/csvs/{sample}_{down_sample}_{replicate}.csv"
	params:
		walltime = "1:00:00",
		ncpus = "1"
	shell:
		"gatk VariantsToTable -V {input.vcf} -O {output} -SMA -F CHROM -F POS -F REF -F ALT -F DP -F TLOD -GF AD -GF AF  -GF GT -GF F1R2 -GF F2R1"

# Collect everything
rule merge:
	input:
		expand("output/csvs/{sample}_{down_sample}_{replicate}.csv", sample=config["original_bams"], down_sample= config["down_samples"], replicate= config["down_sample_replicates"]),
		expand("output/depth/{sample}_{down_sample}_{replicate}.txt", sample=config["original_bams"], down_sample= config["down_samples"], replicate= config["down_sample_replicates"]) 
	output:
		"output/final.txt"
	params:
		walltime = "1:00:00",
		ncpus = "1"
	shell:
		"echo {input} > {output}"



