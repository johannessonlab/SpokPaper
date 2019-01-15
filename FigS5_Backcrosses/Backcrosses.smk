### Backcrosses: Analysis of the SNP distribution of 5x backcrosses of killers to strain S
#############################################################################
# We sequenced the 5xS backcrosses of van der Gaag et al (2000) using Illumina
# HiSeq 2500, as well as our own parental S strains with Illumina HiSeq X. The
# idea is to narrow down the genomic areas that could be responsible for the
# Psk killer types.

# - Pa170m: Wa53 Backcross 5x to S-  in the paper as Psk1xS5-
# - Pa130p: Wa28 Backcross 5x to S+  in the paper as Psk2xS5+
# - Pa200p: Y Backcross 5x to S+     in the paper as Psk5xS5+
# - Pa180p: Wa58 Backcross 5x to S+  in the paper as Psk7xS5+

# To see the whole pipeline do:
# $ snakemake --snakefile Backcrosses.smk --configfile Backcrosses_config.yml --rulegraph | dot -Tpng > rulegraph.png

#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2018/08/14
# ---------------
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 2


# -------------------------------------------------
samples = config["SampleIDs"]

# Illumina reads path:
Illumina = config["Illumina"]

# The reference genome
Podan2file = config["Podan2file"]

TElib = config["TElib"]

rscript = config["rscript"]

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: rawdata_illumina, referencegenome, indexbwa, snpsvcf, snpsvcfnomiss, snpsvcfmaf, bedtoolsTEs, plotinR #snpsvcfthin, indelvcf
# ----------

rule all:
	input:
		# Figure for the paper
		"results/Figure_S6_Backcrosses.pdf"


# ------- PREPARE ALL DATA --------
rule rawdata_illumina:
	""" Prepare a folder with the Illumina data files in decent naming """
	output:
		read1 = "data/Illumina/{sample}_postQC.1.fq.gz",
		read2 = "data/Illumina/{sample}_postQC.2.fq.gz"	
	version: "1"
	params:
		illuminapath = Illumina
	shell:
		"""
		ln -sf {params.illuminapath}/{wildcards.sample}_postQC.1.fq.gz data/Illumina
		ln -sf {params.illuminapath}/{wildcards.sample}_postQC.2.fq.gz data/Illumina
		"""

rule referencegenome:
	""" Prepare a copy for the Podan2 reference genome """
	output:
		"data/Podan2/Podan2_AssemblyScaffoldsmt.fa"
	params:
		refgenome = Podan2file
	shell:
		"""
		ln -s {params.refgenome} data/Podan2
		"""
# ---------------------------------

rule indexbwa:
	""" Index genome with BWA """
	input:
		genome = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa"
	output:
		index = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa.bwt"
	version: "1"
	shell:
		"""
		bwa index {input.genome}
		"""

rule bwa_mem:
	""" Map Illumina reads with BWA """
	input:
		genome = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
		index = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa.bwt",
		read1 = "data/Illumina/{sample}_postQC.1.fq.gz",
		read2 = "data/Illumina/{sample}_postQC.2.fq.gz",
	output:
		bwaoutput = temp("mapping/{sample}/{sample}-to-Podan2.bam.sorted"),
		log = "logs/bwa_mem/{sample}.log"
	params:
		time = "3:30:00",
		threads = 10, 
		refbase = "Podan2",
		rg = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina",
	version: "1"
	shell:
		"""
		(bwa mem {input.genome} {input.read1} {input.read2} -t {params.threads} -R '{params.rg}' -M | samtools view -Su - | samtools sort -l 5 -O bam -T {wildcards.sample}'-to-'{params.refbase} -@ {params.threads} > {output.bwaoutput}) 2> {output.log}
		# -l 5 following Doug
		"""

rule markduplicates:
	""" Mark duplicates in BAM """
	# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
	input:
		bwaoutput = "mapping/{sample}/{sample}-to-Podan2.bam.sorted"
	output:
		mdoutput = temp("mapping/{sample}/{sample}-to-Podan2.sorted.debup.bam"),
		mdmetrics = "mapping/{sample}/{sample}-to-Podan2.sorted.metrics.txt"
	params:
		time = "1:30:00",
		threads = 10,
	version: "1"
	shell:
		"""
		# Using normal Picard
		picard MarkDuplicates I={input.bwaoutput} O={output.mdoutput} M={output.mdmetrics} ASSUME_SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR="temp"

		"""	
		# # VALIDATION_STRINGENCY=ValidationStringency
		# #                               Validation stringency for all SAM files read by this program.  Setting stringency to
		# #                               SILENT can improve performance when processing a BAM file in which variable-length data
		# #                               (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This
		# #                               option can be set to 'null' to clear the default value. Possible values: STRICT,
		# #                               LENIENT, SILENT
		# # CREATE_INDEX=Boolean          Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value:
		# #                               false. This option can be set to 'null' to clear the default value. Possible values:
		# #                               true, false
		# # TMP_DIR (File)  Default value: null. This option may be specified 0 or more times.

rule indexsanddict:
	""" Index reference and dictionary for GATK """ 
	input:
		genome = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
	output:
		indexsamtools = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa.fai",
		diction = "data/Podan2/Podan2_AssemblyScaffoldsmt.dict"
	params:
		time = "1:00:00",
		threads = 1,
	version: "1.1"
	shell:
		"""
		# Make an reference index
		samtools faidx {input.genome}

		# Make a reference dictionary
		picard CreateSequenceDictionary R={input.genome} O={output.diction}
		"""	

rule realigngatk3:
	""" Realign indels in BAM """
	# https://software.broadinstitute.org/gatk/documentation/article.php?id=7156
	input:
		mdoutput = "mapping/{sample}/{sample}-to-Podan2.sorted.debup.bam",
		genome = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
		indexsamtools = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa.fai",
		diction = "data/Podan2/Podan2_AssemblyScaffoldsmt.dict"
	output:
		intervals = "mapping/{sample}/{sample}-to-Podan2.sorted.debup.intervals",
		reoutput = "mapping/{sample}/{sample}-to-Podan2.sorted.debup.realign.bam",
	params:
		time = "3:30:00",
		threads = 8,
		JavaMem = int(8 * 6.8) # A Rackham node contains 128 GB of RAM and 20 compute cores (each core gets at most 6.8 GB).
	version: "1"
	shell:
		"""
		ln -fs /sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar .

		# Create a target list of intervals to be realigned (Identify what regions need to be realigned)
		java -Xmx{params.JavaMem}G -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -I {input.mdoutput} -R {input.genome} -nt {params.threads} -o {output.intervals}

		# Perform realignment of the target intervals, 20 min?
		java -Xmx{params.JavaMem}G -jar GenomeAnalysisTK.jar -T IndelRealigner -I {input.mdoutput} -R {input.genome} --targetIntervals {output.intervals} -o {output.reoutput}
		"""

# ------- General report of BAM files -------

rule qualimap:
	# TODO: STILL NOT WORKING, it produces an error
	input:
		reoutput = "mapping/{sample}/{sample}-to-Podan2.sorted.debup.realign.bam"
	output:
		"mapping/{sample}/{sample}_bamqc/genome_results.txt"
	params:
		time = "2:00:00",
		threads = 4,
	shell:
		"""
		unset DISPLAY # Necessary to work under jobs in Uppmax
		qualimap bamqc -bam {input.reoutput} -outdir mapping/{wildcards.sample}/{wildcards.sample}_bamqc -nt {params.threads} # -outformat pdf

		"""

# ------- SNP calling -------

rule gvcf_gatk3:
	""" Produce a GVCF file from BAM - haploid """
	input:
		nicebam = "mapping/{sample}/{sample}-to-Podan2.sorted.debup.realign.bam",
		ref = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
	output:
		gvcf = "gvcfs/{sample}.g.vcf"
	params:
		time = "10:00:00",
		threads = 8,
		JavaMem = int(8 * 6.8),
		ploidy = 1
	shell:
		"""
		java -Xmx{params.JavaMem}G -jar GenomeAnalysisTK.jar \
		-nct {params.threads} \
		-ploidy {params.ploidy} \
		-I {input.nicebam} -R {input.ref} \
		-T HaplotypeCaller\
		-o {output.gvcf} \
		--bamWriterType CALLED_HAPLOTYPES \
		-newQual \
		-stand_call_conf 20.0 \
		-gt_mode DISCOVERY \
		--emitRefConfidence GVCF
		"""
		# For changes in GATK 3.7.0 see https://software.broadinstitute.org/gatk/blog?id=8692
		# Default of -stand_call_conf is 10 in 3.7.0, which is very low to keep as many variants as possible before filtering
		# -newQual has the new method for calculating QUAL that will be default anyway in 4.0.0

rule genotypeGVCF:
	""" Make a vcf using GenotypeGVCFs """
	input:
		expand("gvcfs/{sample}.g.vcf", sample = samples),
		ref = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
	output:
		rawvcf = "results/Backcrosses.vcf"
	params:
		time = "2:00:00",
		threads = 1,
		JavaMem = int(8 * 6.8),
		ploidy = 1
	run:
		# Create a string in the format --variant path/to/gvcf/sample1 --variant path/to/gvcf/sample2 etc...
		variantlist = ""
		for sample in input[:-1]: # Ignore the reference in the end
			variantlist += "--variant " + sample + " "
		
		gatkcommand = "java -Xmx{}G -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R {} {} -o {} -nt {}".format(params.JavaMem, input.ref, variantlist, output.rawvcf, params.threads)

		shell(gatkcommand)

# ------- RepeatMasking -------

rule repeatmasker:
	""" Use RepeatMasker to find regions that should be filtered out """
	input:
		"data/Podan2/Podan2_AssemblyScaffoldsmt.fa"
	output:
		"RepeatMasker/Podan2_AssemblyScaffoldsmt.fa.out.gff"
	params:
		time = "2:00:00",
		threads = 16,
		TElib = TElib
	shell:
		""" 
		RepeatMasker -pa {params.threads} -a -xsmall -gccalc -gff -excln -lib {params.TElib} -dir RepeatMasker {input}
		"""

# ------- Filtering -------

rule snpsvcf:
	""" Filter resulting vcf file """
	input:
		rawvcf = "results/Backcrosses.vcf"
	output:
		snpsvcf = "results/Backcrosses_snps.vcf",
	shell:
		"vcftools --vcf {input.rawvcf} --remove-indels --recode --recode-INFO-all --stdout | grep -v 'PaMt_NC_001329.3' > {output.snpsvcf}" # Notice I filtered out the mitochondria too

rule snpsvcfnomiss:
	""" Filter resulting vcf file """
	input:
		vcf = "results/Backcrosses_snps.vcf"
	output:
		filteredvcf = "results/Backcrosses_snps.miss1.vcf",
	shell:
		"vcftools --vcf {input.vcf} --max-missing 1 --recode --recode-INFO-all --stdout > {output.filteredvcf}"

rule snpsvcfmaf:
	""" Remove sites that are only different in the reference genome """
	input:
		vcf = "results/Backcrosses_snps.miss1.vcf",
	output:
		filteredvcf = "results/Backcrosses_snps.miss1.maf.vcf",
	shell:
		"bcftools view -q 0.09:minor {input.vcf} > {output.filteredvcf}" # There are 10 samples (so > 0.09)
		# notice however this would remove the centromere too, because they will all be different from S.

rule bedtoolsTEs:
	""" Filter vcf with the repeats from RepeatMasker """
	input:
		vcf = "results/Backcrosses_snps.miss1.maf.vcf",
		# vcf = expand("results/Backcrosses_snps.miss1.thin{thin}.vcf", thin = thining),
		gfffile = "RepeatMasker/Podan2_AssemblyScaffoldsmt.fa.out.gff"
	output:
		filteredvcf = "results/Backcrosses_snps.miss1.maf.NoTEs.vcf",
		# filteredvcf = expand("results/Backcrosses_snps.miss1.thin{thin}.NoTEs.vcf", thin = thining),
	shell:
		"""
		# Get the header
		bcftools view -h {input.vcf} > {output.filteredvcf}
		
		# Filter out the repeats
		bedtools intersect -a {input.vcf} -b {input.gfffile} -v >> {output.filteredvcf}
		"""


# ------- Plotting the final figure for the paper -------

rule plotinR:
	""" Plot the SNPs for all the samples """
	input: 
		"results/Backcrosses_snps.miss1.maf.NoTEs.vcf"
	output:
		"results/Figure_S6_Backcrosses.pdf"
	script:
		rscript

