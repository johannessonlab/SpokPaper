### PoolSeq: Pool sequencing of spores coming from 4- vs 2-spore asci of Y vs Wa87p
#############################################################################
# To confirm that Spok2 is responsible for killing in crosses between Psk-1
# and Psk-5 a pooled sequencing approach was employed. A cross was conducted
# between Wa87 and Y, and spores from 2-spored (spore killing) and 4-spored
# asci (heterozygous for killers) were collected. We obtained 21 progeny from
# 2-spored asci and 63 progeny from 4-spored asci. Progeny from each ascus
# type were grown on the same plate and extracted for sequencing with Illumina
# Hi-seq X.
#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2018/08/14
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1


# -------------------------------------------------
# Illumina reads path:
Illumina = config["Illumina"]

# List of samples to analyze
diploids = config["diploids"]
haploids = config["haploids"] #The parents were probably PaYm and PaWa87p, tho, but they are mostly isogenic

# The reference genome
Podan2file = config["Podan2file"]

# For filtering
thining = config["thining"]

rscript = config["rscript"]

gatk3 = config["gatk3"]

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: rawdata_illumina, referencegenome, indexbwa, snpsvcf, snpsvcfnomiss, snpsvcfthin
# ----------

rule all:
	input:
		# The final figure
		"results/Figure_S8_PoolSeq.pdf",

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
		JavaMem = int(8 * 6.8), # A Rackham node contains 128 GB of RAM and 20 compute cores (each core gets at most 6.8 GB).
		gatk3 = gatk3
	version: "1"
	shell:
		"""
		ln -fs {params.gatk3} .

		# Create a target list of intervals to be realigned (Identify what regions need to be realigned)
		java -Xmx{params.JavaMem}G -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -I {input.mdoutput} -R {input.genome} -nt {params.threads} -o {output.intervals}

		# Perform realignment of the target intervals, 20 min?
		java -Xmx{params.JavaMem}G -jar GenomeAnalysisTK.jar -T IndelRealigner -I {input.mdoutput} -R {input.genome} --targetIntervals {output.intervals} -o {output.reoutput}
		"""

rule gvcf_gatk3_1n:
	""" Produce a GVCF file from BAM - haploid """
	input:
		nicebam = "mapping/{sample}/{sample}-to-Podan2.sorted.debup.realign.bam",
		ref = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
	output:
		gvcf = "gvcfs/{sample}-p1.g.vcf"
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

rule gvcf_gatk3_2n:
	""" Produce a GVCF file from BAM - diploid """
	input:
		nicebam = "mapping/{sample}/{sample}-to-Podan2.sorted.debup.realign.bam",
		ref = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
	output:
		gvcf = "gvcfs/{sample}-p2.g.vcf"
	params:
		time = "10:00:00",
		threads = 8,
		JavaMem = int(8 * 6.8),
		ploidy = 2
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
		expand("gvcfs/{sample}-p1.g.vcf", sample = haploids),
		expand("gvcfs/{sample}-p2.g.vcf", sample = diploids),
		ref = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
	output:
		rawvcf = "results/PoolSeq.vcf"
	params:
		time = "2:00:00",
		threads = 8,
		JavaMem = int(8 * 6.8),
		ploidy = 2
	run:
		# Create a string in the format --variant path/to/gvcf/sample1 --variant path/to/gvcf/sample2 etc...
		variantlist = ""
		for sample in input[:-1]: # Ignore the reference in the end
			variantlist += "--variant " + sample + " "
		
		gatkcommand = "java -Xmx{}G -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R {} {} -o {} -nt {}".format(params.JavaMem, input.ref, variantlist, output.rawvcf, params.threads)

		shell(gatkcommand)

rule snpsvcf:
	""" Filter resulting vcf file """
	input:
		rawvcf = "results/PoolSeq.vcf"
	output:
		snpsvcf = "results/PoolSeq_snps.vcf",
	shell:
		"vcftools --vcf {input.rawvcf} --remove-indels --recode --recode-INFO-all --stdout > {output.snpsvcf}"


rule snpsvcfnomiss:
	""" Filter resulting vcf file """
	input:
		vcf = "results/PoolSeq_snps.vcf"
	output:
		filteredvcf = "results/PoolSeq_snps.miss1.vcf",
	shell:
		"vcftools --vcf {input.vcf} --max-missing 1 --recode --recode-INFO-all --stdout > {output.filteredvcf}"

rule snpsvcfthin:
	""" Filter resulting vcf file """
	input:
		vcf = "results/PoolSeq_snps.miss1.vcf"
	output:
		filteredvcf = expand("results/PoolSeq_snps.miss1.thin{thin}.vcf", thin = thining),
	params:
		thin = thining
	shell:
		"vcftools --vcf {input.vcf} --thin {params.thin} --recode --recode-INFO-all --stdout > {output.filteredvcf}"

# ------- Plotting the final figure for the paper -------

rule plotinR:
	""" Plot the SNPs for all the samples """
	input: 
		expand("results/PoolSeq_snps.miss1.thin{thin}.vcf", thin = thining)
	output:
		"results/Figure_S8_PoolSeq.pdf"
	script:
		rscript

