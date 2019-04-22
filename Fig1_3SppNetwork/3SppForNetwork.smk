# -*- snakemake -*-

import glob
from Bio import SeqIO

### 3SppForNetwork
#############################################################################
# Get 1000 single-copy ortholog groups of Podospora anserina and close sisters
# into a matrix for a network analysis
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/04/14-19
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
samples = config["SampleIDs"] + ["Podan2", "PODCO"]
assembliespath = config["assembliespath"]
podan2 = config["podan2"]
podan2gff = config["podan2gff"]
podan2genes = config["podan2genes"]
PODCO = config["PODCO"]
PODCOgff = config["PODCOgff"]

# Scripts
gff2fasta = config["gff2fasta"]
orthogrs_parser = config["orthogrs_parser"]
query2hitseq = config["query2hitseq"]
fastaconcat = config["fastaconcat"]

# Constants
SAMPLEsize = config["SAMPLEsize"]
# -------------------------------------------------

# ----------
# Rules not submitted to a job # For a cluster
# localrules: 

# ----------

rule all:
	input:
		"concatenated/allonetoones.fa",

# ------- PREPARE ALL DATA --------
## Make symlink of the genomes to work more easily
rule getgenomes:
	""" Make links to the assemblies """
	input: 
		assembliespath + "/{sample}.nice.fa", # I'm expecting a file called eg. PaWa100p.nice.fa
	output:
		"genomes/{sample}.fa"
	shell:
		"ln -sf {input} {output}"

rule getrefgenomes: # Rename the references for convenience
	""" Make links to the assemblies """
	input: 
		podan2 = podan2,
		PODCO = PODCO,
	output:
		podan2 = "genomes/Podan2.fa",
		PODCO = "genomes/PODCO.fa",
	shell:
		"cat {input.podan2} | sed 's;>;>Podan2_;' > {output.podan2};"
		"cat {input.PODCO} | sed 's;>;>PODCO_;' > {output.PODCO};" # I previously already changed the names of the chromosomes


# ----------------------------------
rule getprotsamples:
	""" Get protein sequences """ 
	input:
		genome = assembliespath + "/{sample}.nice.fa",
		gff = assembliespath + "/{sample}.nice.gff3",
	output:
		prots = "proteins/{sample}.fas",
	params:
		gff2fasta = gff2fasta
	shell:
		"""
		# Get CDS translated
		python {params.gff2fasta} {input.genome} {input.gff} --output proteins/{wildcards.sample} -t CDS -p -j --onlynames
	
		# Put the sample ID in the name of the sequences
		sed -i 's;>\\(.*\\);>\\1_{wildcards.sample};' {output.prots} 
		"""

rule getprotpodan:
	""" Get protein sequences of reference genomes""" 
	input:
		genome = podan2,
		gff = podan2gff,
	output:
		prots = "proteins/Podan2.fas",
	params:
		gff2fasta = gff2fasta
	shell:
		"""
		# Get CDS translated
		python {params.gff2fasta} {input.genome} {input.gff} --output proteins/Podan2 -t CDS -p -j --onlynames
	
		# Put the sample ID in the name of the sequences
		sed -i 's;>\\(.*\\);>\\1_Podan2;' {output.prots} 
		"""

rule getprotcomata:
	""" Get protein sequences of reference genomes""" 
	input:
		genome = PODCO,
		gff = PODCOgff,
	output:
		prots = "proteins/PODCO.fas",
	params:
		gff2fasta = gff2fasta
	shell:
		"""
		# Get CDS translated
		python {params.gff2fasta} {input.genome} {input.gff} --output proteins/PODCO -t CDS -p -j --onlynames
	
		# Put the sample ID in the name of the sequences
		sed -i 's;>\\(.*\\);>\\1_PODCO;' {output.prots} 
		"""
# -----------------------------------------

rule orthofinder:
	""" Run OrthoFinder """
	# https://github.com/davidemms/OrthoFinder/blob/master/OrthoFinder-manual.pdf
	input:
		expand("proteins/{sample}.fas", sample = samples)
	output:
		"orthofinder/Orthogroups.csv"
	params:
		threads = 30,
	conda: 
		"envs/orthofinder.yaml"
	shell:
		"""
		# Run Orthofinder
		orthofinder -f proteins -t {params.threads}

		# Move the new folder and change name
		mv proteins/Results_*/* orthofinder/
		"""

rule parseorthogroups1n:
	""" Parse the output of orthogroups """
	input:
		"orthofinder/Orthogroups.csv"
	output:
		"filtering/Podan2_1n.txt",
	params:
		orthogrs_parser = orthogrs_parser,
		threads = 1,
		SAMPLEsize = SAMPLEsize
	shell:
		"""
		# Filter the Orthogroups.csv for groups of one-to-one orthologs present in all samples
		{params.orthogrs_parser} {input} -n1 -b -o filtering -s {params.SAMPLEsize}
		
		# Clean the output a bit
		sed -i 's;_Podan2;;g' filtering/Podan2_1n.txt
		"""

checkpoint query2hitseqfas:
	""" Get fasta files of each ortholog group containing all samples""" 
	input:
		genomes = expand("genomes/{sample}.fa", sample = samples),
		orthologs = "filtering/Podan2_1n.txt",
	output:
		directory("fastahits")
	params: 
		query2hitseq = query2hitseq,
		refgenes = podan2genes,
	run:
		shell("mkdir -p fastahits")
	
		for genome in input.genomes:

			# Make a dir for each sample to put blast results
			sample = genome.rstrip(".fa").split("/")[1]
			shell("mkdir -p fastahits/" + sample)

			# Read the ortholog groups file
			tabopen = open(input.orthologs, 'r')
			tabs = [line.rstrip("\n") for line in tabopen] 			# Read tab file into a list

			# Get fasta
			count = 1
			for ortho in tabs:
				number = "{0:04d}".format(count)
				cmd1 = "python %s %s %s -i %s --tophit --extrabp %d --temp %s >> fastahits/orthologs%s.fas" % (params.query2hitseq, params.refgenes, genome, ortho, 0, "fastahits/" + sample + "/", number)
				shell(cmd1)

				count += 1

rule mafft:
	""" Align ortholog groups with MAFFT """
	input:
		ortholog = "fastahits/orthologs{i}.fas",
	output:
		"alignments/orthologs{i}.fas",
	params:
		threads = 6,
	shell:
		""" 
		mafft --thread {params.threads} --threadit 0 --adjustdirection --anysymbol --maxiterate 1000 --retree 1 --localpair {input.ortholog} > {output}
		"""

def mafftoutput(wildcards):
	""" Make the names of the expected output """
	checkpoint_output = checkpoints.query2hitseqfas.get(**wildcards)#.output[0] # I actually don't need anything from this, just to confirm it's a checkpoint

	# Get the names of all the final alignments
	file = glob.glob("filtering/Podan2_1n.txt")[0]
	num_lines = sum(1 for line in open(file)) # How many lines in the file?
	finalnames = []
	for n in range(1, num_lines + 1):
		number = "{0:04d}".format(n)
		fname = "alignments/orthologs%s.fas" % (number)
		finalnames.append(fname)
		
	return finalnames	

rule concatenation:
	""" Concatenate all the one-to-one ortholog alignments"""
	input:
		mafftoutput,
	output:
		"concatenated/allonetoones.fa"
	params:
		fastaconcat = fastaconcat,
		samples = samples,
	run:
		print(input)
		totalnumsamples = len(params.samples)
		completealignments = [] # In case some alignments didn't get all taxa when BLASTing using nucleotides (as opposed to the protein input of OrthoFinder)
		
		for alignment in input:
			# Get the name back of the alignment
			# alignment = "alignments/" + log.rstrip(".log").split("/")[2]
			
			# Read it
			thisfasta = [seq_record for seq_record in SeqIO.parse(alignment, "fasta")]
			if len(thisfasta) == totalnumsamples:
				completealignments.append(alignment)

		print("Number of alignments with all samples: %d" % len(completealignments))
		alignstring = ' '.join(completealignments)
		
		# Concatenate and rename them
		cmd = params.fastaconcat + " " + alignstring + " | sed 's;_R_;;' | sed 's;\\(>[a-zA-Z0-9\\.]*\\)_\\([a-zA-Z0-9._-]*\\);\\1;' > " + output[0]
		shell(cmd)

