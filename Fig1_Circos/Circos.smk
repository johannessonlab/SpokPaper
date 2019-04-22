# -*- snakemake -*-

### Circos plot of Podospora spp
#############################################################################
#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/04/16
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
anserina = config["anserina"]
pauciseta = config["pauciseta"]
comata = config["comata"]
# -------------------------------------------------

rule all:
	input:
		"circos/circos.png",
		"circos/circos_inv.png"

rule catgenomes:
	input:
		anserina = anserina,
		pauciseta = pauciseta,
		comata = comata
	output:
		"data/all.fa"
	shell:
		"cat {input} > {output}"

rule allvsall:
	input:
		query = "data/all.fa",
		reference = "data/all.fa",
	output:
		delta = "mummer/all.delta",
		deltafilter = "mummer/all.filter",
		coords = "mummer/all.coords",
		coordsfilter = "mummer/all.filter.coords",	
	params:
		time = "1:00:00",
		threads = 8,
	shell:
		"""
		echo
		echo "MUMmer alignment ..."
		nucmer -b 2000 -c 2000 --maxmatch -p mummer/all {input.reference} {input.query} -t {params.threads}

		# Filter the delta
		delta-filter -q {output.delta} > {output.deltafilter}

		# To view a summary of all the alignments produced by NUCmer
		echo "Running show-coords"
		# For Ribbon http://genomeribbon.com/
		echo "...for Ribbon"
		show-coords -r -lTH {output.delta} > {output.coords}
		show-coords -r -lTH {output.deltafilter} > {output.coordsfilter}

		"""
		# --mum  Use anchor matches that are unique in both the reference and query
		# --mumreference  Use anchor matches that are unique in in the reference
        #           but not necessarily unique in the query (default behavior)
        # -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
        # -b|breaklen     Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)


rule makelinks:
	""" Prepare a links with the MUMmer alignments """
	input:
		"mummer/all.coords",
		# "mummer/anspauci.filter.coords", 
		# "mummer/anscomata.filter.coords",
		# "mummer/paucomata.filter.coords",
	output:
		allsites = "links/mummer.txt",
		clean = "links/mummer.clean.txt",
	shell:
		"""
		## Prepare a links file assigning the colors based on the chromosome names
		cat {input} | awk '{{print $10,$1,$2,$11,$3,$4,"color="}}' | awk '{{a=$1; b=$4; c=gensub(/(.+)hromosome_([0-9])([\.0-9]*)/, "chr\\\\2", "g", a); d=gensub(/(.+)hromosome_([0-9])([\.0-9]*)/, "chr\\\\2", "g", b); if( c==d ) {{print $0 c "_a5"}} else {{print $0 c}}}}' > {output.allsites}

		## Remove self alignments
		# Also remove the one region of Ns in chr7 of PODCO that aligns with another N-region in chr6
		cat {output.allsites} | awk '{{if ($1 $2 $3 != $4 $5 $6) {{print}} }}' | grep -v 'comataT_Chromosome_7 3436438 3479306' > {output.clean}
		"""

rule circos:
	""" Run Circos to plot the alignment """
	input:
		"links/mummer.clean.txt",
	output:
		"circos/circos.png"
	shell:
		"cd circos; circos"

rule circosinv:
	""" Run Circos to plot the alignment """
	input:
		"links/mummer.clean.txt",
	output:
		"circos/circos_inv.png"
	shell:
		"cd circos; circos --conf etc/circos_inv.conf -outputfile circos_inv.png"

