# PoolSeq: Pool sequencing of spores coming from 4- vs 2-spore asci of Y vs Wa58p

## Building the environment

This pipeline relies on the AssemblyNanopore environment (but you can call it as you wish).

First, you can start by updating conda.

    $ conda update -n base conda

Now, to create the environment.

    $ conda create -n VariantCalling -c bioconda snakemake-minimal=5.3.0 picard=2.18.11 bwa=0.7.17 samtools=1.9 vcftools=0.1.16

The last part of the pipeline plots the result in R, so I'll need:

    # First I'll activate it to install things inside
    $ conda activate VariantCalling
    $ conda install -c r r-base=3.5.1=h1e0a451_2
    $ conda install -c conda-forge r-ggplot2=3.1.0=r351h6115d3f_1000 r-reshape2=1.4.3=r351h9d2a408_2 r-dplyr=0.7.6=r351h9d2a408_1
    $ conda install -c bioconda r-poppr=2.8.1=r351h470a237_0 r-vcfr=1.8.0=r351h9d2a408_0

And save it in a file.

    $ conda env export > PoolSeq_env.yml

**AN IMPORTANT CAVEAT TO THIS:** My installation of GATK3 is LOCAL. Meaning that I didn't install it along with this environment. The reason is that it requires to download a license ... Instead I have installed locally in the cluster and call that within the snakemake pipeline. It's unfortunate, but this will change with GATK4. The reason why I kept it in GATK3 is because I used the `IndelRealigner` function externally for my Pilon polishing in the genome assemblies. Hence, I decided to stick to a single methodology for producing BAM files, including the SNP calling. In reality using GATK4 would have no effect on our conclusions, since we only need coarse results, not accurate variant calling.

## Prepare your configuration file

This pipeline depends on a given configuration yaml file including the samples, the path to the data, and the reference. Below it's an example, but make sure the paths are correct for you!
    
```yaml
    # Illumina reads path:
    Illumina: "path_to_IlluData/"

    # List of samples to analyze
    diploids: ("Pool2spores", "Pool4spores") # Notice the tuple
    haploids: ("PaYp", "PaWa87m") #The parents were PaYm and PaWa87p, but they are mostly isogenic

    # The reference genome
    Podan2file: "extras/Podan2_AssemblyScaffoldsmt.fa"

    # For filtering
    thining: 200

    # My local installation of GATK3
    gatk3: "localpath/GenomeAnalysisTK.jar"

    rscript: "scripts/FigS8_PoolSeqYvsWa87.R"
```

## Run pipeline in slurm server (like Uppmax)

First, to get an idea of how the pipeline looks like we can make a rulegraph:

In Mac, you need to install graphviz to run the following command. For that you can do `brew install graphviz` using Homebrew, for example. It otherwise works well in Ubuntu.

    $ snakemake --snakefile PoolSeq.smk --configfile PoolSeq_config.yml --rulegraph | dot -Tpng > rulegraph.png

![rulegraph](rulegraph.png "rulegraph of PoolSeq.smk")

To run in local computer:

    $ snakemake --snakefile PoolSeq.smk --configfile PoolSeq_config.yml -p

To run the pipeline in a slurm cluster:

    $ screen -R PoolSeq
    # Important to activate environment!!
    $ conda activate VariantCalling
    $ snakemake --snakefile PoolSeq.smk --configfile PoolSeq_config.yml -p --cluster "sbatch -A projectname -p core -n {params.threads} -t {params.time} --mail-user your@mail.com --mail-type=ALL" -j 10 --keep-going &> PoolSeq_snakemake.log &

