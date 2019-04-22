# 3SppForNetwork: Estimating relationships of *Podospora anserina* and sisters
**S. Lorena Ament-Velasquez**

**2019.04.22**

-----

A general pipeline to create a large concatenated alignment of orthologs from *Podospora anserina*, *P. pauciseta*, and *P. comata*. The output is meant for SplitsTree.

## Building the environment

I ran the pipeline under a [Conda](https://docs.anaconda.com/) environment. Install it first. If you like, you can start by updating it.

    $ conda update -n base conda

To create the environment arbitrarily named `myenv`:

    $ conda create -n myenv

Now, to install software, activate the environment.
    
    $ conda activate myenv

Install the following packages:

    $ conda install -c bioconda snakemake-minimal=5.4.4
    $ conda install -c bioconda biopython=1.72=py37h04863e7_0
    $ conda install -c bioconda gffutils=0.9=py_1
    $ conda install -c bioconda mafft=7.407=1

Unfortunately Snakemake runs in python3 and OrthoFinder requires python 2. So I run the OrthoFinder rule with it's own conda environment. This is achieved by having a separate configuration file in the directory `envs`:

    $ cat envs/orthofinder.yaml
```yaml
    channels:
      - bioconda
      - defaults
      - conda-forge
    dependencies:
      - orthofinder=2.2.6
```

## The configuration file

The pipeline runs with a configuration file containing the sample's IDs, the path to the assemblies, the reference genomes, the scripts, and the number of desired ortholog groups. You should make sure paths work for you.

The repo includes a copy of each script, but notice that you can find the latest versions in [my personal GitHub](https://github.com/SLAment/Genomics).

The reference genomes for *P. anserina*, or Podan2 ([Espagne et al. 2008](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-5-r77)) and *P. comata*, or PODCO ([Silar et al. 2018](https://link.springer.com/article/10.1007/s00438-018-1497-3)) can be downloaded from JGI and NCBI, respectively. For both, I used a modified version of the assemblies and annotations for my convenience that is included in the repo as well. The sequence is the same, I just changed a bit the names. Also notice that there is now a [Podan3](https://genome.jgi.doe.gov/Podan3/Podan3.home.html) in JGI, but as far as I can tell the assembly sequence is the same. Like, the same, same.

So the file looks like: 

    $ cat 3SppForNetwork_config.yaml
```yaml
### 3SppForNetwork configuration file
# ++++++++++++++++++++++++++++++++++++++++++++++
# Get 1000 single-copy ortholog groups of Podospora anserina and close sisters
# into a matrix for a network analysis
# ++++++++++++++++++++++++++++++++++++++++++++++

# Samples names
SampleIDs: ["PaWa100p", "PaWa21m", "PaWa28m", "PaWa46p", "PaWa53m", "PaWa58m", "PaWa63p", "PaWa87p", "PaYp", "CBS237.71m", "PaTgp"]

# Path to alignments in format "sample.nice.fa"
assembliespath: "/path/to/assemblies"

# Data of P. anserina and P. comata reference genomes
podan2: "references/Podan2_AssemblyScaffoldsmt.fa"
podan2genes: "references/Podan2_AssemblyScaffoldsGenesEd_gene.fas"
podan2gff: "references/Podan2_AssemblyScaffoldsmtGenesEd_gh.gff"
PODCO: "references/PODCO_genomic.fas"
PODCOgff: "references/PODCO_genomic.gff3"

## Scripts (included in the repo)
gff2fasta: "scripts/gffutils2fasta.py"
orthogrs_parser: "scripts/orthogrs_parser.py"
query2hitseq: "scripts/query2hitseq.py"
fastaconcat: "scripts/fastaconcat.py"

# Number of sample orthologs
SAMPLEsize: 1000
```

## Run pipeline locally

Get into the folder with this repo's content, for example:

    $ cd /home/lore/1_SpokPaper/1a_3SppForNetwork

For testing without running the pipeline:

    $ snakemake --snakefile 3SppForNetwork.smk --configfile 3SppForNetwork_config.yaml -pn

Now you can run the pipeline. I like to make screen first, then activate the environment, and finally run the pipeline in the background.

    $ screen -R phylo
    $ conda activate myenv
    $ snakemake --snakefile 3SppForNetwork.smk --configfile 3SppForNetwork_config.yaml -p -j 35 --keep-going --use-conda &> 3SppForNetwork.log &

Notice `-j` stands for the number of threads that you want to give to your pipeline. See [Snakemake](https://snakemake.readthedocs.io/en/stable/) documentation for more information. If the thing crashes, let me know :)

The result is in the `concatenated` folder, with a big fasta file ready to be input into SplitsTree. I am providing also the actual list of orthologs used for the paper (selected randomly) in `filtering/Podan2_1n.txt`. If present in the folder, probably the pipeline will start from there instead of running the whole OrthoFinder business, but if ran from scratch then the list of genes will be different.

I should add, the pipeline is far from efficient or fast, but I have it like this so I can run it with a [checkpoint](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?fbclid=IwAR1v29DPDpWqve6yRlnc5vob2uIsxCfZt-NSjxfTtbaOZa4TFRuuqn8VbEk#data-dependent-conditional-execution), given that the actual output of OrthoFinder cannot be predicted, a key aspect of the SnakeMake syntax. Anyway, it should work, just not super fast :P.
