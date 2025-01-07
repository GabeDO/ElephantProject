Structural variant calling snakemake pipeline developed by Dr Gabe O'Reilly at UNCC.

Requires:
snakemake
samtools
bwa
blast


An example directory is provided in this repository.
Setting in 'Example/Config/config.yaml' and 'Example/SnakeMake.sbatch' can be tweeked to suit your avalable HPC resources.

Once the 'setting.smk' file is filled out, you can run "sbatch SnakeMake.sbatch" in the example directory and wait for the pipline to run. This may take some time.

The setting.smk file should be filled out as follows:
##########################################################################################################################################################################
###Leave these values as is
SV = ['TandemDups','Rearrangements']
READ = [1,2]
readtype = ['pair','mate']

##Fill in these
SAMPLES = [] #a list of sample names
ANCESTRAL_SAMPLE_NAME = "" #if you have an ancestral state in your samples (which can be used to mark mutations for polarization) you can put the sample name here. Otherwise, leave blank

CHROM_LIST = [] #a list of chromosome names
CHROM_LIST_NOSEX = []#a list of chromosome names excluding sex chromosomes

REF= "" #Path to your reference

GTFfile =  "" #Path to GTF file - the GTF file must its 9th coloum be an info coloum, with the gene ID being the 3rd elepment when split by a semi colon. e.g.: "X; Y ;GENE_ID"
SampleLocation = "" #Path to where the pipeline will create BAM files and other intermediate files
Source_SampleLocation = "" #Path to your fasta samples. Each sample must be in its own folder named with its sample ID. This can be the same location as 'SampleLocation' if you dont mind everything in the same foilder - although I would not reccomend it. 

scripts = "" #path the the pipelines 'scripts' folder
BlastDB = "../BlastDB/"
##########################################################################################################################################################################

