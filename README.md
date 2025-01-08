Structural variant calling snakemake pipeline developed by Dr Gabe O'Reilly at UNCC.<br/> 
Funded in part by NIH NIGMS MIRA R35 GM133376 to RL Rogers.<br/> 



Requires:<br/> 
snakemake<br/> 
samtools<br/> 
bwa<br/> 
blast<br/> 


An example directory is provided in this repository.<br/> 
Setting in 'Example/Config/config.yaml' and 'Example/SnakeMake.sbatch' can be tweeked to suit your avalable HPC resources.

Once the 'setting.smk' file is filled out, you can run "sbatch SnakeMake.sbatch" in the example directory and wait for the pipline to run. This may take some time.

The setting.smk file should be filled out as follows:<br/> 

###Leave these values as is<br/> 
SV = ['TandemDups','Rearrangements']<br/> 
READ = [1,2]<br/> 
readtype = ['pair','mate']
<br/> 


##Fill in these<br/> 
SAMPLES = [] #a list of sample names<br/> 
ANCESTRAL_SAMPLE_NAME = "" #if you have an ancestral state in your samples (which can be used to mark mutations for polarization) you can put the sample name here. Otherwise, leave blank<br/> 

CHROM_LIST = [] #a list of chromosome names<br/> 
CHROM_LIST_NOSEX = []#a list of chromosome names excluding sex chromosomes<br/> 

REF= "" #Path to your reference<br/> 

GTFfile =  "" #Path to GTF file - the GTF file must its 9th coloum be an info coloum, with the gene ID being the 3rd elepment when split by a semi colon. e.g.: "X; Y ;GENE_ID"<br/> 
SampleLocation = "" #Path to where the pipeline will create BAM files and other intermediate files<br/> 
Source_SampleLocation = "" #Path to your fasta samples. Each sample must be in its own folder named with its sample ID. This can be the same location as 'SampleLocation' if you dont mind everything in the same foilder - although I would not reccomend it. <br/> 

scripts = "" #path the the pipelines 'scripts' folder<br/> 
BlastDB = "../BlastDB/"<br/> 

