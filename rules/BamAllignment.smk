rule bwa_aln:
	input:
		ref=REF,
		r= Source_SampleLocation+"/{sample}/{sample}_{read}.fastq.gz"
	output:
		Sai = SampleLocation+"/{sample}/{sample}_{read}.sai"
	priority:
		10
	threads:
		8
	shell: 
		"bwa aln -t 16 {input.ref} {input.r} > {output.Sai}"

rule bwa_sampe:
	input: 
		ref=REF,
		Sai1 = SampleLocation+"/{sample}/{sample}_1.sai",
		Sai2 = SampleLocation+"/{sample}/{sample}_2.sai",
		r1 =  Source_SampleLocation+"/{sample}/{sample}_1.fastq.gz",
		r2 = Source_SampleLocation+"/{sample}/{sample}_2.fastq.gz"
	output: 
		protected(SampleLocation+"/{sample}/{sample}_sampe_sorted.bam")
	resources:
		cpus_per_task=8
	priority:
		10
	shell: 
		"""
		bwa sampe {input.ref} {input.Sai1} {input.Sai2} {input.r1} {input.r2}| samtools view -Sb| samtools sort -@ 8 -o {output}
		samtools index -@ 8 {output}
		"""

rule samtools_Coverage:
	input:
		bams = SampleLocation+"/{sample}/{sample}_sampe_sorted.bam"
	output:
		SampleLocation+"/{sample}/{sample}_Coverage.out"
	shell:
		"samtools coverage {input.bams} > {output}"