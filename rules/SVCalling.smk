rule samtools_filter:
	input:
		Bam = SampleLocation+"/{sample}/{sample}_sampe_sorted.bam",
	params:
		Elephant = "{sample}"
	output:
		TandemDups = SampleLocation +"/{sample}/{sample}_TandemDupsReads.out",
		Rearrangements = SampleLocation +"/{sample}/{sample}_RearrangementsReads.out"
	shell:
		"""
		samtools view -@ 16 -f 33 -F 284 {input.Bam} | awk '$7==\"=\"' |awk '$4>$8' |awk '$4-$8 > 150'| awk '$4-$8 < 10000'|awk '{{print $3, $4, $7, $8, \"{params.Elephant}\"}}' > {output.TandemDups}
		samtools view -@ 16 -f 33 -F 284 {input.Bam} | awk '{{if($4-$8 > 100000 || $7 != \"=\" ||$8-$4 > 100000)print}}'| awk '{{print $3, $4, $7, $8, \"{params.Elephant}\"}}' > {output.Rearrangements}
		"""

rule pyScript_ClustingingWithinElephant:
	input:
		SampleLocation + "/{sample}/{sample}_{sv}Reads.out"
	params:
		SV_Type = "{sv}"
	output:
		SampleLocation + "/{sample}/{sample}_Clustered{sv}.out"
	script:
		scripts+"WithinElephantClustering.py"

rule CombineWithinElephantClusteringFiles:
	input:
		expand(SampleLocation+"/{sample}/{sample}_Clustered{{sv}}.out", sample = SAMPLES)
	output:
		"Combined_Clustered{sv}.txt"
	shell:
		"sort -t, -k1,1 -k2,2n -k4,3n {input} > {output}"

rule pyScript_ClustingBetweenElephants:
	input:
		"Combined_Clustered{sv}.txt"
	output:
		"{sv}_all.csv"
	script:
		scripts+"ClusteringBetweenElephants.py"

rule Genotype_SVs:
	priority:
		1
	input:
		SVdataIN = "{sv}_all.csv",
		NonTE_Regions = 'NonTE_Regions.bed'
	params:
		SV_Type = "{sv}",
		PolarizationElephant = ANCESTRAL_SAMPLE_NAME,
		AllElephants = SAMPLES,
		sampleloc = SampleLocation
	output:
		"{sv}_all_genotyped.csv"
	script:
		scripts+"SV_GetGenotype.py"

rule Find_Genes_on_SVs:
	priority:
		1
	input:
		GTFfile,
		"{sv}_all_genotyped.csv"
	params:
		SV_Type = "{sv}"
	output:
		"{sv}_all_genotyped_GeneID.csv"
	script:
		scripts+"FindGenesSV.py"



rule pNpS_per_gene:
	input:
		path_vcf = "SNPs/AllChrom_SNPs_NoMissing.vcf.gz",
		GTF = GTFfile,
		DNDS_path = "SNPs_Synonymous_NSSN.csv",
		TD_data = "TandemDups_all_genotyped_GeneID.csv",
		RE_data = "Rearrangements_all_genotyped_GeneID.csv"
	output:
		"pNpS_per_gene.out"
	script:
		scripts+"pNpS_Pergene.py"


rule Find_TEs_td:
	input:
		SV_td = "TandemDups_all_genotyped_GeneID.csv",
		ref = REF,
		DB = "/projects/rogers_research/fnb/GabeO/ElephantProject/BlastDB/Repbase_all.fa"
	params:
		CHR = "{chromosome}"
	output:
		outfile_td = "TE_Files/TandemDups_{chromosome}_genotyped_GeneID_TE.csv",
	shell:
		'python ' + scripts + 'SV_FindTEs.py {input.SV_td} {input.ref} {input.DB} "TandemDups" {params.CHR} {output.outfile_td}'

rule Find_TEs_re:
	input:
		SV_re = "Rearrangements_all_genotyped_GeneID.csv",
		ref = REF,
		DB = "/projects/rogers_research/fnb/GabeO/ElephantProject/BlastDB/Repbase_all.fa"
	params:
		CHR = "{chromosome}"
	output:
		outfile_re = "TE_Files/Rearrangements_{chromosome}_genotyped_GeneID_TE.csv"
	shell:
		'python ' + scripts + 'SV_FindTEs.py {input.SV_re} {input.ref} {input.DB} "Rearrangements" {params.CHR} {output.outfile_re}'

		


rule Merge_TE_output:
	input:
		TD = expand("TE_Files/TandemDups_{chromosome}_genotyped_GeneID_TE.csv", chromosome = CHROM_LIST),
		RE = expand("TE_Files/Rearrangements_{chromosome}_genotyped_GeneID_TE.csv", chromosome = CHROM_LIST)
	output:
		TD = "TandemDups_all_genotyped_GeneID_TE.csv",
		RE = "Rearrangements_all_genotyped_GeneID_TE.csv"
	shell:
		"""
		awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input.TD} > {output.TD}
		awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input.RE} > {output.RE}
		"""