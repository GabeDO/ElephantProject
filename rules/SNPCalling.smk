rule snp_calling_call:
	priority:
		2 
	input:
		bams = SampleLocation+"/{sample}/{sample}_sampe_sorted.bam", 
		ref = REF
	output:
		NonAncest = SampleLocation+"/{sample}/{sample}_SNPs_{chromosome}.vcf.gz"
	shell:
		"""
		bcftools mpileup -r {wildcards.chromosome} -f {input.ref} {input.bams} | bcftools call -m -Oz -o {output.NonAncest}
		"""

rule snp_make_index:
	priority:
		2
	input:
		NonAncest = SampleLocation+"/{sample}/{sample}_SNPs_{chromosome}.vcf.gz"
	output:
		NonAncest = SampleLocation+"/{sample}/{sample}_SNPs_{chromosome}.vcf.gz.tbi"
	shell:
		"""
		tabix -p vcf {input.NonAncest}
		"""

rule snp_vcf_merge:
	priority:
		2
	input:
		VCFFiles = expand(SampleLocation+"/{sample}/{sample}_SNPs_{{chromosome}}.vcf.gz",sample = SAMPLES_Asian),
		index = expand(SampleLocation+"/{sample}/{sample}_SNPs_{{chromosome}}.vcf.gz.tbi",sample = SAMPLES_Asian)
	output:
		MergedFile = "SNPs/ChromosomeVCFs/{chromosome}_Elephant_SNPs.vcf.gz"
	shell:
		"""
		bcftools merge -Oz -o {output.MergedFile} {input.VCFFiles} 
		"""

rule filter_SNPs:
	input:
		NonAncest = "SNPs/ChromosomeVCFs/{chromosome}_Elephant_SNPs.vcf.gz",
	output:
		NonAncest = "SNPs/ChromosomeVCFs/{chromosome}_Elephant_SNPs_NoMissing.vcf.gz",
	shell:
		"""
		bcftools view --types snps -i 'F_MISSING<0.1' -Oz -o  {output.NonAncest} {input.NonAncest}
		"""

rule merge_chrom_VCFs:
	input:
		vcfs = expand("SNPs/ChromosomeVCFs/{chromosome}_Elephant_SNPs_NoMissing.vcf.gz", chromosome = CHROM_LIST)
	output:
		mergedchroms = "SNPs/AllChrom_Elephant_SNPs_NoMissing.vcf.gz"
	shell:
		"""
		bcftools concat {input.vcfs} --threads 16 -Oz -o {output.mergedchroms}
		"""

rule MakeAlleleFreq:
	priority:
		2
	input:
		NonAncest = "SNPs/AllChrom_Elephant_SNPs_NoMissing.vcf.gz",
	output:
		o1 = "SNPs/Elephant_SNPs_NoMissing.geno",
	shell:
		'python ' + scripts + 'SNP_GenoCounts.py {input.NonAncest} {output.o1}'

rule Find_Genes_on_SNPs:
	input:
		GTFfile,
		"SNPs/AllChrom_Elephant_SNPs_NoMissing.vcf.gz"
	output:
		"Genes_Found_Near_SNPs.out"
	script:
		scripts+"FindGenesSNPs.py"

rule Find_Synonymous_SNPs:
	input:
		GFT = GTFfile,
		VCF = "SNPs/ChromosomeVCFs/{chromosome}_Elephant_SNPs_NoMissing.vcf.gz",
		ref = REF
	output:
		all_o = "SNPs/ChromosomeVCFs/{chromosome}_SNPs_Synonymous.csv",
		NSSN_o = "SNPs/ChromosomeVCFs/{chromosome}_SNPs_Synonymous_NSSN.csv"
	script:
		scripts+"SNPs_CheckIfSynonymous.py"

rule Merge_Synonymous_SNPs:
	input:
		all = expand("SNPs/ChromosomeVCFs/{chromosome}_SNPs_Synonymous.csv", chromosome = CHROM_LIST),
		NSSN = expand("SNPs/ChromosomeVCFs/{chromosome}_SNPs_Synonymous_NSSN.csv", chromosome = CHROM_LIST)
	output:
		all_o = "SNPs_Synonymous.csv",
		NSSN_o ="SNPs_Synonymous_NSSN.csv"
	shell:
		"""
		cat {input.all} > {output.all_o}
		cat {input.NSSN} > {output.NSSN_o}
		"""

rule Calculate_Pi:
	input:
		VCF = "SNPs/AllChrom_Elephant_SNPs_NoMissing.vcf.gz"
	params:
		chrom = CHROM_LIST
	output:
		pi_output = "SNPs/Pi_window.out"
	script:
		scripts+"SNP_Pi.py"


#rule TWINS:
#	input:
#		Samples1 = /projects/rogers_research/fnb/GabeO/ElephantProject/Samples/{malaylist}/*.vcf.gz,
#		Samples2 = /projects/rogers_research/fnb/GabeO/ElephantProject/Samples/{malaylist}/*.vcf.gz,
#	output:
#		ele1 = temp(Ele1.vcf.gz)
#		ele2 = temp(Ele2.vcf.gz)
#		ele12 = temp(Ele12.vcf.gz)
#		ele12f = temp(Ele12f.vcf.gz)
#	shell:
#		"""
#		bcftools concat {input.Samples1} -Oz -o {output.ele1}
#		tabix -p vcf -f {output.ele1}
#		bcftools concat {input.Samples2} -Oz -o {output.ele2}
#		tabix -p vcf -f {output.ele2}
#
#
#		bcftools merge -Oz -o {output.ele12} {output.ele2} {output.ele1}
#		bcftools view --types snps -g ^miss -Oz -o {output.ele12f} {output.ele12}
#
#		python twins.py Twins_Filtered.vcf.gz /projects/rogers_research/fnb/GabeO/ElephantProject/SnakeMakePipeline/SNPs/TWINS\!/TwinsOutput.txt
#		"""
#echo "$i $j" >> TwinsOutput.txt

