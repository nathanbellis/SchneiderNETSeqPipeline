configfile: "config.yaml"

#rule all:
#	input:
#		expand("sorted_reads/{sample}sorted.bam.bai", sample=config["samples"]),
#		expand("rDNAcounts/{sample}plus.bed", sample=config["samples"]),
#		expand("rDNAcounts/{sample}minus.bed", sample=config["samples"]),
#		expand("htseq/{sample}counts.txt", sample=config["samples"]),
#		expand("fastqcreports/{sample}_fastqc.html", sample=config["samples"]),
#		#expand("fastqcreports/{sample}fqt_fastqc.html", sample=config["samples"]),
#		expand("fastqcreports/{sample}ca3_fastqc.html", sample=config["samples"]),
#		expand("fastqcreports/{sample}ca5_fastqc.html", sample=config["samples"]),
#		"tables/wt_v_mut.csv"

rule all:
	input:
		"rDNAcounts/SRR13733831plus.bed",
		"rDNAcounts/SRR13733831minus.bed",
		"rDNAcounts/SRR13733832plus.bed",
		"rDNAcounts/SRR13733832minus.bed",
		"htseq/SRR13733831counts.txt",
		"htseq/SRR13733832counts.txt",
		"sorted_reads/SRR13733831sorted.bam.bai",
		"sorted_reads/SRR13733832sorted.bam.bai",
		"fastqcreports/SRR13733831_fastqc.html",
		#"fastqcreports/SRR13733831fqt_fastqc.html",
		"fastqcreports/SRR13733831ca3_fastqc.html",
		"fastqcreports/SRR13733831ca5_fastqc.html",
		"fastqcreports/SRR13733832_fastqc.html",
                #"fastqcreports/SRR13733831fqt_fastqc.html",
                "fastqcreports/SRR13733832ca3_fastqc.html",
		"fastqcreports/SRR13733832ca5_fastqc.html"

#rule fqtrim_collapse:
#	input:
#		"data/samples/{sample}.fastq"
#	output:
#		temp("fqtctemp/{sample}fqt.fastq") 
#	shell:
#		"""
#		module load fqtrim
#              	fqtrim -C -o fastq --outdir fqtctemp {input}
#		mv fqtctemp/{sample}.fastq fqtctemp/{sample}fqt.fastq
#              	"""

rule cutadapt_3:
	input:
		"data/samples/{sample}.fastq"
	output:
		temp("ca3temp/{sample}ca3.fastq")
	shell:
		"""
		module load cutadapt
		cutadapt -g AGNNNNNNNNTG --no-indels -o {output} {input}
		"""

rule cutadapt_5:
	input:
		"ca3temp/{sample}ca3.fastq"
	output:
		temp("ca5temp/{sample}ca5.fastq")
	shell:
		"""
		module load cutadapt
		cutadapt -a CTGTAGGCACCAT --no-indels -e 0.10 -m 12 -o {output} {input}
		"""

rule fastqc_original:
	input:
		"data/samples/{sample}.fastq"
		
	output:	
		zip="fastqcreports/{sample}_fastqc.zip",
		html="fastqcreports/{sample}_fastqc.html"
	params:
		path ="fastqcreports/"
	shell:
		"""
		module load FastQC
		fastqc -o {params.path} {input}
		"""
		
#rule fastqc_fqt:
#	input:
#	    "fqtctemp/{sample}.fastq"
#	output:	
#		zip="fastqcreports/{sample}_fastqc.zip",
#		html="fastqcreports/{sample}_fastqc.html"
#	params:
#		path ="fastqcreports/"
#	shell:
#		"""
#		module load FastQC
#		fastqc -o {params.path} {input}
#		"""

rule fastqc_ca3:
	input:
	    "ca3temp/{sample}.fastq"
		
	output:	
		zip="fastqcreports/{sample}_fastqc.zip",
		html="fastqcreports/{sample}_fastqc.html"
	params:
		path ="fastqcreports/"
	shell:
		"""
		module load FastQC
		fastqc -o {params.path} {input}
		"""
		
rule fastqc_ca5:
	input:
	    "ca5temp/{sample}.fastq"
		
	output:	
		zip="fastqcreports/{sample}_fastqc.zip",
		html="fastqcreports/{sample}_fastqc.html"
	params:
		path ="fastqcreports/"
	shell:
		"""
		module load FastQC
		fastqc -o {params.path} {input}
		"""

rule bowtie2_align:
	input:
		fastq="ca5temp/{sample}ca5.fastq"
	params:
		btext="data/ref/SC1rDNA"
	output:
		bam=temp("mapped_reads/{sample}.temp.sam")
	log:
		"bt2logs/{sample}bowtie2.log"
	shell:	
		"""
		module load Bowtie2
		bowtie2 -x {params.btext} -U {input.fastq} -S {output.bam} 2> {log}
		"""

rule samtools_to_bam:
	input:
		"mapped_reads/{sample}.temp.sam"
	output:
		temp("mapped_reads/{sample}.bam")
	shell:
		"""
		module load SAMtools
		samtools view -bS {input} > {output}
		"""

rule samtools_sort:
	input:
		"mapped_reads/{sample}.bam"
	output:
		"sorted_reads/{sample}sorted.bam"
	shell:
		"""
		module load SAMtools
		samtools sort -o {output} {input}
		"""

rule samtools_index:
	input:
		"sorted_reads/{sample}sorted.bam"
	output:
		"sorted_reads/{sample}sorted.bam.bai"
	shell:
		"""
		module load SAMtools
		samtools index {input} {output}
		"""

rule htseq_gene_counts:
	input:
		"sorted_reads/{sample}sorted.bam"
	output:
		"htseq/{sample}counts.txt"
	shell:
		"""
		module load HTSeq
		htseq-count -s reverse {input} ./data/gtf/genes.gtf > {output}
		"""


rule bedtools_gc_plus:
	input:
		"sorted_reads/{sample}sorted.bam"
	output:
		"counts/{sample}plus.bed"
	shell:
		"""
		module load BEDTools
		bedtools genomecov -d -5 -strand + -ibam {input} > {output}
		"""

rule bedtools_gc_minus:
	input:
		"sorted_reads/{sample}sorted.bam"
	output:
		"counts/{sample}minus.bed"
	shell:
		"""
		module load BEDTools
		bedtools genomecov -d -5 -strand - -ibam {input} > {output}
		"""

rule awk_rDNA_trim:
	input:
		"counts/{sample}.bed"
	output:
		"rDNAcounts/{sample}.bed"
	shell:
		"""
		awk '($1 == "XII" && $2 > 451000 && $2 < 477000)' {input} > {output}
		"""

rule group_stats_wt:
	input:
		expand("rDNAcounts/{sample}{strand}.bed", sample=config["wildtype"], strand=config["strand"]),
	output:
		"tables/wildtype.csv"
	script:
		"scripts/groupstats.r"

rule group_stats_mut:
	input:
		expand("rDNAcounts/{sample}{strand}.bed", sample=config["mutant"], strand=config["strand"]),
	output:
		"tables/mutant.csv"
	script:
		"scripts/groupstats.r"


rule compare_groups:
	input:
		wildtype="tables/wildtype.csv",
		mutant="tables/mutant.csv"
	output:
		"tables/wt_v_mut.csv"
	script:
		"scripts/comparegroups.r"

#rule rDNA_plot:
#	input:
#
#	output:
#
#	script:
#		scripts/rDNAplot.r
