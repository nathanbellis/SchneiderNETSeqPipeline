configfile: "config.yaml"

rule all:
	input:
		expand("sorted_reads/{sample}sorted.bam.bai", sample=config["samples"]),
		expand("rDNAcounts/{sample}plus.bed", sample=config["samples"]),
		expand("rDNAcounts/{sample}minus.bed", sample=config["samples"]),
		expand("htseq/{sample}counts.txt", sample=config["samples"]),
		expand("fastqcreports/{sample}_fastqc.html", sample=config["samples"]),
		expand("fastqcreports/{sample}fqt_fastqc.html", sample=config["samples"]),
		expand("fastqcreports/{sample}ca3_fastqc.html", sample=config["samples"]),
		expand("fastqcreports/{sample}ca5_fastqc.html", sample=config["samples"])

#rule all:
#	input:
#		"rDNAcounts/SRR13733831plus.bed",
#		"rDNAcounts/SRR13733831minus.bed",
#		"rDNAcounts/SRR13733832plus.bed",
#		"rDNAcounts/SRR13733832minus.bed",
#		"htseq/SRR13733831counts.txt",
#		"htseq/SRR13733832counts.txt",
#		"sorted_reads/SRR13733831sorted.bam.bai",
#		"sorted_reads/SRR13733832sorted.bam.bai",
#		"fastqcreports/SRR13733831_fastqc.html",
#		#"fastqcreports/SRR13733831fqt_fastqc.html",
#		"fastqcreports/SRR13733831ca3_fastqc.html",
#		"fastqcreports/SRR13733831ca5_fastqc.html",
#		"fastqcreports/SRR13733832_fastqc.html",
#                #"fastqcreports/SRR13733831fqt_fastqc.html",
#               "fastqcreports/SRR13733832ca3_fastqc.html",
#		"fastqcreports/SRR13733832ca5_fastqc.html"

rule fqtrim_collapse:
	input:
		"data/samples/{sample}.fastq"
	output:
		temp("fqtctemp/{sample}fqt.fastq")
	params:
		prefix="{sample}",
		folder="fqtctemp"
	shell:
		"""
		module load fqtrim
              	fqtrim -C -o fastq --outdir {params.folder} {input}
		mv fqtctemp/{params.prefix}.fastq fqtctemp/{params.prefix}fqt.fastq
              	"""

rule cutadapt_3:
	input:
		"fqtctemp/{sample}fqt.fastq"
	output:
		temp("ca3temp/{sample}ca3.fastq")
	threads:
		8
	shell:
		"""
		module load cutadapt
		cutadapt -g AGNNNNNNNNTG --no-indels -j {threads} -o {output} {input}
		"""

rule cutadapt_5:
	input:
		"ca3temp/{sample}ca3.fastq"
	output:
		temp("ca5temp/{sample}ca5.fastq")
	threads:
		8
	shell:
		"""
		module load cutadapt
		cutadapt -a CTGTAGGCACCAT -j {threads} --no-indels -e 0.10 -m 12 -o {output} {input}
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
		
rule fastqc_fqt:
	input:
	    "fqtctemp/{sample}.fastq"
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

rule star_align:
	input:
		fastq="ca5temp/{sample}ca5.fastq"
	params:
		starfolder="data/ref/",
		prefix="mapped_reads/{sample}"
	output:
		bam="mapped_reads/{sample}Aligned.sortedByCoord.out.bam",
		lfinal="mapped_reads/{sample}Log.final.out",
		lprog="mapped_reads/{sample}Log.progress.out",
		lout="mapped_reads/{sample}Log.out",
		lsj="mapped_reads/{sample}SJ.out.tab"
	threads:
		8
	shell:	
		"""
		module load STAR
		STAR --runThreadN {threads} \
		--genomeDir {params.starfolder} \
		--runMode alignReads \
		--limitBAMsortRAM 40000000000 \
		--alignIntronMax 1 \
		--outFilterMultimapNmax 20 \
		--outFilterMismatchNoverLmax 0.05 \
		--outFilterScoreMinOverLread 0 \
		--outFilterMatchNminOverLread 0 \
		--readFilesIn {input} \
		--outFileNamePrefix {params.prefix} \
		--outSAMtype BAM SortedByCoordinate
		"""

rule move_bam:
	input:
		"mapped_reads/{sample}Aligned.sortedByCoord.out.bam"
	output:
		"sorted_reads/{sample}sorted.bam"
	shell:
		"""
		mv {input} {output}
		"""

#rule samtools_sort:
#	input:
#		"mapped_reads/{sample}Aligned.out.bam"
#	output:
#		"sorted_reads/{sample}sorted.bam"
#	threads:
#		8
#	shell:
#		"""
#		module load SAMtools
#		samtools sort --threads {threads} -o {output} {input}
#		"""

rule samtools_index:
	input:
		"sorted_reads/{sample}sorted.bam"
	output:
		"sorted_reads/{sample}sorted.bam.bai"
	threads:
		8
	shell:
		"""
		module load SAMtools
		samtools index -@ {threads} {input}
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
		awk '($1 == "XII" && $2 >= 451575 && $2 < 458433)' {input} > {output}
		"""

#rule group_stats_wt_plus:
#	input:
#		bedlist=expand("rDNAcounts/{sample}plus.bed", sample=config["wildtype"]),
#	output:
#		"tables/WTplus.csv"
#	params:
#		norm="raw",
#		name="WT"
#	conda:
#        	"envs/SM.NETseq.yaml"
#	script:
#		"scripts/groupstats.r"

#rule group_stats_mut:
#	input:
#		expand("rDNAcounts/{sample}{strand}.bed", sample=config["mutant"], strand=config["strand"]),
#	output:
#		"tables/mutant.csv"
#	script:
#		"scripts/groupstats.r"


#rule compare_groups:
#	input:
#		wildtype="tables/wildtype.csv",
#		mutant="tables/mutant.csv"
#	output:
#		"tables/wt_v_mut.csv"
#	script:
#		"scripts/comparegroups.r"

#rule rDNA_plot:
#	input:
#
#	output:
#
#	script:
#		scripts/rDNAplot.r
