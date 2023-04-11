configfile: "config.yaml"

rule all:
	input:
		expand("output/sorted_reads/{sample}sorted.bam.bai", sample=config["samples"]),
		expand("output/rDNAcounts/{sample}plus.bed", sample=config["samples"]),
		expand("output/rDNAcounts/{sample}minus.bed", sample=config["samples"]),
		expand("output/htseq/{sample}counts.txt", sample=config["samples"]),
		expand("output/fastqcreports/{sample}_fastqc.html", sample=config["samples"]),
		expand("output/fastqcreports/{sample}fqt_fastqc.html", sample=config["samples"]),
		expand("output/fastqcreports/{sample}ca3_fastqc.html", sample=config["samples"]),
		expand("output/fastqcreports/{sample}ca5_fastqc.html", sample=config["samples"])

rule fqtrim_collapse:
	input:
		"data/samples/{sample}.fastq"
	output:
		temp("output/fqtctemp/{sample}fqt.fastq")
	params:
		prefix="{sample}",
		folder="output/fqtctemp"
	shell:
		"""
		module load fqtrim
              	fqtrim -C -o fastq --outdir {params.folder} {input}
		mv output/fqtctemp/{params.prefix}.fastq output/fqtctemp/{params.prefix}fqt.fastq
              	"""

rule cutadapt_5:
	input:
		"output/fqtctemp/{sample}fqt.fastq"
	output:
		temp("output/ca5temp/{sample}ca5.fastq")
	threads:
		8
	shell:
		"""
		module load cutadapt
		cutadapt -g AGNNNNNNNNTG --no-indels -j {threads} -o {output} {input}
		"""

rule cutadapt_3_1:
	input:
		"output/ca5temp/{sample}ca5.fastq"
	output:
		temp("output/ca3temp/{sample}ca31.fastq")
	threads:
		8
	shell:
		"""
		module load cutadapt
		cutadapt -a CTGTAGGCACCAT -j {threads} --no-indels -e 0.10 -m 12 -o {output} {input}
		"""

rule cutadapt_3_2:
        input:
                "output/ca3temp/{sample}ca31.fastq"
        output:
                temp("output/ca3temp/{sample}ca3.fastq")
        threads:
                8
        shell:
                """
                module load cutadapt
                cutadapt -a TAGGCACCATCAA -j {threads} --no-indels -e 0.10 -m 12 -o {output} {input}
                """

rule fastqc_original:
	input:
		"data/samples/{sample}.fastq"
		
	output:	
		zip="output/fastqcreports/{sample}_fastqc.zip",
		html="output/fastqcreports/{sample}_fastqc.html"
	params:
		path ="output/fastqcreports/"
	shell:
		"""
		module load FastQC
		fastqc -o {params.path} {input}
		"""
		
rule fastqc_fqt:
	input:
	    "output/fqtctemp/{sample}.fastq"
	output:	
		zip="output/fastqcreports/{sample}_fastqc.zip",
		html="output/fastqcreports/{sample}_fastqc.html"
	params:
		path ="output/fastqcreports/"
	shell:
		"""
		module load FastQC
		fastqc -o {params.path} {input}
		"""

rule fastqc_ca3:
	input:
		"output/ca3temp/{sample}.fastq"
		
	output:	
		zip="output/fastqcreports/{sample}_fastqc.zip",
		html="output/fastqcreports/{sample}_fastqc.html"
	params:
		path ="output/fastqcreports/"
	shell:
		"""
		module load FastQC
		fastqc -o {params.path} {input}
		"""
		
rule fastqc_ca5:
	input:
		"output/ca5temp/{sample}.fastq"
		
	output:	
		zip="output/fastqcreports/{sample}_fastqc.zip",
		html="output/fastqcreports/{sample}_fastqc.html"
	params:
		path ="output/fastqcreports/"
	shell:
		"""
		module load FastQC
		fastqc -o {params.path} {input}
		"""

rule star_align:
	input:
		fastq="output/ca3temp/{sample}ca3.fastq"
	params:
		starfolder="data/ref/",
		prefix="output/STAR_output/{sample}"
	output:
		bam="output/STAR_output/{sample}Aligned.sortedByCoord.out.bam",
		lfinal="output/STAR_output/{sample}Log.final.out",
		lprog="output/STAR_output/{sample}Log.progress.out",
		lout="output/STAR_output/{sample}Log.out",
		lsj="output/STAR_output/{sample}SJ.out.tab"
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
		"output/STAR_output/{sample}Aligned.sortedByCoord.out.bam"
	output:
		"output/sorted_reads/{sample}sorted.bam"
	shell:
		"""
		mv {input} {output}
		"""

rule samtools_index:
	input:
		"output/sorted_reads/{sample}sorted.bam"
	output:
		"output/sorted_reads/{sample}sorted.bam.bai"
	threads:
		8
	shell:
		"""
		module load SAMtools
		samtools index -@ {threads} {input}
		"""

rule htseq_gene_counts:
	input:
		"output/sorted_reads/{sample}sorted.bam"
	output:
		"output/htseq/{sample}counts.txt"
	shell:
		"""
		module load HTSeq
		htseq-count -s reverse {input} ./data/gtf/genes.gtf > {output}
		"""


rule bedtools_gc_plus:
	input:
		"output/sorted_reads/{sample}sorted.bam"
	output:
		"output/counts/{sample}plus.bed"
	shell:
		"""
		module load BEDTools
		bedtools genomecov -d -5 -strand + -ibam {input} > {output}
		"""

rule bedtools_gc_minus:
	input:
		"output/sorted_reads/{sample}sorted.bam"
	output:
		"output/counts/{sample}minus.bed"
	shell:
		"""
		module load BEDTools
		bedtools genomecov -d -5 -strand - -ibam {input} > {output}
		"""

rule awk_rDNA_trim:
	input:
		"output/counts/{sample}.bed"
	output:
		"output/rDNAcounts/{sample}.bed"
	shell:
		"""
		awk '($1 == "XII" && $2 >= 458433 && $2 <= 467569)' {input} > {output}
		"""
