
# Schneider NET-Seq pipeline

This pipeline is inteded for the use in the Schneider lab on UAB's Cheaha HPC 
cluster. It utilizes snakemake to organize and execute the pipeline, allowing 
for the parallelization.

## File Overview

In the main pipeline folder you will see several files and a folder labeled data.
Files and their uses are as listed:
1. **Snakefile**
    - This file contains the rules that the pipeline will apply to your
    sequencing reads. This should not need to be modified unless you need to change
    steps to the pipeline or maybe if you start to tweak parameters for specific tools.
    It is possible if the sequencing adapters sequences change you may have to
    make those changes in this file. But for the most part this file should stay
    same. 
2. **cluster_config.yaml**
    - This file contains the parameters to run your pipeline on Cheaha. This file
    and config yaml will be the most modified. I will go into detail on what changes
    may be necessary later in this document.
3. **config.yaml**
    - This file contains the names and locations for your fastq files. The default
    config file has names of dummy samples that aren't include. The names of the
    files will have to be changed before running the pipeline.
4. **reset.sh**
    - This is just a handy script to clear the output folder and any .err and 
    .out files from Cheaha. It is not necessary, but it can be helpful if you are
    making significant changes and are debugging as it resets the folder back to 
    default.
5. **runpipeline.sh**
    - This script runs the pipeline. It calls snakemake to use the Snakefile with
    the cluster_config.yaml settings. Only one modification is necessary in this 
    file and that would be to possibly change the number of jobs. 
6. **tempdirclear.sh**
    - This script clears temporary directories for steps ca3 and ca5. Snakemake 
    actually can delete the temp .fastqs itself but can't delete the folders so 
    thats all this does. It is not necessary, but it can help keep your output
    folder cleaner and more organized.

There is one directory called **data** with three directories inside:
1. **gtf**
    - In this folder there is **genes.gtf** which contains the gene annotations 
    for the reference genome this pipeline uses (basically the locations of genes
    in the raw sequence). Additionally there is a **notes.txt** file which contains
    the specific modifications made to this .gtf file.
2. **ref**
    - **SacCer3onerDNArepeat.fna** is the reference genome used for this pipeline.
    Significant edits are made in this genome basically masking out repeated regions
    of the rDNA with Ns. This allows for the rDNA reads to be uniquely mapped 
    according to the aligner. Details on modifications are on the **notes.txt** 
    file in this folder.
    - **STARindex.job** and **STARindex.sh** are the files used to index the 
    genome. You should only have to run these scripts once. I typically split my
    scripts that I plan on submitting through sbatch into a .sh file (which has 
    the base command I want to run) and a .job file (which simply has the slurm
    parameters and calls the .sh file). I'm not sure if this is the best practice
    but it works. Details about how to run these will be in the initial set up
    section.
3. **samples**
    - This folder is empty except for a placeholder file. This is where your 
    .fastq's will go.



## Initial Set Up

To begin, download the folder to your user data storage in Cheaha. 
To do this you can either download the ZIP file from this page or with the 
command <br> `wget https://github.com/nathanbellis/SchneiderNETSeqPipeline/archive/refs/heads/main.zip` <br>
from there you can unzip the zipfile with `unzip main.zip`
Since all of our samples are done in *S. cerevisiae*, you should only need to 
index the genome once. The yeast genome is small enough it should handle copying 
and pasting to and from scratch pretty well.

### Indexing the Genome

Navigate to the data/ref/ folder and run the STARindex.job with the following 
command.<br>
`sbatch STARindex.job` <br>
This should only need to be run once unless you make modifications to the genome 
file or use a different genome.

### Misc. Changes
- In the cluster_config.yaml you will want to go in and change the account parameter
from my username nfbellis to yours.
- You can delete the placeholder file in the data/samples folder
- You may have to change the permissions to all the .sh files to make them 
executable. For this use the command: `chmod +x filename.sh`

## Running the Pipeline

- **STEP 1:** <br>
    Copy the main folder and all files and sub-folders to scratch space.
- **STEP 2:** <br>
    Copy your fastq files into the data/samples folder.
- **STEP 3:** <br>
    **Edit config.yaml**
    ```
    samples:
        SRR13733831: data/samples/SRR13733831.fastq
        SRR13733832: data/samples/SRR13733832.fastq
        SRR13733833: data/samples/SRR13733833.fastq
        SRR13733834: data/samples/SRR13733834.fastq
        SRR13733835: data/samples/SRR13733835.fastq
        SRR13733836: data/samples/SRR13733836.fastq
    ```
    to
    ```
    samples:
        yoursamplename1: data/samples/yoursamplename1.fastq
        yoursamplename2: data/samples/yoursamplename2.fastq
        yoursamplename3: data/samples/yoursamplename3.fastq
        yoursamplename4: data/samples/yoursamplename4.fastq
        yoursamplename5: data/samples/yoursamplename5.fastq
        yoursamplename6: data/samples/yoursamplename6.fastq
    ```
    You can include as many or as few samples as you if you are running it for 
    the first time it may be wise to just run one to make sure everything runs 
    well.
- **STEP 4:** <br>
    **Edit cluster_config.yaml** There shouldn't need to be too many modifications
    made in this file. All of the jobs get submitted to the express partition with one
    CPU unless they specifically are allocated below in this file. All jobs are
    given one hour which is usually more than enough. The htseq count step typically
    runs in 20-30 minutes so it is possible that for larger files you may want to add
    the following addition to cluster_config.yaml
    ```
    htseq_count:
        time: 2:00:00 # Or however long you think it may need. 

    ```
    This can also be applied to any steps where one hour is not enough. <br>
    
    Additionally fqtrim needs additional memory. For most of the steps 1 gigabyte
    (1G in the cluster_config.yaml) is enough, however fqtrim reads the whole 
    file into memory so I have been changing that parameter to the size of your 
    largest file. This pipeline was developed using a max fastq size of 9.75GB, hence
    the 10G in 
    ```
    fqtrim_collapse:
        mem-per-cpu: 10G # Change this to be larger than your largest fastq
    ```
- **STEP 5:** <br>
    **Edit runpipeline.sh** Only one modification really needs to be made here
    and that is changes the --jobs parameter. I developed this pipeline primarily 
    by running 6 fastq's at a time, so I chose 12 jobs as that is the theoretical
    max number of jobs that could be run at one time according to the Snakefile
    rules. Typically I would just set this to double the number of fastq files 
    you plan on running at a time
    ```
    snakemake --snakefile Snakefile \
	--cluster "sbatch -A {cluster.account} \
		-t {cluster.time} \
		-p {cluster.partition} \
		-N {cluster.nodes} \
		-c {cluster.cpus-per-task} \
		-o {cluster.output} \
		-e {cluster.error} \
		-n {cluster.ntasks} \
		--mem-per-cpu {cluster.mem-per-cpu}" \
	--cluster-config cluster_config.yaml \
	--jobs 12 # Change this to double the number of fastqs you want to run
    ```
    
- **STEP 6:** <br>
    **Run the pipeline** <br> 
    First load snakemake with: <br>
    `module load snakemake` <br>
    Make sure you are in the main folder with the Snakemake and yaml files.
    Before I run the pipeline I like to do a dry run with the command: <br>
    `snakemake -np` <br>
    This will use the Snakefile to piece together the pipeline and give you a 
    roadmap of the entire workflow and all the jobs with out actually running the
    commands. This is especially useful if you make any modifications in the 
    Snakefile like adding or removing a rule <br>
    
    If everything looks good you are ready to run the pipeline with the 
    following: <br>
    `./runpipeline.sh` <br>
    The ./ basically just lets the prompt know that this script is in the current 
    directory. If it doesn't run make sure the file is executable. Look at the 
    initial set up section for more information about changing the permissions.
    
- **STEP 7:** <br>
    **Monitor the pipeline** <br>
    While the jobs are submitted to the cluster via SLURM and don't require the 
    session to stay open, snakemake is running in the foreground submitting 
    multiple jobs as they are needed. For this reason the window running the 
    pipeline will need to be kept open. I suggest that you keep the window open 
    the first time you run the pipe line and periodically check on it to make 
    sure the jobs are running and being submitted successfully. You can check 
    on the jobs through the jobs tab on the research computing homepage or by 
    opening a new shell  window and using `squeue -u USERNAME` command. I 
    found the best way to monitor was to use the job tab to monitor jobs as it
    provided more information, and then have another window open where I would 
    check the output files. <br>
    
    If all runs smoothly there are ways to run Snakemake in the background, or 
    another window. If this is something you find you would like to learn. Look 
    into the linux commands `&`, `nohup`, or `tmux`, which all provide ways of 
    running snakemake in the background or even when you are logged out.
    
- **STEP 8**: <br>
    **Copy the output file** <br>
    The pipeline will put all resulting files into a folder labeled output. When
    the full pipeline is run and there are no more jobs running, copy the output 
    folder and move it to your USER_DATA folder for more long term storage. From
    here you can take the files and do further analysis in R. Once data has been
    moved and you are satisfied with the results, you can then delete the 
    main folder from scratch keeping the scratch space clean.

## Step by step breakdown

1. **rule all** <br>
    ```
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
	```
    This rule tells snakemake what files should be created at the end of the 
    pipeline. It then backtracks and works out what rules it needs to run based
    on matching input and output. You want to make sure that all files that 
    all files that have a terminal endpoint are included in this rule. For 
    example the first rule points to the {sample}sorted.bam.bai. No other rule 
    uses this bai file so I need to tell snakemake to make it. However you will
    see I don't need to call the .bam file because it's creation and use is 
    implicit as a necessary step to created the .bed files.
    
2. **rule fqtrim_collapse** <br>
    ```
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
    ```
    This rule works on the original fastq files to collapse duplicate reads due 
    to PCR. It relies on 5' adapters with UMI so it must be done before they are
    removed. -C is the parameter to indicate collapse. This is mentioned in the 
    running the pipeline, but the -C command reads the entire file into memory
    so you need to request slightly more memory than your largest file via 
    SLURM <br>
    [fqtrim documentation](https://ccb.jhu.edu/software/fqtrim/)
4. **rule cutadapt_5** <br>
    ```
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
    ```
    This rule trims the 5' adapter according to the sequence AGNNNNNNNNTG, where 
    N is any nucleotide. -g indicates 5' adapter, and --no-indels means when 
    matching no insertions or deletions will be allow. Unlike fqtrim cutadapt 
    can be multithreaded to speed up performance so I have told cutadapt to use
    8 cores, this must also be requested in the SLURM request via the 
    cluster_config.yaml file or else it will not multithread. This applies to
    all other rules that use multithreading as well. <br>
    [cutadapt documentation](https://cutadapt.readthedocs.io/en/stable/guide.html)
    
5. **rule cutadapt_3_1** <br>
    ```
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
    ```
    This rule is similar to cutadapt_5 except it trims the 3' adapter. The 
    sequence for this is CTGTAGGCACCAT. If the 3' adapter sequence changes, this 
    will need to be changed. Like the previous step --no-indels is used. -e is 
    the error rate. -m 12 is the minimum read length, meaning if a sequence is 
    trimmed to shorter than 12 nts it is discarded.
6. **rule cutadapt_3_2** <br>
    ```
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
    ```
    This is the same rule and settings as the step above however it runs a 
    second 3' trim on the sequence minus the first two nts. The reason for this 
    is that the 5' adapter ends in NTG so occasionally when the 5' adapter is 
    trimmed it will also remove the beginning of the 3' adapter 
    "**CTG**TAGGCACCAT" if the read is just adapters. This means that the 
    subsequent 3' trimming will not recognize the 5' adapter. This step is not 
    necessary but it does help clean up some off target reads.

7. **rule fastqc_original** <br>
    ```
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
    ```
    This runs a fastqc command on the original fastqs. This is a helpful tool 
    for gauging the quality of your initial reads. <br>
    [fastqc documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    
8. **rule fastqc_fqt** <br>
    This runs the same rule as fastqc_original but on the fastq after the fqtrim
    -C rule. It is helpful to have the fastqc's to track the changes being made
    to the fastq before aligning.
    
9. **rule fastqc_ca5** <br>
    This runs the same rule as fastqc_original but on the fastq after the 
    cutadapt_5 rule It is helpful to have the fastqc's to track the changes 
    being made to the fastq before aligning.
10. **rule fastqc_ca3** <br>
    This runs the same rule as fastqc_original but on the fastq after the 
    cutadapt_3_2 rule It is helpful to have the fastqc's to track the changes 
    being made to the fastq before aligning.
11. **rule star_align** <br>
    ```
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
    ```
    This is the alignment step. STAR aligner is used to map reads to the 
    yeast genome. The genome file parameter is the starfolder in params:. This 
    would need to be changed if you planned on using a different genome. There 
    are many parameters all of which are explained in the STAR documentation.
    [STAR documentation](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
    
12. **rule move_bam** <br>
    ```
    rule move_bam:
    	input:
    		"output/STAR_output/{sample}Aligned.sortedByCoord.out.bam"
    	output:
    		"output/sorted_reads/{sample}sorted.bam"
    	shell:
    		"""
    		mv {input} {output}
    		""
    ```
    This step simply moves the bam file output by STAR to a new folder for less
    clutter.
13. **rule samtools_index** <br>
    ```
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
    ```
    This step creates an index file for the bam file created by STAR. It is not
    necessary for any downstream rules in this pipeline, however it is necessary
    for many tools that use .bam files. bam files must be sorted prior to 
    indexing, however STAR can sort as it aligns which is why there is not a 
    separate samtools_sort step. While not a major part of the pipeline, I have 
    included documentation for both samtools and an overview of .bam and .sam 
    file format. This is useful if you want to access the quality of the 
    alignment at a single read level. <br>
    [samtools documentation](http://www.htslib.org/doc/) <br>
    [sam/bam file format](https://samtools.github.io/hts-specs/SAMv1.pdf) <br>
14. **rule htseq_gene_counts** <br>
    ```
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
    ```
    This step creates gene counts for your aligned reads as if you were doing
    mRNA seq. This will not be a heavily used step for analysis but I find it's 
    helpful data for quality control. For good Pol I NET-seq reads, the vast 
    majority of the reads should be coming from rRNA regions (25S, 18S, 5.8S,
    ETS1/2, ITS1/2) if there are alot of protein coding genes or if the ITS and
    ETS regions have low reads this could be an indication of poor IP. A helpful
    command for looking at these reads is `sort -nrk 2 file.counts | less` which
    sorts the counts by highest value.
15. **rule bedtools_gc_plus** <br> 
    ```
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
    ```
    This rule creads a .bed file that counts where the 5' end is at an individual 
    nucleotide leve over the entire genome, set by the -5 parameter. This rule 
    specifically runs for the plus strand. <br>
    [bedtools genomecov documentation](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)

16. **rule bedtools_gc_minus** <br>
    ```
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
    ```
    This rule is the exact same as above except the + is change to a - for the 
    minus strand.
    
17. **rule awk_rDNA_trim** <br>
    ```
    rule awk_rDNA_trim:
    	input:
    		"output/counts/{sample}.bed"
    	output:
    		"output/rDNAcounts/{sample}.bed"
    	shell:
    		"""
    		awk '($1 == "XII" && $2 >= 451575 && $2 < 458433)' {input} > {output}
    		"""
    ```
    This rule trims the bed files to just the rDNA region. The .bed file for the 
    entire genome is relatively large and quite unwieldy in terms of data so 
    this step trims the data to just the region we are interested (Chromosome 
    XII from position 451575 to 458433) making it easier to import for 
    downstream analysis on R, which can be done on a local computer. This step 
    can also be done on R from the full .bed files, which are also included in 
    the final output folder. This step just captures the transcribed region so 
    if you want the NTS reads use the full .bed file and trim in R or change 
    these parameters to the coordinates you want.
       


