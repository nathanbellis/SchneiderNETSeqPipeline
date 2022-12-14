
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

To begin, download the folder to your user data storage in Cheaha. WGET??. Since
all of our samples are done in *S. cerevisiae*, you should only need to align the
genome once. The yeast genome is small enough it should handle copying and pasting
to and from scratch pretty well.

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

    


