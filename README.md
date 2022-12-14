
# Schneider NET-Seq pipeline

This pipeline is inteded for the use in the Schneider lab on UAB's Cheaha HPC 
cluster. It utilizes snakemake to organize and execute the pipeline, allowing 
for the parallelization of jobs.

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



## Initial Set up

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