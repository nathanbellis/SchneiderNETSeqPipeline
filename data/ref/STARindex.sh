#!/bin/bash

module load STAR

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ./ --genomeSAindexNbases 10\
     --genomeFastaFiles ./SacCer3oneAltOneRep.fna --sjdbGTFfile ../gtf/genes.gtf

