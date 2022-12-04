#!/bin/bash

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
	--jobs 12
