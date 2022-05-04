snakemake --use-conda \
--cluster-config ../config/cluster_config.yaml \
--cluster "sbatch --qos=maxjobs500 --time={cluster.duration} --account=scw1641 --job-name=imprinting --export=ALL --no-requeue --signal=2 --mem={cluster.total_mem} --output=smk.{rule}.%J.out --error=smk.{rule}.%J.err --ntasks={cluster.num_cores}" \
--keep-going \
-j 500 $@ 2> smk-"`date +"%d-%m-%Y"`".log 
mail -s "Snakemake has finished" camerond@cardiff.ac.uk < smk-"`date +"%d-%m-%Y"`".log
