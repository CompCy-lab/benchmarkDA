#!/bin/bash

# add slurm module
module purge
module load r/4.1.3

# Set input and output data directory
data_dir=$HOME/Documents/proj/benchmarkDA/data
outdir=$HOME/Documents/proj/benchmarkDA/data
cd $HOME/Documents/proj/benchmarkDA/scripts

## Run real data ##

data_id=pancreas
for p in "delta cell" "alpha cell" "gamma cell" "acinar cell" "beta cell" "ductal cell" "epsilon cell"
    do
    for seed in 43 44 45
        do
        for enr in $(seq 0.75 0.1 0.95)
        do
            sbatch -J "${data_id}-${p}-${enr}-${seed}" \
              --time=0:30:00 \
              --partition=general \
              --mem 8g \
              -o "$HOME/SlurmLog/${data_id}-${p}-${enr}-${seed}.log" \
              --error="$HOME/SlurmLog/${data_id}-${p}-${enr}-${seed}.err" \
              --wrap="Rscript ./make_bm_data.R ${data_dir}/real/${data_id}/pancreas_preprocessed_sce.rds ${seed} '${p}' \
                --pop_enrichment $enr \
                --make_batch_effect no \
                --data_id $data_id \
                --pop_col Factor.Value.inferred.cell.type...authors.labels. \
                --outdir $outdir/real/${data_id}/"
        done
    done
done


data_id=levine32
for p in "pDCs" "CD4 T cells" "CD8 T cells" "Pre B cells" "Mature B cells" "Monocytes" "Basophils"
    do
    for seed in 43 44 45
        do
        for enr in $(seq 0.75 0.1 0.95)
        do
            sbatch -J "${data_id}-${p}-${enr}-${seed}" \
              --time=0:30:00 \
              --partition=general \
              --mem 8g \
              -o "$HOME/SlurmLog/${data_id}-${p}-${enr}-${seed}.log" \
              --error="$HOME/SlurmLog/${data_id}-${p}-${enr}-${seed}.err" \
              --wrap="Rscript ./make_bm_data.R ${data_dir}/real/${data_id}/levine32_preprocessed_sce.rds ${seed} '${p}' \
                --pop_enrichment $enr \
                --make_batch_effect no \
                --data_id $data_id \
                --pop_col cell_type \
                --reduced.dim X \
                --outdir $outdir/real/${data_id}/"
        done
    done
done