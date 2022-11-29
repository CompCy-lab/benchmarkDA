#!/bin/bash

# add slurm module
module purge
module load r/4.1.3

# Set input and output data directory
data_dir=$HOME/Documents/proj/benchmarkDA/data
outdir=$HOME/Documents/proj/benchmarkDA/data
cd $HOME/Documents/proj/benchmarkDA/scripts

## Run synthetic ##

data_id=linear
for p in $(seq 1 1 7)
   do
   for seed in 43 44 45
       do
       for enr in $(seq 0.75 0.1 0.95)
       do
           sbatch -J ${data_id}-M${p}-${enr}-${seed} \
             --time=0:30:00 \
             --partition=general \
             --mem 4g \
             -o $HOME/SlurmLog/${data_id}-M${p}-${enr}-${seed}.log \
             --error=$HOME/SlurmLog/${data_id}-M${p}-${enr}-${seed}.err \
             --wrap="Rscript ./make_bm_data.R ${data_dir}/synthetic/${data_id}/${data_id}_data_bm.RDS ${seed} M${p} \
               --pop_enrichment $enr \
               --make_batch_effect yes \
               --data_id $data_id \
               --outdir $outdir/synthetic/${data_id}/"
       done
   done
done


data_id=branch
for p in $(seq 1 1 8)
   do
   for seed in 43 44 45
       do
       for enr in $(seq 0.75 0.1 0.95)
       do
           sbatch -J ${data_id}-M${p}-${enr}-${seed} \
             --time=0:30:00 \
             --partition=general \
             --mem 4g \
             -o $HOME/SlurmLog/${data_id}-M${p}-${enr}-${seed}.log \
             --error=$HOME/SlurmLog/${data_id}-M${p}-${enr}-${seed}.err \
             --wrap="Rscript ./make_bm_data.R ${data_dir}/synthetic/${data_id}/${data_id}_data_bm.RDS ${seed} M${p} \
               --pop_enrichment $enr \
               --make_batch_effect yes \
               --data_id $data_id \
               --outdir $outdir/synthetic/${data_id}/"
       done
   done
done


data_id=cluster
for p in $(seq 1 1 3)
    do 
    for seed in 43 44 45
        do 
        for enr in $(seq 0.75 0.1 0.95)
        do
            sbatch -J ${data_id}-M${p}-${enr}-${seed} \
              --time=0:30:00 \
              --partition=general \
              --mem 4g \
              -o $HOME/SlurmLog/${data_id}-M${p}-${enr}-${seed}.log \
              --error=$HOME/SlurmLog/${data_id}-M${p}-${enr}-${seed}.err \
              --wrap="Rscript ./make_bm_data.R ${data_dir}/synthetic/${data_id}/${data_id}_data_bm.RDS ${seed} M${p} \
                --pop_enrichment $enr \
                --make_batch_effect yes \
                --data_id $data_id \
                --outdir $outdir/synthetic/${data_id}/"
        done
    done
done


data_id=cluster_balanced
for p in $(seq 1 1 3)
    do 
    for seed in 43 44 45
        do 
        for enr in $(seq 0.75 0.1 0.95)
        do
            sbatch -J ${data_id}-M${p}-${enr}-${seed} \
              --time=0:30:00 \
              --partition=general \
              --mem 4g \
              -o $HOME/SlurmLog/${data_id}-M${p}-${enr}-${seed}.log \
              --error=$HOME/SlurmLog/${data_id}-M${p}-${enr}-${seed}.err \
              --wrap="Rscript ./make_bm_data.R ${data_dir}/synthetic/${data_id}/${data_id}_data_bm.RDS ${seed} M${p} \
                --pop_enrichment $enr \
                --make_batch_effect yes \
                --data_id $data_id \
                --outdir $outdir/synthetic/${data_id}/"
        done
    done
done


## Run real data ##

data_id=covid19-pbmc
for p in RBC B PB CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC
    do
    for seed in 43 44 45
        do
        for enr in $(seq 0.75 0.1 0.95)
        do
            sbatch -J ${data_id}-${p}-${enr}-${seed} \
              --time=0:30:00 \
              --partition=general \
              --mem 8g \
              -o $HOME/SlurmLog/${data_id}-${p}-${enr}-${seed}.log \
              --error=$HOME/SlurmLog/${data_id}-${p}-${enr}-${seed}.err \
              --wrap="Rscript ./make_bm_data.R ${data_dir}/real/${data_id}/single-cell-atlas-pbmc-sars-cov2_sce.rds ${seed} ${p} \
                --pop_enrichment $enr \
                --make_batch_effect no \
                --data_id $data_id \
                --pop_col cell.type.coarse \
                --outdir $outdir/real/${data_id}/"
        done
    done
done


data_id=bcr-xl
for p in CD4_T-cells NK_cells CD8_T-cells B-cells_IgM+ monocytes surface- B-cells_IgM- DC
    do
    for seed in 43 44 45
        do
        for enr in $(seq 0.75 0.1 0.95)
        do
            sbatch -J ${data_id}-${p}-${enr}-${seed} \
              --time=0:30:00 \
              --partition=general \
              --mem 8g \
              -o $HOME/SlurmLog/${data_id}-${p}-${enr}-${seed}.log \
              --error=$HOME/SlurmLog/${data_id}-${p}-${enr}-${seed}.err \
              --wrap="Rscript ./make_bm_data.R ${data_dir}/real/${data_id}/bcr_xl_preprocessed_sce.rds ${seed} ${p} \
                --pop_enrichment $enr \
                --make_batch_effect no \
                --data_id $data_id \
                --pop_col cell_type \
                --reduced.dim X \
                --outdir $outdir/real/${data_id}/"
        done
    done
done