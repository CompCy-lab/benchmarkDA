#!/bin/bash

# add slurm module first
module purge
module load r/4.1.3

# set slurm parameters
time=12:00:00
partition=general

data_id=$1
root=$HOME/Documents/proj/benchmarkDA
cd ${root}/scripts


if [[ "$data_id" == "pancreas" ]]
    then
    data_dir=${root}/data/real/$data_id
    pops=("delta_cell" "alpha_cell" "gamma_cell" "acinar_cell" "beta_cell" "ductal_cell" "epsilon_cell")
    R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=1.2
    beta=80
    downsample=3
    mem=8g
elif [[ "$data_id" == "levine32" ]]
    then
    data_dir=${root}/data/real/$data_id
    pops=("pDCs" "CD4_T cells" "CD8_T cells" "Pre_B cells" "Mature_B cells" "Monocytes" "Basophils")
    R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.6
    beta=36
    downsample=25
    mem=96g
fi

## Run
for pop in "${pops[@]}";
    do
    for pop_enr in $(seq 0.75 0.1 0.95)
        do
        for seed in 43 44 45
            do
            for batch_sd in $batch_vec
                do
                for method in $R_methods
                    do
                    jobid="${data_id}-${pop}-${pop_enr}-${batch_sd}-${method}-${seed}"
                    if [[ "$method" == "cna" ]]; then
                        # enalbe cna env
                        source activate cna
                        cna_bin=$HOME/Documents/proj/benchmarkDA/methods/cna/bm_cna.py
                        sbatch -J "${jobid}" \
                          --time=${time} \
                          --partition=${partition} \
                          --mem ${mem} \
                          -o "$HOME/SlurmLog/${jobid}.out" \
                          --error="$HOME/SlurmLog/${jobid}.err" \
                          --wrap="python $cna_bin \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --k $k \
                            --pop '${pop}' \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --outdir ${root}/benchmark/${data_id}/"
                        # disable cna env
                        conda deactivate
                    elif [[ "$method" == "cna_batch" ]]; then
                        # enalbe cna env
                        source activate cna
                        cna_bin=$HOME/Documents/proj/benchmarkDA/methods/cna/bm_cna.py
                        sbatch -J "${jobid}" \
                          --time=${time} \
                          --partition=${partition} \
                          --mem ${mem} \
                          -o "$HOME/SlurmLog/${jobid}.out" \
                          --error="$HOME/SlurmLog/${jobid}.err" \
                          --wrap="python $cna_bin \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --k $k \
                            --pop '${pop}' \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --outdir ${root}/benchmark/${data_id}/ \
                            --model_batch"
                        # disable cna env
                        conda deactivate
                    elif [[ "$method" == "meld" ]]; then
                        # enable meld env
                        source activate meld
                        meld_bin=$HOME/Documents/proj/benchmarkDA/methods/meld/bm_meld.py
                        sbatch -J "${jobid}" \
                          --time=${time} \
                          --partition=${partition} \
                          --mem ${mem} \
                          -o "$HOME/SlurmLog/${jobid}.out" \
                          --error="$HOME/SlurmLog/${jobid}.err" \
                          --wrap="python $meld_bin \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --k $k \
                            --beta ${beta} \
                            --pop '${pop}' \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --outdir ${root}/benchmark/${data_id}/"
                        conda deactivate
                    else
                        sbatch -J "${jobid}" \
                          --time=${time} \
                          --partition=${partition} \
                          --mem ${mem} \
                          -o "$HOME/SlurmLog/${jobid}.out" \
                          --error="$HOME/SlurmLog/${jobid}.err" \
                          --wrap="Rscript ./run_DA.r \
                            ${data_dir}/${data_id}_data_bm.RDS $method $seed '$pop' \
                            --data_dir ${data_dir}/ \
                            --pop_enrichment $pop_enr \
                            --data_id $data_id \
                            --k $k \
                            --resolution ${resolution} \
                            --downsample ${downsample} \
                            --batchEffect_sd $batch_sd \
                            --outdir ${root}/benchmark/${data_id}/"
                    fi
                done
            done
        done
    done
done
