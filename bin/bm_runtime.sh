#!/bin/bash

# add slurm module first
module purge
module load r/4.1.3

# set slurm parameters
time=6:00:00
partition=general

data_id=$1
root=$HOME/Documents/proj/benchmarkDA
cd ${root}/scripts

if [[ "$data_id" == "test_scale_4000" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=1
    beta=60
    downsample=10
    mem=8g
elif [[ "$data_id" == "test_scale_10000" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=1
    beta=60
    downsample=10
    mem=8g
elif [[ "$data_id" == "test_scale_15000" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=1
    beta=60
    downsample=10
    mem=16g
elif [[ "$data_id" == "test_scale_30000" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=1
    beta=60
    downsample=10
    mem=16g
elif [[ "$data_id" == "test_scale_50000" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=1
    beta=60
    downsample=10
    mem=28g
elif [[ "$data_id" == "test_scale_100000" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=1
    beta=60
    downsample=10
    mem=28g
fi


## Run
for pop in $pops
    do
    for pop_enr in 0.85
        do
        for seed in 43 44 45
            do
            for batch_sd in $batch_vec
                do
                for method in $R_methods
                    do
                    jobid=${data_id}-${pop}-${pop_enr}-${batch_sd}-${method}-${seed}
                    if [[ "$method" == "cna" ]]; then
                        # enalbe cna env
                        source activate cna
                        cna_bin=$HOME/Documents/proj/benchmarkDA/methods/cna/bm_cna.py
                        sbatch -J ${jobid} \
                          --time=${time} \
                          --partition=${partition} \
                          --mem ${mem} \
                          -o $HOME/SlurmLog/${jobid}.out \
                          --error=$HOME/SlurmLog/${jobid}.err \
                          --wrap="python $cna_bin \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --k $k \
                            --pop ${pop} \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --outdir ${root}/benchmark/${data_id}/"
                        # disable cna env
                        conda deactivate
                    elif [[ "$method" == "cna_batch" ]]; then
                        # enalbe cna env
                        source activate cna
                        cna_bin=$HOME/Documents/proj/benchmarkDA/methods/cna/bm_cna.py
                        sbatch -J ${jobid} \
                          --time=${time} \
                          --partition=${partition} \
                          --mem ${mem} \
                          -o $HOME/SlurmLog/${jobid}.out \
                          --error=$HOME/SlurmLog/${jobid}.err \
                          --wrap="python $cna_bin \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --k $k \
                            --pop ${pop} \
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
                        sbatch -J ${jobid} \
                          --time=${time} \
                          --partition=${partition} \
                          --mem ${mem} \
                          -o $HOME/SlurmLog/${jobid}.out \
                          --error=$HOME/SlurmLog/${jobid}.err \
                          --wrap="python $meld_bin \
                            --data_dir ${data_dir} \
                            --data_id ${data_id} \
                            --pop_enr $pop_enr \
                            --k $k \
                            --beta ${beta} \
                            --pop ${pop} \
                            --be_sd $batch_sd \
                            --seed $seed \
                            --outdir ${root}/benchmark/${data_id}/"
                        conda deactivate
                    else
                        sbatch -J ${jobid} \
                          --time=${time} \
                          --partition=${partition} \
                          --mem ${mem} \
                          -o $HOME/SlurmLog/${jobid}.out \
                          --error=$HOME/SlurmLog/${jobid}.err \
                          --wrap="Rscript ./run_DA.r \
                            ${data_dir}/${data_id}_data_bm.RDS $method $seed $pop \
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
