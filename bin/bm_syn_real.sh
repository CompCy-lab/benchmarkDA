#!/bin/bash

# add slurm module first
module purge
module load r/4.1.3

# set slurm parameters
time=4:00:00
partition=general

data_id=$1
root=$HOME/Documents/proj/benchmarkDA
cd ${root}/scripts

if [[ "$data_id" == "cluster" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain milo_batch cydar_batch cna_batch louvain_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=0.2
    beta=33
    downsample=3
    mem=8g
elif [[ "$data_id" == "cluster_balanced" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 3); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain milo_batch cydar_batch cna_batch louvain_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=0.2
    beta=33
    downsample=3
    mem=8g
elif [[ "$data_id" == "linear" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 7); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain milo_batch cydar_batch cna_batch louvain_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=1
    beta=71
    downsample=3
    mem=8g
elif [[ "$data_id" == "branch" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 8); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain milo_batch cydar_batch cna_batch louvain_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.75 1 1.25 1.5; do echo $m; done)
    k=30
    resolution=1
    beta=65
    downsample=3
    mem=8g
elif [[ "$data_id" == "covid19-pbmc" ]]
    then
    data_dir=${root}/data/real/$data_id
    pops=$(for m in RBC B PB CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC; do echo $m; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.5
    beta=25
    downsample=3
    mem=32g
elif [[ "$data_id" == "bcr-xl" ]]
    then
    data_dir=${root}/data/real/$data_id
    pops=$(for m in CD4_T-cells NK_cells CD8_T-cells B-cells_IgM+ monocytes surface- B-cells_IgM- DC; do echo $m; done)
    R_methods=$(for m in milo daseq cydar cna meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolution=0.6
    beta=23
    downsample=10
    mem=32g
fi


## Run
for pop in $pops
    do
    for pop_enr in $(seq 0.75 0.1 0.95)
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
