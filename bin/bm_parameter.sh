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
    R_methods=$(for m in daseq cydar meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolutions=$(for m in 0.1 0.2 1 2 5; do echo $m; done)
    betas=$(for m in 10 20 33 50 60; do echo $m; done)
    steps=$(for s in 30 40 50 60 70; do echo $s; done)
    tols=$(for t in 1.5 2 2.8 4 5; do echo $t; done)
    downsample=3
    mem=8g
elif [[ "$data_id" == "linear" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 7); do echo M$p; done)
    R_methods=$(for m in daseq cydar meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolutions=$(for m in 0.5 1 2 5 10; do echo $m; done)
    betas=$(for m in 30 50 71 90 110; do echo $m; done)
    steps=$(for s in 30 40 50 60 70; do echo $s; done)
    tols=$(for t in 1.5 2.3 3 4 5; do echo $t; done)
    downsample=3
    mem=8g
elif [[ "$data_id" == "branch" ]]
    then
    data_dir=${root}/data/synthetic/$data_id
    pops=$(for p in $(seq 1 1 8); do echo M$p; done)
    R_methods=$(for m in daseq cydar meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolutions=$(for m in 0.5 1 2 5 10; do echo $m; done)
    betas=$(for m in 30 50 65 80 100; do echo $m; done)
    steps=$(for s in 30 40 50 60 70; do echo $s; done)
    tols=$(for t in 1.5 2.4 3 4 5; do echo $t; done)
    downsample=3
    mem=8g
elif [[ "$data_id" == "covid19-pbmc" ]]
    then
    data_dir=${root}/data/real/$data_id
    pops=$(for m in RBC B PB CD14_Monocyte CD8_T CD4_T Platelet NK Granulocyte CD16_Monocyte gd_T pDC DC; do echo $m; done)
    R_methods=$(for m in daseq cydar meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolutions=$(for m in 0.2 0.5 1 2 5; do echo $m; done)
    betas=$(for m in 10 25 40 60; do echo $m; done)
    steps=$(for s in 30 40 50 60 70; do echo $s; done)
    tols=$(for t in 1.2 1.5 2.1 3 4; do echo $t; done)
    downsample=3
    mem=32g
elif [[ "$data_id" == "bcr-xl" ]]
    then
    data_dir=${root}/data/real/$data_id
    pops=$(for m in CD4_T-cells NK_cells CD8_T-cells B-cells_IgM+ monocytes surface- B-cells_IgM- DC; do echo $m; done)
    R_methods=$(for m in daseq cydar meld louvain; do echo $m; done)
    batch_vec=0
    k=30
    resolutions=$(for m in 0.2 0.6 1 2 5; do echo $m; done)
    betas=$(for m in 10 23 40 60; do echo $m; done)
    steps=$(for s in 30 40 50 60 70; do echo $s; done)
    tols=$(for t in 0.25 0.5 0.75 1 1.5; do echo $t; done)
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
                    jobid_prefix=${data_id}-${pop}-${pop_enr}-${batch_sd}-${method}-${seed}
                    if [[ "$method" == "meld" ]]; then
                        # enable meld env
                        source activate meld
                        meld_bin=$HOME/Documents/proj/benchmarkDA/methods/meld/bm_meld.py
                        for beta in $betas
                            do
                            mkdir -p ${root}/benchmark/${method}/${data_id}-beta=${beta}/
                            
                            jobid=${jobid_prefix}-beta=${beta}                           
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
                                --outdir ${root}/benchmark/${method}/${data_id}-beta=${beta}/"
                        done
                        conda deactivate
                    elif [[ "$method" == "louvain" ]]; then
                        for resolution in $resolutions
                            do
                            mkdir -p ${root}/benchmark/${method}/${data_id}-resolution=${resolution}/

                            jobid=${jobid_prefix}-resolution=${resolution}
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
                                --outdir ${root}/benchmark/${method}/${data_id}-resolution=${resolution}/"
                        done
                    elif [[ "$method" == "daseq" ]]; then
                        for step in $steps
                            do
                            mkdir -p ${root}/benchmark/${method}/${data_id}-step=${step}/

                            jobid=${jobid_prefix}-step=${step}
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
                                --step $step \
                                --downsample ${downsample} \
                                --batchEffect_sd $batch_sd \
                                --outdir ${root}/benchmark/${method}/${data_id}-step=${step}/"
                        done
                    elif [[ "$method" == "cydar" ]]; then
                        for tol in $tols
                            do
                            mkdir -p ${root}/benchmark/${method}/${data_id}-tol=${tol}/

                            jobid=${jobid_prefix}-tol=${tol}
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
                                --tol $tol \
                                --downsample ${downsample} \
                                --batchEffect_sd $batch_sd \
                                --outdir ${root}/benchmark/${method}/${data_id}-tol=${tol}/"
                        done
                    fi
                done
            done
        done
    done
done
