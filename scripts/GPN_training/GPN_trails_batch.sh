#!/bin/bash
#SBATCH --account=johan_gpn
#SBATCH --output=gpn_training.out
#SBATCH --error=gpn_training.err
#SBATCH --time=96:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:4
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
# Activate the conda environment
source /home/johanulstrup/miniconda3/etc/profile.d/conda.sh
conda activate GPN

# Run the training script
WANDB_PROJECT=baboon_diversity torchrun --nproc_per_node=$(echo $CUDA_VISIBLE_DEVICES | awk -F',' '{print NF}') -m gpn.ss.run_mlm --do_train --do_eval \
    --report_to wandb --prediction_loss_only True --remove_unused_columns False \
    --dataset_name /home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/steps/dataset_assembly \
    --tokenizer_name /home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/scripts/GPN_training/tokenizer-dna-mlm \
    --soft_masked_loss_weight_train 0.1 --soft_masked_loss_weight_evaluation 0.0 \
    --weight_decay 0.01 --optim adamw_torch \
    --dataloader_num_workers 16 --seed 42 \
    --save_strategy steps --save_steps 10000 --evaluation_strategy steps \
    --eval_steps 10000 --logging_steps 10000 --max_steps 120000 --warmup_steps 1000 \
    --learning_rate 1e-3 --lr_scheduler_type constant_with_warmup \
    --run_name your_run --output_dir /home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/model --model_type GPN \
    --per_device_train_batch_size 64 --per_device_eval_batch_size 64 --gradient_accumulation_steps 1 --total_batch_size 256 \
    --torch_compile \
    --ddp_find_unused_parameters False \
    --bf16 --bf16_full_eval 2>&1 | tee /faststorage/project/johan_gpn/people/johanulsrup/johan_gpn/scripts/GGPN_training/training_log.txt
