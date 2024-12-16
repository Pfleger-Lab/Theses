#!/bin/bash

source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh
unset PYTHONPATH
conda activate /home/glbrc.org/jdietrich3/.conda/envs/biopython
python fastq_library_seq_complete_determ_only_Scarcity_v4.py 2023_0725_unsorted_lib_pass.fastq
conda deactivate