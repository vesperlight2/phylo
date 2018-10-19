#!/bin/bash
project=$1
phylo_dir=$2
cd $phylo_dir
python phylophlan.py -c $project
python phylophlan.py --nproc 10 -u $project >> input/$project/phylo.log 2>&1
