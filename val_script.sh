#!/bin/bash

DIR=$1
NUM=$2
N_DECODE=$3
total_n=$4

for ((i=1; i<NUM; i++)); do
    f=$DIR/model.iter-$i
    if [ -e $f ]; then
        echo $f
        python ~/iclr19-graph2graph/diff_vae/decode.py --num_decode $N_DECODE --test $DIR/../data/test.txt --vocab $DIR/../data/vocab.txt --model $f --use_molatt | python ~/CS6250_project/scripts/bg_score.py $f > $DIR/results.$i
        python ~/CS6250_project/scripts/bg_analyze.py --num_decode $N_DECODE --sim_delta .2 --prop_delta 6 --total_n $total_n < $DIR/results.$i
    fi
done