#!/bin/bash

FILES=$( ls switches/ )

for i in ${FILES[@]}
do
    START=1
    END=10000
    while [ $END -le 60000 ]
    do

        echo $START $END
        DATA=$(sed -n "$START,$END p" switches/$i | tr '\n' ' ')

        OS=$( echo $i | cut -d'-' -f1-2)

        echo $DATA "\"$OS\"" >> switches_dataset3.csv 

        ((END+=10000))
        ((START+=10000))
    done
done
