#!/bin/bash

FILES=$( ls samples/ )

CLASS_LIST="../../analysis/classes_list.csv"
    
# tamanho dos dados por serie
STEP=5000
# STEP=10000

# DATASET_NAME="dataset_class_10k.dat"
DATASET_NAME="../classification/dataset_class_$STEP.dat"

if [ -r $DATASET_NAME ]
then
    rm $DATASET_NAME
fi

for i in ${FILES[@]}
do
    # nome do arquivo
    FILENAME=$(echo $i | sed 's/-60k.csv//')

    echo $FILENAME

    # obtendo a classe a partir do nome
    CLASS=$(cat $CLASS_LIST | grep $FILENAME\" | cut -d, -f2 | tr -d \")

    START=1
    END=$STEP

    while [ $END -le 60000 ]
    do

        echo $START $END
        
        DATA=$(sed -n "$START,$END p" samples/$i | tr '\n' ' ')

        # OS=$( echo $i | cut -d'-' -f1-2)
        # echo $DATA "\"$OS\"" | sed 's/ /,/g' # >> test3.dat
        
        echo $DATA "\"$FILENAME\"" $CLASS | sed 's/ /,/g' >> $DATASET_NAME
        # echo $DATA "\"$FILENAME\"" $CLASS | sed 's/ /,/g' 


        ((END+=$STEP))
        ((START+=$STEP))
    done
done
