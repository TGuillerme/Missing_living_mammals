#!/bin/sh
##########################
#Fix the bionmial names of matrices in a folder according to a list
##########################
#SYNTAX :
#sh fix_binomial.sh <matrices> <changes>
#with
#<matrices> being the path to the folder containing the matrices
#<changes> being the list of elements to change in the matrices - must be a csv file
##########################
#guillert(at)tcd.ie - 03/02/2015
##########################

#INPUT

#Input values
matrices=$1
changes=$2

#number of lines to change
lines=$(wc -l $changes | sed 's/[[:space:]]//g' | sed 's/'"$changes"'//g')

for line in $(seq 1 $lines)
do
    #Splitting the different element of each line
    element_full=$(sed -n ''"$line"'p' $changes)
    matrix=$(echo $element_full | sed 's/[,].*$//')
    original=$(echo $element_full | sed 's/'"$matrix"',//' | sed 's/[,].*$//')
    replace=$(echo $element_full | sed 's/'"$matrix"','"$original"',//')

    #Replacing the $orignal name by the $replace in $matrix
    sed 's/'"$original"'/'"$replace"'/' $matrices/${matrix} > $matrices/${matrix}.tmp
    rm $matrices/${matrix}
    mv $matrices/${matrix}.tmp $matrices/${matrix}
    printf .
done

echo ""