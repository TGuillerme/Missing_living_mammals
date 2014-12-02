#!/bin/sh
##########################
#Selecting files in a folder matching a given list of names
##########################
#SYNTAX :
#sh select.files.sh folder list
#with
#<folder>
#<list>
##########################
#Selecting all the files containing at least one element of the provided list
#version: 0.1
select_files_version="select.files.sh v0.1"
#Update: The tree shape and the fossil/living distribution is fully randomized.
#----
#guillert(at)tcd.ie - 02/12/2014
##########################

#Input
folder=$1
list=$2

#Make the empty folder where to save the results
mkdir selected_files

#Number of elements in the list
numbers=$(cat $list | wc -l | sed 's/[[:space:]]//g')

#Initializing the loop
echo "Searching the $numbers terms:"
for n in $(seq 1 $numbers)
do
    #Selecting the n element of the list
    element=$(sed -n ''"$n"'p' $list)
    #Recursive grep within element within folder
    grep -lr "$element" $folder >> file.to.copy
    #Verbose
    printf .
done
echo "Done"

#keeping only the unique files and remove the spaces
sort file.to.copy | uniq -u | sed 's/ /\\ /g' > files.to.copy

#remove the file with extra names
rm file.to.copy

#Selecting the number of files to copy
files=$(cat files.to.copy | wc -l | sed 's/[[:space:]]//g')

#Looping through the files to copy
echo "Copying the $files files:"
for n in $(seq 1 $files)
do
    #Selecting the file to copy
    file=$(sed -n ''"$n"'p' files.to.copy)
    #Copying the file in the selected_file folder
    cp "$file" selected_files
    #verbose
    printf .
done
echo "Done"

#remove the file with names
rm files.to.copy

#End
