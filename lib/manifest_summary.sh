#!/bin/bash

# Makes a table with columns "#rname" and "genome name".
# This is used to look up genome names to output from
# competitive mapping.

manifest_path=$1
summary_path=$2

#Column headers
printf "\"#rname\"\t\"genome_name\"\n" > $summary_path

#Skip >, then grab everything up to the first space, then everything up to the end of the line
regex='(>)([^( )]+)[ ](.+)'

#Extract and write out rname and genome name as .tsv
while read line
do
    if [[ $line =~ $regex ]]
    then
        printf "\"${BASH_REMATCH[2]}\"\t\"${BASH_REMATCH[3]}\"\n" >> $summary_path
    fi
done < <(cat $manifest_path | grep ">")
