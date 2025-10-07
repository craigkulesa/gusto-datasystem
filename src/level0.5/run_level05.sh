#!/usr/bin/env bash

cpus=1
erase=0
outDir='./build'

while getopts "ej:f:o:" flag
do
    case "${flag}" in
	e) erase=1;;
        j) cpus=${OPTARG};;
	f) list=${OPTARG};;
	o) outDir=${OPTARG};;
    esac
done

echo "# of CPUs: $cpus"
echo "file list: ${list}"
echo "Output base path: ${outDir}"

if [ ! -f $list ]; then
    echo "$list does not exist. Sorry."
    exit
fi

mkdir -p $outDir/B1
mkdir -p $outDir/B2

if [ $erase -eq 1 ]; then
    echo "Cleaning destination B1 and B2 directories..."
    rm -f ${outDir}/B1/*
    rm -f ${outDir}/B2/*
fi
rm -f x??

if [ $cpus -eq 1 ]; then
    echo "Running corrspec on $list in uniprocessor mode"
    time ./src/corrspec $list $outDir
else
    total_lines=$(wc -l <$list)
    ((lines_each = (total_lines + cpus - 1) / cpus))
    split -l $lines_each $list x
    for file in x??; do
	echo "Running corrspec on $file in parallel"
	time ./src/corrspec $file $outDir &
    done
fi
