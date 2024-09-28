#!/usr/bin/env bash

cpus=1
clean=0

while getopts "cj:f:" flag
do
    case "${flag}" in
	c) clean=1;;
        j) cpus=${OPTARG};;
	f) list=${OPTARG};;
    esac
done

echo "# of CPUs: $cpus"
echo "file list: ${list}"

if [ ! -f $list ]; then
    echo "$list does not exist. Sorry."
    exit
fi

if [ $clean -eq 1 ]; then
    echo "Cleaning destination directory..."
    rm -f build/B1/*.fits
    rm -f build/B2/*.fits
fi
rm -f x??
sync

if [ $cpus -eq 1 ]; then
    echo "Running corrspec on $list in uniprocessor mode"
    time ./src/corrspec $list 
else
    total_lines=$(wc -l <$list)
    ((lines_each = (total_lines + cpus - 1) / cpus))
    split -l $lines_each $list x
    for file in x??; do
	echo "Running corrspec on $file in parallel"
	time ./src/corrspec $file &
    done
fi
