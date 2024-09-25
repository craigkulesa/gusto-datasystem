#!/usr/bin/env bash

cpus=1
while getopts j:f: flag
do
    case "${flag}" in
        j) cpus=${OPTARG};;
	f) list=${OPTARG};;
    esac
done

echo "# of CPUs: $cpus"
echo "file list: $list"

if test -f "$list"
then
    echo "$list exists, good."
else
    echo "$list does not exist. Sorry."
    exit
fi

rm -f x??
if [ $cpus -eq 1 ]; then
    echo "Running corrspec on $list in uniprocessor mode"
    time ./corrspec $list 
else
    total_lines=$(wc -l <$list)
    ((lines_each = (total_lines + cpus - 1) / cpus))
    split -l $lines_each $list x
    for file in x??; do
	echo "Running corrspec on $file in parallel"
	time ./corrspec $file &
    done
fi
