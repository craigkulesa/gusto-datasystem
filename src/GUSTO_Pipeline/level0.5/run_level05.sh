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

mkdir -p $outDir

if [ $erase -eq 1 ]; then
    echo "Cleaning destination level 0.5 directory..."
    rm -rf ${outDir}
    mkdir -p ${outDir}
fi
rm -f x??

cd src
make
cd ..

if [ $cpus -eq 1 ]; then
    echo "Running corrspec on $list in uniprocessor mode"
    ./src/corrspec $list $outDir
else
    total_lines=$(wc -l <$list)
    ((lines_each = (total_lines + cpus - 1) / cpus))
    split -l $lines_each $list x
    for file in x??; do
	echo "Running corrspec on $file in parallel"
	./src/corrspec $file $outDir &
    done
    # because this returns immediately, hold in a loop and wait for processes to stop
    # before returning to pipeline or shell prompt
    while true; do
	if pgrep -x "corrspec" > /dev/null; then
            echo "$(date): corrspec is still running..."
	else
            echo "$(date): corrspec is done."
	    break
	fi
	sleep 5
    done
fi
