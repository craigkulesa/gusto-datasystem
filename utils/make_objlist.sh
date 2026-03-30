#!/bin/bash

for i in level0.9/NII*fits; do 
	echo -n $i; ftlist $i[0] K include="OBJECT"; 
	done

