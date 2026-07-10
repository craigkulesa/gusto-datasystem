#!/bin/bash
# Create a list of objects an sequence / scan numbers for either CII or NII data
# the output of this script should be redirected to a file.list
for i in level0.9/CII*fits; do 
	echo -n $i; ftlist $i[0] K include="OBJECT"; 
	done

