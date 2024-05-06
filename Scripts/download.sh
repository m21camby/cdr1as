#!/bin/bash

while IFS= read -r line; do
    echo "Text read from file: $line"
    wget http://gtrd19-10.biouml.org/downloads/19.10/target_genes/Mus%20musculus/genes%20promoter%5b-500,+50%5d/$line
done < "$1"

#wget http://gtrd19-10.biouml.org/downloads/19.10/target_genes/Mus%20musculus/genes%20promoter%5b-500,+50%5d/



