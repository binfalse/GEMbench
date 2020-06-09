#!/bin/bash
#convert -density 150 GSE2109/INIT/clusterFR-NA-Zero.pdf -chop  0x30 GSE2109/INIT/clusterFR-NA-Zero.png


for i in $(find . -type f -name "*.pdf")
do
    echo converting $i to ${i/pdf/png}
    pdftoppm -png -singlefile $i ${i/.pdf/}
    convert ${i/pdf/png} -transparent white -chop 0x30 ${i/pdf/png}
done
