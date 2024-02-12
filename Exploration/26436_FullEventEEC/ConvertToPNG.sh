#!/bin/bash

for i in *pdf
do
   echo $i
   convert -trim -density 300 "$i" \
      -quality 100 "${i/pdf/png}"
done


