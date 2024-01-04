#!/bin/bash

for i in Samples/ALEPH/*1994*
do
   echo "Processing file $i..."

	./Execute --Input $i --Output ${i/ALEPH/ALEPHRecluster} \
		--Reco t --SkipGen true --JetR 0.2,0.4,0.6,0.8,1.0 --GhostSpacing 50
done

