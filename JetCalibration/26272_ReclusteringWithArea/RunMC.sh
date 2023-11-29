#!/bin/bash

for i in Samples/ALEPHMC/*
do
   echo "Processing file $i..."

	./Execute --Input $i --Output ${i/ALEPHMC/ALEPHMCRecluster} \
		--Reco t --Gen tgen --GenBefore --JetR 0.2,0.4,0.6,0.8,1.0 --GhostSpacing 50
	# ./Execute --Input $i --Output ${i/ALEPHMC/ALEPHMCReclusterNoNu} \
	# 	--Reco t --Gen tgen --GenBefore --JetR 0.2,0.4,0.6,0.8,1.0 --SkipNeutrino true
done

