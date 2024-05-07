#!/bin/bash

mkdir -p Samples/ALEPHMCChargeMatchMetricTest

for i in Samples/ALEPHMC/*
do
   echo
   echo Processing file \"$i\"...
   echo

	./Execute --Input "${i}" --Output "${i/ALEPHMC/ALEPHMCChargeMatchMetricTest}" \
		--Gen tgen --Reco t --Fraction 1.0
done
