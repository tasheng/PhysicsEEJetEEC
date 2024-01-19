#!/bin/bash

mkdir -p Samples/ALEPHMCChargeMatch

for i in Samples/ALEPHMC/*
do
   echo
   echo Processing file \"$i\"...
   echo

	./Execute --Input "${i}" --Output "${i/ALEPHMC/ALEPHMCChargeMatch}" \
		--Gen tgen --Reco t --Fraction 1.00
done
