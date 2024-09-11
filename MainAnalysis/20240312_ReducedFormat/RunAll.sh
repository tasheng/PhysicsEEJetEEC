#!/bin/bash

for i in Samples/ALEPHMC/*
do
	./Execute --Input "$i" --Output `basename "$i"`
done

hadd -f LEP1MC1994_Reco.root LEP1MC1994_recons*root

for i in Samples/ALEPHMC/*
do
   rm `basename "$i"`
done

for i in Samples/ALEPHMC/*
do
	./Execute --Input "$i" --Output `basename "$i"` --GenLevel true --Tree tgen
done

hadd -f LEP1MC1994_Gen.root LEP1MC1994_recons*root

for i in Samples/ALEPHMC/*
do
   rm `basename "$i"`
done

for i in Samples/ALEPH/*1994*
do
	./Execute --Input "$i" --Output `basename "$i"`
done

mv LEP1MC1994_Reco.root LEP1MC1994_Gen.root LEP1Data1994_*root Output/
