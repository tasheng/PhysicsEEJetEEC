#!/bin/bash

for i in `seq 1 40`
do
   j=`printf "%02d" "$i"`
	./Execute --Input $ProjectBase/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-0$j.root \
		--Output PlotGen_$j.root \
		--Particle tgen --IsReco false --DoEENormalize true
	./Execute --Input $ProjectBase/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-0$j.root \
		--Output PlotReco_$j.root \
		--Particle t --IsReco true --DoEENormalize true
done

./Execute --Input $ProjectBase/Samples/ALEPH/LEP1Data1994P1_recons_aftercut-MERGED.root \
   --Output PlotData_1.root \
   --Particle t --IsReco true --DoEENormalize true
./Execute --Input $ProjectBase/Samples/ALEPH/LEP1Data1994P2_recons_aftercut-MERGED.root \
   --Output PlotData_2.root \
   --Particle t --IsReco true --DoEENormalize true
./Execute --Input $ProjectBase/Samples/ALEPH/LEP1Data1994P3_recons_aftercut-MERGED.root \
   --Output PlotData_3.root \
   --Particle t --IsReco true --DoEENormalize true
