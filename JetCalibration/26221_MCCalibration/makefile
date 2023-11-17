
default: RunJECClosure

Prepare:
	ln -sf ../../Samples/ALEPHMC Samples

FullRun: Execute
	for j in 2 4 6 8 10; \
	do \
		for i in `ls Samples/ALEPHMCRecluster | cut -b 28-30`; \
		do \
			echo $$i; \
			./Execute --FileName Samples/ALEPHMCRecluster/LEP1MC1994_recons_aftercut-$${i}.root \
				--Output Samples/MatchedR$${j}$${i}.root \
				--RecoTree RecoR$${j} \
				--GenTree GenR$${j} \
				--JEC JECR$${j}.txt \
				--Fraction 1.00; \
			./Execute --FileName Samples/ALEPHMCRecluster/LEP1MC1994_recons_aftercut-$${i}.root \
				--Output Samples/ReverseMatchedR$${j}$${i}.root \
				--RecoTree GenR$${j} \
				--GenTree RecoR$${j} \
				--JEC JECR$${j}.txt \
				--Fraction 1.00; \
		done; \
		hadd -f Samples/AllMatchedR$${j}.root Samples/MatchedR$${j}0*.root; \
		hadd -f Samples/AllReverseMatchedR$${j}.root Samples/ReverseMatchedR$${j}0*.root; \
		rm Samples/MatchedR$${j}0*.root Samples/ReverseMatchedR$${j}0*.root; \
	done

TestRun: Execute
	./Execute --FileName Samples/ALEPHMCRecluster/LEP1MC1994_recons_aftercut-011.root \
		--Output Matched011.root \
		--RecoTree RecoR4 \
		--GenTree GenR4 \
		--JEC JECR4.txt \
		--Fraction 1.00

Execute: JetMatcher.cpp
	CompileRootMacro JetMatcher.cpp Execute

TestRunPlot: ExecutePlot
	./ExecutePlot
	cp JetResponse.pdf ~/WindowsHome/Downloads/

ExecutePlot: JetExploration.cpp
	CompileRootMacro JetExploration.cpp ExecutePlot

ExecuteJEC: JetEnergyCorrection.cpp
	CompileRootMacro JetEnergyCorrection.cpp ExecuteJEC

RunJEC: ExecuteJEC
	time ./ExecuteJEC --File Samples/AllMatchedR2.root  --Tree MatchedTree --Output JECR2.txt --Check JECCheckR2.pdf --R 0.2
	time ./ExecuteJEC --File Samples/AllMatchedR4.root  --Tree MatchedTree --Output JECR4.txt --Check JECCheckR4.pdf --R 0.4
	time ./ExecuteJEC --File Samples/AllMatchedR6.root  --Tree MatchedTree --Output JECR6.txt --Check JECCheckR6.pdf --R 0.6
	time ./ExecuteJEC --File Samples/AllMatchedR8.root  --Tree MatchedTree --Output JECR8.txt --Check JECCheckR8.pdf --R 0.8
	time ./ExecuteJEC --File Samples/AllMatchedR10.root --Tree MatchedTree --Output JECR10.txt --Check JECCheckR10.pdf --R 1.0
	# cp JECCheck*.pdf ~/WindowsHome/Downloads

ExecuteJECClosure: JetEnergyClosure.cpp
	CompileRootMacro JetEnergyClosure.cpp ExecuteJECClosure

RunJECClosure: ExecuteJECClosure
	time ./ExecuteJECClosure --Name Samples/AllMatchedR2.root --Output JECClosureR2.pdf --R 0.2
	time ./ExecuteJECClosure --Name Samples/AllMatchedR4.root --Output JECClosureR4.pdf --R 0.4
	time ./ExecuteJECClosure --Name Samples/AllMatchedR6.root --Output JECClosureR6.pdf --R 0.6
	time ./ExecuteJECClosure --Name Samples/AllMatchedR8.root --Output JECClosureR8.pdf --R 0.8
	time ./ExecuteJECClosure --Name Samples/AllMatchedR10.root --Output JECClosureR10.pdf --R 1.0
	# cp JECClosure*.pdf ~/WindowsHome/Downloads

Deploy:
	cp JEC*txt ../../JetCorrection/23474_MC/

