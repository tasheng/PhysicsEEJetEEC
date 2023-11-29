


default: Prepare FullRun RunComparison

Prepare:
	ln -sf ../../Samples/ALEPHRecluster DataSamples
	ln -sf ../../Samples/ALEPHMCRecluster MCSamples
	ln -sf ../../Samples/DataMCComparison Samples
	ln -sf ../../JetCalibration/26221_MCCalibration/*txt .
	ln -sf ../../JetCalibration/26271_DataRelativeResidual/*txt .
	# ln -sf ../../JetCorrection/23475_DataAbs/*txt .

FullRun: Execute
	./Execute --Input DataSamples/LEP1Data1994P1_recons_aftercut-MERGED.root --Output Samples/DataP1R4.root \
		--Tree RecoR4 --JEC JECR4.txt
	./Execute --Input DataSamples/LEP1Data1994P2_recons_aftercut-MERGED.root --Output Samples/DataP2R4.root \
		--Tree RecoR4 --JEC JECR4.txt
	./Execute --Input DataSamples/LEP1Data1994P3_recons_aftercut-MERGED.root --Output Samples/DataP3R4.root \
		--Tree RecoR4 --JEC JECR4.txt
	./Execute --Input DataSamples/LEP1Data1994P1_recons_aftercut-MERGED.root --Output Samples/DataP1R10.root \
		--Tree RecoR10 --JEC JECR10.txt
	./Execute --Input DataSamples/LEP1Data1994P2_recons_aftercut-MERGED.root --Output Samples/DataP2R10.root \
		--Tree RecoR10 --JEC JECR10.txt
	./Execute --Input DataSamples/LEP1Data1994P3_recons_aftercut-MERGED.root --Output Samples/DataP3R10.root \
		--Tree RecoR10 --JEC JECR10.txt
	./Execute --Input DataSamples/LEP1Data1994P1_recons_aftercut-MERGED.root --Output Samples/DataRP1R4.root \
		--Tree RecoR4 --JEC JECR4.txt,JECRelativeResidualR4.txt
	./Execute --Input DataSamples/LEP1Data1994P2_recons_aftercut-MERGED.root --Output Samples/DataRP2R4.root \
		--Tree RecoR4 --JEC JECR4.txt,JECRelativeResidualR4.txt
	./Execute --Input DataSamples/LEP1Data1994P3_recons_aftercut-MERGED.root --Output Samples/DataRP3R4.root \
		--Tree RecoR4 --JEC JECR4.txt,JECRelativeResidualR4.txt
	./Execute --Input DataSamples/LEP1Data1994P1_recons_aftercut-MERGED.root --Output Samples/DataRP1R10.root \
		--Tree RecoR10 --JEC JECR10.txt,JECRelativeResidualR10.txt
	./Execute --Input DataSamples/LEP1Data1994P2_recons_aftercut-MERGED.root --Output Samples/DataRP2R10.root \
		--Tree RecoR10 --JEC JECR10.txt,JECRelativeResidualR10.txt
	./Execute --Input DataSamples/LEP1Data1994P3_recons_aftercut-MERGED.root --Output Samples/DataRP3R10.root \
		--Tree RecoR10 --JEC JECR10.txt,JECRelativeResidualR10.txt
	./Execute --Input DataSamples/LEP1Data1994P1_recons_aftercut-MERGED.root --Output Samples/DataRRP1R4.root \
		--Tree RecoR4 --JEC JECR4.txt,JECRelativeResidualR4.txt,JECResidualR4.txt
	./Execute --Input DataSamples/LEP1Data1994P2_recons_aftercut-MERGED.root --Output Samples/DataRRP2R4.root \
		--Tree RecoR4 --JEC JECR4.txt,JECRelativeResidualR4.txt,JECResidualR4.txt
	./Execute --Input DataSamples/LEP1Data1994P3_recons_aftercut-MERGED.root --Output Samples/DataRRP3R4.root \
		--Tree RecoR4 --JEC JECR4.txt,JECRelativeResidualR4.txt,JECResidualR4.txt
	./Execute --Input DataSamples/LEP1Data1994P1_recons_aftercut-MERGED.root --Output Samples/DataRRP1R10.root \
		--Tree RecoR10 --JEC JECR10.txt,JECRelativeResidualR10.txt,JECResidualR10.txt
	./Execute --Input DataSamples/LEP1Data1994P2_recons_aftercut-MERGED.root --Output Samples/DataRRP2R10.root \
		--Tree RecoR10 --JEC JECR10.txt,JECRelativeResidualR10.txt,JECResidualR10.txt
	./Execute --Input DataSamples/LEP1Data1994P3_recons_aftercut-MERGED.root --Output Samples/DataRRP3R10.root \
		--Tree RecoR10 --JEC JECR10.txt,JECRelativeResidualR10.txt,JECResidualR10.txt
	for i in `ls MCSamples/ | cut -b 28-30`; do \
		echo $$i; \
		./Execute --Input MCSamples/LEP1MC1994_recons_aftercut-$$i.root --Output Samples/MCP$${i}R4.root \
			--Tree RecoR4 --JEC JECR4.txt; \
		./Execute --Input MCSamples/LEP1MC1994_recons_aftercut-$$i.root --Output Samples/MCP$${i}R10.root \
			--Tree RecoR10 --JEC JECR10.txt; \
	done
	hadd -f Samples/DataRRAllR4.root  Samples/DataRRP[123]R4.root
	hadd -f Samples/DataRRAllR10.root Samples/DataRRP[123]R10.root
	hadd -f Samples/DataRAllR4.root   Samples/DataRP[123]R4.root
	hadd -f Samples/DataRAllR10.root  Samples/DataRP[123]R10.root
	hadd -f Samples/DataAllR4.root    Samples/DataP[123]R4.root
	hadd -f Samples/DataAllR10.root   Samples/DataP[123]R10.root
	hadd -f Samples/MCAllR4.root      Samples/MCP*R4.root
	hadd -f Samples/MCAllR10.root     Samples/MCP*R10.root
	rm Samples/DataP[123]R[48].root Samples/DataRP[123]R[48].root Samples/DataRRP[123]R[48].root
	rm Samples/MCP*R[48].root

TestRun: Execute
	./Execute --Input DataSamples/LEP1Data1994P1_recons_aftercut-MERGED.root \
		--Output DataP1.root \
		--Tree RecoR4 \
		--Fraction 0.05 --JEC JECR4.txt

Execute: ProcessTree.cpp
	CompileRootMacro ProcessTree.cpp Execute -I../../CommonCode/include

ExecuteComparison: CompareSamples.cpp
	CompileRootMacro CompareSamples.cpp ExecuteComparison

RunComparison: ExecuteComparison
	./ExecuteComparison --Data DataAll.root --MC MCAll.root --Output DataMCComparison.pdf
	./ExecuteComparison --Data DataAllFull.root --MC MCAll.root --Output DataMCComparisonFull.pdf
	cp DataMCComparison*.pdf ~/WindowsHome/Downloads/
