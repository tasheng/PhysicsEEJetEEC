default: ProductionRun

Execute: JECFitter.cpp
	CompileRootMacro JECFitter.cpp Execute

Prepare:
	ln -sf ../26271_DataMCComparison/Samples/DataAllR4.root .
	ln -sf ../26271_DataMCComparison/Samples/DataAllR10.root .
	ln -sf ../26271_DataMCComparison/Samples/DataRAllR4.root .
	ln -sf ../26271_DataMCComparison/Samples/DataRAllR10.root .
	ln -sf ../26271_DataMCComparison/Samples/DataRRAllR4.root .
	ln -sf ../26271_DataMCComparison/Samples/DataRRAllR10.root .
	ln -sf ../26271_DataMCComparison/Samples/MCAllR4.root .
	ln -sf ../26271_DataMCComparison/Samples/MCAllR10.root .

ProductionRun: Execute
	# time ./Execute --MC MCAllR4.root --Data DataRAllR4.root --Validation FitValidationR4.pdf \
	# 	--NJet 9 --PCut 3 --State R4_9_3 --Percentage 1.00
	time ./Execute --MC MCAllR10.root --Data DataRAllR10.root --Validation FitValidationR10.pdf \
		--NJet 9 --PCut 3 --State R10_9_3 --Percentage 1.00 --ThetaMin 0.30 --ThetaMax 0.70
	# cp FitValidationR[48].pdf ~/WindowsHome/Downloads/

RunExport: MakeTextFile.cpp
	CompileRootMacro MakeTextFile.cpp ExecuteTextFile
	./ExecuteTextFile --DHFile FitResults.dh --Output JECResidualR4.txt --Order 1 --State R4_9_3
	./ExecuteTextFile --DHFile FitResults.dh --Output JECResidualR10.txt --Order 1 --State R10_9_3

Deploy:
	# mkdir -p ../../JetCorrection/23475_DataAbs/
	# cp JECResidualR4.txt ../../JetCorrection/23475_DataAbs/
	# cp JECResidualR10.txt ../../JetCorrection/23475_DataAbs/
