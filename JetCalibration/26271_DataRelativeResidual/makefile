
DHFile = --DHFile DijetResult.dh
DHFileR = --DHFile DijetResultR.dh

default: Prepare TestRun TestRunSummary TestRunExport

Execute: DijetResidual.cpp
	CompileRootMacro DijetResidual.cpp Execute

Prepare:
	ln -sf ../26271_DataMCComparison/Samples/MCAllR4.root .
	ln -sf ../26271_DataMCComparison/Samples/MCAllR10.root .
	ln -sf ../26271_DataMCComparison/Samples/DataAllR4.root .
	ln -sf ../26271_DataMCComparison/Samples/DataAllR10.root .
	ln -sf ../26271_DataMCComparison/Samples/DataRAllR4.root .
	ln -sf ../26271_DataMCComparison/Samples/DataRAllR10.root .

TestRun: Execute
	./Execute --File MCAllR4.root     --E3Cut 5  --PdfFile ResidualSideMCR4.pdf       $(DHFile)  --State MCR4_5
	./Execute --File MCAllR10.root    --E3Cut 5  --PdfFile ResidualSideMCR10.pdf      $(DHFile)  --State MCR10_5
	./Execute --File DataAllR4.root   --E3Cut 5  --PdfFile ResidualSideDataR4.pdf     $(DHFile)  --State DataR4_5
	./Execute --File DataAllR10.root  --E3Cut 5  --PdfFile ResidualSideDataR10.pdf    $(DHFile)  --State DataR10_5
	./Execute --File DataAllR4.root   --E3Cut 10 --PdfFile ResidualSideDataR4X.pdf    $(DHFile)  --State DataR4_10
	./Execute --File DataAllR10.root  --E3Cut 10 --PdfFile ResidualSideDataR10X.pdf   $(DHFile)  --State DataR10_10
	./Execute --File DataAllR4.root   --E3Cut 3  --PdfFile ResidualSideDataR4X2.pdf   $(DHFile)  --State DataR4_3
	./Execute --File DataAllR10.root  --E3Cut 3  --PdfFile ResidualSideDataR10X2.pdf  $(DHFile)  --State DataR10_3
	./Execute --File DataRAllR4.root  --E3Cut 3  --PdfFile ResidualSideDataRR4X2.pdf  $(DHFileR) --State DataR4_3
	./Execute --File DataRAllR10.root --E3Cut 3  --PdfFile ResidualSideDataRR10X2.pdf $(DHFileR) --State DataR10_3
	# cp *pdf ~/WindowsHome/Downloads/

ExecuteSummary: GrandSummaryPlots.cpp
	CompileRootMacro GrandSummaryPlots.cpp ExecuteSummary

TestRunSummary: ExecuteSummary
	./ExecuteSummary --DHFile DijetResult.dh --Output ResidualSideSummaryR4.pdf \
		--Labels "Data E_{3} < 10 GeV","Data E_{3} < 5 GeV","Data E_{3} < 3 GeV","MC E_{3} < 5 GeV" \
		--States "DataR4_10","DataR4_5","DataR4_3","MCR4_5" \
		--Colors 4,1,2,6 \
		--XMin 0.5 --XMax 0.8 --YMin -0.7 --YMax 0.7
	./ExecuteSummary --DHFile DijetResult.dh --Output ResidualSideSummaryR10.pdf \
		--Labels "Data E_{3} < 10 GeV","Data E_{3} < 5 GeV","Data E_{3} < 3 GeV","MC E_{3} < 5 GeV" \
		--States "DataR10_10","DataR10_5","DataR10_3","MCR10_5" \
		--Colors 4,1,2,6 \
		--XMin 0.5 --XMax 0.8 --YMin -0.7 --YMax 0.7
	# cp *pdf ~/WindowsHome/Downloads/

ExecuteExport: ExportTextFile.cpp
	CompileRootMacro ExportTextFile.cpp ExecuteExport

TestRunExport: ExecuteExport
	./ExecuteExport --DHFile DijetResult.dh --State DataR4_3 --Output JECRelativeResidualR4.txt
	./ExecuteExport --DHFile DijetResult.dh --State DataR10_3 --Output JECRelativeResidualR10.txt

Deploy:
	# mkdir -p ../../JetCorrection/23475_DataRel/
	# cp JECRelativeResidualR4.txt ../../JetCorrection/23475_DataRel/
	# cp JECRelativeResidualR10.txt ../../JetCorrection/23475_DataRel/
