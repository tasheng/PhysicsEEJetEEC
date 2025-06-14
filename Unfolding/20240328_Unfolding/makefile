default: TestRun

debug: Execute
	gdb --args ./Execute --Input $(ProjectBase)/Samples/ALEPHMC/skim_all.root \
		--Output skim_all_Matched.root \
		--Gen tgen --Reco t --Fraction 1.00
TestRun: Execute ExeUnfoldingHist ExeMatchingEffCorr
	./Execute --Input $(ProjectBase)/Samples/ALEPHMC/skim_all.root \
		--Output skim_all_Matched.root \
		--Gen tgen --Reco t --Fraction 1.00
	# ./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-001.root \
	# 	--Output LEP1MC1994_recons_aftercut-001_Matched.root \
	# 	--Gen tgen --Reco t --Fraction 1.00 --MatchingSchemeChoice 1
	# ./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-001.root \
	# 	--Output LEP1MC1994_recons_aftercut-001_Matched.root \
	# 	--Gen tgen --Reco t --Fraction 1.00 --MatchingSchemeChoice 3
	# ./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-001.root \
	# 	--Output LEP1MC1994_recons_aftercut-001_Matched.root \
	# 	--Gen tgen --Reco t --Fraction 1.00 --MatchingSchemeChoice 4
	# ./ExeUnfoldingHist
	# ./ExeMatchingEffCorr --Input matchingScheme2/skim_all_Matched.root \
	# 		     --MatchingEffName matchingScheme2/MatchingEff.root --MakeMatchingEffCorrFactor true
	# ./ExeMatchingEffCorr --Input matchingScheme2/LEP1MC1994_recons_aftercut-001_Matched.root \
	# 		     --MatchingEffName matchingScheme2/MatchingEff.root
	# ./ExeMatchingEffCorr --Input matchingScheme3/LEP1MC1994_recons_aftercut-001_Matched.root \
	# 		     --MatchingEffName matchingScheme3/MatchingEff.root --MakeMatchingEffCorrFactor true
	# ./ExeMatchingEffCorr --Input matchingScheme3/LEP1MC1994_recons_aftercut-001_Matched.root \
	# 		     --MatchingEffName matchingScheme3/MatchingEff.root
	# ./ExeMatchingEffCorr --Input matchingScheme4/LEP1MC1994_recons_aftercut-001_Matched.root \
	# 		     --MatchingEffName matchingScheme4/MatchingEff.root --MakeMatchingEffCorrFactor true
	# ./ExeMatchingEffCorr --Input matchingScheme4/LEP1MC1994_recons_aftercut-001_Matched.root \
	# 		     --MatchingEffName matchingScheme4/MatchingEff.root

FullRun: Execute
	bash MatchSome.sh

Execute: MatchEEC.cpp
	g++ MatchEEC.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

ExeUnfoldingHist: createUnfoldingHistograms.C
	g++ createUnfoldingHistograms.C -o ExeUnfoldingHist \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

ExeMatchingEffCorr: matchingEffCorr.cpp
	g++ matchingEffCorr.cpp -o ExeMatchingEffCorr \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
