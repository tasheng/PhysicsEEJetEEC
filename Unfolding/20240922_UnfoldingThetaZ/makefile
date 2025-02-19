Execute2DDoubleLog: RooUnfoldSimple_DoubleLogBins_ZTheta.cpp
	g++ RooUnfoldSimple_DoubleLogBins_ZTheta.cpp -o Execute2DDoubleLog \
		-I$(RooUnfoldBase)/src -L$(RooUnfoldBase)/ -lRooUnfold  -Wl,-rpath=$(RooUnfoldBase)/ \
		`root-config --cflags --glibs`


ExecuteUnfoldingStab: plotUnfoldingStab.cpp
	g++ plotUnfoldingStab.cpp -o ExecuteUnfoldingStab \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

Execute: fillPairTree.cpp
	g++ fillPairTree.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
