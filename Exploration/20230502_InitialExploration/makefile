

default: Prepare TestRun

Prepare:
	ln -sf $(ProjectBase)/CommonCode/source/Dictionary_rdict.pcm .

TestRun: Execute
	./Execute --Input $(ProjectBase)/Samples/ALEPHMC/LEP1MC1994_recons_aftercut-014.root

Execute: FirstExploration.cpp
	g++ FirstExploration.cpp -o Execute \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
