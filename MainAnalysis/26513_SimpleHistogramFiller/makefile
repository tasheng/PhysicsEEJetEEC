

default: TestRun

TestRun: Execute
	./Execute

Execute: HistogramFiller.cpp
	g++ HistogramFiller.cpp -o Execute `root-config --cflags --libs` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
