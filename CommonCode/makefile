
default: all

all: prepare library/Messenger.o library/BasicUtilities.o library/CATree.o library/Dictionary.o 

prepare:
	mkdir -p library/

library/BasicUtilities.o: source/BasicUtilities.cpp include/BasicUtilities.h
	g++ source/BasicUtilities.cpp -I. -c -o library/BasicUtilities.o -I${RootMacrosBase}/ -std=c++11

library/CATree.o: source/CATree.cpp include/CATree.h
	g++ source/CATree.cpp -I. -c -o library/CATree.o -I${RootMacrosBase}/ -std=c++11

library/Dictionary.o: include/Dictionary.h include/DictionaryObject.h
	rootcint -f source/Dictionary.cxx -c include/DictionaryObject.h include/Dictionary.h
	g++ `root-config --cflags` source/Dictionary.cxx -o library/Dictionary.o -I. -c -fpic

library/Messenger.o: source/Messenger.cpp include/Messenger.h
	g++ source/Messenger.cpp -Iinclude -c -o library/Messenger.o `root-config --cflags` -std=c++17

