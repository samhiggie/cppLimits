all: limit findmodel
limit:
	g++ limit.cc -o limit `root-config --cflags --glibs` -lRooFit -lRooFitCore -lboost_program_options -I ../../../.
	cp limit ~/.local/bin
findmodel:
	g++ findmodel.cc -o findmodel `root-config --cflags --glibs` -lRooFit -lRooFitCore -lboost_program_options -I ../../../.
	cp findmodel ~/.local/bin
export PATH := bin:$(PATH):$(HOME)/.local/bin
clean:
	$(RM) limit 
	$(RM) findmodel 

