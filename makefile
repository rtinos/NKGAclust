NKGAclust : global.o ga_clustering.o util_functions.o file_man.o selection.o statistics.o 
	g++ -Wall global.o ga_clustering.o util_functions.o file_man.o selection.o statistics.o -o NKGAclust

file_man.o : file_man.cpp	
	g++ -Wall -o file_man.o -c file_man.cpp

ga_clustering.o : ga_clustering.cpp	
	g++ -Wall -o ga_clustering.o -c ga_clustering.cpp

global.o : global.cpp	
	g++ -Wall -o global.o -c global.cpp

selection.o : selection.cpp	
	g++ -Wall -o selection.o -c selection.cpp

statistics.o : statistics.cpp	
	g++ -Wall -o statistics.o -c statistics.cpp

util_functions.o : util_functions.cpp	
	g++ -Wall -o util_functions.o -c util_functions.cpp