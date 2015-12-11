cce.o: app/cce.cpp
	g++ -c app/cce.cpp -o ../obj/cce.o  -I. -std=c++11
	g++ -c source/spin/Spin/Spin.cpp -o ../obj/Spin.o  -I. -std=c++11
	g++ -c source/spin/Spin/SpinData.cpp -o ../obj/SpinData.o  -I. -std=c++11
	g++ -c source/spin/SpinCollection/SpinCollection.cpp -o ../obj/SpinCollection.o -I. -std=c++11
	g++ -c source/spin/SpinCollection/SpinSource.cpp -o ../obj/SpinSource.o  -I. -std=c++11
	g++ ../obj/cce.o ../obj/Spin.o ../obj/SpinData.o ../obj/misc.o ../obj/SpinCollection.o ../obj/SpinSource.o -o ../bin/cce -std=c++11
