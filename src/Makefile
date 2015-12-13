../bin/cce: ../obj/Spin.o ../obj/SpinData.o ../obj/SpinCollection.o ../obj/SpinSource.o ../obj/SpinGrouping.o ../obj/SpinCluster.o ../obj/cce.o
	g++ ../obj/cce.o ../obj/Spin.o ../obj/SpinData.o ../obj/SpinCollection.o ../obj/SpinSource.o ../obj/SpinGrouping.o ../obj/SpinCluster.o -o ../bin/cce -std=c++11

../obj/Spin.o: source/spin/Spin/Spin.cpp include/spin/Spin.h
	g++ -c source/spin/Spin/Spin.cpp -o ../obj/Spin.o  -I. -std=c++11

../obj/SpinData.o: source/spin/Spin/SpinData.cpp include/spin/SpinData.h
	g++ -c source/spin/Spin/SpinData.cpp -o ../obj/SpinData.o  -I. -std=c++11

../obj/SpinCollection.o: source/spin/SpinCollection/SpinCollection.cpp include/spin/SpinCollection.h
	g++ -c source/spin/SpinCollection/SpinCollection.cpp -o ../obj/SpinCollection.o -I. -std=c++11

../obj/SpinSource.o: source/spin/SpinCollection/SpinSource.cpp include/spin/SpinSource.h
	g++ -c source/spin/SpinCollection/SpinSource.cpp -o ../obj/SpinSource.o  -I. -std=c++11

../obj/SpinCluster.o: source/spin/SpinCluster/SpinCluster.cpp include/spin/SpinCluster.h
	g++ -c source/spin/SpinCluster/SpinCluster.cpp -o ../obj/SpinCluster.o  -I. -std=c++11

../obj/SpinGrouping.o: source/spin/SpinCluster/SpinGrouping.cpp include/spin/SpinGrouping.h
	g++ -c source/spin/SpinCluster/SpinGrouping.cpp -o ../obj/SpinGrouping.o  -I. -std=c++11

../obj/cce.o: app/cce.cpp
	g++ -c app/cce.cpp -o ../obj/cce.o  -I. -std=c++11

clean: 
	rm ../obj/*.o
