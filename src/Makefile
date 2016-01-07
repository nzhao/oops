../bin/cce: ../obj/Spin.o ../obj/SpinData.o ../obj/SpinCollection.o ../obj/SpinSource.o ../obj/SpinClusterAlgorithm.o ../obj/SpinCluster.o ../obj/SpinInteraction.o ../obj/SpinInteractionComponent.o ../obj/KronProd.o ../obj/SpinSystem.o ../obj/QuantumOperator.o ../obj/HilbertSpaceOperator.o ../obj/LiouvilleSpaceOperator.o ../obj/cce.o
	g++ ../obj/cce.o ../obj/Spin.o ../obj/SpinData.o ../obj/SpinCollection.o ../obj/SpinSource.o ../obj/SpinClusterAlgorithm.o ../obj/SpinCluster.o ../obj/SpinInteraction.o ../obj/SpinInteractionComponent.o ../obj/KronProd.o ../obj/SpinSystem.o ../obj/QuantumOperator.o ../obj/HilbertSpaceOperator.o ../obj/LiouvilleSpaceOperator.o -o  ../bin/cce -std=c++11 -larmadillo

../obj/Spin.o: source/spin/Spin/Spin.cpp include/spin/Spin.h
	g++ -c source/spin/Spin/Spin.cpp -o ../obj/Spin.o  -I. -std=c++11

../obj/SpinData.o: source/spin/Spin/SpinData.cpp include/spin/SpinData.h
	g++ -c source/spin/Spin/SpinData.cpp -o ../obj/SpinData.o  -I. -std=c++11

../obj/SpinCollection.o: source/spin/SpinCollection/SpinCollection.cpp include/spin/SpinCollection.h ../obj/SpinSource.o
	g++ -c source/spin/SpinCollection/SpinCollection.cpp -o ../obj/SpinCollection.o -I. -std=c++11

../obj/SpinSource.o: source/spin/SpinCollection/SpinSource.cpp include/spin/SpinSource.h
	g++ -c source/spin/SpinCollection/SpinSource.cpp -o ../obj/SpinSource.o  -I. -std=c++11

../obj/SpinCluster.o: source/spin/SpinCluster/SpinCluster.cpp include/spin/SpinCluster.h ../obj/SpinClusterAlgorithm.o
	g++ -c source/spin/SpinCluster/SpinCluster.cpp -o ../obj/SpinCluster.o  -I. -std=c++11

../obj/SpinClusterAlgorithm.o: source/spin/SpinCluster/SpinClusterAlgorithm.cpp include/spin/SpinClusterAlgorithm.h
	g++ -c source/spin/SpinCluster/SpinClusterAlgorithm.cpp -o ../obj/SpinClusterAlgorithm.o  -I. -std=c++11

../obj/SpinInteraction.o: source/spin/SpinInteraction/SpinInteraction.cpp include/spin/SpinInteraction.h
	g++ -c source/spin/SpinInteraction/SpinInteraction.cpp -o ../obj/SpinInteraction.o  -I. -std=c++11

../obj/SpinInteractionComponent.o: source/spin/SpinInteraction/SpinInteractionComponent.cpp include/spin/SpinInteractionComponent.h include/misc/misc.h include/spin/SpinInteractionDefine.h
	g++ -c source/spin/SpinInteraction/SpinInteractionComponent.cpp -o ../obj/SpinInteractionComponent.o  -I. -std=c++11

../obj/KronProd.o: source/kron/KronProd.cpp include/kron/KronProd.h include/spin/SpinInteractionDefine.h
	g++ -c source/kron/KronProd.cpp -o ../obj/KronProd.o  -I. -std=c++11

../obj/SpinSystem.o: source/spin/SpinSystem/SpinSystem.cpp include/spin/SpinSystem.h include/spin/SpinInteractionDefine.h
	g++ -c source/spin/SpinSystem/SpinSystem.cpp -o ../obj/SpinSystem.o  -I. -std=c++11

../obj/QuantumOperator.o: source/quantum/QuantumOperator/QuantumOperator.cpp include/quantum/QuantumOperator.h
	g++ -c source/quantum/QuantumOperator/QuantumOperator.cpp -o ../obj/QuantumOperator.o -I. -std=c++11

../obj/HilbertSpaceOperator.o: source/quantum/QuantumOperator/HilbertSpaceOperator.cpp include/quantum/HilbertSpaceOperator.h
	g++ -c source/quantum/QuantumOperator/HilbertSpaceOperator.cpp -o ../obj/HilbertSpaceOperator.o -I. -std=c++11

../obj/LiouvilleSpaceOperator.o: source/quantum/QuantumOperator/LiouvilleSpaceOperator.cpp include/quantum/LiouvilleSpaceOperator.h
	g++ -c source/quantum/QuantumOperator/LiouvilleSpaceOperator.cpp -o ../obj/LiouvilleSpaceOperator.o -I. -std=c++11

../obj/cce.o: app/cce.cpp
	g++ -c app/cce.cpp -o ../obj/cce.o  -I. -std=c++11

clean: 
	rm ../obj/*.o

doc:
	/Applications/Doxygen.app/Contents/Resources/doxygen ../doc/DocGen 

