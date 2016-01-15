../bin/cce: ../obj/misc.o ../obj/Spin.o ../obj/SpinData.o ../obj/SpinCollection.o ../obj/SpinSource.o ../obj/SpinClusterAlgorithm.o ../obj/SpinCluster.o ../obj/SpinInteraction.o ../obj/SpinInteractionComponent.o ../obj/SpinState.o ../obj/KronProd.o ../obj/QuantumOperator.o ../obj/HilbertSpaceOperator.o ../obj/LiouvilleSpaceOperator.o ../obj/QuantumState.o ../obj/PureState.o ../obj/MixedState.o ../obj/QuantumEvolution.o ../obj/QuantumEvolutionAlgorithm.o ../obj/cce.o
	g++ -g ../obj/cce.o ../obj/misc.o ../obj/Spin.o ../obj/SpinData.o ../obj/SpinCollection.o ../obj/SpinSource.o ../obj/SpinClusterAlgorithm.o ../obj/SpinCluster.o ../obj/SpinInteraction.o ../obj/SpinInteractionComponent.o ../obj/SpinState.o ../obj/KronProd.o ../obj/QuantumOperator.o ../obj/HilbertSpaceOperator.o ../obj/LiouvilleSpaceOperator.o ../obj/QuantumState.o ../obj/PureState.o ../obj/MixedState.o ../obj/QuantumEvolution.o ../obj/QuantumEvolutionAlgorithm.o -o  ../bin/cce -larmadillo -L/Applications/MATLAB_R2014b.app/bin/maci64 -lmx -lmat -std=c++11

../obj/Spin.o: source/spin/Spin/Spin.cpp include/spin/Spin.h
	g++ -g -c source/spin/Spin/Spin.cpp -o ../obj/Spin.o  -I.  -std=c++11

../obj/SpinData.o: source/spin/Spin/SpinData.cpp include/spin/SpinData.h
	g++ -g -c source/spin/Spin/SpinData.cpp -o ../obj/SpinData.o  -I.  -std=c++11

../obj/SpinCollection.o: source/spin/SpinCollection/SpinCollection.cpp include/spin/SpinCollection.h ../obj/SpinSource.o
	g++ -g -c source/spin/SpinCollection/SpinCollection.cpp -o ../obj/SpinCollection.o -I.  -std=c++11

../obj/SpinSource.o: source/spin/SpinCollection/SpinSource.cpp include/spin/SpinSource.h
	g++ -g -c source/spin/SpinCollection/SpinSource.cpp -o ../obj/SpinSource.o  -I.  -std=c++11

../obj/SpinCluster.o: source/spin/SpinCluster/SpinCluster.cpp include/spin/SpinCluster.h ../obj/SpinClusterAlgorithm.o
	g++ -g -c source/spin/SpinCluster/SpinCluster.cpp -o ../obj/SpinCluster.o  -I.  -std=c++11

../obj/SpinClusterAlgorithm.o: source/spin/SpinCluster/SpinClusterAlgorithm.cpp include/spin/SpinClusterAlgorithm.h
	g++ -g -c source/spin/SpinCluster/SpinClusterAlgorithm.cpp -o ../obj/SpinClusterAlgorithm.o  -I.  -std=c++11

../obj/SpinInteraction.o: source/spin/SpinInteraction/SpinInteraction.cpp include/spin/SpinInteraction.h
	g++ -g -c source/spin/SpinInteraction/SpinInteraction.cpp -o ../obj/SpinInteraction.o  -I.  -std=c++11

../obj/SpinInteractionComponent.o: source/spin/SpinInteraction/SpinInteractionComponent.cpp include/spin/SpinInteractionComponent.h include/misc/misc.h include/spin/SpinInteractionDefine.h
	g++ -g -c source/spin/SpinInteraction/SpinInteractionComponent.cpp -o ../obj/SpinInteractionComponent.o  -I.  -std=c++11

../obj/SpinState.o: source/spin/SpinInteraction/SpinState.cpp include/spin/SpinState.h
	g++ -g -c source/spin/SpinInteraction/SpinState.cpp -o ../obj/SpinState.o  -I.  -std=c++11

../obj/KronProd.o: source/kron/KronProd.cpp include/kron/KronProd.h include/spin/SpinInteractionDefine.h
	g++ -g -c source/kron/KronProd.cpp -o ../obj/KronProd.o  -I.  -std=c++11

../obj/QuantumOperator.o: source/quantum/QuantumOperator/QuantumOperator.cpp include/quantum/QuantumOperator.h
	g++ -g -c source/quantum/QuantumOperator/QuantumOperator.cpp -o ../obj/QuantumOperator.o -I. -I/Applications/MATLAB_R2014b.app/extern/include/  -std=c++11

../obj/HilbertSpaceOperator.o: source/quantum/QuantumOperator/HilbertSpaceOperator.cpp include/quantum/HilbertSpaceOperator.h
	g++ -g -c source/quantum/QuantumOperator/HilbertSpaceOperator.cpp -o ../obj/HilbertSpaceOperator.o -I. -I/Applications/MATLAB_R2014b.app/extern/include/  -std=c++11

../obj/LiouvilleSpaceOperator.o: source/quantum/QuantumOperator/LiouvilleSpaceOperator.cpp include/quantum/LiouvilleSpaceOperator.h
	g++ -g -c source/quantum/QuantumOperator/LiouvilleSpaceOperator.cpp -o ../obj/LiouvilleSpaceOperator.o -I. -I/Applications/MATLAB_R2014b.app/extern/include/  -std=c++11

../obj/QuantumState.o: source/quantum/QuantumState/QuantumState.cpp include/quantum/QuantumState.h
	g++ -g -c source/quantum/QuantumState/QuantumState.cpp -o ../obj/QuantumState.o -I. -I/Applications/MATLAB_R2014b.app/extern/include/  -std=c++11

../obj/PureState.o: source/quantum/QuantumState/PureState.cpp include/quantum/PureState.h
	g++ -g -c source/quantum/QuantumState/PureState.cpp -o ../obj/PureState.o -I. -I/Applications/MATLAB_R2014b.app/extern/include/  -std=c++11

../obj/MixedState.o: source/quantum/QuantumState/MixedState.cpp include/quantum/MixedState.h
	g++ -g -c source/quantum/QuantumState/MixedState.cpp -o ../obj/MixedState.o -I. -I/Applications/MATLAB_R2014b.app/extern/include/  -std=c++11

../obj/misc.o: source/misc/misc.cpp
	g++ -g -c source/misc/misc.cpp -o ../obj/misc.o  -I.   -std=c++11

../obj/QuantumEvolution.o: source/quantum/QuantumEvolution/QuantumEvolution.cpp include/quantum/QuantumEvolution.h
	g++ -g -c source/quantum/QuantumEvolution/QuantumEvolution.cpp -o ../obj/QuantumEvolution.o  -I. -I/Applications/MATLAB_R2014b.app/extern/include/  -std=c++11

../obj/QuantumEvolutionAlgorithm.o: source/quantum/QuantumEvolution/QuantumEvolutionAlgorithm.cpp include/quantum/QuantumEvolutionAlgorithm.h
	g++ -g -c source/quantum/QuantumEvolution/QuantumEvolutionAlgorithm.cpp -o ../obj/QuantumEvolutionAlgorithm.o  -I. -I/Applications/MATLAB_R2014b.app/extern/include/  -std=c++11

../obj/cce.o: app/cce.cpp
	g++ -g -c app/cce.cpp -o ../obj/cce.o  -I. -I/Applications/MATLAB_R2014b.app/extern/include/  -std=c++11

clean: 
	rm ../obj/*.o

doc:
	/Applications/Doxygen.app/Contents/Resources/doxygen ../doc/DocGen 

