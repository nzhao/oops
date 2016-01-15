../bin/cce: ../obj/misc.o ../obj/Spin.o ../obj/SpinData.o ../obj/SpinCollection.o ../obj/SpinSource.o ../obj/SpinClusterAlgorithm.o ../obj/SpinCluster.o ../obj/SpinInteraction.o ../obj/SpinInteractionComponent.o ../obj/SpinState.o ../obj/KronProd.o ../obj/QuantumOperator.o ../obj/HilbertSpaceOperator.o ../obj/LiouvilleSpaceOperator.o ../obj/QuantumState.o ../obj/PureState.o ../obj/MixedState.o ../obj/QuantumEvolution.o ../obj/QuantumEvolutionAlgorithm.o ../obj/cce.o
	g++ -DARMA_DONT_USE_CXX11 -g ../obj/cce.o ../obj/misc.o ../obj/Spin.o ../obj/SpinData.o ../obj/SpinCollection.o ../obj/SpinSource.o ../obj/SpinClusterAlgorithm.o ../obj/SpinCluster.o ../obj/SpinInteraction.o ../obj/SpinInteractionComponent.o ../obj/SpinState.o ../obj/KronProd.o ../obj/QuantumOperator.o ../obj/HilbertSpaceOperator.o ../obj/LiouvilleSpaceOperator.o ../obj/QuantumState.o ../obj/PureState.o ../obj/MixedState.o ../obj/QuantumEvolution.o ../obj/QuantumEvolutionAlgorithm.o -o  ../bin/cce -DARMA_DONT_USE_WRAPPER -lblas -llapack #-larmadillo #  -L/usr/local/MATLAB/R2014b/bin/glnxa64 -lmx -lmat -Wl,-rpath=/usr/local/MATLAB/R2014b/bin/glnxa64

../obj/Spin.o: source/spin/Spin/Spin.cpp include/spin/Spin.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/spin/Spin/Spin.cpp -o ../obj/Spin.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include #-std=c++0x

../obj/SpinData.o: source/spin/Spin/SpinData.cpp include/spin/SpinData.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/spin/Spin/SpinData.cpp -o ../obj/SpinData.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include #-std=c++0x

../obj/SpinCollection.o: source/spin/SpinCollection/SpinCollection.cpp include/spin/SpinCollection.h ../obj/SpinSource.o
	g++ -DARMA_DONT_USE_CXX11 -g -c source/spin/SpinCollection/SpinCollection.cpp -o ../obj/SpinCollection.o -I.  -I/home/nzhao/code/armadillo-6.400.3/include #-std=c++0x

../obj/SpinSource.o: source/spin/SpinCollection/SpinSource.cpp include/spin/SpinSource.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/spin/SpinCollection/SpinSource.cpp -o ../obj/SpinSource.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include #-std=c++0x

../obj/SpinCluster.o: source/spin/SpinCluster/SpinCluster.cpp include/spin/SpinCluster.h ../obj/SpinClusterAlgorithm.o
	g++ -DARMA_DONT_USE_CXX11 -g -c source/spin/SpinCluster/SpinCluster.cpp -o ../obj/SpinCluster.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include #-std=c++0x

../obj/SpinClusterAlgorithm.o: source/spin/SpinCluster/SpinClusterAlgorithm.cpp include/spin/SpinClusterAlgorithm.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/spin/SpinCluster/SpinClusterAlgorithm.cpp -o ../obj/SpinClusterAlgorithm.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include #-std=c++0x

../obj/SpinInteraction.o: source/spin/SpinInteraction/SpinInteraction.cpp include/spin/SpinInteraction.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/spin/SpinInteraction/SpinInteraction.cpp -o ../obj/SpinInteraction.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include #-std=c++0x

../obj/SpinInteractionComponent.o: source/spin/SpinInteraction/SpinInteractionComponent.cpp include/spin/SpinInteractionComponent.h include/misc/misc.h include/spin/SpinInteractionDefine.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/spin/SpinInteraction/SpinInteractionComponent.cpp -o ../obj/SpinInteractionComponent.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include #-std=c++0x

../obj/SpinState.o: source/spin/SpinInteraction/SpinState.cpp include/spin/SpinState.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/spin/SpinInteraction/SpinState.cpp -o ../obj/SpinState.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include #-std=c++0x

../obj/KronProd.o: source/kron/KronProd.cpp include/kron/KronProd.h include/spin/SpinInteractionDefine.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/kron/KronProd.cpp -o ../obj/KronProd.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include  -I/home/nzhao/code/armadillo-6.400.3/include #-std=c++0x

../obj/QuantumOperator.o: source/quantum/QuantumOperator/QuantumOperator.cpp include/quantum/QuantumOperator.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/quantum/QuantumOperator/QuantumOperator.cpp -o ../obj/QuantumOperator.o -I.  -I/home/nzhao/code/armadillo-6.400.3/include -I/usr/local/MATLAB/R2014b/extern/include/  #-std=c++0x

../obj/HilbertSpaceOperator.o: source/quantum/QuantumOperator/HilbertSpaceOperator.cpp include/quantum/HilbertSpaceOperator.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/quantum/QuantumOperator/HilbertSpaceOperator.cpp -o ../obj/HilbertSpaceOperator.o -I.  -I/home/nzhao/code/armadillo-6.400.3/include -I/usr/local/MATLAB/R2014b/extern/include/  #-std=c++0x

../obj/LiouvilleSpaceOperator.o: source/quantum/QuantumOperator/LiouvilleSpaceOperator.cpp include/quantum/LiouvilleSpaceOperator.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/quantum/QuantumOperator/LiouvilleSpaceOperator.cpp -o ../obj/LiouvilleSpaceOperator.o -I.  -I/home/nzhao/code/armadillo-6.400.3/include -I/usr/local/MATLAB/R2014b/extern/include/  #-std=c++0x

../obj/QuantumState.o: source/quantum/QuantumState/QuantumState.cpp include/quantum/QuantumState.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/quantum/QuantumState/QuantumState.cpp -o ../obj/QuantumState.o -I.  -I/home/nzhao/code/armadillo-6.400.3/include -I/usr/local/MATLAB/R2014b/extern/include/  #-std=c++0x

../obj/PureState.o: source/quantum/QuantumState/PureState.cpp include/quantum/PureState.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/quantum/QuantumState/PureState.cpp -o ../obj/PureState.o -I.  -I/home/nzhao/code/armadillo-6.400.3/include -I/usr/local/MATLAB/R2014b/extern/include/  #-std=c++0x

../obj/MixedState.o: source/quantum/QuantumState/MixedState.cpp include/quantum/MixedState.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/quantum/QuantumState/MixedState.cpp -o ../obj/MixedState.o -I.  -I/home/nzhao/code/armadillo-6.400.3/include -I/usr/local/MATLAB/R2014b/extern/include/  #-std=c++0x

../obj/misc.o: source/misc/misc.cpp
	g++ -DARMA_DONT_USE_CXX11 -g -c source/misc/misc.cpp -o ../obj/misc.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include   #-std=c++0x

../obj/QuantumEvolution.o: source/quantum/QuantumEvolution/QuantumEvolution.cpp include/quantum/QuantumEvolution.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/quantum/QuantumEvolution/QuantumEvolution.cpp -o ../obj/QuantumEvolution.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include -I/usr/local/MATLAB/R2014b/extern/include/  #-std=c++0x

../obj/QuantumEvolutionAlgorithm.o: source/quantum/QuantumEvolution/QuantumEvolutionAlgorithm.cpp include/quantum/QuantumEvolutionAlgorithm.h
	g++ -DARMA_DONT_USE_CXX11 -g -c source/quantum/QuantumEvolution/QuantumEvolutionAlgorithm.cpp -o ../obj/QuantumEvolutionAlgorithm.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include -I/usr/local/MATLAB/R2014b/extern/include/  #-std=c++0x

../obj/cce.o: app/cce.cpp
	g++ -DARMA_DONT_USE_CXX11 -g -c app/cce.cpp -o ../obj/cce.o  -I.  -I/home/nzhao/code/armadillo-6.400.3/include -I/usr/local/MATLAB/R2014b/extern/include/  #-std=c++0x

clean: 
	rm ../obj/*.o

doc:
	usr/Doxygen.app/Contents/Resources/doxygen ../doc/DocGen 

