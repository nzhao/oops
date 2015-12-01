cce.o: src/app/cce.cpp
	g++ -c src/app/cce.cpp -o obj/cce.o -I src/ -std=c++11
	g++ -c src/source/spin/Spin/Spin.cpp -o obj/Spin.o -I src/ -std=c++11
	g++ -c src/source/spin/Spin/SpinData.cpp -o obj/SpinData.o -I src/ -std=c++11
	g++ obj/cce.o obj/Spin.o obj/SpinData.o -o bin/cce -I src/
