(
	mkdir -p build bin
	cd build/
	export CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}:/usr/local/share/eigen3/cmake
	CXX=g++-14 CC=gcc-14 cmake ..
	make
)
