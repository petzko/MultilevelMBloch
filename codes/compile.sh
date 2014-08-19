g++ -c -std=gnu99 PSTDSolver.cpp -I/home/kenny/libs/fftw-3.3.4/fftw-installation/include/
g++ PSTDSolver.o -L/home/kenny/libs/fftw-3.3.4/fftw-installation/lib/ -lfftw3 -o start
