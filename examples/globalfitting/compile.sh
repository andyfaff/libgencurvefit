mpic++ -O3 -funroll-loops -DUSE_MPI -I../../  -fopenmp *.cpp  -o ../build/Debug/motoMPI -lgomp -lpthread -lgencurvefit


On commodore
mpic++ -O3 -funroll-loops -DUSE_MPI -I../../src  -L../../../libgencurvefit -fopenmp *.cpp  -o ~/bin/global_fitter -lgomp -lpthread -lgencurvefit -ldl
