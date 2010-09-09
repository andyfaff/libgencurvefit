mpic++ -O3 -funroll-loops -DUSE_MPI -I../../  -fopenmp *.cpp  -o ../build/Debug/motoMPI -lgomp -lpthread -lgencurvefit


On commodore
mpic++ -O3 -funroll-loops -DUSE_MPI -I../libgencurvefit/  -L../libgencurvefit -fopenmp *.cpp  -o motoMPI -lgomp -lpthread -lgencurvefit
