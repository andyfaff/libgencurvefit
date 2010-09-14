# Define required macros here
SHELL = /bin/sh

CFLAGS = -O3 -funroll-loops -c -fopenmp -std=c99
CC = gcc
STATICLIB=libgencurvefit.a
LIBS=$(STATICLIB)
AR=ar rsc

prefix = /usr/local
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
sharedlibdir = ${libdir}
includedir = ${prefix}/include

OBJS = src/gencurvefit.o src/mt19937p.o src/levenbergMarquardt.o

all: libgencurvefit.a gaussian globalfitting

libgencurvefit.a: $(OBJS)
	$(AR) $@ $(OBJS)
	-@ ($(RANLIB) $@ || true) >/dev/null 2>&1

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $< 

clean:
	rm src/*.o
	rm *.a
	
install-libs: $(LIBS)
	-@if [ ! -d $(DESTDIR)$(exec_prefix)  ]; then mkdir -p $(DESTDIR)$(exec_prefix); fi
	-@if [ ! -d $(DESTDIR)$(libdir)       ]; then mkdir -p $(DESTDIR)$(libdir); fi
	cp $(STATICLIB) $(DESTDIR)$(libdir)
	cd $(DESTDIR)$(libdir); chmod u=rw,go=r $(STATICLIB)
	
install: install-libs
	-@if [ ! -d $(DESTDIR)$(includedir)   ]; then mkdir -p $(DESTDIR)$(includedir); fi
	cp src/gencurvefit.h $(DESTDIR)$(includedir)
	chmod 644 $(DESTDIR)$(includedir)/gencurvefit.h

gaussian: 
	g++ -O3 examples/gaussian/*.cpp -Isrc -L. -fopenmp -o examples/gaussian/gaussian_fitter -lgencurvefit -lgomp -lpthread -lm
	
globalfitting:
	g++ -O3 examples/globalfitting/*.cpp -Isrc -fopenmp -L. -o examples/globalfitting/global_fitter -lgencurvefit -lgomp -lpthread -lm
