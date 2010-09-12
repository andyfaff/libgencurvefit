# Define required macros here
SHELL = /bin/sh

OBJ =  gencurvefit.o mt19937p.o
CFLAG = -O3 -funroll-loops -c -fopenmp
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

all: libgencurvefit.a

libgencurvefit.a: $(OBJS)
	$(AR) $@ $(OBJS)
	-@ ($(RANLIB) $@ || true) >/dev/null 2>&1
	
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
	
