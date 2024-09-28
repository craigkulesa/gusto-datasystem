KERNEL := $(shell uname -s)

CFLAGS := -O2 -Wall #-fsanitize=address

ifeq ($(KERNEL), $(filter $(KERNEL), Linux NetBSD FreeBSD))
      CXX := gcc
      INCLUDES := -I/usr/local/include -I/usr/pkg/include 
      LIBS := -Wl,-rpath,/usr/local/lib -lit -L/usr/pkg/lib -L/usr/local/lib -lm -lfftw3 -lcurl -lcfitsio

else ifeq ($(KERNEL), Darwin)
      CXX := clang
      INCLUDES := -I/usr/local/include -I/opt/pkg/include 
      LIBS := -Wl,-rpath,/usr/local/lib -lit -L/opt/pkg/lib -lm -lfftw3 -lcfitsio -lcurl

else
      $(error unknown OS)

endif

corrspec: corrspec.c corrspec.h callback.c callback.h influx.c influx.h
	$(CXX) -c callback.c $(INCLUDES) $(CFLAGS)
	$(CXX) -c corrspec.c $(INCLUDES) $(CFLAGS)
	$(CXX) -c influx.c $(INCLUDES) $(CFLAGS)
	$(CXX) -o corrspec corrspec.o callback.o influx.o $(INCLUDES) $(CFLAGS) $(LIBS)

all: corrspec

clean:
	rm -f corrspec *.o *~ *.log x?? 

