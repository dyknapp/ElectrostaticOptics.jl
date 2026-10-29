FC=gfortran
FFLAGS=-g
LDFLAGS=-lgfortran -shared -fPIC

qawse.dll: dqawse.f
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $?

clean:
	rm -f *.o *.exe *.dll