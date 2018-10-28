FC=ifort 
srcroot=tetrahedron_all
srcfile=$(srcroot).f90
exefile=$(srcroot).x

%.o : %.f90
	$(FC) -c $?

INCDIR=../src

all: mod
	$(FC) -I$(INCDIR) $(INCDIR)/*.o  $(srcfile) -o $(exefile)

mod:
	cd ../src/ && make clean; make
run: 
	./$(exefile) 

runall: all run

clean:
	rm -f *.o *.mod $(exefile) 

cleanAll: clean
	cd ../src && make clean
