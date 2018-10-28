FC=ifort 

%.o : %.f90
	$(FC) -c $?

module: parmkind.o  tetrahedron_geom.o


clean:
	rm -f *.o *.mod 
