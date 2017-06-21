FC = ifort

FLAGS = -mcmodel=large

FFLAGS = -O3 -qopenmp#-axMIC-AVX512 

OPTIONS = -DFFTFINE

INC = -I$SCINET_FFTW_INC -L$SCINET_FFTW_LIB -lfftw3f
all: rad3d.x

rad3d.x: params.o pencil_fft.o calc.o write.o rad3d.o fourn.for
	$(FC) $(FLAGS) $(FFLAGS) -o $@ $^ $(INC)

pencil_fft.o: pencil_fft.f90
	$(FC) $(FLAGS) $(FFLAGS) -c $< $(INC)

params.o: params.f90
	$(FC) $(FLAGS) $(FFLAGS) -c $< 

calc.o: calc.f90 
	$(FC) $(FLAGS) $(FFLAGS) $(INC) -c $<

write.o: write.f90
	$(FC) $(FLAGS) $(FFLAGS) -c $<

rad3d.o: rad3d.f90
	$(FC) $(FLAGS) $(FFLAGS) -c $<

clean: 
	rm *.o
	rm *.mod

