#	Makefile

SRCS = ../Modules/COLLISION_INPUTS_MODULE.f90 \
       ../Modules/PHYSICS_MODULE.f90 \
       ../Modules/GRID_MODULE.f90 \
       ../Modules/CURRENT_ENERGY_MODULE.f90 \
       ../Modules/MPI_INFO_MODULE.f90 \
       ../Modules/DCS_MODULE.f90 \
       ../Modules/EIKONAL_MODULE.f90 \
			 ../Modules/LM_MODULE.f90 \
       ../Potentials/He+H/HeH_CUBIC_SPLINE_MODULE.f90 \
       ../Potentials/potential.f90 \
       ../Potentials/He+He/HeHe_AbInitCore.f90 \
       ../Potentials/He+H/HeH_Uni.f90 \
       ../Potentials/He+H/HeH_Spline.f90 \
       ../Potentials/He+O/HeO_Uni.f90 \
       ../Potentials/He+O/HeO_Pi_Core.f90 \
       ../Potentials/He+O/HeO_Sg_Core.f90 \
       ../Potentials/H+H/Pot_Singlet.f90 \
       ../Potentials/H+H/Pot_Triplet.f90 \
       ../Potentials/H+H/pot_mjj.f90 \
       ../Potentials/H+H/spline.f90 \
       ../Potentials/H+H/lj.f90 \
       ../Routines/input_reader.f90 \
			 ../Routines/angle_fit.f90 \
       ../Routines/universal_amplitude.f90 \
       ../Routines/phaseshift.f90 \
       ../Routines/grid.f90 \
       ../Routines/bessel.f90 \
       ../Routines/total_cross_section.f90 \
       ../Routines/diffusion_cross_section.f90 \
       ../Routines/diff_cross_section.f90 \
       ../Routines/legendre.f90 \
       ../Routines/pw_amp.f90 \
       ../Routines/hard_sphere_test.f90 \
       ../Routines/timing_test.f90 \
       ../Routines/eikonal.f90 \
       ../Routines/eik_phase.f90 \
       ../Routines/eik_amp.f90 \
       ../Routines/eik_diff_cross_section.f90 \
       ../Routines/eik_total_cross_section.f90 \
       ../Routines/frame_converter.f90 \
       ../Routines/phase_writer.f90 \
       ../Routines/prob_density.f90 \
       ../Routines/critical_angle.f90 \
       ../Routines/depth_finder.f90 \
       ../Routines/pot_tester.f90 \
       ../Routines/reduced_coordinates.f90 \
       ../Routines/diff_cross3d.f90 \
       ../Routines/ave_energy_loss.f90 \
       ../Routines/call_tests.f90 \
       ../Routines/real_imag_amp.f90 \
       ../Routines/randy.f90 \
       ../Routines/phase_JWKB.f90 \
       ../Routines/He_tcs.f90 \
       ../Routines/HH_tcs.f90 \
       ../Routines/collide.f90

OBJS = $(SRCS:.f90=.o) 

#Compile and link flags
CC	    = mpif90 
CFLAGS  = -O3 
DEBUG   =
LDFLAGS =  
EXE     = COLLIDE 

#Link and executable create
$(EXE) : $(OBJS) 
	$(CC) $(DEBUG) $(LDFLAGS) -o $@ $(OBJS) 

# Compile using a pattern rule 
%.o : %.f90 $(HDRS)   
	$(CC) -c $(CFLAGS) $(DEBUG) -o $@ $<

clean:
	rm -f ../Potentials/H+H/*.o ../Potentials/He+He/*.o ../Potentials/He+H/*.o ../Potentials/He+O/*.o ../Potentials/H+H/*.o ../Potentials/*.o ../Routines/*.o ../Modules/*.o core *.mod

clobber:
	rm -f ../Potentials/H+H/*.o ../Potentials/He+He/*.o ../Potentials/He+H/*.o ../Potentials/He+O/*.o ../Potentials/H+H/*.o ../Potentials/*.o ../Modules/*.o ../Routines/*.o core *.mod ../Data/*.dat $(EXE) 

