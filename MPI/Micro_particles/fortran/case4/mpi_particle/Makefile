.SUFFIXES:
FC = mpiifort
# FC = mpif90
OBJDIR=objs
SRCDIR=src

SRC=commondata.F90 mpi_starts.F90 main.F90 initial.F90 fluid.F90 output.F90 freeall.F90\
	particle_force.F90 particle_update.F90  particle_bounceback.F90 \
	message_exc.F90 message_send_all.F90 message_particle.F90
objects=$(addprefix $(OBJDIR)/, $(subst .F90,.o,$(SRC)))

macros =

# define movingFrame / stationaryFrame
# define linear	/	quadratic

default: mpi_test

mpi_test : dirs $(objects)
	$(FC) $(macros) -o $@ $(objects)

$(OBJDIR)/%.o: $(SRCDIR)/%.F90
	$(FC) $< $(macros) -c -o $@


.PHONY : clean  clean_all dirs
clean :
	rm -f *.mod *.plt *.dat mpi_test nohup $(objects)

# clean_all: 
# 	rm -f *.plt *.dat mpi_test  $(objects) 

dirs:
	mkdir -p $(OBJDIR)/