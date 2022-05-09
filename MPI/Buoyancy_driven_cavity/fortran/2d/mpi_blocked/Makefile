.SUFFIXES:
FC = mpiifort
OBJDIR=objs
SRCDIR=src

SRC=module.F90 main.F90 initial.F90 output.F90 evolution_f.F90 evolution_g.F90 check.F90\
	message_exchange.F90

objects=$(addprefix $(OBJDIR)/, $(subst .F90,.o,$(SRC)))

macros = -DsteadyFlow -DHorizontalWallsNoslip -DVerticalWallsNoslip \
	-DSideHeatedCell -DHorizontalWallsAdiabatic -DVerticalWallsConstT

#define steadyFlow   
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT

default: mpi_test

mpi_test : dirs $(objects)
	$(FC) $(macros) -o $@ $(objects)

$(OBJDIR)/%.o: $(SRCDIR)/%.F90
	$(FC) $< $(macros) -c -o $@


.PHONY : clean  clean_all dirs
clean :
	rm -f mpi_test $(objects) 

clean_all: 
	rm -f *.plt *.dat mpi_test  $(objects) 

dirs:
	mkdir -p $(OBJDIR)/