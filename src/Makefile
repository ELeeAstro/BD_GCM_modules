## To compile, select the compile options below and enter 'make' into the terminal
## the command 'make clean' removes all compiled files and executable files for recompilation

# for gfortran Compiler
#======================
F90          = gfortran
F90LINKER    = gfortran
DEFS      =
#FFLAGS   	= -Og -pipe -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -Wconversion -fbounds-check
FFLAGS   = -O3 -pipe
INCLUDES  =
LFLAGS    = $(FFLAGS)

# for nvfortran Compiler
#======================
# F90          = nvfortran
# F90LINKER    = nvfortran
# DEFS      =
# FFLAGS   = -fast
# INCLUDES  =
# LFLAGS    = $(FFLAGS)


# for ifort Compiler
#======================
# F90          = ifort
# F90LINKER    = ifort
# DEFS      =
# FFLAGS   = -O3
# INCLUDES  =
# LFLAGS    = $(FFLAGS)


OBJECTS = \
    dopri5.o \
	  k_Rosseland_mod.o \
	  IC_mod.o \
		BD_kappa_poly_mod.o \
		BD_MLT_mod.o \
		BD_clouds_mod.o \
    BD_ver_diff_mod.o \
		BD_clouds_adv_mod.o \
		dry_conv_adj_mod.o \
		CE_mod.o \
		dc_decsol.o \
		decsol.o \
		seulex.o \
		mini_ch_precision.o \
		mini_ch_class.o \
		mini_ch_read_reac_list.o \
		mini_ch_chem.o \
		mini_ch_i_seulex.o \
	  ts_Toon_scatter_mod.o \
    FMS_BD_RC.o

# executable statement
EXECS  = ../FMS_BD_RC

.SUFFIXES: .o .f90 .F .f

default: FMS_BD_RC

FMS_BD_RC:  $(OBJECTS)
	 $(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

clean:
	rm -f *.o *.mod *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.F.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<