PROG = crest

OBJDIR = build
#--------------------------------------------------------------------------
 OSTYPE=LINUXI
 #MKLROOT=/opt/intel/mkl
#--------------------------------------------------------------------------

MODULES += classes.o 
MODULES += setdg.o
MODULES += atmasses.o
MODULES += ls_rmsd.o
MODULES += iomod.o
MODULES += filemod.o
MODULES += syscheck.o
MODULES += zdata.o
MODULES += strucreader.o
MODULES += msmod.o
MODULES += thermo.o

CFILES += symmetry2.o
CFILES += symmetry_i.o
CFILES += signal.o
CFILES += sigterm.o
CFILES += 

OBJS1 += data.o
OBJS1 += io.o
OBJS1 +=

OBJS2 += drsp.o
OBJS2 += flexi.o
OBJS2 += trialmd.o
OBJS2 += neighbor.o
OBJS2 += prmat.o
OBJS2 += readl2.o
OBJS2 += readl.o
OBJS2 += scratch.o
OBJS2 += sdfio.o
OBJS2 += shortpath.o
OBJS2 += solve.o

OBJS3 += axis3.o
OBJS3 += cleanup.o
OBJS3 += confcross.o
OBJS3 += confopt.o
OBJS3 += confparse.o
OBJS3 += constraining.o
OBJS3 += cregen.o
OBJS3 += ccegen.o
OBJS3 += deprotonate.o
OBJS3 += dft_propcalc.o
OBJS3 += ensemblecomp.o
OBJS3 += eval_timer.o
OBJS3 += freqmasses.o
OBJS3 += identifiers.o
OBJS3 += intmodes.o
OBJS3 += ligands.o
OBJS3 += minitools.o
OBJS3 += ncigeo.o
OBJS3 += newcregen.o
OBJS3 += normmd.o
OBJS3 += ompmklset.o
OBJS3 += optinplace.o
OBJS3 += printouts.o
OBJS3 += progressbar.o
OBJS3 += propcalc.o
OBJS3 += protonate.o
OBJS3 += quicksort.o
OBJS3 += rotcompare.o
OBJS3 += select.o
OBJS3 += sortens.o
OBJS3 += spline.o
OBJS3 += tautomerize.o
OBJS3 += thermocalc.o
OBJS3 += trackorigin.o
OBJS3 += utilities.o
OBJS3 += xtbmodef.o
OBJS3 += zsort.o
OBJS3 +=

OBJS4 += geo.o
OBJS4 += ztopology.o
OBJS4 += floydwarshall.o
OBJS4 += bondconstraint.o
OBJS4 += stereoisomers.o
OBJS4 += marqfit.o
OBJS4 += entropy.o
OBJS4 += entropic.o
OBJS4 += reactor.o
OBJS4 += compress.o
OBJS4 += relaxensemble.o
OBJS4 += acidbase.o
OBJS4 += biasmerge.o
OBJS4 += arvo.o
OBJS4 += solvtool.o
OBJS4 += solvtool_misc.o

#mainfile for crest binary
MAINFILES += confscript1.o
MAINFILES += confscript2_main.o
MAINFILES += confscript2_misc.o
MAINFILES += confscript3.o
MAINFILES += msreact.o
MAINFILES += pka.o
MAINFILES += crest_main.o 

OBJS := $(MODULES) $(CFILES) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(MAINFILES)
#--------------------------------------------------------------------------

vpath % .
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))
#--------------------------------------------------------------------------

ifeq ($(OSTYPE),LINUXI)
  PREOPTS :=
  FC := ifort
  CC := icc
  LINKER = ifort -g -O -static -fopenmp  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
#  LINKER = ifort -static -fopenmp  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
#  LIBS = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
  LIBS := -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group 
  PREFLAG := -E -P
  FFLAGS := -g -O -qopenmp -align array64byte -module $(OBJDIR)
#  FFLAGS := -g -O0
#  FFLAGS := -O -axAVX -qopenmp -align array64byte 
#  FFLAGS += -check all -g -debug all -traceback -fpe0
  CCFLAGS += -D_Float128=__float128
  CCFLAGS := -O -DLINUX
endif

ifeq ($(OSTYPE),LINUXG)
  PREOPTS :=
  FC := gfortran
  CC := gcc
  LINKER := gfortran -static-libgfortran -fopenmp 
  LIBS :=
  PREFLAG :=
  FFLAGS :=  -J$(MODDIR)
  CCFLAGS := -O -DLINUX
endif


ifeq ($(OSTYPE),SOMETHINGELSE)
  PREOPTS :=
  FC := 
  CC := 
  LINKER := 
  LIBS :=
  PREFLAG := 
  FFLAGS := 
  CCFLAGS := 
endif


.PHONY: all
.PHONY: setup
.PHONY: clean

all: setup
all: $(PROG)
setup: $(OBJDIR)

#--------------------------------------------------------------------------
# Setting the rules to make the object files:
#--------------------------------------------------------------------------
# implizite Regel zur Erzeugung von *.o aus *.F ausschalten
%.o: %.F
fortran.o : $(CUDA)/src/fortran.c
	gcc $(CCFLAGS) -I$(CUDA)/src/ $(INC) -c  $(CUDA)/src/fortran.c

# aus *.F mache ein *.f
%.f: %.F
	@echo "making $@ from $<"
	$(FC) $(PREFLAG) $(PREOPTS) $< -o $@

# aus *.f mache ein *.o
$(OBJDIR)/%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

$(OBJDIR)/%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.F90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# aus *.c mache ein *.o
$(OBJDIR)/%.o: %.c
	@echo "making $@ from $<"
	$(CC) $(CCFLAGS) -c $< -o $@

# linken
$(PROG): $(OBJS)
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)

$(OBJDIR):
	mkdir -p $@

#aufraeumen
clean:
	rm -f $(filter %.o,$(OBJS)) $(PROG) 
	rm -f $(OBJDIR)/*.mod *.mod $(PROG) 
	rm -f $(patsubst %.F, %.f, $(wildcard *.F))

