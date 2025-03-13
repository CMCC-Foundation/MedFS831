SHELL = /bin/sh

############################################################################
#
#        Makefile for input4NEMO 
#
############################################################################
#
#    Copyright 2020 Antonio Mariani, CMCC, Bologna
#
############################################################################
#
# INTEL 

F90=ifort
#
RM  = rm -f

EXTINC = $(FFLAGS)
EXTLIB = $(LDFLAGS)


SRC =  prep_SLA_3dvar_C7_EUR_DT_7km.f90 prep_SLA_3dvar_C7_EUR_7km.f90

EXEC = $(SRC:.f90=.exe)


all: $(EXEC)
	@echo $(EXEC) are compiled

$(EXEC) :
	$(F90)  $(@:.exe=.f90)  -o $@  $(EXTINC) $(EXTLIB)
	
clean:
	$(RM) *.o *.mod $(EXEC)
