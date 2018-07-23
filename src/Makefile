# Makefile for DNABEND

#HOME	= /media/3AAB264357D77635/Programming/NOMeSeq_analysis_v2/source/
#BASE	= ..
#INCLS	= $(BASE)/include
#SOURCE	= $(BASE)/source
#BIN	= $(BASE)/bin
#YAMLINCL = /work2/gpeters/Appz/linux/usr/local/include/
#GSLINCL = /usr/include/gsl
#GMPINCL = /import/bc2/soft/app/gmp/4.2/Linux/include
#GSLINCL = $(HOME)/gsl-1.7/include
#CFLAGS	= -I$(INCLS) -I$(GSLINCL) -O3 -fopenmp -std=c++11
#CFLAGS	= -I$(INCLS) -I$(GSLINCL) -O3
#CFLAGS	= -I$(INCLS) -pg -Wall -I$(GSLINCL) -fopenmp
RCPPINCL = /work/gbioinfo/linux/lib64/R/library/Rcpp/include/
RINCL = /usr/share/R/include/
#CFLAGS = -O3 -std=c++11 -I$(YAMLINCL) -I$(RCPPINCL) -I$(RINCL)
CFLAGS = -O3 -std=c++11 -I$(RCPPINCL) -I$(RINCL)
#LIBS    =  /work2/gpeters/Appz/linux/usr/local/lib/libyaml-cpp.a
#/import/bc2/soft/app/gmp/4.2/Linux/lib/libgmp.a 
#/import/bc2/soft/app/gmp/4.2/Linux/lib/libgmp.la


# OBJS lists the files that do *not* contain a "main" and are used in the current build
# cf: grep main *.c -L | sed 's/\.c/.o/g'
OBJS	= parameters.o sequence.o usefullfunctions.o fastafile.o DNAbinding_object.o background.o binding_object_model.o DNAbindobj_vector.o forward_backward.o

# VPATH allows make to consider files in all our directories
VPATH   = $(CFLAGS)
CC      = g++


# utility functions #######################################################

#.PHONY: print clean
#print:
#	@echo BASE=$(BASE)
#	@echo BIN=$(BIN)
#	@echo SEDBIN=$(SEDBIN)
#	@echo CFLAGS=$(CFLAGS)
#	@echo LIBS=$(LIBS)
#	@echo OBJS=$(OBJS)
#	@echo INCLS=$(INCLS)
#	@echo GSLINCL=$(GSLINCL)

#clean:
#	rm $(SOURCE)/*.o $(SOURCE)/.*.d
#	rm -i $(SOURCE)/MotEvo_OOP



# define new pattern rules for compilation #################################
# (includes path info & dependency compatibility)

# executables
NOMeSeq_analysis: NOMeSeq_analysis.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)


# object files
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $< 




# automatic dependency ###########################################

CCM = g++
CCMFLAG = -M  # CCMFLAG is the flag needed to generate dependencies

# DEPS variable contains a list of all dependency makefiles
DEPS = $(patsubst %.c,.%.d,$(wildcard *.c))

# pattern rule for creating dependency files
SEDBIN  = $(subst /,\\/,$(BIN))# escapes the path separator, "/" for sed
.%.d: %.c
	@echo -n making $@ from $<...
	@$(SHELL) -ec "$(CCM) $(CCMFLAG) $(CFLAGS) $<  > $@"
# 
	@echo done.

include $(DEPS)

#@$(SHELL) -ec "$(CCM) $(CCMFLAG) $(CFLAGS) $< | sed -e 's/$*.o/$(SEDBIN)\/$*.o $@/g' > $@"
# shortcuts #############################################################


# etags for emacs jumps
tags:	TAGS
TAGS:   $(wildcard *.c) $(wildcard $(INCLS)/*.h)
	etags -C $(INCLS)/*.h $(SOURCE)/*.c 








