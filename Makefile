#Some Makefile for CLASS_SZ.
#Julien Lesgourgues, 28.11.2011 for the CLASS part.
#Boris Bolliet 2015- for the class_sz part.

# To install class_sz all you need to do is to set the environment variable
# PATH_TO_COSMOPOWER_ORGANIZATION to the path to the cosmopower-organization directory
# and then run the makefile.


MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
OUTDIR = $(MDIR)/output
file_path=$(MDIR)/python/classy_szfast/classy_szfast/config.py
PATH_TO_COSMOPOWER_ORGANIZATION=



.base:
	if ! [ -e $(WRKDIR) ]; then mkdir -p $(WRKDIR)/lib; fi;
	if ! [ -e $(OUTDIR) ]; then mkdir $(OUTDIR); fi;
	touch build/.base

	DEPTH=3; \
	if [ -z "$(PATH_TO_COSMOPOWER_ORGANIZATION)" ]; then \
		COSMOPOWER_DIR=$$(find .. -maxdepth $$DEPTH -type d -name "cosmopower-organization" -print -quit); \
		if [ -z "$$COSMOPOWER_DIR" ]; then \
			COSMOPOWER_DIR=$$(find . -maxdepth $$DEPTH -type d -name "cosmopower-organization" -print -quit); \
		fi; \
		if [ -n "$$COSMOPOWER_DIR" ]; then \
			echo "Found cosmopower-organization directory at: $$(realpath "$$COSMOPOWER_DIR")"; \
			PATH_TO_COSMOPOWER_ORGANIZATION=$$(realpath "$$COSMOPOWER_DIR"); \
		else \
			echo "--> cosmopower-organization directory not found within $$DEPTH levels up or down.\n\
--> We will install it one level up!"; \
			cd ..;\
			mkdir cosmopower-organization;\
			cd cosmopower-organization;\
			git clone https://github.com/cosmopower-organization/lcdm.git;\
			git clone https://github.com/cosmopower-organization/mnu.git;\
			git clone https://github.com/cosmopower-organization/mnu-3states.git;\
			git clone https://github.com/cosmopower-organization/ede.git;\
			git clone https://github.com/cosmopower-organization/neff.git;\
			git clone https://github.com/cosmopower-organization/wcdm.git;\
			cd ..;\
		fi; \
	else \
		echo "Using provided PATH_TO_COSMOPOWER_ORGANIZATION: $(PATH_TO_COSMOPOWER_ORGANIZATION)"; \
		PATH_TO_COSMOPOWER_ORGANIZATION=$$(realpath "$(PATH_TO_COSMOPOWER_ORGANIZATION)"); \
	fi; \
	export PATH_TO_COSMOPOWER_ORGANIZATION; \
	rm -f $(file_path); \
	touch $(file_path); \
	echo "path_to_cosmopower_organization = '$$PATH_TO_COSMOPOWER_ORGANIZATION'" >> $(file_path)


vpath %.c source:tools:main:test
vpath %.o build
vpath .base build

########################################################
###### LINES TO ADAPT TO YOUR PLATEFORM ################
########################################################

# your C compiler:
# CC       = gcc-12
# on Mac M1:
CC       = /usr/bin/clang


# your tool for creating static libraries:
AR        = ar rv


PYTHON ?= python


# your optimization flag
#OPTFLAG = -O4 -ffast-math #-march=native
# on Mac M1

OPTFLAG = -O3 #-ffast-math #-ffast-math #-arch x86_64
#OPTFLAG = -O3 
#OPTFLAG = -Ofast -ffast-math #-march=native
#OPTFLAG = -fast

# your openmp flag (comment for compiling without openmp)
#OMPFLAG   = -fopenmp
# on Mac M1
OMPFLAG   = -Xclang -fopenmp
#OMPFLAG   = -mp -mp=nonuma -mp=allcores -g
#OMPFLAG   = -openmp

# all other compilation flags
CCFLAG = -g -fPIC
LDFLAG = -g -fPIC

#on Mac M1
LDFLAG += -lomp

# leave blank to compile without HyRec, or put path to HyRec directory
# (with no slash at the end: e.g. hyrec or ../hyrec)
HYREC = hyrec


########################################################
###### IN PRINCIPLE THE REST SHOULD BE LEFT UNCHANGED ##
########################################################

# pass current working directory to the code
CCFLAG += -D__CLASSDIR__='"$(MDIR)"'

# where to find include files *.h
INCLUDES =  -I../include -I/usr/local/include/ 


# automatically add external programs if needed. First, initialize to blank.
EXTERNAL =

# eventually update flags for including HyRec
ifneq ($(HYREC),)
vpath %.c $(HYREC)
CCFLAG += -DHYREC
#LDFLAGS += -DHYREC
INCLUDES += -I../hyrec
EXTERNAL += hyrectools.o helium.o hydrogen.o history.o
endif


# if you would like to print the warnings, make sure you dont set the following flags
CFLAGS = -Wno-parentheses-equality -Wno-format -Wno-return-type -Wno-comment -Wno-incompatible-pointer-types -Wno-deprecated-non-prototype -Wno-return-type -Wno-visibility


%.o:  %.c .base
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(CFLAGS) $(INCLUDES) -c ../$< -o $*.o

TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o  evolver_ndf15.o arrays.o parser.o quadrature.o hyperspherical.o common.o trigonometric_integrals.o r8lib.o class_sz_tools.o class_sz_custom_profiles.o class_sz_custom_bias.o Patterson.o fft.o

SOURCE = input.o background.o thermodynamics.o perturbations.o primordial.o nonlinear.o transfer.o spectra.o lensing.o class_sz.o class_sz_clustercounts.o

INPUT = input.o

PRECISION = precision.o

BACKGROUND = background.o

THERMO = thermodynamics.o

PERTURBATIONS = perturbations.o

TRANSFER = transfer.o

PRIMORDIAL = primordial.o

SPECTRA = spectra.o

NONLINEAR = nonlinear.o

LENSING = lensing.o

OUTPUT = output.o

CLASS_SZ = class_sz_driver.o


C_TOOLS =  $(addprefix tools/, $(addsuffix .c,$(basename $(TOOLS))))
C_SOURCE = $(addprefix source/, $(addsuffix .c,$(basename $(SOURCE) $(OUTPUT))))
C_TEST = $(addprefix test/, $(addsuffix .c,$(basename $(TEST_DEGENERACY) $(TEST_LOOPS) $(TEST_TRANSFER) $(TEST_NONLINEAR) $(TEST_PERTURBATIONS) $(TEST_THERMODYNAMICS))))
C_MAIN = $(addprefix main/, $(addsuffix .c,$(basename $(CLASS_SZ))))
C_ALL = $(C_MAIN) $(C_TOOLS) $(C_SOURCE)
H_ALL = $(addprefix include/, common.h svnversion.h $(addsuffix .h, $(basename $(notdir $(C_ALL)))))
PRE_ALL = cl_ref.pre clt_permille.pre
INI_ALL = explanatory.ini lcdm.ini
MISC_FILES = Makefile CPU psd_FD_single.dat myselection.dat myevolution.dat README bbn/sBBN.dat external_Pk/* cpp
PYTHON_FILES = python/classy.pyx python/setup.py python/cclassy.pxd python/test_class.py


all: class_sz libclass.a classy_sz

libclass.a: $(TOOLS) $(SOURCE) $(EXTERNAL)
	$(AR)  $@ $(addprefix build/, $(TOOLS) $(SOURCE) $(EXTERNAL))

class_sz: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS_SZ) 
	$(CC) $(CFLAGS) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -g -o class_sz $(addprefix build/,$(notdir $^)) -lgsl -lgslcblas -lfftw3 -lm


tar: $(C_ALL) $(C_TEST) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)
	tar czvf class_sz.tar.gz $(C_ALL) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)

classy_sz: libclass.a python/classy.pyx python/cclassy.pxd
	cd python/classy_szfast; pip install -r requirements.txt -e . 
ifdef OMPFLAG
	cp python/setup.py python/autosetup.py
else
	grep -v "lgomp" python/setup.py > python/autosetup.py
endif
	cd python; pip install -e .
	


clean: 
	rm -rf $(WRKDIR);
	rm -f libclass.a
	rm -f $(MDIR)/python/classy.c
	rm -rf $(MDIR)/python/build

