# Please define environmental variables $(ARMA_LIB) and $(ARMA_INC_PATH) to enable armadillo 
# If you have Malab installed, pleas set $(HAS_MATLAB) = true, 
# and define $(MATLAB_LIB_PATH) and $(MATLAB_INC_PATH) in your .bashrc file.
# If you want to link the program by INTEL MKL libaray, set $(USE_MKL)=true in your .bashrc file.

DESTINATION := ../bin/cce_solution
LIBPATH := .  
INCLUDES := .  $(ARMA_INC_PATH) 

RM := rm -f
PS = cpp

CC = mpic++
CPPFLAGS = $(addprefix -I,$(INCLUDES))
CPPFLAGS += -g -DARMA_DONT_USE_WRAPPER 

ifeq ($(USE_MKL), true)
	LIBS := 
	#CPPFLAGS  +=  -DARMA_BLAS_LONG -DMKL_ILP64 -qopenmp -I${MKLROOT}/include 
	#CXXLINKS  :=  -parallel -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_ilp64 -lpthread -lm 
	CPPFLAGS  +=  -DARMA_BLAS_LONG -DMKL_ILP64 -I${MKLROOT}/include 
	CXXLINKS  :=  -parallel -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_ilp64 -lpthread -lm
else
	LIBS :=  $(ARMA_LIB) 
	CXXLINKS  := 
endif


ifeq ($(HAS_MATLAB), true)
	LIBS += mx mat eng
	CPPFLAGS += -DHAS_MATLAB
	LIBPATH := $(MATLAB_LIB_PATH)
	INCLUDES += $(MATLAB_INC_PATH)
endif

SOURCE := $(wildcard application/*.$(PS)) $(wildcard source/*/*.$(PS)) $(wildcard source/*/*/*.$(PS))
OBJS := $(patsubst %.$(PS),%.o,$(notdir $(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))

OBJPATH := ../obj
vpath %.cpp $(sort $(dir $(SOURCE)))
BINPATH := ../bin

.PHONY : all clean

all : $(DESTINATION)

clean :
	@$(RM) $(OBJPATH)/*.o
	@$(RM) $(OBJPATH)/*.d
	@$(RM) $(DESTINATION)

doc :
	/Applications/Doxygen.app/Contents/Resources/doxygen ../doc/DocGen 

$(addprefix $(OBJPATH)/, $(DEPS)) : $(OBJPATH)/%.d : %.cpp | $(OBJPATH)
	@$(CC) -MM $< -I. -I$(MATLAB_INC_PATH) > $@
	@sed -i="" '1s/^/..\/obj\//g' $@
	@rm $@=

-include $(addprefix $(OBJPATH)/, $(DEPS))

$(DESTINATION) : $(addprefix $(OBJPATH)/, $(OBJS)) | $(BINPATH)
	$(CC) -o $(DESTINATION) $(addprefix $(OBJPATH)/, $(OBJS)) -L$(LIBPATH) $(addprefix -l,$(LIBS)) $(CXXLINKS)

$(addprefix $(OBJPATH)/, %.o) : %.cpp | $(OBJPATH)
	$(CC) -c $< -o $@ $(CPPFLAGS)

$(OBJPATH) :
	@mkdir -p $(OBJPATH)

$(BINPATH) :
	@mkdir -p $(BINPATH)

