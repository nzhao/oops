# Please define environmental variables $(ARMA_LIB) and $(ARMA_INC_PATH) to enable armadillo 
# If you have Malab installed, pleas set $(HAS_MATLAB) = true, 
# and define $(MATLAB_LIB_PATH) and $(MATLAB_INC_PATH) in your .bashrc file.
# If you want to link the program by INTEL MKL libaray, set $(USE_MKL)=true in your .bashrc file.
# If you are using INTEL MPI, set $(MPI_INTEL)=true in your .bashrc file.
ROOT_DIR  :=  $(shell pwd)

BINPATH := ../bin
OBJPATH := ../obj

APP := EnsNVeBath EnsNVnBath EnsNVLattice SingleSampleCCE_NV_C13 MatExpTest test magR_new

DESTINATION := $(addprefix $(BINPATH)/, $(APP))
LIBPATH := .  
INCLUDES := .  $(ARMA_INC_PATH) 

RM := rm -f
PS = cpp

CC = mpic++ -g
CPPFLAGS = $(addprefix -I,$(INCLUDES)) -qopenmp

ifeq ($(USE_11), true)
	CPP11FLAG := -std=gnu++11
else
	CPP11FLAG := 
endif

ifeq ($(USE_MKL), true)
	LIBS := 
	CXXLINKS  := -parallel -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_ilp64 -lpthread -lm -ldl
	CPPFLAGS  +=  -DARMA_BLAS_LONG -DMKL_ILP64 -m64 -I${MKLROOT}/include -DARMA_DONT_USE_WRAPPER -DUSE_MKL
else
	LIBS :=  $(ARMA_LIB) 
	CXXLINKS  := 
	CPPFLAGS += -g -DARMA_DONT_USE_WRAPPER 
endif


ifeq ($(HAS_MATLAB), true)
	LIBS += mx mat eng
	CPPFLAGS += -DHAS_MATLAB
	LIBPATH := $(MATLAB_LIB_PATH) -Wl,-rpath=$(MATLAB_LIB_PATH)
	INCLUDES += $(MATLAB_INC_PATH)
endif

SOURCE := $(wildcard application/$(APP).$(PS)) $(wildcard source/*/*.$(PS)) $(wildcard source/*/*/*.$(PS))
OBJS := $(patsubst %.$(PS),%.o,$(notdir $(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))
APP_OBJ := $(addsuffix .o, $(addprefix $(OBJPATH)/, $(APP)) )

vpath %.cpp $(sort $(dir $(SOURCE)))

FC          := mpif77
ifeq ($(MPI_INTEL), true)
  FFLAGS      := -i8 -I${MKLROOT}/include
  FLINKS      := -nofor-main -cxxlib -parallel -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
else
  FFLAGS      := -fdefault-integer-8 -m64 -I${MKLROOT}/include
#  CXXLINKS    :=  -L${GFORTRAN_PATH} -lpthread -lm -ldl -lgfortran
  CXXLINKS    := -lgfortran -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
endif

NCC         := nvcc
NCFLAGS     := -arch=sm_30 -O3 -DCUBLAS_GFORTRAN -ccbin g++ -Xcompiler "$(CPPFLAGS)"
NCLINKER    := -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusparse -lboost_program_options -L/usr/local/lib #-lboost_system
CXXLINKS    += $(NCLINKER)
FLINKS      += $(NCLINKER)

.PHONY : all clean

all : $(DESTINATION)

clean :
	@$(RM) $(OBJPATH)/*.o
	@$(RM) $(OBJPATH)/*.d
	@$(RM) $(DESTINATION)
	@$(RM) $(OBJPATH)/*.so

doc :
	/Applications/Doxygen.app/Contents/Resources/doxygen ../doc/DocGen 

$(addprefix $(OBJPATH)/, $(DEPS)) : $(OBJPATH)/%.d : %.cpp | $(OBJPATH)
	@g++ -MM $(CPP11FLAG) $< $(addprefix -I, $(INCLUDES))> $@
	@sed -i="" '1s/^/..\/obj\//g' $@
	@rm $@=

-include $(addprefix $(OBJPATH)/, $(DEPS))

EXPV_LIB      :=  expv
EXPV_DIR      :=  $(ROOT_DIR)/$(OBJPATH)
EXPV_LINKS    :=  -L$(ROOT_DIR)/$(OBJPATH) -Wl,-rpath,$(EXPV_DIR) -l$(EXPV_LIB)
FLINKS        +=  $(EXPV_LINKS)
CXXLINKS      +=  $(EXPV_LINKS)

ifeq ($(MPI_INTEL), true)
$(DESTINATION) : $(addprefix $(OBJPATH)/, $(OBJS)) $(OBJPATH)/$(EXPV_LIB) | $(BINPATH)
	$(FC) -o $@ $(OBJPATH)/$(notdir $@).o $(filter-out $(APP_OBJ), $(addprefix $(OBJPATH)/, $(OBJS))) -L$(LIBPATH) $(addprefix -l,$(LIBS)) $(FLINKS)
else
$(DESTINATION) : $(addprefix $(OBJPATH)/, $(OBJS)) $(OBJPATH)/$(EXPV_LIB) | $(BINPATH)
	$(CC) -o $@ $(OBJPATH)/$(notdir $@).o $(filter-out $(APP_OBJ), $(addprefix $(OBJPATH)/, $(OBJS))) -L$(LIBPATH) $(addprefix -l,$(LIBS)) $(CXXLINKS)
endif

$(addprefix $(OBJPATH)/, %.o) : %.cpp | $(OBJPATH)
	$(CC) -c $< -o $@ $(CPPFLAGS)

$(OBJPATH) :
	@mkdir -p $(OBJPATH)

$(BINPATH) :
	@mkdir -p $(BINPATH)

$(OBJPATH)/$(EXPV_LIB):
	@make -C expv
