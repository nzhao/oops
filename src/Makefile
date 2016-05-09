# Please define environmental variables $(ARMA_LIB) and $(ARMA_INC_PATH) to enable armadillo 
# If you have Malab installed, pleas set $(HAS_MATLAB) = true, 
# and define $(MATLAB_LIB_PATH) and $(MATLAB_INC_PATH) in your .bashrc file.
# If you want to link the program by INTEL MKL libaray, set $(USE_MKL)=true in your .bashrc file.
# If you are using INTEL MPI, set $(MPI_INTEL)=true in your .bashrc file.
BINPATH := ../bin
OBJPATH := ../obj

#APP := EnsembleCCE SingleSampleCCE MatExpTest test
APP := EnsNVeBath EnsNVnBath SingleSampleCCE_NV_C13 MatExpTest test 

DESTINATION := $(addprefix $(BINPATH)/, $(APP))
LIBPATH := .  
INCLUDES := .  $(ARMA_INC_PATH) 

RM := rm -f
PS = cpp

CC = mpic++
CPPFLAGS = $(addprefix -I,$(INCLUDES))

ifeq ($(USE_MKL), true)
	LIBS := 
	CXXLINKS  :=  -parallel -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_ilp64 -lpthread -lm -ldl
	CPPFLAGS  +=  -DARMA_BLAS_LONG -DMKL_ILP64 -m64 -I${MKLROOT}/include -DARMA_DONT_USE_WRAPPER 
else
	LIBS :=  $(ARMA_LIB) 
	CXXLINKS  := 
	CPPFLAGS += -g -DARMA_DONT_USE_WRAPPER 
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
APP_OBJ := $(addsuffix .o, $(addprefix $(OBJPATH)/, $(APP)) )

vpath %.cpp $(sort $(dir $(SOURCE)))

EXPOKIT_SRC := $(wildcard expokit/*.f)
EXPOKIT_OBJ := $(OBJPATH)/expokit.o $(OBJPATH)/mataid.o
KRYLOV_OBJ  := $(OBJPATH)/krylov_expv_f.o $(OBJPATH)/krylov_expv.o $(OBJPATH)/DLARAN.o
vpath %.f $(sort $(dir $(EXPOKIT_SRC)))
vpath %.c $(sort $(dir wildcard source/math/*.c))
FC          := mpif77
ifeq ($(MPI_INTEL), true)
  FFLAGS      := -I${MKLROOT}/include
  FLINKS      := -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl -nofor-main -cxxlib
else
  FFLAGS      := -m64 
#  CXXLINKS    :=  -L${GFORTRAN_PATH} -lpthread -lm -ldl -lgfortran
  CXXLINKS    := -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl -lgfortran
endif

ZKMV_MKL_SRC  := $(wildcard source/math/*.f)
ZKMV_MKL_OBJ  := $(OBJPATH)/main_mkl.o $(OBJPATH)/hamvec_zgexpv_w_mkl_profile.o
vpath %.f $(sort $(dir $(ZKMV_MKL_SRC)))
ZKMV_GPU_SRC  := $(wildcard source/math/*.cu)
ZKMV_GPU_OBJ  := $(OBJPATH)/main_cache.o $(OBJPATH)/zgexpv_cache.o $(OBJPATH)/hamvec_cuda3.o $(OBJPATH)/get_grid.o
CUDAAPI_DIR   := /usr/local/cuda-6.0/src
CUDAAPI_SRC   := fortran.c cusparse_fortran.c
CUDAAPI_OBJ   := $(OBJPATH)/fortran.o $(OBJPATH)/cusparse_fortran.o
NCC         := nvcc
NCFLAGS     := -arch=sm_30 -O3 -DCUBLAS_GFORTRAN -ccbin g++ -Xcompiler "$(CPPFLAGS)"
NCLINKER    := -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusparse -lboost_program_options -L/usr/local/lib -lboost_system
vpath %.cu $(sort $(dir $(ZKMV_GPU_SRC)))
vpath %.c  $(CUDAAPI_DIR)
CXXLINKS    += $(NCLINKER)

.PHONY : all clean

all : $(DESTINATION)

clean :
	@$(RM) $(OBJPATH)/*.o
	@$(RM) $(OBJPATH)/*.d
	@$(RM) $(DESTINATION)

doc :
	/Applications/Doxygen.app/Contents/Resources/doxygen ../doc/DocGen 

$(addprefix $(OBJPATH)/, $(DEPS)) : $(OBJPATH)/%.d : %.cpp | $(OBJPATH)
	@g++ -MM $< $(addprefix -I, $(INCLUDES))> $@
	@sed -i="" '1s/^/..\/obj\//g' $@
	@rm $@=

-include $(addprefix $(OBJPATH)/, $(DEPS))

ifeq ($(MPI_INTEL), true)
$(DESTINATION) : $(addprefix $(OBJPATH)/, $(OBJS)) $(KRYLOV_OBJ) $(EXPOKIT_OBJ) $(ZKMV_MKL_OBJ) $(ZKMV_GPU_OBJ) $(CUDAAPI_OBJ) | $(BINPATH)
	$(FC) -o $@ $(OBJPATH)/$(notdir $@).o $(filter-out $(APP_OBJ), $(addprefix $(OBJPATH)/, $(OBJS))) $(KRYLOV_OBJ) $(EXPOKIT_OBJ) $(ZKMV_MKL_OBJ) $(ZKMV_GPU_OBJ) $(CUDAAPI_OBJ) -L$(LIBPATH) $(addprefix -l,$(LIBS)) $(FLINKS)
else
$(DESTINATION) : $(addprefix $(OBJPATH)/, $(OBJS)) $(KRYLOV_OBJ) $(EXPOKIT_OBJ) $(ZKMV_MKL_OBJ) $(ZKMV_GPU_OBJ) $(CUDAAPI_OBJ) | $(BINPATH)
	$(CC) -o $@ $(OBJPATH)/$(notdir $@).o $(filter-out $(APP_OBJ), $(addprefix $(OBJPATH)/, $(OBJS))) $(KRYLOV_OBJ) $(EXPOKIT_OBJ) $(ZKMV_MKL_OBJ) $(ZKMV_GPU_OBJ) $(CUDAAPI_OBJ) -L$(LIBPATH) $(addprefix -l,$(LIBS)) $(CXXLINKS)
endif

$(addprefix $(OBJPATH)/, %.o) : %.cpp | $(OBJPATH)
	$(CC) -c $< -o $@ $(CPPFLAGS)

$(OBJPATH) :
	@mkdir -p $(OBJPATH)

$(BINPATH) :
	@mkdir -p $(BINPATH)

$(OBJPATH)/expokit.o : expokit.f | $(OBJPATH)
	$(FC) -c $< -o $@ $(FFLAGS)

$(OBJPATH)/mataid.o : mataid.f | $(OBJPATH)
	$(FC) -c $< -o $@ $(FFLAGS)

$(OBJPATH)/main_mkl.o : main_mkl.f | $(OBJPATH)
	$(FC) $< -c -o $@ $(FFLAGS)

$(OBJPATH)/hamvec_zgexpv_w_mkl_profile.o : hamvec_zgexpv_w_mkl_profile.f | $(OBJPATH)
	$(FC) $< -c -o $@ $(FFLAGS)

$(OBJPATH)/main_cache.o : main_cache.f | $(OBJPATH)
	$(FC) $< -c -o $@ $(FFLAGS)

$(OBJPATH)/zgexpv_cache.o : zgexpv_cache.f | $(OBJPATH)
	$(FC) $< -c -o $@ $(FFLAGS)

$(OBJPATH)/hamvec_cuda3.o : hamvec_cuda3.cu | $(OBJPATH)
	$(NCC) $< -c -o $@ $(NCFLAGS)

$(OBJPATH)/get_grid.o : get_grid.cu | $(OBJPATH)
	$(NCC) $< -c -o $@ $(NCFLAGS)

$(OBJPATH)/fortran.o : $(CUDAAPI_DIR)/fortran.c | $(OBJPATH)
	$(NCC) $< -c -o $@ $(NCFLAGS)

$(OBJPATH)/cusparse_fortran.o : $(CUDAAPI_DIR)/cusparse_fortran.c | $(OBJPATH)
	$(NCC) $< -c -o $@ $(NCFLAGS)

$(OBJPATH)/krylov_expv_f.o : krylov_expv_f.f | $(OBJPATH)
	$(FC) -c $< -o $@ $(FFLAGS)

$(OBJPATH)/DLARAN.o : DLARAN.f | $(OBJPATH)
	$(FC) -c $< -o $@ $(FFLAGS)

$(OBJPATH)/krylov_expv.o : krylov_expv.c | $(OBJPATH)
	$(CC) -c $< -I./include/math -o $@ $(FFLAGS)

