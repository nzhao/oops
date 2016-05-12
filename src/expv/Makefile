OBJPATH       :=  ../../obj

EXPV_LIB      :=  libexpv.so

FC            := mpif77
ifeq ($(MPI_INTEL), true)
  FFLAGS      := -I${MKLROOT}/include
  FLINKS      := -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl -nofor-main -cxxlib
else
  FFLAGS      := -m64 
  CXXLINKS    := -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl -lgfortran
endif

CC          :=  mpic++
CFLAGS      :=  -m64 -I./include -parallel -L${MKLROOT}/lib/intel64

CFLAGS        +=  -fPIC
FFLAGS        +=  -fPIC

NCC         :=  nvcc
NCFLAGS     :=  -arch=sm_30 -O3 -DCUBLAS_GFORTRAN -ccbin g++ -Xcompiler "$(CFLAGS)" --compiler-options '-fPIC'
NCLINKER    :=  -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusparse -lboost_program_options -L/usr/local/lib -lboost_system

EXPOKIT_SRC   :=  expokit.f mataid.f DLARAN.f
KRYLOV_F_SRC  :=  krylov_expv_f.f
KRYLOV_C_SRC  :=  krylov_expv.c
ZKMV_MKL_F_SRC:=  hamvec_zgexpv_w_mkl_profile.f main_mkl.f
ZKMV_MKL_C_SRC:=  hamvec_func3.cpp
ZKMV_GPU_F_SRC:=  zgexpv_cache.f main_cache.f
ZKMV_GPU_C_SRC:=  hamvec_cuda3.cu get_grid.cu
CUDAAPI_DIR   :=  /usr/local/cuda-6.0/src
CUDAAPI_SRC   :=  fortran.c cusparse_fortran.c

vpath %.f   source source/expokit
vpath %.c   source $(CUDAAPI_DIR)
vpath %.cpp source
vpath %.cu  source

EXPOKIT_OBJ   :=  $(EXPOKIT_SRC:%.f=$(OBJPATH)/%.o)
KRYLOV_F_OBJ  :=  $(KRYLOV_F_SRC:%.f=$(OBJPATH)/%.o)
KRYLOV_C_OBJ  :=  $(KRYLOV_C_SRC:%.c=$(OBJPATH)/%.o)
ZKMV_MKL_F_OBJ:=  $(ZKMV_MKL_F_SRC:%.f=$(OBJPATH)/%.o)
ZKMV_MKL_C_OBJ:=  $(ZKMV_MKL_C_SRC:%.cpp=$(OBJPATH)/%.o)
ZKMV_GPU_F_OBJ:=  $(ZKMV_GPU_F_SRC:%.f=$(OBJPATH)/%.o)
ZKMV_GPU_C_OBJ:=  $(ZKMV_GPU_C_SRC:%.cu=$(OBJPATH)/%.o)
CUDAAPI_OBJ   :=  $(CUDAAPI_SRC:%.c=$(OBJPATH)/%.o)

KRYLOV_OBJ    :=  $(KRYLOV_F_OBJ) $(KRYLOV_C_OBJ)
ZKMV_MKL_OBJ  :=  $(ZKMV_MKL_F_OBJ) $(ZKMV_MKL_C_OBJ)
ZKMV_GPU_OBJ  :=  $(ZKMV_GPU_F_OBJ) $(ZKMV_GPU_C_OBJ)

$(EXPV_LIB):$(EXPOKIT_OBJ) $(KRYLOV_OBJ) $(ZKMV_MKL_OBJ) $(ZKMV_GPU_OBJ) $(CUDAAPI_OBJ) | $(OBJPATH)
	@$(NCC) --shared -o $(OBJPATH)/$(EXPV_LIB) $(EXPOKIT_OBJ) $(KRYLOV_OBJ) $(ZKMV_MKL_OBJ) $(ZKMV_GPU_OBJ) $(CUDAAPI_OBJ)
#	@ar -rcs $(EXPV_LIB) $(EXPOKIT_OBJ) $(KRYLOV_OBJ) $(ZKMV_MKL_OBJ) $(ZKMV_GPU_OBJ) $(CUDAAPI_OBJ)

$(EXPOKIT_OBJ):$(OBJPATH)/%.o:%.f | $(OBJPATH)
	$(FC) -c $< -o $@ $(FFLAGS)

$(KRYLOV_F_OBJ):$(OBJPATH)/%.o:%.f | $(OBJPATH)
	$(FC) -c $< -o $@ $(FFLAGS)

$(KRYLOV_C_OBJ):$(OBJPATH)/%.o:%.c | $(OBJPATH)
	$(CC) -c $< -o $@ $(CFLAGS)

$(ZKMV_MKL_F_OBJ):$(OBJPATH)/%.o:%.f | $(OBJPATH)
	$(FC) -c $< -o $@ $(FFLAGS)

$(ZKMV_MKL_C_OBJ):$(OBJPATH)/%.o:%.cpp | $(OBJPATH)
	$(CC) -c $< -o $@ $(CFLAGS)

$(ZKMV_GPU_F_OBJ):$(OBJPATH)/%.o:%.f | $(OBJPATH)
	$(FC) -c $< -o $@ $(FFLAGS)

$(ZKMV_GPU_C_OBJ):$(OBJPATH)/%.o:%.cu | $(OBJPATH)
	$(NCC) -c $< -o $@ $(NCFLAGS)

$(CUDAAPI_OBJ):$(OBJPATH)/%.o:%.c | $(OBJPATH)
	$(NCC) -c $< -o $@ $(NCFLAGS)

