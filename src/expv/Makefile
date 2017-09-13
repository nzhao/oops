OBJPATH       :=  ../../obj

EXPV_LIB      :=  libexpv.so

FC            := mpif77
ifeq ($(MPI_INTEL), true)
  FFLAGS      := -I${MKLROOT}/include -i8 -qopenmp
  FLINKS      := -nofor-main -cxxlib -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl
else
  FFLAGS      := -m64 -fdefault-integer-8 -fopenmp
  CXXLINKS    := -lgfortran -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl
endif

CC          :=  mpic++
CFLAGS      :=  -fopenmp -m64 -I./include -L${MKLROOT}/lib/intel64 -DUSE_MKL -I${MKLROOT}/include
#CFLAGS      :=  -m64 -I./include -parallel -L${MKLROOT}/lib/intel64

CFLAGS        +=  -fPIC
FFLAGS        +=  -fPIC

NCC         :=  nvcc
NCFLAGS     :=  -arch=sm_30 -O3 -DCUBLAS_GFORTRAN -ccbin g++ -Xcompiler "$(CFLAGS)" --compiler-options '-fPIC'
NCLINKER    :=  -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusparse -lboost_program_options -L/usr/local/lib -lboost_system

MPI_ZKEXPV_SRC:=  mpi_zkexpv.f mpi_zcsrkexpv.f mpi_hzcsrkexpv.f
MPI_HAMVEC_SRC:=  mpi_hamvec.c mpi_zcsrkhmv.c ham.c
EXPOKIT_SRC   :=  expokit.f mataid.f DLARAN.f
KRYLOV_F_SRC  :=  krylov_expv_f.f
KRYLOV_C_SRC  :=  krylov_expv.c
ZKMV_MKL_F_SRC:=  hamvec_zgexpv_w_mkl_profile.f main_mkl.f
ZKMV_MKL_C_SRC:=  hamvec_func3.cpp
ZKMV_GPU_F_SRC:=  zgexpv_cache.f main_cache.f
ZKMV_GPU_C_SRC:=  hamvec_cuda3.cu get_grid.cu
CUDAAPI_DIR   :=  /usr/local/cuda/src
CUDAAPI_SRC   :=  fortran.c cusparse_fortran.c

vpath %.f   source source/expokit
vpath %.c   source $(CUDAAPI_DIR)
vpath %.cpp source
vpath %.cu  source

MPI_ZKEXPV_OBJ:=  $(MPI_ZKEXPV_SRC:%.f=$(OBJPATH)/%.o)
MPI_HAMVEC_OBJ:=  $(MPI_HAMVEC_SRC:%.c=$(OBJPATH)/%.o)
EXPOKIT_OBJ   :=  $(EXPOKIT_SRC:%.f=$(OBJPATH)/%.o)
KRYLOV_F_OBJ  :=  $(KRYLOV_F_SRC:%.f=$(OBJPATH)/%.o)
KRYLOV_C_OBJ  :=  $(KRYLOV_C_SRC:%.c=$(OBJPATH)/%.o)
ZKMV_MKL_F_OBJ:=  $(ZKMV_MKL_F_SRC:%.f=$(OBJPATH)/%.o)
ZKMV_MKL_C_OBJ:=  $(ZKMV_MKL_C_SRC:%.cpp=$(OBJPATH)/%.o)
ZKMV_GPU_F_OBJ:=  $(ZKMV_GPU_F_SRC:%.f=$(OBJPATH)/%.o)
ZKMV_GPU_C_OBJ:=  $(ZKMV_GPU_C_SRC:%.cu=$(OBJPATH)/%.o)
CUDAAPI_OBJ   :=  $(CUDAAPI_SRC:%.c=$(OBJPATH)/%.o)

MPI_OBJ       :=  $(MPI_ZKEXPV_OBJ) $(MPI_HAMVEC_OBJ)
KRYLOV_OBJ    :=  $(KRYLOV_F_OBJ) $(KRYLOV_C_OBJ)
ZKMV_MKL_OBJ  :=  $(ZKMV_MKL_F_OBJ) $(ZKMV_MKL_C_OBJ)
ZKMV_GPU_OBJ  :=  $(ZKMV_GPU_F_OBJ) $(ZKMV_GPU_C_OBJ)

$(EXPV_LIB):$(MPI_OBJ) $(EXPOKIT_OBJ) $(KRYLOV_OBJ) $(ZKMV_MKL_OBJ) $(ZKMV_GPU_OBJ) $(CUDAAPI_OBJ) | $(OBJPATH)
	$(NCC) --shared -o $(OBJPATH)/$(EXPV_LIB) $(MPI_OBJ) $(EXPOKIT_OBJ) $(KRYLOV_OBJ) $(ZKMV_MKL_OBJ) $(ZKMV_GPU_OBJ) $(CUDAAPI_OBJ)
#	@ar -rcs $(EXPV_LIB) $(EXPOKIT_OBJ) $(KRYLOV_OBJ) $(ZKMV_MKL_OBJ) $(ZKMV_GPU_OBJ) $(CUDAAPI_OBJ)

$(MPI_HAMVEC_OBJ):$(OBJPATH)/%.o:%.c | $(OBJPATH)
	$(CC) -c $< -o $@ $(CFLAGS)

$(MPI_ZKEXPV_OBJ):$(OBJPATH)/%.o:%.f | $(OBJPATH)
	$(FC) -c $< -o $@ $(FFLAGS)

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

$(OBJPATH) :
	@mkdir -p $(OBJPATH)

