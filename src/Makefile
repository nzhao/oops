ARMA_INC_PATH   := ${HOME}/code/armadillo/include
matlabroot      := /usr/local/MATLAB/R2013a

ARMA_LIB        := blas lapack
MATLAB_INC_PATH := $(matlabroot)/extern/include
MATLAB_LIB_PATH := $(matlabroot)/bin/glnxa64

DESTINATION := ../bin/cce
LIBS :=  $(ARMA_LIB) mx mat eng
LIBPATH := $(MATLAB_LIB_PATH) 
INCLUDES := . $(ARMA_INC_PATH) $(MATLAB_INC_PATH) 

RM := rm -f
PS = cpp

#CC=g++
CC = mpic++
CPPFLAGS = -g -DARMA_DONT_USE_WRAPPER
CPPFLAGS += $(addprefix -I,$(INCLUDES))

SOURCE := $(wildcard app/*.$(PS)) $(wildcard source/*/*.$(PS)) $(wildcard source/*/*/*.$(PS))
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
#	@rm -rf $(OBJPATH)
#	@rm -rf $(BINPATH)

doc :
	/Applications/Doxygen.app/Contents/Resources/doxygen ../doc/DocGen 

$(addprefix $(OBJPATH)/, $(DEPS)) : $(OBJPATH)/%.d : %.cpp | $(OBJPATH)
	@$(CC) -MM $< -I. -I$(MATLAB_INC_PATH) > $@
	@sed -i="" '1s/^/..\/obj\//g' $@
	@rm $@=

-include $(addprefix $(OBJPATH)/, $(DEPS))

$(DESTINATION) : $(addprefix $(OBJPATH)/, $(OBJS)) | $(BINPATH)
	$(CC) -o $(DESTINATION) $(addprefix $(OBJPATH)/, $(OBJS)) -L$(LIBPATH) $(addprefix -l,$(LIBS))

$(addprefix $(OBJPATH)/, %.o) : %.cpp | $(OBJPATH)
	$(CC) -c $< -o $@ $(CPPFLAGS)

$(OBJPATH) :
	@mkdir -p $(OBJPATH)

$(BINPATH) :
	@mkdir -p $(BINPATH)

