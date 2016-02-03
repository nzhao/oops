# Please define environmental variables $(ARMA_LIB) and $(ARMA_INC_PATH) to enable armadillo 
# If you have Malab installed, pleas set $(HAS_MATLAB) = true, 
# and define $(MATLAB_LIB_PATH) and $(MATLAB_INC_PATH) in your .bashrc file.

DESTINATION := ../bin/cce
LIBS :=  $(ARMA_LIB) 
LIBPATH := .  
INCLUDES := .  $(ARMA_INC_PATH) 

RM := rm -f
PS = cpp

#CC=g++
CC = mpic++
CPPFLAGS = $(addprefix -I,$(INCLUDES))
CPPFLAGS += -g -DARMA_DONT_USE_WRAPPER 

ifeq ($(HAS_MATLAB), true)
	LIBS += mx mat eng
	CPPFLAGS += -DHAS_MATLAB
	LIBPATH := $(MATLAB_LIB_PATH)
	INCLUDES += $(MATLAB_INC_PATH)
endif

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

