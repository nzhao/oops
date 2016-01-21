DESTINATION := ../bin/cce
LIBS :=  $(ARMA_LIB) mx mat eng
LIBPATH := $(MATLAB_LIB_PATH) 
INCLUDES := . $(ARMA_INC_PATH) $(MATLAB_INC_PATH) 

RM := rm -f
PS = cpp

CC=g++
CPPFLAGS = -g
CPPFLAGS += $(addprefix -I,$(INCLUDES))

SOURCE := $(wildcard app/*.$(PS)) $(wildcard source/*/*.$(PS)) $(wildcard source/*/*/*.$(PS))
OBJS := $(patsubst %.$(PS),%.o,$(notdir $(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))

OBJPATH := ../obj
vpath %.cpp $(sort $(dir $(SOURCE)))

.PHONY : all clean

all : $(DESTINATION)

clean :
	@$(RM) $(OBJPATH)/*.o
	@$(RM) $(OBJPATH)/*.d
	@$(RM) $(DESTINATION)

doc :
	/Applications/Doxygen.app/Contents/Resources/doxygen ../doc/DocGen 

$(addprefix $(OBJPATH)/, $(DEPS)) : $(OBJPATH)/%.d : %.cpp 
	@$(CC) -MM $< -I. -I$(MATLAB_INC_PATH) > $@
	@sed -i="" '1s/^/..\/obj\//g' $@
	@rm $@=

-include $(addprefix $(OBJPATH)/, $(DEPS))

$(DESTINATION) : $(addprefix $(OBJPATH)/, $(OBJS))
	$(CC) -o $(DESTINATION) $(addprefix $(OBJPATH)/, $(OBJS)) -L$(LIBPATH) $(addprefix -l,$(LIBS))

$(addprefix $(OBJPATH)/, %.o) : %.cpp
	$(CC) -c $< -o $@ $(CPPFLAGS)
