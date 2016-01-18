DESTINATION := ../bin/cce
LIBS :=  blas lapack
INCLUDES := . $(ARMA_PATH)

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

.PHONY : all deps clean rebuild

all : $(DESTINATION)

clean :
	@$(RM) $(OBJPATH)/*.o
	@$(RM) $(OBJPATH)/*.d
	@$(RM) $(DESTINATION)

rebuild: clean all 

$(addprefix $(OBJPATH)/, $(DEPS)) : $(OBJPATH)/%.d : %.cpp 
	@$(CC) -MM -MMD $< -o $(OBJPATH)/$@ -I.
	@sed -i '1s/^/..\/obj\//g' $(OBJPATH)/$@ 

-include $(addprefix $(OBJPATH)/, $(DEPS))

$(DESTINATION) : $(addprefix $(OBJPATH)/, $(OBJS))
	$(CC) -o $(DESTINATION) $(addprefix $(OBJPATH)/, $(OBJS)) $(addprefix -l,$(LIBS))

$(addprefix $(OBJPATH)/, %.o) : %.cpp
	$(CC) -c $< -o $@ $(CPPFLAGS)
