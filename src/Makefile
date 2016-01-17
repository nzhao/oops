DESTINATION := ../bin/cce
LIBS :=  blas lapack
INCLUDES := . /home/nzhao/code/armadillo-6.400.3/include

RM := rm -f
#C,CC或CPP文件的后缀
PS=cpp

# GNU Make的隐含变量定义
CC=g++
CPPFLAGS = -g# -Wall -O3 -march=i486
CPPFLAGS += $(addprefix -I,$(INCLUDES))
#CPPFLAGS += -MMD

# #以下部分无需修改
#SOURCE := $(wildcard *.$(PS))
#OBJS := $(patsubst %.$(PS),%.o,$(SOURCE))
SOURCE := $(wildcard app/*.$(PS)) $(wildcard source/*/*.$(PS)) $(wildcard source/*/*/*.$(PS))
OBJS := $(patsubst %.$(PS),%.o,$(notdir $(SOURCE)))

vpath %.cpp $(sort $(dir $(SOURCE)))

DEPS := $(patsubst %.o,%.d,$(OBJS))
MISSING_DEPS := $(filter-out $(wildcard $(DEPS)),$(DEPS))
MISSING_DEPS_SOURCES := $(wildcard $(patsubst %.d,%.$(PS),$(MISSING_DEPS)))

.PHONY : all deps objs clean rebuild

all : $(DESTINATION)

#deps : $(DEPS)
	#$(CC) -MM -MMD $(SOURCE) -I.

objs : $(OBJS)

clean :
	@$(RM) *.o
	@$(RM) *.d
	@$(RM) $(DESTINATION)

rebuild: clean all 

#ifneq ($(MISSING_DEPS),) 
#$(MISSING_DEPS) : 
#	@$(RM) $(patsubst %.d,%.o,$@)
#endif

$(DEPS) : %.d : %.cpp 
	$(CC) -MM -MMD $< -o ../obj/$@ -I.

include $(DEPS)

$(DESTINATION) : $(OBJS)
	$(CC) -o $(DESTINATION) $(OBJS) $(addprefix -l,$(LIBS))

$(OBJS) : %.o : %.cpp
	$(CC) -c $< -o $@ $(CPPFLAGS)
