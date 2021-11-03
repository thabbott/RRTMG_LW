# Set search path and source file list
DIRS := $(shell cat Filepath)
VPATH := $(foreach dir,$(DIRS),$(wildcard $(dir)))
SOURCES := $(shell cat Srcfiles)
OBJS := $(addsuffix .o, $(basename $(SOURCES)))
.SUFFIXES: .f90 .o

all: $(BASE)/rrtmg_lw_wrapper.f90 $(OBJS)
	$(F2PY) -I. $(OBJS) -m rrtmg_lw_wrapper $(BASE)/rrtmg_lw_wrapper.f90

.f90.o:
	$(FC) $(FFLAGS) $<

include Depends
