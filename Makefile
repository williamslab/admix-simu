CPPSRCS= main.cc marker.cc
CSRCS= # admutils.c mcio.c # files for Nick's I/O
OBJS= $(patsubst %.cc,%.o,$(CPPSRCS)) $(patsubst %.c,%.o,$(CSRCS))
EXEC= mixer

GPP = g++
GCC = gcc
DEFINES= 
#CFLAGS = -g -Wall $(DEFINES)
# optimized; remove asserts
#CFLAGS = -O2 -Wall -DNDEBUG $(DEFINES)
CFLAGS = -O2 -Wall $(DEFINES)
# profiling: run program normally then do:
# (note: I haven't read about the options below, I just found them in a
#  Dr. Dobb's article.)
#     `gprof -b -z hapi gmon.out`
#CFLAGS = -pg -O2 -Wall $(DEFINES)

LIBS = 

# dependency variables / commands
DEPDIR = .deps
df = $(DEPDIR)/$(*F)

all: $(EXEC)

$(EXEC): $(OBJS) $(HEADERS)
	$(GPP) -o $(EXEC) $(OBJS) $(CFLAGS) $(LIBS)

# This way of building dependencies (per-file) described at
# http://make.paulandlesley.org/autodep.html

.c.o:
	@mkdir -p $(DEPDIR)
	$(GCC) -MMD $(CFLAGS) -o $@ -c $<
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d

.cc.o:
	@mkdir -p $(DEPDIR)
	$(GPP) -MMD $(CFLAGS) -o $@ -c $<
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d


# include the .P dependency files, but don't warn if they don't exist (the -)
-include $(CPPSRCS:%.cc=$(DEPDIR)/%.P)
-include $(CSRCS:%.c=$(DEPDIR)/%.P)
# The following applies if we don't use a dependency directory:
#-include $(SRCS:.cc=.P)

tags: $(SRCS) *.h
	ctags --language-force=c++ --extra=+q --fields=+i --excmd=n *.c *.cc *.h

clean:
	rm -f $(EXEC) $(OBJS)

clean-deps:
	rm -f $(DEPDIR)/*.P
