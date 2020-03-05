NAME=constantq~
cflags=-DFLEXT_INLINE -DFLEXT_ATTRIBUTES=1

# source files
$(NAME).class.sources = $(wildcard *.cpp)

# help files
datafiles = $(wildcard *-help.pd)

# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
PDLIBBUILDER_DIR=./pd-lib-builder/
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
