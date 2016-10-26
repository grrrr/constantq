

### Makefile.pdlibbuilder ###

Helper makefile for Pure Data external libraries.
Written by Katja Vetter March-June 2015 for the public domain. No warranties.
Inspired by Hans Christoph Steiner's Makefile Template and Stephan Beal's
ShakeNMake.

GNU make version >= 3.81 required.


### characteristics ###


* defines build settings based on autodetected OS and architecture
* defines rules to build Pd class- or lib executables from C or C++ sources
* defines rules for libdir installation
* defines convenience targets for developer and user
* evaluates implicit dependencies for non-clean builds


### basic usage ###


In your Makefile, define your Pd lib name and class files, and include
Makefile.pdlibbuilder at the end of the Makefile. Like so:


      # Makefile for mylib

      lib.name = mylib

      class.sources = myclass1.c myclass2.c

      datafiles = myclass1-help.pd myclass2-help.pd README.txt LICENSE.txt

      include Makefile.pdlibbuilder


For files in class.sources it is assumed that class name == source file
basename. The default target builds all classes as individual executables
with Pd's default extension for the platform. For anything more than the
most basic usage, read the documentation sections in Makefile.pdlibbuilder.


### paths ###


Makefile.pdlibbuilder >= v0.4.0 supports pd path variables which can be
defined not only as make command argument but also in the environment, to
override platform-dependent defaults:

PDDIR:
Root directory of 'portable' pd package. When defined, PDINCLUDEDIR and 
PDBINDIR will be evaluated as $(PDDIR)/src and $(PDDIR)/bin.

PDINCLUDEDIR:
Directory where Pd API m_pd.h should be found, and other Pd header files.
Overrides the default search path.

PDBINDIR:
Directory where pd.dll should be found for linking (Windows only). Overrides
the default search path.

PDLIBDIR:
Root directory for installation of Pd library directories. Overrides the
default install location.


### documentation ###


This README.md provides only basic information. A large comment section inside
Makefile.pdlibbuilder lists and explains the available user variables, default
paths, and targets. The internal documentation reflects the exact functionality
of the particular version. A tips&tricks page is in the works. 


### examples ###


Here are a few projects using the Makefile.pdlibbuilder approach:

https://github.com/pure-data/helloworld

https://github.com/electrickery/pd-cyclone (stable)

https://github.com/porres/pd-cyclone (experimental)

https://git.iem.at/pd/iemguts

https://git.iem.at/pd/iemnet

https://git.iem.at/pd/iem_ambi

https://git.iem.at/pd/mediasettings

https://git.iem.at/pd-gui/punish

https://github.com/residuum/PuRestJson

More examples will be referenced here when they are available. 
