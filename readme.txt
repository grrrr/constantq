constantq~ - spectral analysis for Pure Data and Max

Copyright (c) 2010-2012 Thomas Grill (gr@grrrr.org)
For information on usage and redistribution, and for a DISCLAIMER OF ALL
WARRANTIES, see the file, "license.txt," in this distribution.  

Donations for further development of the package are highly appreciated.
Visit https://www.paypal.com/xclick/business=gr%40grrrr.org&item_name=constantq&no_note=1&tax=0&currency_code=EUR

----------------------------------------------------------------------------

Goals/features of the package:

The constantq~ object spectrally analyzes an audio stream for defined frequency bands. 
It can be seen as an optimized implementation of a band pass filterbank.

----------------------------------------------------------------------------

IMPORTANT INFORMATION for all PD users:

Put the pd-msvc/constantq~.dll, pd-linux/constantq~.pd_linux or pd-darwin/constantq~.pd_darwin file
into the extra folder of the PD installation.

Put the constantq~-help.pd file into the doc\5.reference subfolder of your PD installation.

----------------------------------------------------------------------------

IMPORTANT INFORMATION for all Max/MSP users:

For Mac OSX put the max-darwin/constantq~.mxo file into the folder 
/Applications/Max*/Cycling '74/externals

For Windows put the max-msvc\ constantq~.mxe file into the folder
C:\program files\common files\Cycling '74\externals (english version)

Put the constantq~.maxhelp file into the max-help folder.

============================================================================


BUILDING from source
--------------------

You will need the flext C++ layer for PD and Max/MSP externals to compile this.
See http://grrrr.org/ext/flext
Download, install and compile the package.
Afterwards you can proceed with building this external.

Further dependencies are
- blitz++ (http://blitz.sourceforge.net)
- FFTW3 (http://www.fftw.org)

Please note that in terms of compilation blitz++ is a very demanding piece of software.
Older compilers are regularly choking on it. Reported to work flawlessly is e.g. clang >= 3.1.


pd/Max - Windows - Microsoft Visual C, Borland C++, MinGW:
----------------------------------------------------------
Start a command shell with your eventual build environment
(e.g. run vcvars32.bat for Microsoft Visual Studio)

then run
 ..\flext\build.bat
(you would have to substitute ..\flext with the respective path to the flext package)


pd/Max - OSX/Linux - GCC:
-------------------------
From a shell run
bash ../flext/build.sh
(you would have to substitute ../flext with the respective path to the flext package)


============================================================================

Version history:

0.2 - changed methods for frequency scales: actually INCOMPATIBLE to version 0.1 

0.1 - first release

