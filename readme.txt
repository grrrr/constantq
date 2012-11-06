constantq~ - Constant-Q spectral analysis

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

Put the help-constantq~.pd file into the doc\5.reference subfolder of your PD installation.

----------------------------------------------------------------------------

IMPORTANT INFORMATION for all Max/MSP users:

For Mac OSX put the max-osx/constantq~.mxd file into the folder 
/Library/Application Support/Cycling '74/externals

For Windows put the max-msvc\ constantq~.mxe file into the folder
C:\program files\common files\Cycling '74\externals (english version)

Put the constantq~.help file into the max-help folder.

============================================================================


BUILDING from source
--------------------

You will need the flext C++ layer for PD and Max/MSP externals to compile this.
See http://grrrr.org/ext/flext
Download, install and compile the package.
Afterwards you can proceed with building this external.


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

