# Solar system ephemeris access library

## Repository contents

* `libephaccess` is a C library for calculating planetary positions and other
  data from [SPICE](https://naif.jpl.nasa.gov/naif/toolkit.html) files
* `ephcalculator` is a command-line utility based on `libephaccess`
* `python` directory contains Python wrapper for `libephaccess`
* `eph.sln` is a MS Visual Studio solution file for `libephaccess` and
  `ephcalculator`
* `setup.py` is a Distutils setup script for the Python wrapper

## Build

### libephaccess

On Windows, build the solution `eph.sln` in Visual Studio.

On Linux or OS X, `cd` to the `libephaccess` folder and run

```
make CONF=<configuration>
```

where `<configuration>` is one of:

- `ReleaseStatic32` (Linux 32-bit static library)
- `ReleaseShared32` (Linux 32-bit shared library)
- `ReleaseStatic64` (Linux 64-bit static library)
- `ReleaseShared64` (Linux 64-bit shared library)
- `ReleaseStatic64OSX` (OS X static library)
- `ReleaseShared64OSX` (OS X shared library)

### ephcalculator

On Windows, build the solution `eph.sln` in Visual Studio.

On Linux or OS X, `cd` to the `ephcalculator` folder and run

```
make CONF=<configuration>
```

where `<configuration>` is one of:

- `Release32` (Linux 32-bit executable)
- `Release64` (Linux 64-bit executable)
- `Release64OSX` (OS X executable)

### Python wrapper for libephaccess

#### Direct installation from repository

Direct installation from Gitlab repository is supported for both
Python 2.7 and Python 3:

```
pip install git+ssh://git@gitlab.iaaras.ru:10022/iaaras/ephemeris-access.git
```
(use `pip2` or `pip3` to distinguish between Python 2 and Python 3)

or

```
python -m pip install git+ssh://git@gitlab.iaaras.ru:10022/iaaras/ephemeris-access.git
```
(use `python2` or `python3` to distinguish between Python 2 and Python 3)

C compiler and linker are needed in the system to compile and link the C modules
into a shared library. On Linux or OS X, `gcc` or `clang` should suffice. On
Windows, different versions of Microsoft Visual Studio are required depending on
Python major version. For Python 3, MSVS 2015 or later (including “Express”
versions) is required. For Python 2, a
[special distribution](https://www.microsoft.com/en-us/download/details.aspx?id=44266)
of MSVS is required.

#### Manual installation

From parent folder of the downloaded `ephemeris-access` folder, run

```
pip install ephemeris-access/
```

or 

```
python -m pip install ephemeris-access/
```

(The above notes regarding Python 2/3 and C compiler apply here, too)


#### Usage without installation

`ephaccess` folder can be used as is, as long as it resides in the same
directory as user's Python program. But the shared library still must be built
beforehand, in this case with the following command:

```
python setup.py build_ext --inplace
```

(The above notes regarding C compiler apply here, too)

## Usage

### Ephemeris files

Different ephemeris files, as long as they are in SPICE formats (BSP and binary
PCK) can be loaded into `libephaccess`. EPM2017 ephemeris files, currently recommended
for usage, can be found on [IAA RAS FTP server](ftp://ftp.iaaras.ru/pub/epm/EPM2017/SPICE/).
[`epm2017.bsp`](ftp://ftp.iaaras.ru/pub/epm/EPM2017/SPICE/epm2017.bsp) contains coordinates
and velocities of the Sun, planets, and the Moon; and also TT-TDB.
[`moonlibr_epm2017.bpc`](ftp://ftp.iaaras.ru/pub/epm/EPM2017/SPICE/moonlibr_epm2017.bpc)
contains lunar libration angles. 

### API

* For `libephaccess`, see the comments in the [`ephaccess.h`](libephaccess/ephaccess.h) file. Also,
  [`ephcalculator.c`](ephcalculator/ephcalculator.c) is a good example of `libephaccess` usage.
* For `ephcalculator`, a build-in help message is provided.
* For Python wrapper, see [`python/test.py`](python/test.py) example script.

