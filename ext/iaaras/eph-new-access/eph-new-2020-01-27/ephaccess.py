#------------------------------------------------------------------------------

import sys
import os
import platform

from ctypes import (c_int, c_char_p, c_double, c_void_p,
                    POINTER, CDLL, RTLD_GLOBAL)

#------------------------------------------------------------------------------

EPH_AU = 1
EPH_KM = 2
EPH_SEC = 3
EPH_DAY = 4

EPH_SSB        = 0
EPH_MERCURY    = 1
EPH_VENUS      = 2
EPH_EMB        = 3
EPH_MARS_BC    = 4
EPH_JUPITER_BC = 5
EPH_SATURN_BC  = 6
EPH_URANUS_BC  = 7
EPH_NEPTUNE_BC = 8
EPH_PLUTO_BC   = 9
EPH_SUN       = 10
EPH_MOON     = 301
EPH_EARTH    = 399

EPH_TT_MINUS_TDB = 1000000001

_ENCODING = 'utf-8'

#------------------------------------------------------------------------------

class EphAccessException(Exception):
    def __init__(self, value):
        if sys.version_info > (3, 0) and type(value) is bytes:
            self.value = value.decode(_ENCODING)
        else:
            self.value = value

    def __str__(self):
        return self.value

class EphAccess(object):

    _crt = None
    _lib = None

    @staticmethod
    def _init_static(path=None):
        
        if not path:
            cur_file = os.path.abspath(__file__)
            dirname = os.path.dirname(cur_file)
            if not dirname:
                dirname = '.'
            path = dirname + '/lib'
            
        if os.name == 'posix' and not platform.mac_ver()[0] and not platform.system().startswith("CYGWIN"):
            arch = platform.architecture()[0]
            path += "/linux"
            if arch == '32bit':
                path += "/x86"
            elif arch == '64bit':
                path += "/x64"
            else:
                raise EphAccessException("unknown platform " + arch)
            
            EphAccess._lib = CDLL(path + "/libephaccess.so", mode=RTLD_GLOBAL)
            
        elif os.name == 'nt' or platform.system().startswith("CYGWIN"):
            arch = platform.architecture()[0]
            path += "/windows"
            if arch == '32bit':
                path += "/x86"
            elif arch == '64bit':
                path += "/x64"
            else:
                raise EphAccessException("unknown platform " + arch)
            if os.path.exists(path + "/msvcr100.dll"):
                EphAccess._crt = CDLL(path + "/msvcr100.dll")
            if os.path.exists(path + "/msvcr110.dll"):
                EphAccess._crt = CDLL(path + "/msvcr110.dll")
            if os.path.exists(path + "/msvcr120.dll"):
                EphAccess._crt = CDLL(path + "/msvcr120.dll")
            if os.path.exists(path + "/vcruntime140.dll"):
                EphAccess._crt = CDLL(path + "/vcruntime140.dll")
            EphAccess._lib = CDLL(path + "/ephaccess.dll")
        elif platform.mac_ver()[0]:
            path += "/osx"
            EphAccess._lib = CDLL(path + "/libephaccess.dylib", mode=RTLD_GLOBAL)
        else:
            raise EphAccessException('unsupported OS: ' + os.name)
        _lib = EphAccess._lib
        _lib.ephObjectByName.restype = c_int
        _lib.ephObjectByName.argtypes = [c_char_p]
        _lib.ephCreate.restype = c_void_p
        _lib.ephCreate.argtypes = None
        _lib.ephLastError.restype = c_char_p
        _lib.ephLastError.argtypes = [c_void_p]
        _lib.ephLoadFile.restype = c_int
        _lib.ephLoadFile.argtypes = [c_void_p, c_char_p]
        _lib.ephSetTimeUnits.restype = c_int
        _lib.ephSetTimeUnits.argtypes = [c_void_p, c_int]
        _lib.ephSetDistanceUnits.restype = c_int
        _lib.ephSetDistanceUnits.argtypes = [c_void_p, c_int]
        _lib.ephGetLeftmostJD.restype = c_double
        _lib.ephGetLeftmostJD.argtypes = [c_void_p]
        _lib.ephGetRightmostJD.restype = c_double
        _lib.ephGetRightmostJD.argtypes = [c_void_p]
        _lib.ephCalculateRectangular.restype = c_int
        _lib.ephCalculateRectangular.argtypes = [c_void_p, c_int, c_int,
                                                 c_double, c_double,
                                                 POINTER(c_double),
                                                 POINTER(c_double)]
        _lib.ephCalculateEulerAngles.restype = c_int
        _lib.ephCalculateEulerAngles.argtypes = [c_void_p, c_int, c_double,
                                                 c_double, POINTER(c_double),
                                                 POINTER(c_double)]
        _lib.ephCalculateTimeDiff.restype = c_int
        _lib.ephCalculateTimeDiff.argtypes = [c_void_p, c_int, c_double,
                                              c_double, POINTER(c_double)]
        _lib.ephDestroy.restype = None
        _lib.ephDestroy.argtypes = [c_void_p]

    def __init__(self, path=None):
        if EphAccess._lib is None:
            self._init_static(path)
        # Capture a reference to the _lib to access it in the __del__ method
        # because at interpreter shutdown, the module's global variables
        # are set to None
        self._lib = EphAccess._lib
        self._eph = self._lib.ephCreate()
        if self._eph is None:
            raise EphAccessException('ephCreate() error')

    def __del__(self):
        if hasattr(self, '_lib'):
            self._lib.ephDestroy(self._eph)

    def _check_result(self, result):
        if result < 0:
            raise EphAccessException(EphAccess._lib.ephLastError(self._eph))
        return result

    def object_by_name(self, name):
        code = self._lib.ephObjectByName(name.encode(_ENCODING))
        if code < 0:
            raise EphAccessException('unknown object \'%s\'' % name)
        return code

    def set_distance_units(self, units):
        return self._check_result(self._lib.ephSetDistanceUnits(self._eph,
                                                                units))

    def set_time_units(self, units):
        return self._check_result(self._lib.ephSetTimeUnits(self._eph, units))

    def load_file (self, filename):
        self._check_result(self._lib.ephLoadFile(self._eph,
                                                 filename.encode(_ENCODING)))

    def get_leftmost_jd(self):
        return self._check_result(self._lib.ephGetLeftmostJD(self._eph))

    def get_rightmost_jd(self):
        return self._check_result(self._lib.ephGetRightmostJD(self._eph))

    def calculate_rectangular (self, body, reference, jd0, jd1):
        pos = (c_double * 3)()
        vel = (c_double * 3)()
        self._check_result(self._lib.ephCalculateRectangular(
                           self._eph, body, reference, jd0, jd1, pos, vel))
        pos = [pos[i] for i in range(3)]
        vel = [vel[i] for i in range(3)]
        return [pos, vel]
    
    def calculate_euler_angles (self, frame, jd0, jd1):
        if frame is None:
            frame = 0
        pos = (c_double * 3)()
        vel = (c_double * 3)()
        self._check_result(self._lib.ephCalculateEulerAngles(
                           self._eph, frame, jd0, jd1, pos, vel))
        pos = [pos[i] for i in range(3)]
        vel = [vel[i] for i in range(3)]
        return [pos, vel]

    def calculate_time_diff(self, code, jd0, jd1):
        arr = (c_double * 1)()
        self._check_result(self._lib.ephCalculateTimeDiff(self._eph, code, jd0,
                                                          jd1, arr))
        return arr[0]

#------------------------------------------------------------------------------
   
