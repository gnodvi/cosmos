#! /usr/bin/env python
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

from ephaccess import *

eph = EphAccess (".")


eph.load_file ("epm2017.bsp")
eph.load_file ("moonlibr_epm2017.bpc")

eph.set_distance_units(EPH_KM)
eph.set_time_units(EPH_SEC)


print(eph.get_leftmost_jd())
print(eph.get_rightmost_jd())

print (eph.calculate_rectangular(eph.object_by_name("moon"),
                                 eph.object_by_name("sun"), 2451545.0, 0.0))

print (eph.calculate_euler_angles (None, 2451545.0, 0.0))

print (eph.calculate_time_diff (EPH_TT_MINUS_TDB, 2451545.0, 0.0))

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
