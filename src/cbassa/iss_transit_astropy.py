#!/usr/bin/env python3
# 
# ---------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Circle
import astropy.units as u

from astropy.time import Time
from astropy.coordinates import get_body, solar_system_ephemeris

from astropy.coordinates import EarthLocation, GCRS, FK5
from astropy.coordinates import SkyCoord, PrecessedGeocentric, CartesianRepresentation
from astropy.wcs import wcs
from astropy.coordinates.earth_orientation import precession_matrix_Capitaine, nutation_components2000B
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product, matrix_transpose
from astropy.coordinates import cartesian_to_spherical

#quit ()

# ---------------------------------------------------------------------------------------
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
# ---------------------------------------------------------------------------------------

#quit ()

# From astropy, but with units
def nutation_matrix(epoch):
    """
    Nutation matrix generated from nutation components.

    Matrix converts from mean coordinate to true coordinate as
    r_true = M * r_mean
    """
    # TODO: implement higher precision 2006/2000A model if requested/needed
    epsa, dpsi, deps = nutation_components2000B(epoch.jd)  # all in radians

    epsa *= u.rad
    dpsi *= u.rad
    deps *= u.rad
    
    return matrix_product(rotation_matrix(-(epsa + deps), 'x'),
                          rotation_matrix(-dpsi, 'z'),
                          rotation_matrix(epsa, 'x'))

# ---------------------------------------------------------------------------------------
if __name__ == "__main__":
    
    # Set time
    t = Time("2004-06-08T10:09:17", scale="utc", format="isot")
    dts = np.arange(5)*u.s
    
    # Set location
    loc = EarthLocation (lat=48.2579*u.deg, lon=17.0272*u.deg, height=202.0*u.m)

    # Set tle
    tle = ["ISSd",
           "1 25544U 98067A   04160.42390752  .00014992  00000-0  13290-3 0  9491",
           "2 25544  51.6329  10.4117 0005395 206.7073 225.7658 15.68815833316945"]


    #quit ()

    # ---------------------------------------------------------------------------------------
    # Set satellite
    satellite = twoline2rv(tle[1], tle[2], wgs84)
    
    p = np.zeros(len(dts)*3).reshape(3, len(dts)).astype("float32")
    v = np.zeros(len(dts)*3).reshape(3, len(dts)).astype("float32")
    
    for i, dt in enumerate(dts):
        tm = t+dt
        p[:, i], v[:, i] = satellite.propagate(tm.datetime.year, tm.datetime.month, tm.datetime.day, tm.datetime.hour, tm.datetime.minute, tm.datetime.second)

    # ---------------------------------------------------------------------------------------

    #quit ()
    
    # Nutation and precession
    epsa, dpsi, deps = nutation_components2000B(t.jd)
    rp = precession_matrix_Capitaine(t, Time("J2000"))
    rn = nutation_matrix(t)
    re = rotation_matrix(-dpsi*np.cos(epsa), "z")
    r = matrix_product(rp, rn, re)

    
    pj2000 = np.dot(r, np.array(p))
    vj2000 = np.dot(r, np.array(v))
    satpos = CartesianRepresentation(x=pj2000[0], y=pj2000[1], z=pj2000[2], unit=u.km)
    obspos, obsvel = loc.get_gcrs_posvel(obstime=t)

    #quit ()

    dr = satpos-obspos
    dist, dec, ra = cartesian_to_spherical(dr.x, dr.y, dr.z)
    
    pos = SkyCoord (ra=ra, dec=dec, frame="gcrs")
    print (pos)
    
    #quit ()

    # Load solary system ephemeris
    #
    #solar_system_ephemeris.set ("de432s")  # [SSL: CERTIFICATE_VERIFY_FAILED]
    solar_system_ephemeris.set ("de421.bsp")

    #quit ()

    # Get position of Venus
    pvenus = get_body ("venus", t, loc)
    avenus = np.arcsin(6051.8*u.km/pvenus.distance.to(u.km))
    
    # Get position of sun
    psun = get_body ("sun", t, loc)
    asun = np.arcsin(u.solRad/psun.distance.to(u.km))

    #quit ()  ERROR
    
    
    # Set up wcs
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [psun.ra.degree, psun.dec.degree]
    w.wcs.crpix = np.array([0.0, 0.0])
    w.wcs.cd = np.array([[1.0, 0.0], [0.0, 1.0]])
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    # Offset position of venus
    rxv, ryv = w.wcs_world2pix(pvenus.ra.degree, pvenus.dec.degree, 1)

    # Offset position of ISS
    rxs, rys = w.wcs_world2pix(pos.ra.degree, pos.dec.degree, 1)
    
    # Create figure
    fig, ax = plt.subplots()
    ax.set_aspect(1.0)

    # Plot Sun
    ax.add_artist(Circle((0.0, 0.0), asun.to(u.deg).value, color="y"))

    # Plot Venus
    ax.add_artist(Circle((rxv, ryv), avenus.to(u.deg).value, color="k"))

    # Plot ISS
    ax.plot(rxs, rys, "gray")
    
    ax.set_xlim(0.3, -0.3)
    ax.set_ylim(-0.3, 0.3)
    ax.grid(alpha=0.5)
    ax.set_xlabel("R.A. offset (deg)")
    ax.set_ylabel("Decl. offset (deg)")
    
    plt.tight_layout()
    plt.savefig("transit_astropy.png", bbox_inches="tight")

# ---------------------------------------------------------------------------------------

#urllib.error.URLError: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate (_ssl.c:1108)>

# ---------------------------------------------------------------------------------------
