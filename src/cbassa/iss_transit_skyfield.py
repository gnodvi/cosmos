#!/usr/bin/env python
#
# -----------------------------------------------------------------------------

# from author /github/cbassa
# 

import numpy as np

from skyfield.api import Topos, load
from skyfield.api import EarthSatellite

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import astropy.units as u
from astropy.wcs import wcs

if __name__ == "__main__":
    
    tle = ["ISSd",
           "1 25544U 98067A   04160.42390752  .00014992  00000-0  13290-3 0  9491",
           "2 25544  51.6329  10.4117 0005395 206.7073 225.7658 15.68815833316945"]

    satellite = EarthSatellite (tle[1], tle[2], tle[0])

    #ts = load.timescale () # not work without internet
    ts = load.timescale (builtin=True) 
    
    t    = ts.utc (2004, 6, 8, 10, 9, 17)
    tsat = ts.utc (2004, 6, 8, 10, 9, range(16, 20))
    
    location = Topos (latitude_degrees=48.2579, longitude_degrees=17.0272, elevation_m=208.0)
    
    ra, dec, dist = (satellite-location).at(tsat).radec() #(epoch='date')

    print ("")
    print ("Satellite ra/dec/dist: ", ra, dec, dist.km)

    planets = load ('de421.bsp')
    
    sun    = planets['sun']
    earth  = planets['earth']
    venus  = planets['venus']
    
    psun   = earth.at(t).observe(sun).apparent().radec()
    asun   = np.arcsin(u.solRad.to(u.km)/psun[2].km)*u.rad
    pvenus = earth.at(t).observe(venus).apparent().radec()
    avenus = np.arcsin(6051.8/pvenus[2].km)*u.rad

    # Set up wcs
    w = wcs.WCS (naxis=2)
    
    w.wcs.crval = [psun[0].to(u.deg).value, psun[1].to(u.deg).value]
    w.wcs.crpix = np.array([0.0, 0.0])
    w.wcs.cd    = np.array([[1.0, 0.0], [0.0, 1.0]])
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    # Offset position of venus
    rxv, ryv = w.wcs_world2pix (pvenus[0].to(u.deg).value, pvenus[1].to(u.deg).value, 1)

    # Offset position of ISS
    rxs, rys = w.wcs_world2pix (ra.to(u.deg).value, dec.to(u.deg).value, 1)
    
    # Create figure
    fig, ax = plt.subplots ()
    ax.set_aspect (1.0)

    # Plot Sun
    ax.add_artist (Circle((0.0, 0.0), asun.to(u.deg).value, color="y"))

    # Plot Venus
    ax.add_artist (Circle((rxv, ryv), avenus.to(u.deg).value, color="k"))

    # Plot ISS
    #ax.plot(rxs, rys, "gray")
    ax.plot (rxs, rys, "red")
    
    ax.set_xlim (0.3, -0.3)
    ax.set_ylim (-0.3, 0.3)
    ax.grid (alpha = 0.5)
    ax.set_xlabel ("R.A. offset (deg)")
    ax.set_ylabel ("Decl. offset (deg)")
    
    plt.tight_layout ()
    plt.savefig ("transit_skyfield.png", bbox_inches="tight")

    print ("")
    
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
    
