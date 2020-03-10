=======================
The Skyfield To-Do List
=======================

This list includes both a section on immediately-reachable goals that
are appropriate for sprints and collaboration, and longer-term goals
that the code base is not quite ready for yet but that we do not want to
forget.

Promises for next version
=========================

* For #145, skip deflections of planets that can’t affect an observation.

* For #145, create a good syntax for combining two ephemerides.

* Expand a bit on the documentation for stars, now that they can have an
  epoch for their position.

Sprint Possibilities
====================

* Skyfield rounds times to the nearest millisecond when turning them
  into datetimes, to avoid times like 16:18:00 from being displayed as
  16:18:00.000512.  But we should be able to be more clever, and adjust
  the number of digits retained to match the precision of the Julian
  date 64-bit float — because, after all, dates closer to zero (long
  before the modern era) will have far more digits of precision, while
  dates in the far future might have noise even in the millisecond
  digit.  We should learn what Python did in the same situation with the
  repr() of its floats, and see if something similar could let us
  choose, for a given date, the minimum number of nonzero digits that
  would let the exact floating point value be reconstituted if the
  datetime were turned back into a Skyfield Time.

* If we are going to allow times like TT to be submitted using
  calendar dates, then we should probably provide methods that would
  fetch them back as calendar dates.

* After running a line like:

    dates = Time(utc=(1980, 1, range(1000)))

  displaying ``dates`` on the screen shows way too much output.  The
  ``repr()`` should be re-crafted so that the length of the array is
  stated, but the actual entries excerpted so that only the first few
  and the last few are shown.  Maybe the ``repr()`` should use calendar
  dates too instead of just showing the JD numbers?

* We should rename ``ecliptic_position()`` to ``ecliptic_xyz()``.

* Make JD give a sensible error if the first argument is a number or
  something.

* Load and use the various offsets between UTC and TAI that were in
  effect before 1972.

* When I wrote `add_deflection()` and needed to know whether Jupiter
  itself is available in an ephemeris, or whether the Jupiter Barycenter
  should be used in its place, I tried writing the test `if name not in
  ephemeris:`.  But this sent the code into an infinite loop!  Why does
  `in` not work on an ephemeris object?  This should be fixed, and a
  test written to keep it fixed.

* The deflection code should really use integer identifiers instead of
  using names, which are slower because they need decoding.

* The deflection code should have a quick way to reach in and ask an
  ephemeris for a raw position in au, without having to spin up a body
  object and have it spin up a `Distance` object.

* We currently download most SPICE kernels from NAIF, but have to use
  FTP for fetching DE422.  Are the files from the two sites equivalent
  and do they have the same data?  Should we prefer one or the other?

* We should have an illustration of Earth satellite heights above the
  surface, plotted against a blue atmosphere fading out into the black
  of space as the plot goes upwards towards the top.

* Iterating across an ephemeris should probably give you object names.

* Update the earth satellite IPython Notebook to the new API, or remove
  it in favor of keeping notebooks in the separate astronomy notebook
  repository.

* The ``Timescale`` object does not currently know where its data comes
  from, so its ``repr()`` is pretty uninformative.  Should it someday be
  put in charge of loading its data from the data files specified, so
  that its ``repr()`` can print out which leap second file and delta T
  file it is using?  Should it also display how up to date the files
  are, and what leap seconds it knows about?  Also: it should be told
  the expiration date of all of its data, so that it can print it out as
  part of its ``repr()``.



* In `stars.rst`, document the other alternatives for how to set the RA
  and dec of a new Star object.

* What happens if an angle that's a vector of values has .dstr() or
  .hstr() called on it?  (And, those might not be those good method
  names.)  Make it return a list of strings.

* Solar and lunar eclipses.

* We should implement comets and asteroids using the standard formulae
  (can we find a good vector version, that will match the rest of our
  approach?) for a Keplerian orbit.

* There should be routines for downloading current astronomical data.
  Each routine should take an optional filename but should also have a
  default that the corresponding "load" routine uses by default.  All
  such routines, because they are likely be called later with the same
  filename as last time, should by default just skip the download if the
  file is present (and perhaps recent enough, if a ``days_old=``
  parameter is provided?) unless an ``overwrite=False`` standard
  parameter is overridden to give it the value ``True``.  (Or does the
  possibility of ``days_old=0`` give us "overwrite" for free?)  Nice
  data would be:

  * Comet orbital elements
  * Asteroid orbital elements
  * Earth-satellite orbital elements
  * (What others can we think of?)

* If anyone is a NumPy expert, I would love comments on whether my code
  is at all idiomatic, or whether I'm doing things inefficiently and in
  such a quirky way that no one else will ever understand it.

* If one of the parameters of nutationlib or precessionlib was
  incorrect, would our tests detect it?  It would be nice to have an
  automated method of determining which particular constants in the code
  base our tests are not sensitive to.

* Should a simple "print position.radec()" print something prettier?

* We should maybe support light seconds or minutes (or both?) as ways to
  express and to set a distance.  As we keep adding them, should we
  really keep adding more named parameters to the constructor?  Maybe we
  should instead start using class method constructors instead?

Adding more smarts to ephemeris handling
========================================

* When we add or subtract vectors, the new `VectorSum` needs to inherit
  an `ephemeris` from one of the segments being combined.  Right now it
  just grabs the first one.  What if, of the choices of ephemeris among
  the segments, grabbing the first one gives us an ephemeris that is
  missing several key large bodies, and so we run into an exception when
  `.apparent()` tries to compute gravitational deflection?  Should we be
  more intelligent in our choice?  Or should we even combine the various
  ephemerides our segments might offer?  Or should we specifically go
  ahead and look for the deflectors we need and try to find a segment
  for them each?

* And additionally: the error when not enough bodies are available for
  deflection maybe someday needs to be more helpful.

Reading List for Brandon
========================

* https://naif.jpl.nasa.gov/pub/naif/FIDO/misc/njb/src/geom.c

For 2.0
=======

* Remove old deprecation warnings for pre-1.0 behaviors.

* Remove support and tests for old ephemeris Python packages.

Longer-term goals
=================

* Make all objects that are `.observe()`’d from Earth include a
  sublatitude and sublongitude coordinate stating the position on Earth
  from which they appear directly overhead.  When complete, make a note
  at the PyEphem GitHub issue:

  https://github.com/brandon-rhodes/pyephem/issues/16

* When Earth Satellites are implemented, include the orbit number of a
  satellite's current position in the public attributes that are set on
  the resulting position object, as promised in PyEphem GitHub issue:

  https://github.com/brandon-rhodes/pyephem/issues/15


.. documentation
     writing up SkyField solutions to PyEphem questions on Stack Overflow
     section on accuracy of each algorithm involved
     logo?
   performance
     Is all this vectorization worth it?
       Run a loop to compute N planet positions.
       Do the same computation using a vector of N jd's.
       Compare the runtimes under both C Python and PyPy.
       Might have to do numpypy thing; do it in skyfield/__init__.py.
       If they both show a difference, then YES it is worth it.
       Could Star() become a whole catalog of stars processed in parallel?
     What routines are taking the most time when the tests are run?
     Try to take advantage of jplephem's ability to use bundles
