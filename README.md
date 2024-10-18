# py_observationplanner

Produce a list of observable objects on a given altitude and azimuth window for a given location.  The use case here is specifically so that people with restricted observing locations (e.g. set up on a balcony with a limited view of the sky) can see which DSO's will cross their visible window and at what time.

The sample object database is from "All Splendours, No Fuzzies." https://www.ocrasc.ca/All%20Splendor.html

Just run the Python script in the current directory: requires astropy and astroplan packages

`code`
pip install astropy
pip install astroplan
`code`

Latitude/longitude and the definition of the observing window are all hard-coded.
