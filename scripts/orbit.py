#!/usr/bin/env python3

from glide.common_components.view_geometry import carruthers_orbit
from glide.science.model import default_vol
from glide.science.plotting import orbit_svg
from glide.science.orbit import viewgeom2ts


vol = default_vol()
vg = carruthers_orbit(duration=30*6, num_obs=40)
orbit_svg(vol, viewgeom2ts(vg)).save('orbit.gif')