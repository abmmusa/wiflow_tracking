The map matching class hierarchy has grown a bit unwieldy, so here is a quick summary.

osmmatcher.py

OSMMatcher.py is a finished map matching implementation subclassing GPSMatcher, for matching against OSM maps. It reads the map from an osmdb file, and splits any long edges into 20 meter segments. Use osmmatcher_run.py to feed OSMMatcher trip files and output matched segments. 

livebusmatcher.py

LiveBusMatcher is a finished map matching implementation subclassing GPSMatcher, customized to handle multiple overlapping but non-connected routes, and multiple agencies.

gpsmatcher.py

GPSMatcher is a superclass that knows about geometry. In particular, it maintains a spatial index for states, and uses this to produce candidate states for the viterbi algorithm. Subclasses should override geometry_of_observation() and geometry_of_state() to make use of the spatial index.

viterbi.py

All matchers build upon the Viterbi algorithm, which is implemented in viterbi.py. This is a plain viterbi, which knows only about states, and nothing about geography or geometry. This is customizable for all sorts of things, and can be sped up by supplying a candidate_states method.


