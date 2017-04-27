import unittest
from OSMMatcher import OSMMatcher
import time

class OSMMatcherTest(unittest.TestCase):
    def test_reading(self):
        
        trip_file = open("../moscow_trips_large/trip_343_1.txt", 'r')
        
        m = OSMMatcher("../planet-120328_moscow_large.osmdb", constraint_length=300, MAX_DIST=100)
        V=None
        p={}
        
        for row in trip_file:
            (location_lat, location_lon, _) = row.strip("\n").split(" ")
            start_time = time.time()
            (V, p) = m.step(( float(location_lat), float(location_lon) ), V, p)
            print "step time: " + str(time.time() - start_time)
            
            #max_prob_state = max(V,key=lambda x:V[x])
            #print "state: " + str(max_prob_state) + ", " + str(V[max_prob_state])
            
            #if max_prob_state in p:
            #    print "p: " + str(p[max_prob_state])
            
            # check for accuracy here somehow
            #assert(m.current_shape(V)=="101")

if __name__ == '__main__':
    unittest.main()
