import unittest
from location import TripLoader
from busmatcher import BusMatcher

class BusMatcherTest(unittest.TestCase):
    def test_reading(self):

        all_trips = TripLoader.get_all_trips("test/trips/0/")        
        for trip in all_trips[:10]:
            m = BusMatcher('test/shapes_up.txt')
            V=None
            p={}
            for location in trip.locations:                
                V,p = m.step((location.latitude,location.longitude),V,p)
            assert(m.current_shape(V)=="101")

if __name__ == '__main__':
    unittest.main()
