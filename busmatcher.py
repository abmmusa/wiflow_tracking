from gpsmatcher import GPSMatcher
from viterbi import Viterbi
import itertools
from mathfunclib import complementary_normal_distribution_cdf
from spatialfunclib import *

#transition probabilities
FROM_SELF_TO_NEXT = 0.95 #0.98
FROM_SELF_TO_SELF = 0.98
FROM_SELF_TO_UNKNOWN = 0.0000025
FROM_UNKNOWN_TO_STATES = 0.0000025  #0.00000025
FROM_UNKNOWN_TO_UNKNOWN = 0.98 #0.7

EMISSION_SIGMA = 75 # meters of GPS error standard deviation
EMISSION_MEAN = 0  #mean for gaussian distribution used to calculate emission probability
EMISSION_UNKNOWN = 0.0001
# jakob: used to be the below...
#EMISSION_UNKNOWN = 1 - ((stdNormCDFLib.standardNormalCDF[emission_unknown_index] - 0.5)*2)

class BusMatcher(GPSMatcher):

    def __init__(self, shapefile, EMISSION_SIGMA=75, STEP=0.1):        
        hmm = self.shapefile_to_hmm(shapefile)

        # precompute probability table
        emission_probabilities = map(lambda x: complementary_normal_distribution_cdf(x,0,EMISSION_SIGMA), 
                                     range(0,int(3.0*EMISSION_SIGMA)))
        
        def emission_probability(state, coord):            
            if(state=='unknown'):
                return EMISSION_UNKNOWN
            (shape,edge) = state
            projected_point = project_onto_segment(edge,coord)
            distance = fast_distance(projected_point,coord)
            if(int(distance) >= 3 * EMISSION_SIGMA):
                return 0            
            return emission_probabilities[int(distance)]

        GPSMatcher.__init__(self,hmm,emission_probability,constraint_length=1200)
        

    def current_shape(self,V):
        (shape, edge) = max(V,key=lambda x:V[x])
        return shape

    def geometry_of_state(self, state):
        """ Overridden to support spatial indexing in GPSMatcher"""
        if state == 'unknown': return None
        else:
            (route, edge) = state
            return edge

    def shapefile_to_hmm(self, busshapefile):
        """ turns a file on the form
        shape,lat,lon,seq,len
        ...
        to an HMM
        """
        shapefile = open(busshapefile,'r')
        shapefile.readline() # skip the first line
        shapes = {}
        coords = [s.strip().split(',') for s in shapefile]

        for shape,lat,lon,_,_ in coords:
            if not shape in shapes: shapes[shape]=[]
            shapes[shape].append((float(lat),float(lon)))

        outgoing = []
        # add transitions between all consecutive edges
        for shape in shapes:
            # first build a list of states (edges)
            from_c = shapes[shape]            
            to_c = from_c[1:]
            to_c.append(from_c[0])
            from_edges = list(itertools.izip(from_c,to_c))

            # then build our hmm, a sequential mapping of edge->edge
            to_edges = from_edges[1:]
            to_edges.append(from_edges[0])
            edge_pairs = itertools.izip(from_edges,to_edges)
                    
            # build a dict which has all the location data in it, seems 
            outgoing+=[((shape, from_edge), 
                    [((shape,to_edge),FROM_SELF_TO_NEXT),
                     ((shape,from_edge),FROM_SELF_TO_SELF),
                     ('unknown',FROM_SELF_TO_UNKNOWN)]) 
                   for from_edge,to_edge in edge_pairs]
            
            # add edges from unknown to each state
            outgoing+=[('unknown',
                       [((shape,edge),0.1/len(from_c)) for edge in from_edges])]

        outgoing+=[('unknown',[('unknown',0.9)])]
        return dict(outgoing)

## everything below here is stale leftovers

    def directionProb(self,ewmaAngle):
        """ calculates additional emission probability in a Hidden Markov Model based on
	direction of travel and a long term exponentially weighted moving average
	of direction of travel (EWMA) """

        if(ewmaAngle == -1000):
            return 1
        segX = self.lat2-self.lat1
        segY = self.lon2-self.lon1
        segLen = vectorLen(segX,segY)
        if(segLen == 0):
            return 1
        segBearing = spatialfunclib.path_bearing(self.lat1,self.lon1,self.lat2,self.lon2)
        angle = math.fabs(ewmaAngle-segBearing)
        if(angle > 180):
            angle = 360 - angle
            standardDirProb = (angle - DIRECTION_MEAN) / float(DIRECTION_SIGMA)
            index = int(standardDirProb * (1/stdNormCDFLib.STDNORMCDF_STEP_SIZE)) + int((3 / stdNormCDFLib.STDNORMCDF_STEP_SIZE))
            if (index >= len(stdNormCDFLib.standardNormalCDF)):
                return 0.0		
            directionProbability = 1 - ((stdNormCDFLib.standardNormalCDF[index] - 0.5)*2)
            return directionProbability


    def old_emissionProb(self,state, coord, emission_probabilities):
            """Takes into account distance from a shapeSegment, direction alignment to the segment, and time
            of day to account for bus schedules.
            
            input: observation (timePoint), id of unknown state,
            #and ewmaAngle (float)--long term direction average
            #output: probability (float)"""

            if(state=='unknown'):
                return EMISSION_UNKNOWN
#            if(self.shapeId != checkTimeSchedule(self.shapeId,observation.timestamp)):
#                return 0
            (shape,segment) = state
            distance = distance_to_segment(segment,coord)
            if(int(distance) >= 3 * EMISSION_SIGMA):
                return 0            
            return emission_probabilities[distance/STEP]
#normal_distribution_cdf(distance,0,EMISSION_SIGMA)

#            standardDist = (distance - EMISSION_MEAN) / float(EMISSION_SIGMA)
#            index = int(standardDist * (1/stdNormCDFLib.STDNORMCDF_STEP_SIZE)) + int((3 / stdNormCDFLib.STDNORMCDF_STEP_SIZE))
#
#            distanceProbability = 1 - ((stdNormCDFLib.standardNormalCDF[index] - 0.5)*2)
#            directionProbability = self.directionProb(ewmaAngle)
#            return distanceProbability * directionProbability



