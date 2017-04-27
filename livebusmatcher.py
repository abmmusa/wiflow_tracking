from gpsmatcher import GPSMatcher
from viterbi import Viterbi
import itertools
from mathfunclib import complementary_normal_distribution_cdf
from spatialfunclib import *
import zmq
import json

#transition probabilities
FROM_SELF_TO_NEXT = 0.98 #0.98
FROM_SELF_TO_SELF = 0.98
FROM_SELF_TO_UNKNOWN = 0.0000025
#FROM_UNKNOWN_TO_STATES = 0.0000025  #0.00000025
FROM_UNKNOWN_TO_UNKNOWN = 0.98 #0.7

EMISSION_MAXDIST = 150
#EMISSION_SIGMA = 75 # meters of GPS error standard deviation
#EMISSION_MEAN = 0  #mean for gaussian distribution used to calculate emission probability
EMISSION_UNKNOWN = 0.01
# jakob: used to be the below...
#EMISSION_UNKNOWN = 1 - ((stdNormCDFLib.standardNormalCDF[emission_unknown_index] - 0.5)*2)

import couchdb

class LiveBusMatcher(GPSMatcher):

    def __init__(self, EMISSION_SIGMA=15, STEP=0.1):        
        couch = couchdb.Server("http://tables.cs.uic.edu:5984/")
        livebus_db = couch['livebus']
        
        hmm = self.couch_to_hmm(livebus_db)

        # these hold separate V and p dictionaries for each unique mac address
        # jakob - it's going to add up: should probably garbage collect these after some time has passed
        self.Vs={}
        self.Ps={}

        # precompute probability table
        #emission_probabilities = map(lambda x: complementary_normal_distribution_cdf(x,0,EMISSION_SIGMA), 
        #                             range(0,int(3.0*EMISSION_SIGMA)))
        
        def emission_probability(state, emission):            
            (_,rest)=state
            (coord_agency,coord)=emission

            if(rest=='unknown'):
                (state_agency,_)=state
                if state_agency != coord_agency: return 0
                else:
                    return EMISSION_UNKNOWN

            ((state_agency,state_shape),edge)=state

            # one agency never ever drives another's routes
            if state_agency != coord_agency: 
                return 0

            projected_point = project_onto_segment(edge,coord)
            dist = distance(projected_point,coord)
            if(dist < EMISSION_MAXDIST):
                return 1
            else:
                return 0
#            elif(int(dist) >= 3 * EMISSION_SIGMA):
#                return 0            
#            return emission_probabilities[int(dist)]

        priors = {}
        agencies = livebus_db.view("_design/main/_view/agency_route_count",None,group=True)
        rowcount = len(agencies.rows)
        for row in agencies.rows:
            priors[(row.key,"unknown")]=1.0/rowcount

        GPSMatcher.__init__(self,hmm,emission_probability,constraint_length=120,priors=priors)
        
    def match_message(self, msg):
        mac = msg["mac"]
        if not mac in self.Vs:
            self.Vs[mac]=None
            self.Ps[mac]={}

        V = self.Vs[mac]
        p = self.Ps[mac]

        obs = (msg["agency"],(msg["lat"],msg["lon"]))
        self.Vs[mac],self.Ps[mac] = self.step(obs,V,p)

    def current_shape(self,mac):
        V=self.Vs[mac]
        sortedV = sorted(V,key=lambda x:V[x],reverse=True)

        sumV = {}
        current_shape = None
        for (agency_shape,edge) in sortedV:
            if edge == "unknown": 
                shape=edge
            else:
                (agency,shape)=agency_shape

            if current_shape!=shape and not (shape in sumV):
                sumV[shape]=V[(agency_shape,edge)]
#                print "adding",mac,sumV[shape],shape
            elif shape==current_shape:
                sumV[shape]+=V[(agency_shape,edge)]
#                print "incrementing",mac,sumV[shape],shape
            current_shape = shape
            
        sortedV = sorted(sumV,key=lambda x:sumV[x],reverse=True)
        shape = sortedV[0]

#        print sumV[sortedV[0]]

        if(len(sortedV)==1 or
           sumV[sortedV[0]]>0.7): 
            return ("dummy",shape)
        else:
            return (agency,"undecided")


    def geometry_of_state(self, state):
        """ Overridden to support spatial indexing in GPSMatcher"""
        (agency,rest)=state
        if rest == 'unknown': return None
        else:
            (route, edge) = state
            return edge

    def geometry_of_observation(self, obs):
        (agency,rest)=obs
        return rest

    def couch_to_hmm(self, livebus_db):
        """ turns a file on the form
        shape,lat,lon,seq,len
        ...
        to an HMM
        """
        shapes = {}
        for row in livebus_db.view('_design/main/_view/approved_routes',None,include_docs=True):
            coords = row.doc["coordinates"]
            shape_id = row.doc["_id"]
            agency = row.doc["agency"]
            shape = (agency,shape_id)
            print shape,row.doc["name"]
            for lat,lon in coords:
                if not shape in shapes: shapes[shape]=[]
                shapes[shape].append((float(lat),float(lon)))

        outgoing = []
        hmm = {}
        # add transitions between all consecutive edges
        for shape in shapes:
            (agency,_)=shape

            # first build a list of states (edges)
            from_c = shapes[shape]            
            to_c = from_c[1:]
            to_c.append(from_c[0])
            from_edges = list(itertools.izip(from_c,to_c))

            # then build our hmm, a sequential mapping of edge->edge
            to_edges = from_edges[1:]
            to_edges.append(from_edges[0])
            edge_pairs = itertools.izip(from_edges,to_edges)
                    
            # build a dict which has all the location data in it, and add this to the hmm
            hmm.update(dict([((shape, from_edge), 
                    [((shape,to_edge),FROM_SELF_TO_NEXT),
                     ((shape,from_edge),FROM_SELF_TO_SELF),
                     ((agency,'unknown'),FROM_SELF_TO_UNKNOWN)]) 
                   for from_edge,to_edge in edge_pairs]))
            
            # add edges from unknown to each state, and back to unknown
            if not (agency,"unknown") in hmm:
                hmm[(agency,"unknown")]=[((agency,'unknown'),FROM_UNKNOWN_TO_UNKNOWN)]
            # jakob: this is a little buggy, should depend on the number of out edges
            hmm[(agency,"unknown")]+=[((shape,edge),(1-FROM_UNKNOWN_TO_UNKNOWN)) for edge in from_edges]

        return hmm

if __name__ == '__main__':
    matcher = LiveBusMatcher()

    context = zmq.Context()
    incoming = context.socket(zmq.SUB)
    incoming.connect("tcp://131.193.34.157:5001");
    incoming.setsockopt(zmq.SUBSCRIBE,"")

    outgoing = context.socket(zmq.PUB)
    outgoing.connect("tcp://131.193.34.157:5000");
    
    while True:
	msg = json.loads(incoming.recv())
        if "lat" in msg: 
            matcher.match_message(msg)
            (match_agency,match_shape)=matcher.current_shape(msg["mac"])            
            outgoing.send(json.dumps({"type":"classification","mac":msg["mac"],"agency":msg["agency"],"route":match_shape}))
