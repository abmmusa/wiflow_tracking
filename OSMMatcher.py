from streetmap import *
from gpsmatcher import GPSMatcher
from spatialfunclib import *
from mathfunclib import *
import time
import sys

EMISSION_SIGMA = 50.0 #25.0
EMISSION_UNKNOWN = 0.01

TRANSITION_UNKNOWN = 0.01
TRANSITION_UNKNOWN_UNKNOWN = 0.9
TRANSITION_SELF = .5

MAX_SPEED_M_PER_S = 20.0

class OSMMatcher(GPSMatcher):
    def __init__(self, mapdb, constraint_length=300, MAX_DIST=100):
        hmm = self.mapdb_to_hmm(mapdb)
        
        # precompute probability table
        emission_probabilities = map(lambda x: complementary_normal_distribution_cdf(x,0,EMISSION_SIGMA), 
                                     range(0,int(3.0*EMISSION_SIGMA)))
        
        def emission_probability(state, coord):
            if(state=='unknown'):
                return EMISSION_UNKNOWN
            edge = state
            projected_point = project_onto_segment(edge,coord)
            distance = haversine_distance(projected_point,coord)
            if(int(distance) >= 3 * EMISSION_SIGMA):
                return 0            
            return emission_probabilities[int(distance)]
        
        sys.stdout.write("Initing GPS matcher... ")
        sys.stdout.flush()
        GPSMatcher.__init__(self,hmm,emission_probability,constraint_length,MAX_DIST,priors={'unknown': 1.0})
        sys.stdout.write("done.\n")
        sys.stdout.flush()
    
    #    
    # add nodes between original map nodes at 20m apart        
    #
    def recursive_map_subdivide(self, themap, node):
        
        #counter is used to assign different values to multiple edges from a node while subdividing
        counter=0
        
        #subdivide edges between this node and all of its outnodes
        node_outnodes = list(node.out_nodes)
        for nextnode in node_outnodes:
            dist = haversine_distance(node.coords(), nextnode.coords()) 
            
            if( dist > MAX_SPEED_M_PER_S):
                (node_lat, node_lon) = node.coords()
                (nextnode_lat, nextnode_lon) = nextnode.coords()
                #get new node MAX_SPEED_M_PER_S apart from this node
                (newlat, newlon)=point_along_line(node_lat, node_lon,  nextnode_lat, nextnode_lon, MAX_SPEED_M_PER_S/dist)
                id = str(node.id)+"."+str(counter)
                newnode = Node(id, newlat, newlon)
                #make nextnode of this node as out_nodes of new node
                newnode.out_nodes = [nextnode]
                
                #remove this node's nextnode from out_nodes and put new node in the this node's out_nodes
                node.out_nodes.remove(nextnode)
                node.out_nodes.insert(0,newnode)
                
                #put new node in the map
                themap.nodes[id] = newnode
                
                #recursively subdivide for the newly created node
                self.recursive_map_subdivide(themap, newnode)
                counter+=1
            else:
                continue
    
    #
    # for all original map nodes call recursive_map_subdevide() to
    # to devide the map edges in 20m segements
    #
    def map_subdivide(self,themap):
        orig_map = themap
        for node in orig_map.nodes.values():
            self.recursive_map_subdivide(themap, node)
    
    def mapdb_to_hmm(self,mapdb):
        themap = StreetMap()
        themap.load_osmdb(mapdb)
        
        sys.stdout.write("Subdividing map... ")
        sys.stdout.flush()
        
        #subdivide the map in 20 m segemetns
        self.map_subdivide(themap)
        
        sys.stdout.write("into " + str(len(themap.nodes)) + " nodes.\n")
        sys.stdout.flush()
        
        sys.stdout.write("Creating HMM... ")
        sys.stdout.flush()
        
        hmm = {}
        all_edges=[]
        all_edges_with_id=[]
        for node in themap.nodes.values():
            for nextnode in node.out_nodes:
                from_edge=(node.coords(), nextnode.coords())
                all_edges.append(from_edge)

                from_edge_with_id=((node.id, nextnode.id), node.coords(), nextnode.coords())
                all_edges_with_id.append(from_edge_with_id)

                next_edges=[] # list of all next hops
                for nextnextnode in nextnode.out_nodes:
                    to_edge=(nextnode.coords(),nextnextnode.coords())
                    next_edges.append(to_edge)

                hmm[from_edge]=[(edge,(1-TRANSITION_SELF-TRANSITION_UNKNOWN)/len(next_edges)) for edge in next_edges]+\
                    [('unknown',TRANSITION_UNKNOWN),(from_edge,TRANSITION_SELF)]
        
        hmm['unknown']=[('unknown',TRANSITION_UNKNOWN_UNKNOWN)]#,(edge,(1-TRANSITION_UNKNOWN_UNKNOWN)/len(all_edges)) for edge in all_edges]
        hmm['unknown']+=[(edge,(1-TRANSITION_UNKNOWN_UNKNOWN)/len(all_edges)) for edge in all_edges]
        #hmm['unknown']=[(edge,(1-TRANSITION_UNKNOWN_UNKNOWN)/len(all_edges)) for edge in all_edges]
        
        sys.stdout.write("done.\n")
        sys.stdout.flush()
        
        return hmm
