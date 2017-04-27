from streetmap import *
from gpsmatcher import GPSMatcher
from spatialfunclib import *
from mathfunclib import *
from viterbi import Viterbi
import time

EMISSION_SIGMA = 100.0
DISTANCE_SIGMA = 300.0

EMISSION_UNKNOWN = 0.00001

TRANSITION_UNKNOWN = 0.00001
TRANSITION_UNKNOWN_UNKNOWN = 0.9
TRANSITION_SELF = .5

MAX_SPEED_M_PER_S = 20.0

TRANSMISSION_PROB=1.0/300.0


# read kde values for different ssi binning from files
kde_file_40 = open('emission_prob_distribution/kde-output-40+').readlines()
kde_file_40_50 = open('emission_prob_distribution/kde-output-40-50').readlines()
kde_file_50_60 = open('emission_prob_distribution/kde-output-50-60').readlines()
kde_file_60_70 = open('emission_prob_distribution/kde-output-60-70').readlines()
kde_file_70 = open('emission_prob_distribution/kde-output-70-').readlines()
kde_file_all = open('emission_prob_distribution/kde-output-all').readlines()

# make list from above data. each list will contain kde values for 1m to 300m distances
kde_data_40 = map(lambda x:float(x),  kde_file_40)
kde_data_40_50 = map(lambda x:float(x),  kde_file_40_50)
kde_data_50_60 = map(lambda x:float(x),  kde_file_50_60)
kde_data_60_70 = map(lambda x:float(x),  kde_file_60_70)
kde_data_70 = map(lambda x:float(x),  kde_file_70)
kde_data_all= map(lambda x:float(x),  kde_file_all) # for no obs calculation


class WiFiMatcher:
    def __init__(self, mapdb):
        self.previous_obs = None

        hmm = self.mapdb_to_hmm(mapdb)
        
        #emission_probabilities = map(lambda x: complementary_normal_distribution_cdf(x,0,EMISSION_SIGMA),range(0,int(3.0*EMISSION_SIGMA)))
            
        priors=dict([(state,1.0/len(hmm)) for state in hmm])
        
        self.viterbi = Viterbi(hmm,self.emission_probability,
                          constraint_length=2500, # BE CAREFUL with it. walking may take long time and higher value may be needed here
                          priors=priors)


    def emission_probability(self, edge, obs):            
        if(edge=='unknown'):
            return EMISSION_UNKNOWN

        (edge_end1,edge_end2)=edge
            
        emission_prob=1.0;
        for obs_single_mon in obs.values():
            (lat,lon,ssi)=map(lambda x:float(x), obs_single_mon)
                
            dist1=fast_distance(edge_end1, (lat,lon))
            dist2=fast_distance(edge_end2, (lat,lon))
                
            prob_for_mon = self.prob_mon_given_segment(dist1, dist2, ssi)
            #if prob_for_mon != 0 and ssi !=0:
            #    print ssi,prob_for_mon
            emission_prob *= prob_for_mon

        #multiply by transmission prob    
        emission_prob *= TRANSMISSION_PROB    
            
        #if emission_prob!=0.0:
        #    print emission_prob
        return emission_prob




    #TODO: now taking mean of kde value for two endpoints only. need to take sum of all in between. 
    #not doing now because of complexity that arises for edges parallel to monitor like following scenario
    # edge: --------
    #
    #
    # edge: --------
    # mon:     .
    def prob_mon_given_segment(self,dist1, dist2, ssi):
            #kde values are for 1m to 300m and the indices are from 0 to 299
        dist1_int=int(dist1-1)
        dist2_int=int(dist2-1)
            
        if(dist1_int<0 or dist1_int>299 or dist2_int<0 or dist2_int>299):
            if ssi==0:
                return 1
            else:
                return 0
            
        kde_table_to_use=[]
            
        if ssi==0: #no observation
            kde_table_to_use = list(kde_data_all)
            prob_detection=((kde_table_to_use[dist1_int]+kde_table_to_use[dist2_int])/2)
            return 1-prob_detection

        else:
            if ssi>=-40:
                kde_table_to_use = list(kde_data_40)
            elif ssi<-40 and ssi>=-50:
                kde_table_to_use = list(kde_data_40_50)
            elif ssi<-50 and ssi>=-60:
                kde_table_to_use = list(kde_data_50_60)
            elif ssi<-60 and ssi>=-70:
                kde_table_to_use = list(kde_data_60_70)
            else:
                kde_table_to_use = list(kde_data_70)
            
            prob_detection=((kde_table_to_use[dist1_int]+kde_table_to_use[dist2_int])/2)

            return prob_detection
        




    def step(self,obs,V,p):    
        #if V is not None:
        #    print "step",len(V.keys()),time.time()
        if self.previous_obs != None:
            for int_obs in self.interpolated_obs(self.previous_obs, obs):
                V,p = self.viterbi.step(int_obs,V,p)        
                
        V,p = self.viterbi.step(obs,V,p)
        self.previous_obs = obs
        return V,p

    def interpolated_obs(self,prev,obs):
        return []
    



    #    
    # add nodes between original map nodes at 20m apart        
    #
    def recursive_map_subdivide(self, themap, node):
        #counter is used to assign different values to multiple edges from a node while subdividing
        counter=0

        #subdivide edges between this node and all of its outnodes
        node_outnodes = list(node.out_nodes)
        for nextnode in node_outnodes:
            dist = fast_distance(node.coords(), nextnode.coords()) 
        
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
    # for all original map nodes call recursive_map_subdevide()
    # to devide the map edges in 20m segements
    #
    def map_subdivide(self,themap):
        orig_map = themap
        for node in orig_map.nodes.values():
            self.recursive_map_subdivide(themap, node)



    def mapdb_to_hmm(self,mapdb):
        themap = StreetMap()
        themap.load_osmdb(mapdb)

        #subdivide the map in 20 m segemetns
        self.map_subdivide(themap)

        #self.print_map_nodes(themap)
        self.print_nodes_outnodes_id(themap)
        self.write_map_nodes(themap)    
            

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
                    [(from_edge,TRANSITION_SELF)]
#                    [('unknown',TRANSITION_UNKNOWN),(from_edge,TRANSITION_SELF)]


        hmm['unknown']=[('unknown',TRANSITION_UNKNOWN_UNKNOWN)]
        hmm['unknown']+=[(edge,(1-TRANSITION_UNKNOWN_UNKNOWN)/len(all_edges)) for edge in all_edges]
        #hmm['unknown']=[(edge,(1-TRANSITION_UNKNOWN_UNKNOWN)/len(all_edges)) for edge in all_edges]

        self.write_hmm_edges_endpoints_with_id(all_edges_with_id)
        
        return hmm

    #
    # print map nodes. debug function
    #
    def print_map_nodes(self,themap):
        for node in themap.nodes.values():
            printstr="node="+str(node.id)+" out_nodes="
            for outnode in node.out_nodes:
                printstr += str(outnode.id)+" "
            print printstr
        

    #
    # print lat,lon of map nodes and all of its outnodes
    #
    def print_nodes_outnodes_latlon(self,themap):
        file=open('test/wifi/map_vis/nodes_outnodes_latlon.txt', 'w')
        for node in themap.nodes.values():
            (lat,lon)=node.coords()
            printstr=str("%.6f"%lat)+","+str("%.6f"%lon)
            for outnode in node.out_nodes:
                (lat,lon)=node.coords()
                printstr +=" "+str("%.6f"%lat)+","+str("%.6f"%lon)
                
            #print printstr
            file.write(printstr)
            file.write("\n")
            
        file.close()



    def print_nodes_outnodes_id(self,themap):
        file=open('test/wifi/map_vis/map_nodes_outnodes_id.txt', 'w')
        for node in themap.nodes.values():
            
            printstr=str(node.id)
            for outnode in node.out_nodes:
                printstr +=" "+str(outnode.id)
                
            #print printstr
            file.write(printstr)
            file.write("\n")
            
        file.close()



    #
    # write map nodes to a file, format: <nodeid lat lon>
    #
    def write_map_nodes(self,themap):
        f_nodes=open('test/wifi/map_vis/map_nodes.txt', 'w')
        for node in themap.nodes.values():
            (lat,lon)=node.coords()
            printstr=str(node.id)+" "+str("%.6f"%lat)+" "+str("%.6f"%lon)+"\n"
            f_nodes.write(printstr)

    #
    # write midpoint of all edges in the hmm for visualizing using google earth (without nodeid)
    #    
    def write_hmm_edges_midpoint(self,hmm_edges):
        f_edges=open('test/wifi/map_vis/map_edges_midpoint.txt', 'w')
        
        for edge in hmm_edges:
            ((lat1,lon1),(lat2,lon2))=edge
            write_str = str( (lat1+lat2)/2 ) + " " + str( (lon1+lon2)/2 )+ "\n"
            f_edges.write(write_str)

        f_edges.close()

    #
    # write midpoint of all edges in the hmm for visualizing using google earth (with edgeid)
    # edgeid is node1id-node2id
    #    
    def write_hmm_edges_midpoint_with_id(self,hmm_edges):
        f_edges=open('test/wifi/map_vis/map_edges_midpoint.txt', 'w')
        
        for edge in hmm_edges:
            ((idnode1,idnode2),(lat1,lon1),(lat2,lon2))=edge
            write_str = str(idnode1)+"-"+str(idnode2)+" "+ str( (lat1+lat2)/2 ) + " " + str( (lon1+lon2)/2 )+ "\n"
            f_edges.write(write_str)

        f_edges.close()

    #
    # write endpoints of edge with id instead of midpoint (with nodeid)
    # edgeid is node1id-node2id
    #    
    def write_hmm_edges_endpoints_with_id(self,hmm_edges):
        f_edges=open('test/wifi/map_vis/map_edges_endpoints.txt', 'w')
        
        for edge in hmm_edges:
            ((idnode1,idnode2),(lat1,lon1),(lat2,lon2))=edge
            write_str = str(idnode1)+"-"+str(idnode2)+" "+str(lat1)+" "+str(lon1)+" "+str(lat2)+" "+str(lon2)+"\n"
            f_edges.write(write_str)

        f_edges.close()
        

#WiFiMatcher("../mapgenerate/ground_truth.osmdb")
