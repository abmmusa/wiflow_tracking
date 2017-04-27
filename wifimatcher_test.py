import unittest
from location import TripLoader
from WiFiMatcher import WiFiMatcher
from spatialfunclib import fast_distance
import sys
import os


MIN_ALLOWED_SPEED_M_S = 1 #3.6 km/hr 
TOTAL_MONITORS = 7

MAX_TIMEGAP_BEFORE_NEW_VITERBI_COMPUTAION = 600 #seconds

#read_path = "test/wifi/input_data"
#write_path = "test/wifi/output_data"

#read_path = "test/wifi/input_data_all_ourmacs"
#write_path = "test/wifi/output_data_all_ourmacs"

#read_path = "test/wifi/input_data_allmacs"
#write_path = "test/wifi/output_data_allmacs"


#read_path = "test/wifi/input_data_mon_removed_one"
#write_path = "test/wifi/output_data_mon_removed_one"

#read_path = "test/wifi/input_data_mon_removed_two"
#write_path = "test/wifi/output_data_mon_removed_two"

read_path = "test/wifi/input_data_mon_removed_three"
write_path = "test/wifi/output_data_mon_removed_three"




log_file=open('test/wifi/log_viterbi.txt', 'w');

m = WiFiMatcher("../mapinference/trr12/roosevelt_small.osmdb")  


class WiFiMatcherTest(unittest.TestCase):
    def write_viterbi_path(self, V, p, fname, starttime):
        #print V.values()
        #print V[max(V,key=lambda x:V[x])]
        
        try:
            edges = p[max(V,key=lambda x:V[x])]
            #print edges
        except: #TODO: investigate later (empty V, so no edges)
            log_file.write("no edges for "+fname)
            print "no edges for "+fname
            return

        f = open(write_path+"/"+fname,'w')
        counter = 0
        for edge in edges:
            if edge=='unknown': #TODO: investigate later (somtime there is no edge, example 00-11-20-f8-ab-6b)
                print "closing file due to unknown edge..."
                f.close()
                return
            
            ((lat1,lon1),(lat2,lon2))=edge
            time=starttime+counter
            write_str = str(time)+","+str( (lat1+lat2)/2 ) + "," + str( (lon1+lon2)/2 )+ "\n"
            #print write_str
            f.write(write_str)
            counter += 1

        f.close()        

    #
    # produce interpolated viterbi path. Interpolates for no observation
    #    
    def write_interpolated_viterbi_path(self,V,p,fname,inittime,input_data):
        input_detection_time_lastnext={}

        # fill up next_detection dictionary from input data
        for i in range(1, len(input_data)):
            timelast=int(input_data[i-1].rstrip().split(' ')[0])
            timenext=int(input_data[i].rstrip().split(' ')[0])

            input_detection_time_lastnext[i-1] = (timelast,timenext)


        #create a noobs from input data (any entry is fine as all entries contains all monitors)
        noobs={}
        data = input_data[0].rstrip().split(' ')
        for k in range(0,TOTAL_MONITORS):
            noobs[k]=[ data[3*k+1], data[3*k+2], 0 ]
            noobs[k]=[ data[3*k+1], data[3*k+2], 0 ]



        try:
            edges = p[max(V,key=lambda x:V[x])]
        except: #TODO: investigate later (empty V, so no edges)
            log_file.write("no edges for "+fname)
            print "no edges for "+fname
            return


        # fill up the time_edge dictionary with edge output from viterbi at each epoch
        starttime = inittime        
        time_edge={}
        for edge in edges:
            if edge=='unknown': 
                print "contains unknown edge... can't procede"
                return

            time_edge[starttime]=edge
            starttime += 1

        
        f = open(write_path+"/"+fname,'w')

        for i in range(0,len(input_detection_time_lastnext)):
            (detection_last, detection_next) = input_detection_time_lastnext[i]
            timegap = detection_next - detection_last

            # write viterbi output to file for detection observations
            time = detection_last
            #print time
            ((lat1,lon1),(lat2,lon2))= time_edge[time]
            write_str = str(time)+","+str( (lat1+lat2)/2 ) + "," + str( (lon1+lon2)/2 )+ "\n"
            f.write(write_str)
            

            if timegap > 1: #edges onwards are for no observation

                # fill up all edges during no observation in a list
                all_edges_during_noobs = []
                for i in range(1,timegap):
                    current_time = detection_last + i
                    edge_at_currenttime = time_edge[current_time]
                    
                    if edge_at_currenttime not in all_edges_during_noobs:
                        all_edges_during_noobs.append(edge_at_currenttime)


                
                # find the total probability for all edges during no observation
                prob_sum_all_states = 0.0
                for edge in all_edges_during_noobs:
                    prob_sum_all_states += m.emission_probability(edge,noobs)

                # now find time at each state for
                time_after_rounding = 0.0   # time extra (+ or -) because of rounding
                counter = 1
                for edge in all_edges_during_noobs:
                    time_at_thisedge = (m.emission_probability(edge_at_currenttime,noobs)/prob_sum_all_states)*(timegap-1) + time_after_rounding
                    rounded_time_at_thisedge = round(time_at_thisedge)
                    time_after_rounding = time_at_thisedge - rounded_time_at_thisedge
                    #print "time at this edge:"+str(time_at_thisedge) +" "+str(rounded_time_at_thisedge)+" "+str(time_after_rounding)
                    
                    # for each edges write that edge to file for all epochs during no observation
                    for k in range(0, int(rounded_time_at_thisedge)):
                        currenttime = detection_last + counter
                        #print "writing for time"+ str(currenttime)
                        ((lat1,lon1),(lat2,lon2)) = edge
                        write_str = str(currenttime)+","+str( (lat1+lat2)/2 ) + "," + str( (lon1+lon2)/2 )+ "\n"
                        f.write(write_str)
                        counter += 1

        f.close()
                    
                

    def run_viterbi(self, fname):
        counter = 0

        V=None
        p={}

        f = open(read_path+"/"+fname, 'r')
        data = f.readlines();
        f.close()

        #starttime is the obs start time for the current segment of input for viterbi computation
        #it is used in write_viterbi_path() to assign timestamp for each edge in computed viterbi 
        #path corresponding to each observation 
        starttime=int(data[0].rstrip().split(' ')[0]) 


        #input data format: <epoch mon1_lat mon1_lon mon1_ssi mon2_lat mon2_lon mon2_ssi ...>
        for i in range(1, len(data)):
            datalast=data[i-1].rstrip().split(' ')
            datanext=data[i].rstrip().split(' ')

            # save time of last and next observation
            timelast=datalast[0]
            timenext=datanext[0]

            # create a dictionary with observation from all monitors where obs for
            # each monitors is <lat,lon,ssi>
            obslast={}
            obsnext={}
            for k in range(0,TOTAL_MONITORS):
                obslast[k]=[ datalast[3*k+1], datalast[3*k+2], datalast[3*k+3] ]
                obsnext[k]=[ datanext[3*k+1], datanext[3*k+2], datanext[3*k+3] ]
                
            #print i

            V,p = m.step( obslast, V, p )

            timegap = int(timenext)-int(timelast)
            if(timegap>1):
                #create no observation with 0 ssi for all monitors
                noobs={}
                for k in range(0,TOTAL_MONITORS):
                    noobs[k]=[ datalast[3*k+1], datalast[3*k+2], 0 ]
                    noobs[k]=[ datanext[3*k+1], datanext[3*k+2], 0 ]


                for i in range(1, timegap):
                    V,p=m.step( noobs, V, p )



                ############################
                ## Old format viterbi running
                ###########################
                # distance_last_next = fast_distance( (float(latlast), float(lonlast)), (float(latnext), float(lonnext)) )
                # if( ((distance_last_next != 0) and (distance_last_next/timegap < MIN_ALLOWED_SPEED_M_S)) or ((distance_last_next) == 0 and (timegap>420) ) ): 
                #     self.write_viterbi_path(V,p,fname+"_"+str(counter), starttime)
                #     print "terminated currrent viterbi computation and started new one"
                #     counter += 1
                #     starttime=int(timenext) #update starttime for next segment of input
                #     m = WiFiMatcher("../mapinference/trr12/roosevelt_small.osmdb")
                #     V=None
                #     p={}
                
                # else:
                #     for i in range(1, timegap):
                #         #print "lastseen"+" "+latlast+" "+lonlast+" "+latnext+" "+lonnext+" "+str(i)+" "+str(timegap-1)
                #         #V,p = m.step( ("lastseen", ( float(latlast), float(lonlast) ), ( float(latnext), float(lonnext) ), i, timegap-1), V, p )
                #     V,p=m.step( ("no_obs",i), V, p )

                
                    

        self.write_viterbi_path(V,p,fname+"_"+str(counter),starttime)

        #self.write_interpolated_viterbi_path(V,p,fname+"_interpolated_"+str(counter),starttime,data)
        

    def test_reading(self):
        files = os.listdir(read_path)
        #files=open('test/wifi/macs_rest')

        #for f in ['james-drive4']:

        for f in files:
            fname = f.rstrip()
            print fname
            if not (fname.startswith(".")):
                self.run_viterbi(fname)

                    
if __name__ == '__main__':
    unittest.main()
