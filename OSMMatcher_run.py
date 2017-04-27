from OSMMatcher import OSMMatcher
import spatialfunclib
import math

MAX_SPEED_M_PER_S = 20.0

class MatchOSM:
    def __init__(self, osmdb_filename, constraint_length, max_dist):
        self.matcher = OSMMatcher(osmdb_filename, constraint_length, max_dist)
        self.constraint_length = constraint_length
    
    def process_trip(self, trip_directory, trip_filename):
        
        trip_file = open(trip_directory + "/" + trip_filename, 'r')
        raw_observations = map(lambda x: x.strip("\n").split(" "), trip_file.readlines())
        raw_observations_count = len(raw_observations)
        trip_file.close()
        
        V=None
        p={}
        
        obs = []
        obs_states = []
        max_prob_p = None
        
        for i in range(1, raw_observations_count):
            (prev_lat, prev_lon, prev_time) = raw_observations[i - 1]
            (curr_lat, curr_lon, curr_time) = raw_observations[i]
            
            elapsed_time = int(float(curr_time) - float(prev_time))
            
            if (i > 1):
                start_time = 1
            else:
                start_time = 0
            
            (prev_int_lat, prev_int_lon, prev_int_time) = (prev_lat, prev_lon, prev_time)
            
            for j in range(start_time, (elapsed_time + 1)):
                
                fraction_along = (float(j) / float(elapsed_time))
                (int_lat, int_lon) = spatialfunclib.point_along_line(float(prev_lat), float(prev_lon), float(curr_lat), float(curr_lon), fraction_along)
                
                int_distance = spatialfunclib.haversine_distance((float(prev_int_lat), float(prev_int_lon)), (float(int_lat), float(int_lon)))
                
                if (int_distance > MAX_SPEED_M_PER_S):
                    
                    int_steps = int(math.ceil(int_distance / MAX_SPEED_M_PER_S))
                    int_step_distance = (int_distance / float(int_steps))
                    int_step_time = (float(elapsed_time) / float(int_steps))
                    
                    for k in range(1, int_steps):
                        step_fraction_along = ((k * int_step_distance) / int_distance)
                        (step_int_lat, step_int_lon) = spatialfunclib.point_along_line(float(prev_int_lat), float(prev_int_lon), float(int_lat), float(int_lon), step_fraction_along)
                        
                        (V, p) = self.matcher.step((float(step_int_lat), float(step_int_lon)), V, p)
                        
                        max_prob_state = max(V, key=lambda x: V[x])
                        #max_prob_p = p[max_prob_state]
                        #
                        #if (len(max_prob_p) == self.constraint_length):
                        #    obs_states.append(max_prob_p[0])
                        obs_states.append(max_prob_state)
                        
                        obs.append((step_int_lat, step_int_lon, (float(prev_int_time) + (float(k) * int_step_time))))
                
                (V, p) = self.matcher.step((float(int_lat), float(int_lon)), V, p)
                
                max_prob_state = max(V, key=lambda x: V[x])
                #max_prob_p = p[max_prob_state]
                #
                #if (len(max_prob_p) == self.constraint_length):
                #    obs_states.append(max_prob_p[0])
                obs_states.append(max_prob_state)
                
                obs.append((int_lat, int_lon, (int(prev_time) + j)))
                (prev_int_lat, prev_int_lon, prev_int_time) = (int_lat, int_lon, (int(prev_time) + j))
        
        #if (len(max_prob_p) < self.constraint_length):
        #    obs_states.extend(max_prob_p)
        #else:
        #    obs_states.extend(max_prob_p[1:])
        
        #print "obs: " + str(len(obs))
        #print "obs states: " + str(len(obs_states))
        assert(len(obs_states) == len(obs))
        
        out_file = open(trip_directory + "/matched_" + trip_filename, 'w')
        
        for i in range(0, len(obs)):
            (obs_lat, obs_lon, obs_time) = obs[i]
            out_file.write(str(obs_lat) + " " + str(obs_lon) + " " + str(obs_time) + " ")
            
            if (obs_states[i] == "unknown"):
                out_file.write(str(obs_states[i]) + "\n")
            else:
                (in_node_coords, out_node_coords) = obs_states[i]
                out_file.write(str(in_node_coords[0]) + " " + str(in_node_coords[1]) + " " + str(out_node_coords[0]) + " " + str(out_node_coords[1]) + "\n")
        
        out_file.close()

import sys, getopt
import os
if __name__ == '__main__':
    
    osmdb_filename = "planet-120328_moscow_large.osmdb"
    constraint_length = 300
    max_dist = 350
    trip_directory = "moscow_trips_large/0/"
    
    (opts, args) = getopt.getopt(sys.argv[1:],"c:m:o:t:h")
    
    for o,a in opts:
        if o == "-c":
            constraint_length = int(a)
        elif o == "-m":
            max_dist = int(a)
        elif o == "-o":
            osmdb_filename = str(a)
        elif o == "-t":
            trip_directory = str(a)
        elif o == "-h":
            print "Usage: python OSMMatcher_run.py [-c <constraint_length>] [-m <max_dist>] [-o <osmdb_filename>] [-t <trip_directory>] [-h]"
            exit()
    
    #print "constraint length: " + str(constraint_length)
    #print "max dist: " + str(max_dist)
    #print "osmdb filename: " + str(osmdb_filename)
    #print "trip directory: " + str(trip_directory)
    
    match_osm = MatchOSM(osmdb_filename, constraint_length, max_dist)
    
    all_trip_files = filter(lambda x: x.startswith("trip_") and x.endswith(".txt"), os.listdir(trip_directory))
    
    for i in range(0, len(all_trip_files)):
        sys.stdout.write("\rProcessing trip " + str(i + 1) + "/" + str(len(all_trip_files)) + "... ")
        sys.stdout.flush()
        
        match_osm.process_trip(trip_directory, all_trip_files[i])
    
    sys.stdout.write("done.\n")
    sys.stdout.flush()
