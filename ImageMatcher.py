from streetmap import *
from gpsmatcher import GPSMatcher
from spatialfunclib import *
from pylibs import spatialfunclib
from mathfunclib import *
from viterbi import Viterbi
import time

import pickle
import sys
sys.path.insert(0, 'image_processing')
from image_processor import ImageProcessor


MAX_DISTANCE_SIGMA = 20.0

EMISSION_SIGMA = 100.0
DISTANCE_SIGMA = 300.0

EMISSION_UNKNOWN = 0.00001

TRANSITION_UNKNOWN = 0.00001
TRANSITION_UNKNOWN_UNKNOWN = 0.9
TRANSITION_SELF = .5

MAX_SPEED_M_PER_S = 20.0

TRANSMISSION_PROB=1.0/300.0

debug=False

class ImageMatcher:
    def __init__(self, ref_frames_data_filename, ref_pickle_filename, test_pickle_filename):
        print "init..."
        self.previous_obs = None
        self.image_processor = ImageProcessor()

        self.descriptors_ref = self.image_processor.load_sift(ref_pickle_filename)
        self.descriptors_test = self.image_processor.load_sift(test_pickle_filename)

        hmm = self.ref_frames_data_to_hmm(ref_frames_data_filename)
        
        #emission_probabilities = map(lambda x: complementary_normal_distribution_cdf(x,0,EMISSION_SIGMA),range(0,int(3.0*EMISSION_SIGMA)))
            
        priors=dict([(state,1.0/len(hmm)) for state in hmm])
        
        self.viterbi = Viterbi(hmm,self.emission_probability,
                          constraint_length=2500, # BE CAREFUL with it. walking may take long time and higher value may be needed here
                          priors=priors)


    def emission_probability(self, frame_ref, frame_test): #frame_ref is the state and frame_test is the obs
        des_test = descriptors_test[frame_test]
        des_ref = descriptors_ref[frame_ref]

        match_count=self.find_match_from_des(des_ref, des_test) 

        return match_count/1000.0 


    # def emission_probability(self, edge, obs):            
    #     if(edge=='unknown'):
    #         return EMISSION_UNKNOWN

    #     (edge_end1,edge_end2)=edge
            
    #     emission_prob=1.0;
    #     for obs_single_mon in obs.values():
    #         (lat,lon,ssi)=map(lambda x:float(x), obs_single_mon)
                
    #         dist1=fast_distance(edge_end1, (lat,lon))
    #         dist2=fast_distance(edge_end2, (lat,lon))
                
    #         prob_for_mon = self.prob_mon_given_segment(dist1, dist2, ssi)
    #         #if prob_for_mon != 0 and ssi !=0:
    #         #    print ssi,prob_for_mon
    #         emission_prob *= prob_for_mon

    #     #multiply by transmission prob    
    #     emission_prob *= TRANSMISSION_PROB    
            
    #     #if emission_prob!=0.0:
    #     #    print emission_prob
    #     return emission_prob




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
    



    def ref_frames_data_to_hmm(self,ref_frames_data_filename):
        # read ref data from file and put in ref_frames_data list
        ref_frames_data = []
        data_ref = open(ref_frames_data_filename, 'r').readlines()
        for line in data_ref:
            data_r = line.strip().split(' ')
            time_ref, lat_ref, lon_ref, frame_ref = float(data_r[0]), float(data_r[1]), float(data_r[2]), int(data_r[3]) 
            ref_frames_data.append((time_ref, lat_ref, lon_ref, frame_ref))


        # create the hmm
        hmm = {}    
        all_ref_frames = []
        # transition_probabilities = map(lambda x: complementary_normal_distribution_cdf(x,0,MAX_DISTANCE_SIGMA),range(0,int(3.0*MAX_DISTANCE_SIGMA)))

        for i in range(len(ref_frames_data)):
            #print i, ref_frames_data[i][3]
            next_frames_w_prob = []
            for j in range(i, len(ref_frames_data)):
                #print i, j, ref_frames_data[i][3], ref_frames_data[j][3], len(ref_frames_data)

                time_from, lat_from, lon_from, frame_from = ref_frames_data[i][0], ref_frames_data[i][1], ref_frames_data[i][2], ref_frames_data[i][3]
                time_to, lat_to, lon_to, frame_to = ref_frames_data[j][0], ref_frames_data[j][1], ref_frames_data[j][2], ref_frames_data[j][3]
                all_ref_frames.append(frame_from)
                distance=spatialfunclib.distance4(lat_from, lon_from, lat_to, lon_to)

                if distance < 3*MAX_DISTANCE_SIGMA:
                    transition_prob = complementary_normal_distribution_cdf(distance,0,MAX_DISTANCE_SIGMA)
                    next_frames_w_prob.append((frame_to, transition_prob))

                    if j == len(ref_frames_data)-1:
                        hmm[frame_from] = next_frames_w_prob
                        continue
                else:
                    hmm[frame_from] = next_frames_w_prob
                                
        # TODO: deal with unknown
        #hmm['unknown']=[('unknown',TRANSITION_UNKNOWN_UNKNOWN)]
        #hmm['unknown']+=[(frame,(1-TRANSITION_UNKNOWN_UNKNOWN)/len(all_ref_frames)) for frame in all_ref_frames]
            
        if debug:
            for key in hmm.keys():
                print key, hmm[key]
                print "\n"


        return hmm
        


#WiFiMatcher("../mapgenerate/ground_truth.osmdb")
