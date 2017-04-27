import sys, getopt
import os
from ImageMatcher import ImageMatcher

data_ref='data_uic2/time_loc_frame_interpolated_smaller_round1.txt'
data_test='data_uic2/time_loc_frame_interpolated_smaller_round2.txt'
pickle_ref='pickles_uic2/keypoints_descriptors_smaller_round_1.pkl'
pickle_test='pickles_uic2/keypoints_descriptors_smaller_round_2.pkl'

print data_ref

if __name__ == '__main__':
    print "starting..."
    ImageMatcher(data_ref, pickle_ref, pickle_test)
