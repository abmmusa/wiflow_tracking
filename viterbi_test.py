import unittest
from viterbi import Viterbi

class ViterbiTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_wikipedia(self):
        hmm = {
            'Rainy' : [('Rainy', 0.7),('Sunny', 0.3)],
            'Sunny' : [('Rainy', 0.4),('Sunny', 0.6)]}

        start_probabilities={'Rainy' : 0.6, 'Sunny': 0.4}
        
        emission_probabilities = {'Rainy' : {'walk':0.1,'shop':0.4,'clean':0.5},
                                  'Sunny' : {'walk':0.6,'shop':0.3,'clean':0.1}}

        vit = Viterbi(hmm, 
                      lambda state, obs: emission_probabilities[state][obs])

        (v,p) = vit.step('walk',start_probabilities)        
        (v,p) = vit.step('shop',v,p)
        (v,p) = vit.step('clean',v,p)

        max_state = max(v,key=lambda x:v[x])
        assert(p[max_state]==['Sunny','Rainy','Rainy'])

if __name__ == '__main__':
    unittest.main()
