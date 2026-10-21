import unittest
from ecosystem import experiment,features,dataset
class TestEcosystem(unittest.TestCase):
    def test_reproducible(self):self.assertEqual(experiment(9,'need',None),experiment(9,'need',None))
    def test_expansion(self):
        fixed=experiment(9,'fixed',None);adaptive=experiment(9,'need',None)
        self.assertEqual(adaptive['winner'],'hybrid')
        self.assertLess(adaptive['heldout_mse'],fixed['heldout_mse']*.2)
    def test_features(self):self.assertEqual(features(dataset()[0]).shape,(240,4))
if __name__=='__main__':unittest.main()
