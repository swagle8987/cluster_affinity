from ca_calculator.rooted_affinity import cluster_affinity,cluster_affinity_cluster
import dendropy
from dendropy.datamodel.taxonmodel import Taxon
import unittest
import sys
import pdb


class TestStringMethods(unittest.TestCase):
    
    def test_no_cluster_match(self):
        t1_string = "((A,B),(C,D));"
        t2_string = "(A,C);"
        t1 = dendropy.Tree.get(data = t1_string, schema="newick", rooting="default-rooted")
        cluster = dendropy.Tree.get(data = t2_string, schema="newick", rooting="default-rooted")
        self.assertEqual(cluster_affinity(cluster,t1),1)

    def test_full_cluster_match(self):
        t1_string = "((A,B),(C,D));"
        t1 = dendropy.Tree.get(data = t1_string, schema="newick", rooting="default-rooted")
        t2_string = "(A,B);"
        cluster = dendropy.Tree.get(data = t2_string, schema="newick", rooting="default-rooted")
        self.assertEqual(cluster_affinity(cluster,t1),0)

    def test_to_diameter(self):
        t1_string = "(A,(B,(C,D)));"
        t2_string = "(D,(C,(A,B)));"
        t1 = dendropy.Tree.get(data = t1_string, schema="newick", rooting="default-rooted")
        t2 = dendropy.Tree.get(data = t2_string, schema="newick", rooting="default-rooted")
        self.assertEqual(cluster_affinity(t1,t2),2)

    def test_no_distance(self):
        t1_string = "(A,(B,(C,D)));"
        t1 = dendropy.Tree.get(data = t1_string, schema="newick", rooting="default-rooted")
        t2 = dendropy.Tree.get(data = t1_string, schema="newick", rooting="default-rooted")
        self.assertEqual(cluster_affinity(t1,t1),0)

    def test_polytomies_but_equal(self):
        t1_string = "(A,(B,C,D));"
        t2_string = "(A,(B,(C,D)));"
        t1 = dendropy.Tree.get(data = t1_string, schema="newick", rooting="default-rooted")
        t2 = dendropy.Tree.get(data = t2_string, schema="newick", rooting="default-rooted")
        self.assertEqual(cluster_affinity(t1,t2),0)

    def test_polytomies_but_unequal(self):
        t1_string = "(D,(B,C,A));"
        t2_string = "(A,(B,(C,D)));"
        t1 = dendropy.Tree.get(data = t1_string, schema="newick", rooting="default-rooted")
        t2 = dendropy.Tree.get(data = t2_string, schema="newick", rooting="default-rooted")
        self.assertEqual(cluster_affinity(t1,t2),1)

if __name__=='__main__':
    unittest.main()
