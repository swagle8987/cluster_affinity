import unittest
from ..cluster_computation import *
from dendropy import Tree,TaxonNamespace
from dendropy.simulate import treesim
import math
import numpy as np

class TestClusterComputation:
    t1 = Tree.get(data="((A,B),(C,D));",schema="newick",rooting="default-rooted")
    t2 = Tree.get(data="((A,C),(B,D));",schema="newick",rooting="default-rooted")

    def test_cluster_affinity_zero(self):
        dist = cluster_affinity.rooted_cluster_affinity(self.t1,self.t1)
        assert dist == 0

    def test_cluster_affinity(self):
        dist = cluster_affinity.rooted_cluster_affinity(self.t1,self.t2)
        assert dist == 2
    
    def test_cluster_affinity_tau(self):
        ntax = 100
        taxon_ns = TaxonNamespace(["l{}".format(i) for i in range(ntax)])
        for i in range(1000):
            t1 = treesim.birth_death_tree(birth_rate=1.0,death_rate=0,num_extant_tips=len(taxon_ns),taxon_namespace=taxon_ns)
            t2 = treesim.birth_death_tree(birth_rate=1.0,death_rate=0,num_extant_tips=len(taxon_ns),taxon_namespace=taxon_ns)
            dist = cluster_affinity.rooted_cluster_affinity(t1,t2)
            assert dist >= 0,"{} {} {} {}".format(t1.as_string(schema="newick"),
                                               t2.as_string(schema="newick"), 
                                               (np.min(cluster_affinity.rooted_cluster_affinity_matrix(t1,t2),axis=0)),
                                               cluster_affinity.rooted_cluster_affinity(t1,t2))
            assert dist <= math.ceil(ntax*ntax - 2*ntax)/4
