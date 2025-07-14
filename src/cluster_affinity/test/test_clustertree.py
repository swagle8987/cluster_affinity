from ..cluster_computation import *
from ete4 import Tree
import math
import pytest


class TestClusterComputation:
    t1 = Tree("((A,B),(C,D));")
    t2 = Tree("((A,C),(B,D));")

    @pytest.mark.prop
    def test_unrooted_cdist(self):
        dist = cluster_affinity.unrooted_cdist(set(["A", "C"]), self.t2, 4)
        assert dist == 0

    @pytest.mark.prop
    def test_rooted_cdist(self):
        dist = cluster_affinity.rooted_cdist(set(["A", "C"]), self.t2)
        assert dist == 0

    @pytest.mark.prop
    def test_cluster_affinity_zero(self):
        dist = cluster_affinity.rooted_cluster_affinity(self.t1, self.t1)
        assert dist == 0

    @pytest.mark.prop
    def test_cluster_affinity(self):
        dist = cluster_affinity.rooted_cluster_affinity(self.t1, self.t2)
        assert dist == 2

    @pytest.mark.prop
    def test_cluster_support_zero(self):
        dist = cluster_affinity.rooted_cluster_support(self.t1, self.t1)
        assert dist == 0

    @pytest.mark.prop
    def test_cluster_support(self):
        dist = cluster_affinity.rooted_cluster_support(self.t1, self.t2)
        assert dist == 1

    @pytest.mark.prop
    def test_unrooted_cluster_affinity_zero(self):
        dist = cluster_affinity.unrooted_cluster_affinity(self.t1, self.t1)
        assert dist == 0

    @pytest.mark.prop
    def test_unrooted_cluster_affinity(self):
        dist = cluster_affinity.unrooted_cluster_affinity(self.t1, self.t2)
        assert dist == 2

    @pytest.mark.slow
    def test_cluster_affinity_tau(self):
        ntax = 100
        labels = ["l{}".format(i) for i in range(ntax)]
        for i in range(1000):
            t1 = Tree()
            t1.populate(ntax, names=labels)
            t2 = Tree()
            t2.populate(ntax, names=labels)
            dist = cluster_affinity.rooted_cluster_affinity(t1, t2)
            tau = cluster_affinity.calculate_rooted_tau(t1)
            assert dist >= 0, "{} {} {}".format(
                t1.write(), t2.write(), cluster_affinity.rooted_cluster_affinity(t1, t2)
            )
            assert dist <= math.ceil(ntax * ntax - 2 * ntax) / 4
            assert dist <= tau

    @pytest.mark.slow
    def test_cluster_support_phi(self):
        ntax = 100
        labels = ["l{}".format(i) for i in range(ntax)]
        for i in range(1000):
            t1 = Tree()
            t1.populate(ntax, names=labels)
            t2 = Tree()
            t2.populate(ntax, names=labels)
            dist = cluster_affinity.rooted_cluster_support(t1, t2)
            phi = cluster_affinity.calculate_rooted_phi(t1)
            assert dist >= 0, "{} {} {}".format(
                t1.write(), t2.write(), cluster_affinity.rooted_cluster_affinity(t1, t2)
            )
            assert dist <= phi

    @pytest.mark.prop
    def test_transfer_index_zero(self):
        dist = transfer_index.compute_transfer_index(self.t1, self.t1,cost="Rooted Cluster Affinity")
        assert dist == 0

    @pytest.mark.prop
    def test_transfer_index_nonzero(self):
        dist = transfer_index.compute_transfer_index(self.t1, self.t2,cost="Rooted Cluster Affinity")
        assert dist == 2

    @pytest.mark.slow
    def test_transfer_index_mock_value(self):
        ntax = 100
        labels = ["l{}".format(i) for i in range(ntax)]
        for i in range(1000):
            t1 = Tree()
            t1.populate(ntax, names=labels)
            t2 = Tree()
            t2.populate(ntax, names=labels)
            dist = cluster_affinity.rooted_cluster_affinity(t1, t2)
            tau = cluster_affinity.calculate_rooted_tau(t1)
            rdist = transfer_index.compute_transfer_index(
                t1, t2, cost="Rooted Cluster Affinity"
            )
            print(dist, rdist)
            assert dist == rdist, "{} {} {}".format(i, dist, rdist)
            assert rdist <= tau

    @pytest.mark.slow
    @pytest.mark.optional
    def test_transfer_index_unrooted_mock_value(self):
        ntax = 100
        labels = ["l{}".format(i) for i in range(ntax)]
        for i in range(1000):
            t1 = Tree()
            t1.populate(ntax, names=labels)
            t2 = Tree()
            t2.populate(ntax, names=labels)
            dist = cluster_affinity.unrooted_cluster_affinity(t1, t2)
            tau = cluster_affinity.calculate_unrooted_tau(t1)
            rdist = transfer_index.compute_transfer_index(
                t1, t2, cost="Unrooted Cluster Affinity"
            )
            print(dist,rdist)
            assert dist == rdist, "{} {} {}".format(i,dist,rdist)
            assert dist <= tau
            assert rdist<=tau
