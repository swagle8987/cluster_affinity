#!/usr/bin/env python3
#
#

from ..cluster_computation.heavy_path_decompositions import *
from ..cluster_computation.transfer_index import *
import pytest
from ete4 import Tree


class TestHeavyTreeCreation:
    t1 = Tree("(((A,B),C),((D,E),F));")

    @pytest.mark.prop
    def test_get_maximal_paths(self):
        annotate_heavy_nodes(self.t1)
        paths = get_maximal_heavy_paths(self.t1)
        assert len(paths) == len(self.t1)
        len_paths = [len(i) for i in paths]
        assert min(len_paths) ==1
        assert max(len_paths) ==4
        for i in paths:
            assert i[-1].is_leaf
            assert not i[0].get_prop("is_heavy")

    @pytest.mark.prop
    def test_PathTree(self):
        p = get_maximal_heavy_paths(self.t1)
        ptrees = []
        for i in p:
            ptrees.append(PathTree(i))
        l = list(self.t1.leaves())[0]
        path_search_tree = l.get_prop("path_tree")
        assert path_search_tree.get_first() == path_search_tree.get_last()
        assert path_search_tree.get_last() == l
        assert path_search_tree.D == 1
        v = l.parent
        path_parent_tree = v.get_prop("path_tree")
        assert path_parent_tree.get_first() == path_parent_tree.get_last()
        assert path_parent_tree.get_last() == v
        assert path_parent_tree.D == len(v)

    @pytest.mark.prop
    def test_update_path(self):
        p = get_maximal_heavy_paths(self.t1)
        ptrees = []
        for i in p:
            ptrees.append(PathTree(i))
        l = list(self.t1.search_leaves_by_name("A"))[0]
        path_search_tree = l.get_prop("path_tree")
        update_path(l,-2)
        assert path_search_tree.D == -1
        v = l.parent
        path_parent_tree = v.get_prop("path_tree")
        assert path_parent_tree.D == len(v)
        parent_search_tree = path_search_tree.get_parent()


    @pytest.mark.slow
    def test_get_maximal_paths_exhaustive(self):
        ntax=100
        labels = ["l{}".format(i) for i in range(ntax)]
        for i in range(1000):
            t1 = Tree()
            t1.populate(ntax,names=labels)
            annotate_heavy_nodes(t1)
            paths = get_maximal_heavy_paths(t1)
            assert len(paths) == len(t1)
            len_paths = [len(i) for i in paths]
            for i in paths:
                assert i.get_last().is_leaf
                assert not i.get_first().get_prop("is_heavy")
