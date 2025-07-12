#!/usr/bin/env python3

import pytest
from ..reader import *


class TestReader:
    rooted_base_newick_string = "((A,B),(C,D));"
    rooted_extended_newick_string = "[&R] ((A,B),(C,D));"
    unrooted_extended_newick_string = "[&U] ((A,B),(C,D));"

    rooted_base_nexus_string = """#NEXUS\nBegin TREES;\nTree testtree1 = ((A,B),(C,D));\nEND;"""
    rooted_extended_nexus_string = """#NEXUS\nBegin TREES;\nTree testtree1 = [&R] ((A,B),(C,D));\nEND;"""
    unrooted_extended_nexus_string = """#NEXUS\nBegin TREES;\ntree testtree1 = [&U] ((A,B),(C,D));\nEND;"""

    @pytest.mark.prop
    def test_read_rooted_base_newick(self):
        t1, rooting = get_tree(self.rooted_base_newick_string)
        assert rooting == True
        assert len(t1) == 4

    @pytest.mark.prop
    def test_read_rooted_extended_newick(self):
        t1, rooting = get_tree(self.rooted_extended_newick_string)
        assert rooting == True
        assert len(t1) == 4

    @pytest.mark.prop
    def test_read_unrooted_extended_newick(self):
        t1, rooting = get_tree(self.unrooted_extended_newick_string)
        assert rooting == False
        assert len(t1) == 4

    @pytest.mark.prop
    def test_read_rooted_base_nexus(self):
        t1, rooting = get_tree(self.rooted_base_nexus_string)
        assert t1.get_prop("tree_name") == "testtree1"
        assert rooting == True
        assert len(t1) == 4

    @pytest.mark.prop
    def test_read_rooted_extended_nexus(self):
        t1, rooting = get_tree(self.rooted_extended_nexus_string)
        assert t1.get_prop("tree_name") == "testtree1"
        assert rooting == True
        assert len(t1) == 4

    @pytest.mark.prop
    def test_read_unrooted_extended_nexus(self):
        t1, rooting = get_tree(self.unrooted_extended_nexus_string)
        assert t1.get_prop("tree_name") == "testtree1"
        assert rooting == False
        assert len(t1) == 4
