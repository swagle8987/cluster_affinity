import math
import heapq
import numpy as np

'''
    cluster_affinity: Tree -> Tree -> int
'''
def rooted_cluster_affinity(t1,t2):
    t1_cmap = convert_tree_to_cmap(t1)
    tree_dist = 0
    for i in t1_cmap.values():
        tree_dist += cluster_tree_dist(i,t2)
    return tree_dist

'''
    cluster_tree_dist: Cluster -> Tree -> Int
'''
def cluster_tree_dist(c,t2):
    mindist = math.inf
    intersection_lookup = dict()
    size_lookup = dict()
    for i in t2.postorder_node_iter():
        intersection = 0
        size = 0
        if i.is_leaf():
            size = 1
            if i.taxon.label in c:
                intersection = 1
            else:
                intersection = 0
                size_lookup[i] = 1
        else:
            for ch in i.child_node_iter():
                intersection += intersection_lookup[ch]
                size += size_lookup[ch]
        intersection_lookup[i] = intersection
        size_lookup[i] = size
        newdist = len(c) + size - 2*intersection
        if mindist > newdist:
            mindist = newdist
    return mindist

def convert_tree_to_cmap(t):
    cluster_map = dict()
    for i in t.postorder_node_iter():
        if i.is_leaf():
            c = {i.taxon.label}
        else:
            c = set()
            for ch in i.child_node_iter():
                c = c | cluster_map[ch]
        cluster_map[i] = c
    return cluster_map
