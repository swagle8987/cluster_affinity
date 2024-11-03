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


def convert_tree_size_array(tree):
    tree_lookup = dict()
    tree_array = np.full(len(tree.nodes()),0)
    ind = 0
    for i in tree.leaf_node_iter():
        tree_lookup[i] = ind
        tree_array[ind] = 1
        ind += 1
    ## We do internal nodes seperately to ensure that the leaf nodes form a matrix at the top
    for j in tree.postorder_internal_node_iter():
        tree_lookup[i] = ind
        for ch in i.child_node_iter():
                tree_array[ind] += tree_array[tree_lookup[ch]]
    return tree_array,tree_lookup

def construct_intersection_matrix(a1,a2):
    return np.full((len(a1),len(a2)),0)

def construct_distance_matrix(a1,a2):
    arr = np.full((len(a1),len(a2)),0)
    arr += a1
    arr += a2[:,np.newaxis]
    return arr

def rooted_cluster_affinity_matrix(t1,t2):
    t1_arr,t1_lookup = convert_tree_size_array(t1)
    t2_arr,t2_lookup = convert_tree_size_array(t2)
    intersection_mat = construct_intersection_matrix(t1_arr,t2_arr)
    distance_mat = construct_distance_matrix(t1_arr,t2_arr)
    for i in t1.leaf_node_iter():
        for j in t2.leaf_node_iter():
            if i.taxon.label == j.taxon.label:
                intersection_mat[t1_lookup[i],t2_lookup[j]] = 1
            else:
                intersection_mat[t1_lookup[i],t2_lookup[j]] = 0
    assert np.all(np.min(intersection_mat,axis=0) < 1)
    breakpoint()
    for u in t2.postorder_internal_node_iter():
        for ch in u.child_node_iter():
            intersection_mat[:,t2_lookup[u]] += intersection_mat[:,t2_lookup[ch]]
    assert np.all(np.min(intersection_mat,axis=0) >= 0)
    breakpoint()
    for v in t1.postorder_internal_node_iter():
        for ch in v.child_node_iter():
            intersection_mat[t1_lookup[v],:] += intersection_mat[t1_lookup[ch],:]
    assert np.all(np.min(intersection_mat,axis=0) >= 0)
    breakpoint()
    assert np.all(2*intersection_mat <= distance_mat), breakpoint()
    distance_mat -= 2*intersection_mat
    return distance_mat

#def rooted_cluster_affinity(t1,t2):
#    return np.sum(np.min(rooted_cluster_affinity_matrix(t1,t2),axis=0))




