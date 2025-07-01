#!/usr/bin/env python3
#
import ete4
from .heavy_path_decompositions import *

def update_path(x,d):
    path_interval = x.get_prop("path_tree")
    path_interval.D += d
    path_interval.minval += d
    path_interval.maxval += d
    while path_interval.parent:
        parent_interval = path_interval.parent
        sibling_interval = path_interval.get_sibling()
        if parent_interval.get_last() == path_interval.get_last():
            sibling_interval.D += d
            sibling_interval.minval += d
            sibling_interval.maxval += d
        parent_interval.minval = min(path_interval.minval,sibling_interval.minval) + parent_interval.D
        parent_interval.maxval = max(path_interval.maxval,sibling_interval.maxval) + parent_interval.D
        path_interval = parent_interval


def add_leaf_general(x,counter):
    update_path(x,-2)
    path = x.get_prop("path")
    while not path.get_first().is_root():
        x = path.get_first().parent
        update_path(x,-2)
    counter = counter+1
    return counter


def remove_leaf_general(x,counter):
    update_path(x,2)
    path = x.get_prop("path")
    while not path.get_first().is_root():
        x = path.get_first().parent
        update_path(x,2)
    counter = counter-1
    return counter


def construct_leaf_map(t):
    map = dict()
    for i in t:
        map[i.name] = i
    return map

def compute_transfer_index(t1,t2):
    annotate_heavy_nodes(t1)
    paths = get_maximal_heavy_paths(t2)
    path_search_trees = [PathTree(i) for i in paths]
    lmap = construct_leaf_map(t2)
    counter = 0
    udist = 0
    rdist = 0
    for i in t1.leaf_names():
        curNode = i
        x = lmap[i]
        counter = add_leaf_general(x,counter)
        while not curNode.is_root:
            if curNode.get_prop("is_heavy"):
                siblingNode = curNode.get_sisters()[0]
                for l in siblingNode:
                    counter=add_leaf_general(lmap[l],counter)
                min1 = min([i.minval for i in path_search_trees]) + counter
                min2 = n- max([i.maxval for i in path_search_trees]) -counter
                rdist += min1
                udist += min(min1,min2)
                curNode = curNode.parent
            else:
                for l in curNode:
                    counter = remove_leaf_general(lmap[l],counter)

    return rdist,udist
