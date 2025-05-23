import math
from alive_progress import alive_bar
from sys import stderr

'''
    cluster_affinity: Tree -> Tree -> int
'''
def rooted_cluster_affinity(t1,t2):
    t1_cmap = convert_tree_to_cmap(t1)
    tree_dist = 0
    with alive_bar(len(t1_cmap), bar="bubbles", spinner="radioactive", file=stderr) as bar:
        for i in t1_cmap.values():
            tree_dist += cluster_to_tree_affinity(i,t2)
            bar()
    return tree_dist

def rooted_cluster_support(t1,t2):
    t1_cmap = convert_tree_to_cmap(t1)
    tree_dist = 0
    with alive_bar(len(t1_cmap), bar="bubbles", spinner="radioactive", file=stderr) as bar:
        for i in t1_cmap.values():
            tree_dist += cluster_to_tree_affinity(i,t2)/len(i)
            bar()
    return tree_dist

def unrooted_cluster_affinity(t1,t2):
    t1_cmap = convert_tree_to_cmap(t1)
    n = len(t1.leaf_nodes())
    tree_dist = 0
    flag = 0
    for i in t1_cmap:
        if i._parent_node != t1.seed_node:
            tree_dist += unrooted_cdist(t1_cmap[i],t2,n)
        elif flag == 0:
            tree_dist += unrooted_cdist(t1_cmap[i],t2,n)
            flag = 1
    return tree_dist

def unrooted_cdist(c,t2,n):
    mindist = math.inf
    maxdist = 0
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
        if maxdist < newdist and len(c) != n:
            maxdist = maxdist
    return min(mindist,n-maxdist)

'''
    cluster_to_tree_affinity: Cluster -> Tree -> Int
'''
def cluster_to_tree_affinity(c,t2):
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

def calculate_rooted_tau(t):
    tau = 0
    n = len(t.leaf_nodes())
    sizemap = dict()
    for i in t.postorder_node_iter():
        if i.is_leaf():
            s = 1
        else:
            s = 0
            for ch in i.child_node_iter():
                s += sizemap[ch]
        sizemap[i] = s
        tau += min(s-1,n-s)
    return tau

def calculate_harmonic_number(i):
    gamma = 0.57721566490153286060651209008240243
    return gamma + math.log(i) + (1/(2*i)) + (1/(12*i**2)) + (1/(120*i**4))

def calculate_rooted_phi(t):
    phi = 0
    n = len(t.leaf_nodes())
    sizemap = dict()
    for i in t.postorder_node_iter():
        if i.is_leaf():
            s = 1
        else:
            s = 0
            for ch in i.child_node_iter():
                s += sizemap[ch]
        sizemap[i] = s
        phi += min(s-1,n-s)/s
    return phi

def calculate_unrooted_tau(t):
    tau = 0
    n = len(t.leaf_nodes())
    sizemap = dict()
    flag = 0
    for i in t.postorder_node_iter():
        if i.is_leaf():
            s = 1
        else:
            s = 0
            for ch in i.child_node_iter():
                s += sizemap[ch]
        sizemap[i] = s
        if i._parent_node != t.seed_node:
            tau += min(s-1,n-s-1)
        elif flag == 0:
            tau += min(s-1,n-s-1)
            flag = 1
    return tau


