import math
from alive_progress import alive_bar,config_handler
from sys import stderr


## config for progress bar
config_handler.set_global(bar="bubbles",spinner="radioactive",file=stderr)

'''
    cluster_affinity: Tree -> Tree -> int
'''
def rooted_cluster_affinity(t1,t2):
    t1_cmap = t1.get_cached_content(prop="name")
    tree_dist = 0
    with alive_bar(len(t1_cmap)) as bar:
        for i in t1.traverse("postorder"):
            dist = rooted_cdist(t1_cmap[i],t2)
            i.add_prop("ca_dist",dist/max(1,min(len(t1_cmap[i])-1,len(t1)-len(t1_cmap[i]))))
            tree_dist += dist
            bar()
    return tree_dist

def rooted_cluster_support(t1,t2):
    t1_cmap = t1.get_cached_content(prop="name")
    tree_dist = 0
    with alive_bar(len(t1_cmap)) as bar:
        for i in t1.traverse("postorder"):
            dist = rooted_cdist(t1_cmap[i],t2)/len(t1_cmap[i])
            i.add_prop("cs_dist",dist)
            tree_dist += dist
            bar()
    return tree_dist

def unrooted_cluster_affinity(t1,t2):
    t1_cmap = t1.get_cached_content(prop="name")
    n = len(t1)
    tree_dist = 0
    flag = 0
    with alive_bar(n) as bar:
        for i in t1.traverse("postorder"):
            if i.parent and not i.parent.is_root:
                tree_dist += unrooted_cdist(t1_cmap[i],t2,n)
            elif flag == 0:
                tree_dist += unrooted_cdist(t1_cmap[i],t2,n)
                flag = 1
            bar()
    return tree_dist

def unrooted_cdist(c,t2,n):
    mindist = math.inf
    maxdist = 0
    intersection_lookup = dict()
    t2lookup = t2.get_cached_content(prop="name")
    for i in t2.traverse("postorder"):
        intersection = 0
        if i.is_leaf:
            if i.name in c:
                intersection = 1
            else: 
                intersection = 0
        else:
            for ch in i.children:
                intersection += intersection_lookup[ch.id]
        intersection_lookup[i.id] = intersection
        newdist = len(c) + len(t2lookup[i]) - 2*intersection
        if mindist > newdist:
            mindist = newdist
        if maxdist < newdist and len(c) != n:
            maxdist = maxdist
    return min(mindist,n-maxdist)

'''
    rooted_cdist: Cluster -> Tree -> Int
'''
def rooted_cdist(c,t2):
    mindist = math.inf
    intersection_lookup = dict()
    t2lookup = t2.get_cached_content(prop="name")
    for i in t2.traverse("postorder"):
        intersection = 0
        if i.is_leaf:
            if i.name in c:
                intersection = 1
            else:
                intersection = 0
        else:
            for ch in i.children:
                intersection += intersection_lookup[ch.id]
        intersection_lookup[i.id] = intersection
        newdist = len(c) + len(t2lookup[i]) - 2*intersection
        if mindist > newdist:
            mindist=newdist
    return mindist

def calculate_rooted_tau(t):
    tau = 0
    n = len(t)
    sizemap = dict()
    for i in t.traverse("postorder"):
        if i.is_leaf:
            s = 1
        else:
            s = 0
            for ch in i.children:
                s += sizemap[ch.id]
        sizemap[i.id] = s
        tau += min(s-1,n-s)
    return tau 

def calculate_rooted_phi(t):
    phi = 0
    n = len(t)
    sizemap = dict()
    for i in t.traverse("postorder"):
        if i.is_leaf:
            s = 1
        else:
            s = 0
            for ch in i.children:
                s += sizemap[ch.id]
        sizemap[i.id] = s
        phi += min(s-1,n-s)/s
    return phi

def calculate_unrooted_tau(t):
    tau = 0
    n = len(t)
    sizemap = dict()
    flag = 0
    for i in t.traverse("postorder"):
        if i.is_leaf:
            s = 1
        else:
            s = 0
            for ch in i.children:
                s += sizemap[ch.id]
        sizemap[i.id] = s
        if i.parent and not i.parent.is_root:
            tau += min(s-1,n-s-1)
        elif flag == 0:
            tau += min(s-1,n-s-1)
            flag = 1
    return tau


