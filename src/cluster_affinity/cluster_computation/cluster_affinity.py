import math
from alive_progress import alive_bar, config_handler
from sys import stderr
from typing import TypeAlias
import ete4


cluster: TypeAlias = set[str]

## config for progress bar
config_handler.set_global(bar="bubbles", spinner="radioactive", file=stderr)


def rooted_cluster_affinity(t1: ete4.Tree, t2: ete4.Tree, disable_bar=False) -> float:
    t1_cmap = t1.get_cached_content(prop="name")
    tree_dist = 0
    with alive_bar(len(t1_cmap), disable=disable_bar) as bar:
        for i in t1.traverse("postorder"):
            dist = rooted_cdist(t1_cmap[i], t2)
            i.add_prop(
                "c_dist",
                dist / max(1, min(len(t1_cmap[i]) - 1, len(t1) - len(t1_cmap[i]))),
            )
            tree_dist += dist
            bar()
    return tree_dist


def rooted_cluster_support(t1: ete4.Tree, t2: ete4.Tree) -> float:
    t1_cmap = t1.get_cached_content(prop="name")
    tree_dist = 0
    with alive_bar(len(t1_cmap)) as bar:
        for i in t1.traverse("postorder"):
            dist = rooted_cdist(t1_cmap[i], t2) / len(t1_cmap[i])
            i.add_prop("c_dist", dist)
            tree_dist += dist
            bar()
    return tree_dist


def unrooted_cluster_affinity(t1, t2):
    t1_bipartitions = t1.edges()
    t2_bipartitions = list(t2.edges())
    n = len(t1)
    tree_dist = 0
    with alive_bar(n) as bar:
        for i in t1_bipartitions:
            if len(i[1]) > 0 and len(i[1])<n:
                mindist = math.inf
                for j in t2_bipartitions:
                    if len(j[1])>0 and len(j[1])<n:
                        c = set([l.name for l in i[1]])
                        x = set([l.name for l in j[1]])
                        cdx=len(c^x)
                        cdist = min(cdx, n-cdx)
                        if cdist<0:
                            raise RuntimeError()
                        if cdist < mindist:
                            mindist = cdist
                tree_dist += mindist
    return tree_dist

def unrooted_cdist(c, t2, n):
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
        newdist = len(c) + len(t2lookup[i]) - 2 * intersection
        if mindist > newdist:
            mindist = newdist
        if maxdist < newdist and len(c) != n:
            maxdist = maxdist
    return min(mindist, n - maxdist)


"""
    rooted_cdist: Cluster -> Tree -> Int
"""


def rooted_cdist(c: cluster, t2: ete4.Tree) -> int:
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
        newdist = len(c) + len(t2lookup[i]) - 2 * intersection
        if mindist > newdist:
            mindist = newdist
    return mindist


def calculate_rooted_tau(t: ete4.Tree) -> int:
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
        tau += min(s - 1, n - s)
    return tau


def calculate_rooted_phi(t: ete4.Tree) -> float:
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
        phi += min(s - 1, n - s) / s
    return phi


def calculate_unrooted_tau(t):
    tau = 0
    n = len(t)
    sizemap = dict()
    flag = 0
    for i in t.edges():
        if len(i[0]) > 0 and len(i[0]) < n:
            tau += min(len(i[0])-1,n-len(i[0])-1)
    return tau


def calculate_unrooted_phi(t):
    phi = 0
    n = len(t)
    sizemap = dict()
    flag = 0
    for i in t.edges():
        if len(i[0]) > 0 and len(i[0]) < n:
            phi += min(len(i[0])-1,n-len(i[0])-1)/len(i[0])
    return phi
