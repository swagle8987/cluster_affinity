#!/usr/bin/env python3

import ete4
import math

def annotate_heavy_nodes(t):
    count = 0
    t1_cmap = t.get_cached_content(prop="name")
    for n in t.traverse("postorder"):
        n.add_prop("path_id",count)
        if n.is_leaf:
            pass
        else:
            u,v = n.children
            if len(t1_cmap[u]) <= len(t1_cmap[v]):
                u.add_prop("is_heavy",False)
                v.add_prop("is_heavy",True)
            else:
                u.add_prop("is_heavy",True)
                v.add_prop("is_heavy",False)


def get_maximal_heavy_paths(t):
    paths = []
    for l in t:
        v = l
        p = [l]
        while v.is_leaf or v.get_prop("is_heavy",default=False):
            v = v.parent
            p.append(v)
        paths.append(p)
    return paths

class PathTree:
    def __init__(self,p,parent=None):
        self.rootpath = HeavyPath(p)
        self.parent = parent
        if len(p)>1:
            self.left_child = PathTree(p[:(math.floor(len(p)/2))],self)
            self.right_child = PathTree(p[(math.floor(len(p)/2)):],self)
            self.is_leaf = False
            self.D = 0
        else:
            self.is_leaf = True
            self.node = p[0]
            self.node.add_prop("path_tree",self)
            self.D = len(self.node)
        self.minval = len(self.rootpath.get_last())
        self.maxval = len(self.rootpath.get_first())


    def get_first(self):
        return self.rootpath.get_first()

    def get_last(self):
        return self.rootpath.get_last()

    def get_parent(self):
        return self.parent

    def get_sibling(self):
        if self.parent.get_first() == self.get_first():
            return self.parent.right_child
        else:
            return self.parent.left_child

class HeavyPath:
    def __init__(self,path=[]):
        if path:
            self.first = path[0]
            self.last = path[-1]
        else:
            self.first = None
            self.last=None
        self.nodes = set(path)


    def get_first(self):
        return self.first

    def get_last(self):
        return self.last

    def contains(self,x):
        return x in self.nodes
