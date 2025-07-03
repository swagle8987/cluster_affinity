#!/usr/bin/env python3

import math
from collections import deque
from itertools import islice


def annotate_heavy_nodes(t):
    count = 0
    t1_cmap = t.get_cached_content(prop="name")
    for n in t.traverse("postorder"):
        n.add_prop("path_id",count)
        count += 1
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



def get_heavy_child(node):
    children = node.children
    light_children = [i for i in children if not i.get_prop("is_heavy")]
    for i in children:
        if i.get_prop("is_heavy"):
            heavy_child = i
            break
    return heavy_child,light_children


def get_maximal_heavy_paths(t):
    annotate_heavy_nodes(t)
    nodes = deque()
    paths = []
    nodes.append(t)
    while nodes:
        v = nodes.popleft()
        p = HeavyPath([v])
        while not v.is_leaf:
            heavy_child,light_children = get_heavy_child(v)
            nodes.extend(light_children)
            v = heavy_child
            p.append(heavy_child)
        paths.append(p)
    return paths

class PathTree:
    def __init__(self,p,parent=None):
        self.rootpath = p
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
        self.minval = len(self.rootpath[-1])
        self.maxval = len(self.rootpath[0])


    def get_first(self):
        return self.rootpath[0]

    def get_last(self):
        return self.rootpath[-1]

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
        for v in path:
            v.add_prop("path",self)
        self.nodes = set(path)
        self.path = path


    def get_first(self):
        return self.first

    def get_last(self):
        return self.last

    def __contains__(self,x):
        return x in self.nodes

    def __getitem__(self,key):
        return self.path[key]

    def slice(self,stop):
        return self.path[:stop]

    def slice(self,start,stop,step=None):
        return self.path[start:stop:step]

    def __setitem__(self,key,value):
        self.path[key] = value

    def __len__(self):
        return len(self.path)

    def __str__(self):
        return " -> ".join([str(i.get_prop("path_id")) for i in self.path])

    def __repr__(self):
        return "Heavy Path {}".format(self.get_first().name)

    def append(self,x):
        self.last = x
        self.nodes.add(x)
        self.path.append(x)
        x.add_prop("path",self)
