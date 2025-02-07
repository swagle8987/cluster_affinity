from extendedtree import ExtendedTree, NewNode
from math import floor

class IntervalNode:

    def __init__(self, a, b) -> None:
       self.start = a
       self.end = b
       self.left_child = None
       self.right_child = None
       self.parent = None

    def set_left_child(self,lnode):
        self.left_child = lnode
        lnode.set_parent(self)

    def set_right_child(self,rnode):
        self.right_child = rnode
        rnode.set_parent(self)

    def set_parent(self, parent):
        self.parent = parent

    def get_interval(self):
        return (self.start,self.end)

class PathSearchTree:
    def __init__(self,path):
        self.nodes = set(path)
        self.path = path
        self.root = IntervalNode(0,len(path)-1)
        self.D = dict()
        self.minval = dict()
        self.maxval = dict()
        self.interval_lookup = {(self.root.start,self.root.end):self.root}
        nnodes = [self.root]
        while nnodes:
            n = nnodes.pop()
            self.minval[n.start,n.end] = self.path[n.end].size
            self.maxval[n.start,n.end] = self.path[n.start].size
            if n.start != n.end:
                self.D[n.start,n.end] = 0
                l = n.end - n.start
                lnode = IntervalNode(n.start,n.start + floor(l/2))
                rnode = IntervalNode(n.start+floor(l/2)+1,n.end)
                n.set_left_child(lnode)
                n.set_right_child(rnode)
                self.interval_lookup[(lnode.start,lnode.end)] = lnode
                self.interval_lookup[(rnode.start,rnode.end)] = rnode
                nnodes.extend([lnode,rnode])
            else:
                self.D[n.start,n.end] = self.path[n.start].size


    def update_path(self,l,d):
        x = self.path.index(l)
        self.D[x,x] = self.D[x,x] + d
        self.minval[x,x] = self.minval[x,x] + d
        self.maxval[x,x] = self.maxval[x,x] +  d
        a,b = x,x
        while self.root.start != a and self.root.end != b and a and b:
            a_p,b_p = self.get_parent_interval((a,b)).get_interval()
            a_s,b_s = self.get_sibling_interval((a,b)).get_interval()
            if b == b_p:
                self.D[a_s,b_s] = self.D[a_s,b_s] + d
                self.minval[a_s,b_s] = self.minval[a_s,b_s] + d
                self.maxval[a_s,b_s] = self.maxval[a_s,b_s] + d
            self.minval[a_p,b_p] = min(self.minval[a_s,b_s],self.minval[a,b])+self.D[a_p,b_p]
            self.maxval[a_p,b_p] = max(self.maxval[a_s,b_s],self.maxval[a,b])+self.D[a_p,b_p]
            a,b = a_p,b_p

    def get_parent_interval(self,interval):
        return self.interval_lookup[interval].parent

    def get_sibling_interval(self, interval):
        n = self.interval_lookup[interval]
        if n.parent == None:
            ValueError("Root node has no parent")
        elif n.start == n.parent.start:
            return n.parent.right_child
        else:
            return n.parent.left_child

    def contains(self, x):
        if x in self.nodes:
            return True
        else:
            return False

    def __str__(self) -> str:
        return ",".join([str(i) for i in self.path])


class HeavyPathDecomposition:
    def __init__(self,tree: ExtendedTree):
        self.paths = []
        self.tree = tree
        visited = set()
        for l in tree.leaf_node_iter():
            path = []
            next_node = l
            while next_node not in visited:
                path.append(next_node)
                visited.add(next_node)
                if next_node.parent_node and next_node.is_heavy():
                    next_node = next_node._parent_node
                else:
                    break
            self.paths.append(PathSearchTree(path[::-1]))


    def get_path(self,x):
        for i in self.paths:
            if i.contains(x):
                return i

    def __str__(self) -> str:
        return "\n".join([str(i) for i in self.paths])


