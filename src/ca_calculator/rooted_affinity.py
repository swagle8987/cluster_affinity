import numpy as np
import dendropy
import pdb


def cluster_affinity_matrix(tree1, tree2):
    n1 = len(tree1.leaf_nodes())
    n2 = len(tree2.leaf_nodes())
    c1 = np.full(len(tree1.nodes()),0)
    c2 = np.full(len(tree2.nodes()),0)
    h = np.full((len(c1),len(c2)),133223)
    for i,v in enumerate(tree1.postorder_node_iter()):
        v.label = i
        if v.is_leaf():
            c1[i] = 1
        else:
            for ck in v.child_node_iter():
                c1[i] += c1[ck.label]
    for i,v in enumerate(tree2.postorder_node_iter()):
        v.label = i
        if v.is_leaf():
            c2[i] = 1
        else:
            for ck in v.child_node_iter():
                c2[i] += c2[ck.label]
    for c in tree1.leaf_node_iter():
        for x in tree2.leaf_node_iter():
            if c.taxon.label == x.taxon.label:
                h[c.label][x.label] = 0
            else:
                h[c.label][x.label] = 2
    for c in tree1.leaf_node_iter():
        for x in tree2.postorder_internal_node_iter():
            for xk in x.child_node_iter():
                h[c.label][x.label] = 1
                if h[c.label][xk.label] < c2[xk.label]:
                    h[c.label][x.label] = -1
                    break
            h[c.label][x.label] += c2[x.label]
    for c in tree1.postorder_internal_node_iter():
        for x in tree2.leaf_node_iter():
            for ck in c.child_node_iter():
                h[c.label][x.label] =1
                if h[ck.label][x.label] < c1[ck.label]:
                    h[c.label][x.label] = -1
                    break
            h[c.label][x.label] += c1[c.label]
    for c in tree1.postorder_internal_node_iter():
        for x in tree2.postorder_internal_node_iter():
            h[c.label][x.label] =0
            for ck in c.child_node_iter():
                h[c.label][x.label] += h[ck.label][x.label] - c2[x.label]
            h[c.label][x.label] += c2[x.label]
    return h,c1,c2


if __name__=="__main__":
    t1 = dendropy.Tree.get(path="1.tre",schema="newick")
    t2 = dendropy.Tree.get(path="2.tre",schema="newick")
    h = cluster_affinity_matrix(t1,t2)
    print(h)
    print(np.min(h,axis=0))
    print(np.sum(np.min(h,axis=0)))
    print(np.min(h,axis=1))
    print(np.sum(np.min(h,axis=1)))

