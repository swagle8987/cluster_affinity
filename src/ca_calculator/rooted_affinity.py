import numpy as np
import dendropy
import pdb

### Trying to speed things up
### First map each leaf to another leaf or set leaf to 1
### For a node v with x1,x2 as children
### c_v = c1 n c2
### c_v n y = c1 n y u c2 n y
### 

### for v with some n children x1,x2...xn
### mappings contains the value |v| - 2|v n c| for all v
### mappings[v] = |v| - 2|v n c|

## Make it simpler
def cluster_affinity_cluster(cluster,tree,mapping= lambda x:x):
    mappings = np.full(len(tree.nodes()),0) ## Mappings of the cluster to each node
    for i,v in enumerate(tree.postorder_node_iter()):
        v.label = i 
        if v.is_leaf():
            if mapping(v.taxon.label) in cluster:
                mappings[i] = -1
            else:
                mappings[i] = 1
        else:
            for child_v in v.child_node_iter():
                mappings[i] += mappings[child_v.label]
    dist =  len(cluster) + np.min(mappings) ## Add the cluster size missing from the value
    return dist


def cluster_affinity(t1,t2,mapping = lambda x:x):
    dist = 0
    cluster_lookup = dict()
    for i,v in enumerate(t1.postorder_node_iter()):
        if v.is_leaf():
            cluster_lookup[v] = {v.taxon.label}
            if t2.taxon_namespace.has_taxon_label(v.taxon.label):
                pass
            else:
                dist += 1
        else:
            c_v = set()
            for child_v in v.child_node_iter():
                c_v = c_v | cluster_lookup[child_v]
            dist += cluster_affinity_cluster(c_v,t2,mapping)
            cluster_lookup[v] = c_v
    return dist


if __name__=="__main__":
    t1 = dendropy.Tree.get(path="1.tre",schema="newick")
    t2 = dendropy.Tree.get(path="2.tre",schema="newick")
    h = cluster_affinity_matrix(t1,t2)
    print(h)
    print(np.min(h,axis=0))
    print(np.sum(np.min(h,axis=0)))
    print(np.min(h,axis=1))
    print(np.sum(np.min(h,axis=1)))

