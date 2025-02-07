from heavy_path import *
import dendropy
from extendedtree import ExtendedTree
from cluster_affinity import rooted_cluster_affinity

t = ExtendedTree.get(data="((A,B),((C,D),(F,E)));", schema="newick",rooting="default-rooted")
t2 = ExtendedTree.get(data="((A,D),((C,B),(F,E)));", schema="newick",rooting="default-rooted")
print(rooted_cluster_affinity(t,t))
print(rooted_cluster_affinity(t,t2)) ## d((A,B)) + d((C,D)) + d((C,D,E,F)) = 1 + 1 + 2 = 4
