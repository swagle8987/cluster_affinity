The Asymmetric Cluster Affinity cost is a phylogenetic cost based on calculating the symmetric difference between the cluster representations of trees. Currently the CLI tool supports calculating the cluster affinity distance from the source tree to the target tree.


### Installation
Cluster Affinity is available in PyPi and can be installed as ``pip install cluster-affinity``. Note that the package is built for Python 3.10 or higher.


### Tutorial
---
Currently the CLI tool supports comparing two trees and outputting the cluster affinity cost and cluster support cost.

The following command computes the cluster affinity cost between two trees "t1.tre" and "t2.tre" and returns the normalized cluster affinity cost. The cost is normalized on a scale of 0-1 where 1 is the maximum possible cluster affinity cost for t1. 
``
cluster_affinity t1 t2 --filetype input_filetype
``
where t1 and t2 are paths to the trees.

There is a similar command for the cluster support cost as well 
``
cluster_support t1 t2 --filetype input_filetype
``
where t1 and t2 are paths to the trees.

There is also a command that generates a matrix showing the pairwise comparison between the trees using the cluster affinity cost.

``
cluster_matrix t --filetype --outfile matrix.png
``
where t is the path to the treefile containing a list of input trees.

