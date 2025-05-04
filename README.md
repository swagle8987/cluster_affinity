The Asymmetric Cluster Affinity cost is a phylogenetic cost based on calculating the symmetric difference between the cluster representations of trees. Currently the CLI tool supports calculating the cluster affinity distance from the source tree to the target tree.


### Installation
Cluster Affinity is available in PyPi and can be installed as pip install cluster-affinity. Note that the package is built for Python 3.10 or higher. Cluster Affinity depends on dendropy, numpy and pytest. 


### Tutorial
---
Currently the CLI tool supports comparing two trees and outputting the cluster affinity cost. The CLI command for the same is
``
cluster_affinity t1 t2 --filetype input_filetype
``
where t1 and t2 are paths to the trees.

There is also a command to generate a matrix representation for the trees:

``
cluster_matrix t --filetype --outfile matrix.png
``
where t is the path to the treefile containing a list of input trees.

