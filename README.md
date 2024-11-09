The Asymmetric Cluster Affinity cost is a phylogenetic cost based on calculating the symmetric difference between the cluster representations of trees. Currently the CLI tool supports calculating the cluster affinity distance from the source tree to the target tree.


### Installation
Cluster Affinity is available in PyPi and can be installed as pip install cluster_affinity. Note that the package is built for Python 3.10 or higher. Cluster Affinity depends on dendropy, numpy and pytest. 


### Tutorial
---
Currently the CLI tool supports comparing two trees and outputting the cluster affinity cost. The CLI command for the same is
``
cluster_affinity t1 t2
``
where t1 and t2 are paths to newick representations of the trees.
