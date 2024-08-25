# General issues
- Code lacks clarity and does too many things. Need defining of a specific feature list and target audience
- Code lacks comments

## main.py
Should be the logic that handles i/o and runs the code.
Should have parser logic as well
* contains matrix i/o
* contains other logic like finding core clusters and the like
* contains tau calculation
* contains cost calculation from trees


## rooted_affinity.py
Should contain only the logic for rooted affinity calculation
The process itself is three steps: 
1. Decompose trees into list of nodes
2. For each node in t1:
    2.1. For each node in t2:
        Calculate the set difference for t1
3. Return sum of minimum values
All can be within same function ideally
However, the idea of reporting all calculations makes things slightly different
Hence, here is the proposal
1. Convert trees into Reportable trees with each node acting as a publisher
2. For each tree calculation, the publisher emits a calculation event. In this case, each node emits an event when a cluster is compared to it
3. The sum and minimum do not change. 

## reportable.py
Technically could be it's own library. Expose a modified version of Dendropy tree that lets functions be registered as emit events and keeps a subscriber log on hand to save



## Needed files/modules:
1. data_extract -- to handle extracting the matrix from the reportable logger

