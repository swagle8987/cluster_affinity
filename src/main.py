import dendropy
from utr_affinity import UtrAffinityCalculator
from rtu_affinity import RtuAffinityCalculator
from rooted_affinity import RootedAffinityCalculator
from extendedtree import ExtendedTree
import argparse

def main(args):
    """
    Main function that initializes the trees and computes the tree affinity distance
    :arg: arguments from the system
    """
    t1 = ExtendedTree.get(path=args.tree_1,schema="newick")
    t2 = ExtendedTree.get(path=args.tree_2,schema="newick")
    if args.stat:
        logging_level = 1
    else:
        logging_level = 0 
    if t1.is_rooted and t2.is_rooted:
        affinity_calculator = RootedAffinityCalculator(logging_level=logging_level)
    elif not t1.is_rooted and t2.is_rooted:
        affinity_calculator = UtrAffinityCalculator(logging_level=logging_level)
    elif t1.is_rooted and not t2.is_rooted:
        affinity_calculator = RtuAffinityCalculator(logging_level=logging_level)
    else:
        raise TypeError("No tree is rooted")
    
    if (not t1.is_rooted or not t2.is_rooted):
        dist,edge = affinity_calculator.calc_affinity(t1,t2)
        t1.reroot_at_edge(edge)
    else:
        dist = affinity_calculator.calc_affinity(t1,t2)
    print(dist)
    t1.writeTreeToFile(args.output_tree_1,schema="newick")
    t2.writeTreeToFile(args.output_tree_2,schema="newick")
    if args.stat:
        df = affinity_calculator.get_event_dataframe("node_computation")
        df.to_csv(args.stat)
    return

if __name__=="__main__":
    program_description = """
    Computes the tree affinity distance between two trees t1 and t2 and outputs the annotated resultant trees.
    Note that the tree affinity distance is asymmetric and is computed from tree 1 to tree 2 
    The resultant tree 1 (left tree) has every node annotated with the distance to tree 2 and the node/s that it maps to in tree2.
    The resultant tree 2 (right tree) has every node annotated with a support (sum of the intersections that map to the node divided by the size of the cluster located at the node) and a list of all nodes that map to it
    """
    parser = argparse.ArgumentParser(description=program_description,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('tree_1', metavar="t1", help="The tree from which the distance is calculated")
    parser.add_argument('tree_2', metavar="t2", help="The tree to which the distance is calculated")
    parser.add_argument('--output_tree_1','-o1', help="The output file for the resultant tree 1, defaults to tree1_result.tre",default="tree1_result.tre")
    parser.add_argument('--output_tree_2','-o2', help="The output file for the resultant tree 2, defaults to tree2_result.tre",default="tree2_result.tre")
    parser.add_argument("--stat","-s", help="The output file for distributions for each cluster, defaults to none", default="")

    args = parser.parse_args()
    main(args)


