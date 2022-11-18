import dendropy
from affinity import tree_affinity
from extendedtree import ExtendedTree
import argparse

def main(args):
    """
    Main function that initializes the trees and computes the tree affinity distance
    :arg: arguments from the system
    """
    t1 = ExtendedTree.get(path=args.tree_1,schema="newick")
    t2 = ExtendedTree.get(path=args.tree_2,schema="newick")
    dist = tree_affinity(t1,t2)
    print(dist)
    t1.writeTreeToFile(args.output_tree_1,schema="newick")
    t2.writeTreeToFile(args.output_tree_2,schema="newick")
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

    args = parser.parse_args()
    main(args)


