from extendedtree import ExtendedTree
from ca_calculator import RootedAffinityCalculator
import argparse
import random
import pdb

def main(args):
    t1 = ExtendedTree.get(path=args.tree_1,schema="newick")
    orig_tree = t1.clone(1)
    for i in range(args.nops):
            c = RootedAffinityCalculator(logging_level=1)
            if args.operation == "spr":
                n = random.choice(t1.nodes(filter_fn=lambda x:True if x.parent_node else False))
                subtree = [i for i in n.postorder_iter()]
                k = random.choice(t1.nodes(filter_fn=lambda x:True if x not in subtree and not x.is_leaf() else False))
            elif args.operation=="nni":
                n = random.choice(t1.nodes(filter_fn=lambda x:True if x.parent_node.parent_node else False))
                k = n.parent_node.parent_node
            t1.spr(n, k)
            dist = c.calc_affinity(orig_tree,t1,True)
            orig_tree.reset_tree()
            orig_tree.enrich_tree()
            print(args.tree_1,i,dist)

     


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="spr test suite",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('tree_1',metavar='t1',help="starting tree")
    parser.add_argument('operation',metavar='opr',help="Operation to be performed")
    parser.add_argument('nops',metavar='nops',help="number of operations to be done",type=int)
    args = parser.parse_args()
    main(args)
