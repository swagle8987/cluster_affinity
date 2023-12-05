import dendropy
import numpy as np
from ca_calculator import cluster_affinity_matrix
import argparse

def process_trees(args):
    """
    Main function that initializes the trees and computes the tree affinity distance
    :arg: arguments from the system
    """
    t1 = dendropy.Tree.get(path=args.tree_1,schema="newick")
    t2 = dendropy.Tree.get(path=args.tree_2,schema="newick")

    mapping = dict()
    if None and args.match_file:
        with open(args.match_file,"r") as match_file:
            for line in match_file.readlines():
                match = [i.trim() for i in line.trim().split(":")]
                mapping[match[0]] = match[1]

    distance_matrix,c1,c2  = cluster_affinity_matrix(t1,t2)



    dt1_t2 = np.min(distance_matrix,axis=0).astype(float)
    dt2_t1 = np.min(distance_matrix,axis=1).astype(float)
    if args.normalize or args.threshold > 0.0:
        n1 = len(t1.leaf_nodes())
        n2 = len(t1.leaf_nodes())
        c1[c1<n1/2] -= 1
        c1[c1>=n1/2] *= -1
        c1[c1<=-n1/2] += n1
        c1[c1==0] += 1
        c2[c2<n2/2] -= 1
        c2[c2>=n2/2] *= -1
        c2[c2<=-n2/2] += n2
        c2[c2==0] += 1
    if args.normalize:
        ct1_t2 = np.sum(dt1_t2/c1)
        ct2_t1 = np.sum(dt2_t1/c2)
    else:
        ct1_t2 = np.sum(dt1_t2)
        ct2_t1 = np.sum(dt2_t1)
    print(f"Distance from t1 to t2 {ct1_t2}")
    print(f"Distance from t2 to t1 {ct2_t1}")
    for i in t1.leaf_nodes():
        if np.min(distance_matrix[i.label])<dt1_t2[i.label]:
            print(i)
    if args.threshold > 0.0:
        core_clusters_t1_t2 = find_core_clusters(t1,np.argwhere(dt1_t2/c1 > max(1,args.threshold/100)))
        core_clusters_t2_t1 = find_core_clusters(t2,np.argwhere(dt2_t1/c2 > max(1,args.threshold/100)))
        print(f"High performing clusters for t1{core_clusters_t1_t2}")
        print(f"High performing clusters for t1{core_clusters_t2_t1}")

    if args.output_tree_1:
        annotate_nodes(t1,lambda x:True,"cost",lambda x:dt1_t2[x.label])
        t1.write(path=args.output_tree_1,schema=args.schema,suppress_internal_node_labels=False,suppress_annotations=False)
    if args.output_tree_2:
        annotate_nodes(t2,lambda x:True,"cost",lambda x:dt2_t1[x.label])
        t2.write(path=args.output_tree_2,schema=args.schema,suppress_internal_node_labels=False,suppress_annotations=False)



def find_core_clusters(tree,idx):
    core_nodes = [i.leaf_nodes() for i in tree.find_nodes(lambda x:np.isin(x.label,idx))]
    core_clusters = []
    for i in core_nodes:
        cluster_string = "{ "
        for l in i:
            cluster_string+=l.taxon.label+","
        cluster_string += " }"
        core_clusters.append(cluster_string)
    return core_clusters


def annotate_nodes(tree,filter_func,annotate_label,annotate_func):
    nodes = tree.find_nodes(filter_func)
    for n in nodes:
        n.annotations.add_new(annotate_label,annotate_func(n))
    
def main():
    program_description = """
    Computes the tree affinity distance between two trees t1 and t2 and outputs the annotated resultant trees.
    Note that the tree affinity distance is asymmetric and is computed from tree 1 to tree 2 
    The resultant tree 1 (left tree) has every node annotated with the distance to tree 2 and the node/s that it maps to in tree2.
    The resultant tree 2 (right tree) has every node annotated with a support (sum of the intersections that map to the node divided by the size of the cluster located at the node) and a list of all nodes that map to it
    """
    parser = argparse.ArgumentParser(description=program_description,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('tree_1', metavar="t1", help="The tree from which the distance is calculated")
    parser.add_argument('tree_2', metavar="t2", help="The tree to which the distance is calculated")
    parser.add_argument('--output_tree_1','-o1', metavar="o1",help="The output file for the resultant tree 1, defaults to tree1_result.tre",default="")
    parser.add_argument('--output_tree_2','-o2',metavar="o2", help="The output file for the resultant tree 2, defaults to tree2_result.tre",default="")
    parser.add_argument("--stat","-s", help="The output file for distributions for each cluster, defaults to none", default="")
    parser.add_argument("--schema",help="The output file format",default="newick") 
    parser.add_argument("--relative",help="If cluster affinity support should be used instead",action="store_true") 
    parser.add_argument("--normalize",help="If distance should be normalized",action="store_true") 
    parser.add_argument("--threshold",help="The threshold as a percentage of the cost above which a cluster is a core cluster", type=float,default=0)

    args = parser.parse_args()
    process_trees(args)

if __name__=="__main__":
    main()


