import dendropy
import numpy as np
from ca_calculator import cluster_affinity_matrix
import argparse
import pandas as pd

def get_tau(n,c):
    l = c
    l[l<n/2] -= 1
    l[l>=n/2] *= -1
    l[l<=-n/2] += n
    return l

def save_matrix(distance_matrix,t1,t2,filepath):
    ar = np.dtype([('c1','object'),('c2','object'),('dist','float')])
    d = np.zeros(distance_matrix.shape,ar)
    for idx,c in enumerate(t1.nodes()):
        c_cluster = ",".join([str(i.taxon.label) for i in c.leaf_nodes()])
        for idx2,x in enumerate(t2.nodes()):
            x_cluster = ",".join([str(i.taxon.label) for i in x.leaf_nodes()])
            d[idx][idx2] = (c_cluster,x_cluster,distance_matrix[c.label][x.label])
    np.savetxt(filepath,d.flatten(),fmt=['{%s}',"{%s}","%f"],delimiter="\t")

def process_trees(args):
    """
    Main function that initializes the trees and computes the tree affinity distance
    :arg: arguments from the system
    """
    t1 = dendropy.Tree.get(path=args.tree_1,schema="newick")
    t2 = dendropy.Tree.get(path=args.tree_2,schema="newick")

    mapping = lambda x:x
    if args.mapping_file:
        mapping = process_mapping_file(args.mapping_file)
    distance_matrix,c1,c2  = cluster_affinity_matrix(t1,t2,mapping)
    distance_matrix = distance_matrix.astype(float)
    dt1_t2 = np.min(distance_matrix,axis=1).astype(float)
    dt2_t1 = np.min(distance_matrix,axis=0).astype(float)
    if args.normalize or args.threshold > 0.0:
        n1 = len(t1.leaf_nodes())
        n2 = len(t2.leaf_nodes())
        tau1 = get_tau(n1,c1)
        tau2 = get_tau(n2,c2)
    if args.normalize:
        ct1_t2 = np.sum(dt1_t2)/np.sum(tau1)
        ct2_t1 = np.sum(dt2_t1)/np.sum(tau2)
        tau1[tau1==0] = 1
        tau2[tau2==0] = 1
        dt1_t2 /= tau1
        dt2_t1 /= tau2
    else:
        ct1_t2 = np.sum(dt1_t2)
        ct2_t1 = np.sum(dt2_t1)
    if args.stat:
        save_matrix(distance_matrix,t1,t2,args.stat)
    print(f"{ct1_t2} {ct2_t1}")
    if args.threshold > 0.0:
        core_clusters_t1_t2 = find_core_clusters(t1,np.argwhere(dt1_t2/tau1 > args.threshold/100),args.min_cluster_size)
        print("\n".join(core_clusters_t1_t2))
    elif args.threshold < 0.0:
        core_clusters_t1_t2 = find_core_clusters(t1,np.argwhere(dt1_t2/tau1 < abs(args.threshold)/100),args.min_cluster_size)
        print("\n".join(core_clusters_t1_t2))
        

    if args.output_tree_1:
        for n in t1.nodes():
            n.edge.annotations.add_new("rel_cost",dt1_t2[n.label])
        t1.write(path=args.output_tree_1,schema=args.schema,suppress_internal_node_labels=False,suppress_annotations=False)
    if args.output_tree_2:
        for n in t2.nodes():
            n.edge.annotations.add_new("rel_cost",dt2_t1[n.label])
        t2.write(path=args.output_tree_2,schema=args.schema,suppress_internal_node_labels=False,suppress_annotations=False)

def process_mapping_file(file_path):
    matching = {}
    with open(file_path,"r") as file:
        for line in file.readlines():
            l,r = [i.strip() for i in line.strip().split(":")]
            matching[l] = r
    return lambda x:matching[x] if x in matching else x


def find_core_clusters(tree,idx,min_size):
    core_nodes = [i.leaf_nodes() for i in tree.find_nodes(lambda x:np.isin(x.label,idx)) if len(i.leaf_nodes()) > min_size]
    core_clusters = []
    for i in core_nodes:
        cluster_string = "{ "
        for l in i:
            cluster_string+=l.taxon.label+","
        cluster_string += " }"
        core_clusters.append(cluster_string)
    return core_clusters


    
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
    parser.add_argument("--normalize",help="If distance should be normalized",action="store_true") 
    parser.add_argument("--threshold",help="The threshold as a percentage of the cost above which a cluster is a core cluster", type=float,default=0)
    parser.add_argument("--min_cluster_size",help="The minimum cluster size for the threshold to count", type=float,default=0)
    parser.add_argument("--mapping_file",help="The mapping file to replace taxa in one tree with another")
    args = parser.parse_args()
    process_trees(args)

if __name__=="__main__":
    main()


