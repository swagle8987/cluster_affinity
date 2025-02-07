import argparse

from cluster_computation import unrooted_cluster_affinity
from cluster_computation import rooted_cluster_affinity,calculate_rooted_tau
import dendropy

from utils import peek_line,convert_dict_to_2d_array, make_matrix_image, check_input_trees

import os

def cluster_matrix_script():
    parser = argparse.ArgumentParser(
            prog='Cluster Matrix',
            description='Calculates the Asymmetric Cluster Affinity matrix for a list of trees',
    )

    parser.add_argument('t', help='The treefile containing a list of trees')
    parser.add_argument('--filetype', help="The input file format")
    parser.add_argument('--outfile', help="The output file",default="matrix.png")

    args = parser.parse_args()

    ftype = args.filetype
    if not ftype:
        ftype = "nexus" if peek_line(args.t) == "#NEXUS" else "newick"

    tlist = dendropy.TreeList.get(path=args.t,schema=ftype,rooting="default-rooted")
    if check_input_trees(tlist):
        matrix_encoding = dict()
        for indx,i in enumerate(tlist):
            matrix_encoding[indx] = dict()
            tau = calculate_rooted_tau(i)
            for jndx,j in enumerate(tlist):
                matrix_encoding[indx][jndx] = rooted_cluster_affinity(i,j)/tau
        matrix = convert_dict_to_2d_array(matrix_encoding)
        make_matrix_image(matrix,output_path=args.outfile)
        print(matrix)

def cluster_affinity_script():

    parser = argparse.ArgumentParser(
            prog='Cluster Affinity',
            description='Calculates the Asymmetric Cluster Affinity cost from t1 to t2',
    )

    parser.add_argument('t1', help='The source tree from which the cost is to be calculated')
    parser.add_argument('t2', help='The target tree to which is to be calculated')
    parser.add_argument('-t','--filetype', help="The input file format")


    args = parser.parse_args()

    ftype = args.filetype
    if not ftype:
        ftype = "nexus" if peek_line(args.t1) == "#NEXUS" else "newick"

    tns = dendropy.TaxonNamespace(label="taxa")
    t1 = dendropy.Tree.get(path=args.t1,taxon_namespace=tns,schema=ftype,rooting="default-rooted")
    t2 = dendropy.Tree.get(path=args.t2,taxon_namespace=tns,schema=ftype,rooting="default-rooted")

    dist = -1
    if check_input_trees([t1,t2]):
            dist = rooted_cluster_affinity(t1,t2)/calculate_rooted_tau(t1)
    print(dist)


if __name__=="__main__":
    cluster_affinity_script()
