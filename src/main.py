import argparse

from cluster_computation import rooted_cluster_affinity,calculate_rooted_tau,rooted_cluster_support,calculate_rooted_phi
from ete4 import Tree,nexus
from ete4.smartview.explorer import add_tree,explore

from ete4.smartview import CircleFace, TextFace,Layout,PropFace,BASIC_LAYOUT,HeatmapFace,LegendFace
from utils import peek_line,convert_dict_to_2d_array, make_matrix_image, check_input_trees

import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap,rgb2hex,to_rgb


import os


COLOR_RANGE=[rgb2hex(to_rgb("xkcd:evergreen")),rgb2hex("#ffad01"),rgb2hex("#ff0000")]
cmap = LinearSegmentedColormap.from_list("custom_tree_map",colors=COLOR_RANGE)

def generate_layout(name,color_only):

    def draw_tree(tree):
        yield LegendFace("Range of the cost:\n 0 indicates a perfect match while 1 indicates the maximum cost.\n The relative width of the edges is also an indicator of cost", "continuous",value_range=(0,1),color_range=COLOR_RANGE[::-1],position="aligned",anchor=(0,0))
        yield { 'node-height-min':0,
                'collapsed':{"shape":"outline"}}

    def draw_node(node):
        if not color_only:
            if node.is_leaf:
                yield PropFace("name",fs_min=11,position="right")
            else:
                yield PropFace("c_dist",fmt="%.3f",fs_min=11,position="right")
        color = rgb2hex(cmap(node.props["c_dist"])[:3])
        yield {
                "hz-line": {
                        "stroke":color,
                        "stroke-width":max(1,node.props["c_dist"]*5)
                    },
                "vt-line":{
                    "stroke":color,
                    "stroke-width":max(1,node.props["c_dist"]*5)
                    }
                }
    return Layout(name=name,draw_node=draw_node,draw_tree=draw_tree)


def get_tree(path,ftype):
    if ftype=="newick":
        return Tree(open(path),parser=1)
    elif ftype=="nexus":
        return list(nexus.loads(open(path).read()).values())[0]
    else:
        raise RuntimeError("Incorrect input file format for {}".format(path))

def run_script(cost,args):
    ftype = args.filetype
    if not ftype:
        ftype = "nexus" if peek_line(args.t1) == "#NEXUS" else "newick"
    t1 = get_tree(args.t1,ftype)
    t2 = get_tree(args.t2,ftype)
    if check_input_trees([t1,t2]):
        if cost=="cluster_affinity":
            dist = rooted_cluster_affinity(t1,t2)/calculate_rooted_tau(t1)
            name="Cluster Affinity Layout"
        else:
            dist = rooted_cluster_support(t1,t2)/calculate_rooted_phi(t1)
            name="Cluster Support Layout"
        if not args.cli:
            layout=generate_layout(name,args.color_only)
            add_tree(t1,name="Source tree ({})".format(os.path.basename(args.t1)),layouts=[layout])
            explore(t2,name="Target tree ({})".format(os.path.basename(args.t2)),layouts=[BASIC_LAYOUT])
            print("Press any key to stop the server and finish")
            input()
        else:
            print(dist)


def cluster_matrix():
    parser = argparse.ArgumentParser(
            prog='Cluster Matrix',
            description='Computes the cluster matrix for pairwise costs',
    )
    parser.add_argument('t',help='The treefile containing the multiple trees')
    parser.add_argument('outfile',help='The output image path')
    parser.add_argument('--title',help='The title for the matrix')
    parser.add_argument('--cost',help='The cost to use for the pairwise comparison')
    parser.add_argument('--filetype',help='The input filetype')
    parser.add_argument('--average',help="Adds a column displaying average cost to the matrix",action="store_true")

    args = parser.parse_args()
    ftype = args.filetype
    if not ftype:
        ftype = "nexus" if peek_line(args.t) == "#NEXUS" else "newick"

    if ftype=="newick":
        with open(args.t) as treefile:
            trees = dict()
            for ind,t in enumerate(treefile.read().strip().split("\n")):
                print(ind)
                trees["tree_{}".format(ind+1)] = Tree(t.strip(),parser=1)
    elif ftype=="nexus":
        trees = nexus.loads(open(args.t).read())
    else:
        raise RuntimeWarning("Unsupported file format")

    matrix = []
    for i in trees:
        row = []
        if args.cost == "cluster_support":
            cluster_phi = calculate_rooted_phi(trees[i])
        else:
            cluster_tau = calculate_rooted_tau(trees[i])
        for j in trees:
            if args.cost == "cluster_support":
                row.append(rooted_cluster_support(trees[i],trees[j])/cluster_phi)
            else:
                row.append(rooted_cluster_affinity(trees[i],trees[j])/cluster_tau)
        if args.average:
            row.append(sum(row)/len(row))
        matrix.append(row)
    if args.average:
        xlabels = list(trees.keys()) + ["average"]
    else:
        xlabels=list(trees.keys())
    make_matrix_image(matrix,args.outfile,xlabels=xlabels,ylabels=trees.keys(),cmap=cmap)

    


def cluster_affinity_script():

    parser = argparse.ArgumentParser(
            prog='Cluster Affinity',
            description='Calculates the Asymmetric Cluster Affinity cost from t1 to t2',
    )

    parser.add_argument('t1', help='The source tree from which the cost is to be calculated')
    parser.add_argument('t2', help='The target tree to which the cost is to be calculated')
    parser.add_argument('-t','--filetype', help="The input file format")
    parser.add_argument('--cli',help="Disables interactive browser", action='store_true')
    parser.add_argument('--color_only',help="Disables node annotations(useful when visualizing dense trees)",action='store_true')


    args = parser.parse_args()

    run_script("cluster_affinity",args)


def cluster_support_script():

    parser = argparse.ArgumentParser(
            prog='Cluster Support',
            description='Calculates the Asymmetric Cluster Support cost from t1 to t2',
    )

    parser.add_argument('t1', help='The source tree from which the cost is to be calculated')
    parser.add_argument('t2', help='The target tree to which the cost is to be calculated')
    parser.add_argument('-t','--filetype', help="The input file format")
    parser.add_argument('--cli',help="Disables interactive browser", action='store_true')
    parser.add_argument('--color_only',help="Disables node annotations(useful when visualizing dense trees)",action='store_true')

    args = parser.parse_args()

    run_script("cluster_support",args)

if __name__=="__main__":
    cluster_affinity_script()
