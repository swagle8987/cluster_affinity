#!/usr/bin/env python3
#
#
import importlib.metadata
import argparse

from .reader import get_tree
from .cluster_computation import (
    compute_transfer_index,
    calculate_rooted_phi,
    calculate_rooted_tau,
    calculate_unrooted_phi,
    calculate_unrooted_tau,
    rooted_cluster_affinity,
    rooted_cluster_support,
    unrooted_cluster_affinity,
)

from .visualization import start_web_server
import ete4
import pathlib


def check_input_trees(tlist) -> bool:
    listtaxa = set(tlist[0].leaf_names())
    for i in tlist:
        for l in i.leaf_names():
            if l not in listtaxa:
                return False
    return True


def check_strictly_binary(t)->bool:
    for i in t.traverse():
        if len(i.children)>2 or len(i.children)==1:
            return False
    return True


def get_default_args(name, description):
    parser = argparse.ArgumentParser(
        prog=name,
        description=description,
    )

    parser.add_argument(
        "t1",
        help="The source tree from which the cost is to be calculated",
        type=pathlib.Path,
    )
    parser.add_argument(
        "t2",
        help="The target tree to which the cost is to be calculated",
        type=pathlib.Path,
    )
    parser.add_argument("--filetype", help="The input file format")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--cli", help="Disables interactive browser", action="store_true"
    )
    group.add_argument(
        "--color_only",
        help="Disables node annotations(useful when visualizing dense trees)",
        action="store_true",
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Version number",
        action="version",
        version=importlib.metadata.version("cluster_affinity"),
    )
    parser.add_argument(
            "--unrooted",
            action="store_true"
            )
    return parser


def cluster_affinity():
    cluster_cost_script(support=False)


def cluster_support():
    cluster_cost_script(support=True)

def get_dist(cost,t1,t2,t1_is_rooted,t2_is_rooted,**kwargs):
    dist = -1
    if (not ("cli" in kwargs and kwargs["cli"])) or not (check_strictly_binary(t1) and check_strictly_binary(t2)):
        if t1_is_rooted and t2_is_rooted:
            cost = "Rooted "+cost
            if "Support" in cost:
                dist = rooted_cluster_support(t1,t2)/calculate_rooted_phi(t1)
            else:
                dist = rooted_cluster_affinity(t1,t2)/calculate_rooted_tau(t1)
        elif not t1_is_rooted:
            cost = "Unrooted "+cost
            if "Support" in cost:
                raise RuntimeError("Unrooted Cluster Support is not supported")
            else:
                dist = unrooted_cluster_affinity(t1,t2)/calculate_unrooted_tau(t1)
        else:
            raise RuntimeError("Rooted to unrooted comparisons are not supported")
    else:
        if t1_is_rooted and t2_is_rooted:
            cost = "Rooted " + cost
            dist = compute_transfer_index(
                t1, t2, cost=cost, annotate_cost=False
            ) / calculate_rooted_tau(t1)
        elif not t1_is_rooted:
            cost = "Unrooted " + cost
            if "Support" in cost:
                raise RuntimeError("Unrooted Cluster Support is not supported")
            dist = compute_transfer_index(
                t1, t2, cost=cost, annotate_cost=False
            ) / calculate_unrooted_tau(t1)
        else:
            raise RuntimeError("Rooted to unrooted comparisons are not supported")
    return dist

def cluster_cost_script(support=False):
    if support:
        cost = "Cluster Support"
    else:
        cost = "Cluster Affinity"
    parser = get_default_args(
        name=cost,
        description="Calculates the Asymmetric {} cost from t1 to t2".format(cost),
    )
    args = parser.parse_args()
    t1, t1_is_rooted = get_tree(args.t1, args.filetype)
    t2, t2_is_rooted = get_tree(args.t2, args.filetype)
    if args.unrooted:
        t1_is_rooted = False
        t2_is_rooted = False
    if not t1_is_rooted and len(t1.children)> 2:
        newnode = ete4.Tree()
        newnode.add_child(child=t1)
        newnode.add_child(child=t1.children[0].detach())
        t1 = newnode
    if not t2_is_rooted and len(t2.children)> 2:
        newnode = ete4.Tree()
        newnode.add_child(child=t2)
        newnode.add_child(child=t2.children[0].detach())
        t2 = newnode
    if check_input_trees([t1, t2]):
        dist = get_dist(cost,t1,t2,t1_is_rooted,t2_is_rooted,cli=args.cli)
    else:
        raise RuntimeError(
            "Tree {} and {} have different taxa set".format(args.t1, args.t2)
        )
    if not args.cli:
        start_web_server(t1, t2, cost, args.color_only, args.t1.name, args.t2.name)
    else:
        print(dist)
