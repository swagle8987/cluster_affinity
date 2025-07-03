#!/usr/bin/env python3
#
import ete4
from .heavy_path_decompositions import *
from alive_progress import alive_bar
from sys import stderr


def update_path(x, d):
    path_interval = x.get_prop("path_tree")
    path_interval.D += d
    path_interval.minval += d
    path_interval.maxval += d
    while path_interval.parent:
        parent_interval = path_interval.parent
        sibling_interval = path_interval.get_sibling()
        if parent_interval.get_last() == path_interval.get_last():
            sibling_interval.D += d
            sibling_interval.minval += d
            sibling_interval.maxval += d
        parent_interval.minval = (
            min(path_interval.minval, sibling_interval.minval) + parent_interval.D
        )
        parent_interval.maxval = (
            max(path_interval.maxval, sibling_interval.maxval) + parent_interval.D
        )
        path_interval = parent_interval
    return path_interval.get_first(), path_interval.get_last()


def add_leaf_general(x, counter):
    first, last = update_path(x, -2)
    while not first.is_root:
        x = first.parent
        first, last = update_path(x, -2)
    counter = counter + 1
    return counter


def remove_leaf_general(x, counter):
    first, last = update_path(x, 2)
    while not first.is_root:
        x = first.parent
        first, last = update_path(x, 2)
    counter = counter - 1
    return counter


def construct_leaf_map(t):
    map = dict()
    for i in t:
        map[i.name] = i
    return map


def compute_transfer_index(t1, t2, cost, annotate_cost=True, disable_bar=False):
    annotate_heavy_nodes(t1)
    paths = get_maximal_heavy_paths(t2)
    path_search_trees = [PathTree(i) for i in paths]
    lmap = construct_leaf_map(t2)
    visited = set()
    counter = 0
    distances = {
        "Rooted Cluster Affinity": 0,
        "Unrooted Cluster Affinity": 0,
        "Rooted Cluster Support": 0,
        "Unrooted Cluster Support": 0,
    }
    if cost not in distances:
        raise RuntimeError("invalid cost provided")
    n = len(t1)
    with alive_bar(
        len(t1),
        bar="bubbles",
        spinner="radioactive",
        file=stderr,
        disable=True,
    ) as bar:
        for i in t1:
            curNode = i
            x = lmap[i.name]
            counter = add_leaf_general(x, counter)
            while not curNode.is_root:
                if (
                    len(curNode) * 2 >= (len(curNode.parent))
                    and curNode.parent not in visited
                ):
                    siblingNode = curNode.get_sisters()[0]
                    for l in siblingNode:
                        counter = add_leaf_general(lmap[l.name], counter)
                    min1 = min(i.minval for i in path_search_trees) + counter
                    min2 = n - max(i.maxval for i in path_search_trees) - counter
                    distances["Rooted Cluster Affinity"] += min1
                    distances["Rooted Cluster Support"] += min1 / len(curNode)
                    distances["Unrooted Cluster Affinity"] += min(min1, min2)
                    distances["Unrooted Cluster Support"] += min(min1, min2) / len(
                        curNode
                    )
                    annotate_node_with_cost(curNode, cost, min1, min2, len(t1))
                    visited.add(curNode.parent)
                    curNode = curNode.parent
                else:
                    for l in curNode:
                        counter = remove_leaf_general(lmap[l.name], counter)
                    break
            if curNode.is_root:
                for l in curNode:
                    counter = remove_leaf_general(lmap[l.name], counter)
    return distances[cost]


def annotate_node_with_cost(v, label, min1, min2, tree_size):
    if label == "Rooted Cluster Affinity" or label == "Rooted Cluster Support":
        v.add_prop(label, (min1 / max(1, min(len(v) - 1, tree_size - len(v)))))
    elif label == "Unrooted Cluster Affinity" or label == "Unrooted Cluster Support":
        v.add_prop(
            label, (min(min1, min2) / max(1, min(len(v) - 1, tree_size - len(v) - 1)))
        )
