def tree_affinity(tree1,tree2):
    """TODO: Docstring for tree_affinity.
    :returns: TODO

    """
    tree_affinity_cost = 0
    for node in tree1.postorder_node_iter():
        ca_cost,outnodes = cluster_affinity(tree1.cluster_lookup[node.label],tree2)
        tree_affinity_cost += ca_cost
        node.annotations.add_new("out_nodes"," || ".join([str(i) for  i in outnodes]))
        node.annotations.add_new("out_cost",ca_cost)
        for out in outnodes:
            outnode = tree2.node_lookup[out]
            val = outnode.annotations.get_value("in_nodes","")
            val = val + " || {}".format(node.label)
            outnode.annotations.drop(name="in_nodes")
            outnode.annotations.add_new("in_nodes",val)
    return tree_affinity_cost

def cluster_affinity(cluster,tree2):
    """TODO: Docstring for cluster_affinity.
    :returns: TODO

    """
    intersection_lookup= {}
    min_distance = 10000000000000000
    outnodes = []
    for node in tree2.postorder_node_iter():
        intersection_size = 0
        if node.is_leaf():
            if str(node.taxon) in cluster:
                intersection_size = 1
            else:
                intersection_size = 0
        else:
            for child in node.child_node_iter():
                intersection_size += intersection_lookup[child.label]
        distance = len(cluster) + len(tree2.get_cluster(node.label)) - (2 * intersection_size)
        intersection_lookup[node.label] = intersection_size
        if distance < min_distance:
            min_distance = distance
            outnodes = [node.label]
        elif distance == min_distance:
            outnodes.append(node.label)
    return min_distance,outnodes
                


