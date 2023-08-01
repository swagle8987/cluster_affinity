from .reportable import Reportable


class RootedAffinityCalculator(Reportable):
    """
    Calculates the rooted cluster affinity cost
    """

    def __init__(self, *args, **kwargs):
        Reportable.__init__(self, *args, **kwargs)

    def calc_affinity(self, tree1, tree2):
        """ Calculates the rooted affinity cost between two rooted trees
        @param tree1
        @param tree2
        :returns: affinity cost from tree1 to tree2
        """
        tree_affinity_cost = 0
        for node in tree1.postorder_node_iter():
            intersection_lookup = {}
            ca_cost, outnodes = self.cluster_affinity(tree1.cluster_lookup[node.label], tree2,intersection_lookup)
            tree_affinity_cost += ca_cost
            node.annotations.add_new("maps_to", " || ".join([str(i) for i in outnodes]))
            node.annotations.add_new("mapping_cost",ca_cost)
            for out in outnodes:
                outnode = tree2.node_lookup[out]
                val = outnode.annotations.get_value("mapped_by","")
                val = val + " || {}".format(node.label)
                outnode.annotations.drop(name="mapped_by")
                outnode.annotations.add_new("mapped_by", val)
                val_num = outnode.annotations.get_value("mapping_frequency","0")
                val = int(val_num) + 1
                outnode.annotations.drop(name="mapping_frequency")
                outnode.annotations.add_new("mapping_frequency", val)
        return tree_affinity_cost

    def cluster_affinity(self, cluster, tree2, intersection_lookup):
        """Calculates the rooted affinity cost between the cluster and a tree
        @param cluster
        @param tree
        :returns: affinity cost from tree1 to tree2
        """
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
            self.log_event("node_computation",
                           {"left_cluster": str(set(cluster)),
                            "right_cluster": str(set(tree2.get_cluster(node.label))),
                            "distance": distance,
                            "intersection": intersection_size})
            intersection_lookup[node.label] = intersection_size
            if distance < min_distance:
                min_distance = distance
                outnodes = [node.label]
            elif distance == min_distance:
                outnodes.append(node.label)
        return min_distance, outnodes
