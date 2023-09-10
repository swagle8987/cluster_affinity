from .reportable import Reportable
import pdb


class RootedAffinityCalculator(Reportable):
    """
    Calculates the rooted cluster affinity cost
    """

    def __init__(self, *args, **kwargs):
        Reportable.__init__(self, *args, **kwargs)

    def calc_affinity(self, tree1, tree2, relative=False):
        """ Calculates the rooted affinity cost between two rooted trees
        @param tree1
        @param tree2
        :returns: affinity cost from tree1 to tree2
        """
        tree_affinity_cost = 0
        node_lookups = dict()
        for node in tree1.postorder_node_iter():
            intersection_lookup = {}
            ca_cost, outnodes = self.cluster_affinity(tree1.cluster_lookup[node.label], tree2,intersection_lookup, relative)
            tree_affinity_cost += ca_cost
            node.annotations.add_new("maps_to", " || ".join([str(i.label) for i in outnodes]))
            node.annotations.add_new("relative_mapping_cost",ca_cost)
            node_lookups[node] = set(outnodes)
            for out in outnodes:
                outnode = tree2.node_lookup[out.label]
                val = outnode.annotations.get_value("mapped_by","")
                val = val + " || {}".format(node.label)
                outnode.annotations.drop(name="mapped_by")
                outnode.annotations.add_new("mapped_by", val)
                val_num = outnode.annotations.get_value("mapping_frequency","0")
                val = int(val_num) + 1
                outnode.annotations.drop(name="mapping_frequency")
                outnode.annotations.add_new("mapping_frequency", val)
        def is_node_antichain(n):
            if n.parent_node == None:
                return False
            parent_nodes_antichains = {i.parent_node for i in node_lookups[n]}
            if len(parent_nodes_antichains & node_lookups[n.parent_node]) > 0:
                return False
            else:
                return True

        tree1.color_nodes(is_node_antichain)
        return tree_affinity_cost


    def cluster_affinity(self, cluster, tree2, intersection_lookup, relative):
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
            if distance<0:
                pdb.set_trace()
            relative_distance = distance/len(cluster) if relative else distance
            self.log_event("node_computation",
                           {"left_cluster": str(set(cluster)),
                            "right_cluster": str(set(tree2.get_cluster(node.label))),
                            "distance": distance,
                            "relative_distance":relative_distance,
                            "intersection": intersection_size})
            intersection_lookup[node.label] = intersection_size
            if relative_distance< min_distance:
                min_distance = relative_distance
                outnodes = [node]
            elif relative_distance== min_distance:
                outnodes.append(node)
        return min_distance, outnodes
