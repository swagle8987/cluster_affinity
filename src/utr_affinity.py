from reportable import Reportable
from pprint import pprint

class UtrAffinityCalculator(Reportable):
    """
    Calculates the unrooted cluster affinity cost
    """

    def __init__(self, *args, **kwargs):
        Reportable.__init__(self, *args, **kwargs)

    def calc_affinity(self, t1, t2):
        t2.generate_clusters()
        tree_affinity_cost = 0
        tree1 = TreeInfo(t1,t2)
        clusters = []
        for u in tree1.dpmap:
            for v in tree1.dpmap[u]:
                clusters.append(tree1.dpmap[u][v])
        sortedc = sorted(clusters, key=lambda x:len(x.cluster))
        for c in sortedc:
            d=0
            if len(c.cluster)>1:
                intersection_lookup = {}
                d,on = self.cluster_affinity(c.cluster, t2, intersection_lookup)
                for k in c.sub_clusters:
                    d += k.distance
                c.distance = d
        mindist = 10000000000000000
        edge = ()
        for u in tree1.dpmap:
            for v in tree1.dpmap[u]:
                dist = tree1.dpmap[u][v].distance + tree1.dpmap[v][u].distance
                if dist < mindist:
                    mindist = dist
                    edge =[i for i in u.get_incident_edges() if i in v.get_incident_edges()] 

        return mindist, edge[0]

    def tree_affinity(self, tree1, tree2):
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

class TreeInfo():
    def __init__(self, t1, t2):
        self.pmap = dict()
        self.dpmap = dict()
        self.emap = dict()
        self.tree = t1
        self.leaf_set = set([str(i) for i in t1.poll_taxa()])
        for i in self.tree.postorder_node_iter():
            self.pmap[i] = dict()
            self.dpmap[i] = dict()
            for v in i.adjacent_nodes():
                self.dpmap[i][v] = ClusterInfo()
        for l in self.tree.leaf_nodes():
            v = l.parent_node
            self.pmap[v][l] = frozenset({str(l.taxon)})
            self.pmap[l][v] = frozenset(self.leaf_set - self.pmap[v][l])
            self.dpmap[v][l].cluster = self.pmap[v][l]
            self.dpmap[v][l].distance = 0
        self.calculate_clusters()
        self.clusters = dict()

    def calculate_clusters(self):
        l = self.tree.leaf_nodes()[0]
        self.tree.reroot_at_node(l,suppress_unifurcations=False)
        for v in self.tree.postorder_node_iter():
            if v.is_internal() and v.parent_node:
                u,w = v.child_nodes()
                x = v.parent_node
                self.pmap[x][v] = frozenset(self.pmap[v][u] | self.pmap[v][w])
                self.pmap[v][x] = frozenset(self.leaf_set - self.pmap[x][v])
                self.dpmap[x][v].sub_clusters.extend([self.dpmap[v][u],self.dpmap[v][w]])
                self.dpmap[x][v].cluster = self.pmap[v][u] | self.pmap[v][w]
                self.dpmap[w][v].sub_clusters.extend([self.dpmap[v][x],self.dpmap[v][u]])
                self.dpmap[w][v].cluster = self.pmap[v][u] | self.pmap[v][x]
                self.dpmap[u][v].sub_clusters.extend([self.dpmap[v][w],self.dpmap[v][x]])
                self.dpmap[u][v].cluster = self.pmap[v][w] | self.pmap[v][x]
    


    def cluster_affinity(self, cluster, tree2):
        """Calculates the rooted affinity cost between the cluster and a tree
        @param cluster
        @param tree
        :returns: affinity cost from tree1 to tree2
        """
        intersection_lookup = dict()
        min_distance = 10000000000000000
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
        return min_distance

class ClusterInfo():
    def __init__(self):
        self.cluster = {}
        self.sub_clusters = []
        self.distance = -10
