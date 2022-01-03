# -*- coding: utf-8 -*-
"""
This module implements the quad-edge data structure used by Gubias and Stolfi,
as well as defining classes for storing multiple edge objects.
"""

# ---------------------------------- Imports ----------------------------------

# Repo module imports
import utils
import numpy as np

# --------------------------------- Edge class --------------------------------

class Edge():
    """
    This class represents a single edge as defined in the Guibas and Stolfi
    edge algebra formalism.

    Attributes
    ----------
    index : int
        unique edge index
    org : int
        index of edge origin point
    dest : int
        index of edge destination point
    sym : int
        index of symetric edge
    onext : int
        index of next ccw edge connected to the origin point
    oprev : int
        index of previous ccw edge connected to the origin point
    deactivate : bool
        status of edge in triangulation. False if the edge is still part of
        the triangulation.
    """
    def __init__(self, idx, org, dst, s, onxt, oprv):
        self.index = idx
        self.org = org
        self.dest = dst
        self.sym = s
        self.onext = onxt
        self.oprev = oprv
        self.deactivate = False

    def __str__(self):
        return (
                f"EDGE| Index: {self.index}, Org: {self.org}, "
                f"Dest: {self.dest}, Sym: {self.sym}, "
                f"Onext: {self.onext}, Oprev: {self.oprev}, "
                f"Deactivate: {self.deactivate} |"
                )

    def shift_indices(self, shift_edges, shift_points):
        """
        This function is used to shift the indices of an edge. This is used
        when merging sets of edges to ensure each edge has a unique index.

        Parameters
        ----------
        shift_edges : int
            integer to shift edge indices by
        shift_points : int
            integer to shift orgin and destination points by
        """
        self.index += shift_edges
        self.org += shift_points
        self.dest += shift_points
        self.sym += shift_edges
        self.onext += shift_edges
        self.oprev += shift_edges

def setup_edge(origin, dest, edge_idx):
    """
    This function takes in the index of two points and creates an edge array,
    as well as an edge for the symetric edge.

    Parameters
    ----------
    origin : int
        index of origin point
    dest : int
        index of destination point
    edge_idx : int
        The index of the new edge

    Returns
    -------
    edge : Edge class
        Edge connecting org to dest
    edge_sym : Edge class
        Edge connecting dest to org
    """
    e1_idx = edge_idx
    e2_idx = edge_idx + 1

    edge = Edge(e1_idx, origin, dest, e2_idx, e1_idx, e1_idx)
    edge_sym = Edge(e2_idx, dest, origin, e1_idx, e2_idx, e2_idx)

    return edge, edge_sym

# -------------------------------- Edges class --------------------------------

class Edges():
    def __init__(self):
        self.edges = []
        self.num_edges = 0
        self.inner = None
        self.outer = None

    def push_back(self, new_edge):
        self.edges.append(new_edge)
        self.num_edges += 1

    def set_extreme_edges(self, left_most_edge, right_most_edge):
        self.inner = left_most_edge
        self.outer = right_most_edge

    def __str__(self):
        return f'EDGES | Edges: {[str(edge) for edge in self.edges]}, Num_edges: {self.num_edges}, Inner: {self.inner}, Outer: {self.outer} |'

    def onext(self, edge):
        return self.edges[edge].onext

    def oprev(self, edge):
        return self.edges[edge].oprev

    def set_onext(self, dest_edge, source_edge):
        self.edges[dest_edge].onext = source_edge

    def set_oprev(self, dest_edge, source_edge):
        self.edges[dest_edge].oprev = source_edge

    def dest(self, edge):
        return self.edges[edge].dest

    def set_dest(self, dest_edge, source_edge):
        self.edges[dest_edge].dest = source_edge

    def org(self, edge):
        return self.edges[edge].org

    def set_org(self, dest_edge, source_edge):
        self.edges[dest_edge].org = source_edge

    def sym(self, edge):
        return self.edges[edge].sym

    def set_sym(self, dest_edge, source_edge):
        self.edges[dest_edge].sym = source_edge

    def set_lprev(self, dest_edge, source_edge):
        self.edges[self.edges[dest_edge].onext].sym = source_edge

    def lprev(self, edge):
        return self.edges[self.edges[edge].onext].sym

    def set_rprev(self, dest_edge, source_edge):
        self.edges[self.edges[dest_edge].sym].onext = source_edge

    def rprev(self, edge):
        return self.edges[self.edges[edge].sym].onext

    def set_lnext(self, dest_edge, source_edge):
        self.edges[self.edges[dest_edge].sym].oprev = source_edge

    def lnext(self, edge):
        return self.edges[self.edges[edge].sym].oprev

    def set_rnext(self, dest_edge, source_edge):
        self.edges[self.edges[dest_edge].oprev].sym = source_edge

    def rnext(self, edge):
        return self.edges[self.edges[edge].oprev].sym

    def splice(self, edge1, edge2):
        """
        This function is used when adding another edge to a point in the
        triangulation which has existing edges connected to it. The next and
        previous ccw edges are updated accordingly for each of the input edges.

        Parameters
        ----------
        edge1 : Edge class
        edge2 : Edge class
        """
        self.set_oprev(self.onext(edge1), edge2)
        self.set_oprev(self.onext(edge2), edge1)

        edge_1_holder = self.onext(edge1)
        edge_2_holder = self.onext(edge2)

        self.set_onext(edge1, edge_2_holder)
        self.set_onext(edge2, edge_1_holder)

    def connect(self, edge1, edge2):
        """
        This function takes two seperated edges and creates a new edge
        connecting the two.

        Parameters
        ----------
        edge1 : Edge class
        edge2 : Edge class

        Returns
        -------
        out : int
            index of the created edge
        """

        last_index = self.num_edges
        edge, edge_sym = setup_edge(self.dest(edge1), self.org(edge2), last_index)

        self.push_back(edge)
        self.push_back(edge_sym)

        edge1_sym_oprev = self.oprev(self.sym(edge1))
        new_edge_index = edge.index

        self.splice(new_edge_index, edge1_sym_oprev)
        self.splice(self.sym(new_edge_index), edge2)

        return new_edge_index

    def kill_edge(self, e):
        """
        This function removes an edge from the triangulation by setting the
        status of edge.deactivate to True. The function also fixed the
        connecting edges too.

        Parameters
        ----------
        e : Edge class
            edge to remove from the triangulation
        """
        # Fix the local triangulation
        self.splice(e, self.edges[e].oprev)
        self.splice(self.edges[e].sym, self.edges[self.edges[e].sym].oprev)

        # Set the status of the edge and it's symetric edge to kill
        self.edges[e].deactivate = True
        self.edges[self.edges[e].sym].deactivate = True

# ---------------------------- Triangulation class ----------------------------

class TriangulationEdges(Edges):
    def __init__(self, points_subset):
        super().__init__()
        self.points = points_subset
        if len(points_subset) == 2:
            p1, p2 = 0, 1
            edge, edge_sym = setup_edge(p1, p2, 0)
            self.push_back(edge)
            self.push_back(edge_sym)

            left_most_edge = edge.index
            right_most_edge = self.edges[edge.index].sym

            self.set_extreme_edges(left_most_edge, right_most_edge)
            return
        elif len(points_subset) == 3:
            p1, p2, p3 = 0, 1, 2
            # Create the first two edges of the triangle
            edge1, edge1_sym = setup_edge(p1, p2, 0)
            self.push_back(edge1)
            self.push_back(edge1_sym)

            edge2, edge2_sym = setup_edge(p2, p3, 2)
            self.push_back(edge2)
            self.push_back(edge2_sym)

            self.splice(edge1_sym.index, edge2.index)

            pt1 = points_subset[self.edges[edge1.index].org]
            pt2 = points_subset[self.edges[edge1.index].dest]
            pt3 = points_subset[p3]

            if utils.is_right_turn(pt1, pt2, pt3):
                # Points are in CCW orientiaton
                c = self.connect(edge2.index, edge1.index)
                self.set_extreme_edges(edge1.index, edge2_sym.index)
                return

            if utils.is_left_turn(pt1, pt2, pt3):
                # Points are in CW orientiaton
                c = self.connect(edge2.index, edge1.index)
                self.set_extreme_edges(self.edges[c].sym, c)
                return

            # Points are collinear
            self.set_extreme_edges(edge1.index, edge2_sym.index)
        else:
            raise ValueError("Initializing subset of size 3 not supported")

    def __str__(self):
        return f"TEDGES | {self.points}, " + super().__str__()

    def get_min_x(self):
        return min(self.points, key = lambda t: t[0])[0]

    def get_max_x(self):
        return max(self.points, key = lambda t: t[0])[0]

    def get_min_y(self):
        return min(self.points, key = lambda t: t[1])[1]

    def get_max_y(self):
        return max(self.points, key = lambda t: t[1])[1]


    def shift_indices(self, shift_edges, shift_points):
        for edge in self.edges:
            edge.shift_indices(shift_edges, shift_points)

    def combine_data(self, appendage_triangulation):
        # Calculate the capacity of the new edges array
        base_triangulation_num_edges = self.num_edges
        base_triangulation_num_points = len(self.points)
        appendage_triangulation_num_edges = appendage_triangulation.num_edges

        appendage_triangulation.shift_indices(base_triangulation_num_edges, base_triangulation_num_points)

        self.edges.extend(appendage_triangulation.edges)
        self.num_edges += appendage_triangulation_num_edges
        self.points.extend(appendage_triangulation.points)

        return self
