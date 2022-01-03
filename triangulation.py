import math
import utils
import edge_topology
import matplotlib.pyplot as plt
import numpy as np

def find_correct_edge(triangulation, left_candidate, right_candidate, left_candidate_valid):
    a = triangulation.points[triangulation.edges[right_candidate].dest]
    b = triangulation.points[triangulation.edges[right_candidate].org]
    c = triangulation.points[triangulation.edges[left_candidate].org]
    deciding_point = triangulation.points[triangulation.dest(left_candidate)]
    result = not utils.in_circle(a, b, c, deciding_point) and left_candidate_valid
    return result

def generalized_candidate_finder(hull, candidate, reference_point_1, reference_point_2, left=True):
    if left:
        candidate_neighbor = hull.oprev(candidate)
    else:
        candidate_neighbor = hull.onext(candidate)

    candidate_neighbor_dest = hull.dest(candidate_neighbor)
    candidate_dest = hull.dest(candidate)
    is_ccw = utils.is_right_turn(reference_point_1, reference_point_2, hull.points[candidate_neighbor_dest])
    next_candidate_valid = utils.in_circle(reference_point_2, reference_point_1, hull.points[candidate_dest], hull.points[candidate_neighbor_dest])
    while is_ccw and not next_candidate_valid:

        if left:
            next_candidate = hull.oprev(candidate)
        else:
            next_candidate = hull.onext(candidate)

        hull.kill_edge(candidate)
        candidate = next_candidate

        if left:
            next_candidate = hull.oprev(candidate)
        else:
            next_candidate = hull.onext(candidate)

        candidate_neighbor_dest = hull.dest(next_candidate)
        candidate_dest = hull.dest(candidate)
        is_ccw = utils.is_right_turn(reference_point_1, reference_point_2, hull.points[candidate_neighbor_dest])
        next_candidate_valid = utils.in_circle(reference_point_2, reference_point_1, hull.points[candidate_dest], hull.points[candidate_neighbor_dest])

    return hull, candidate

def lowest_common_tangent(left_hull, right_hull):
    left_hull_outer_edge = left_hull.outer
    right_hull_inner_edge = right_hull.inner

    left_reference_edge_origin = left_hull.points[left_hull.org(left_hull_outer_edge)]
    left_reference_edge_dest = left_hull.points[left_hull.dest(left_hull_outer_edge)]

    right_reference_edge_origin = right_hull.points[right_hull.org(right_hull_inner_edge)]
    right_reference_edge_dest = right_hull.points[right_hull.dest(right_hull_inner_edge)]

    is_right_incorrect = utils.is_left_turn(right_reference_edge_origin, right_reference_edge_dest, left_hull.points[left_hull.edges[left_hull_outer_edge].org])
    is_left_incorrect = utils.is_right_turn(left_reference_edge_origin, left_reference_edge_dest, right_hull.points[right_hull.edges[right_hull_inner_edge].org])
    while is_left_incorrect or is_right_incorrect:
        if is_right_incorrect:
            right_hull_inner_edge = right_hull.lnext(right_hull_inner_edge)
            right_reference_edge_origin = right_hull.points[right_hull.org(right_hull_inner_edge)]
            right_reference_edge_dest = right_hull.points[right_hull.dest(right_hull_inner_edge)]
        elif is_left_incorrect:
            left_hull_outer_edge = left_hull.rprev(left_hull_outer_edge)
            left_reference_edge_origin = left_hull.points[left_hull.org(left_hull_outer_edge)]
            left_reference_edge_dest = left_hull.points[left_hull.dest(left_hull_outer_edge)]
        is_right_incorrect = utils.is_left_turn(right_reference_edge_origin, right_reference_edge_dest, left_hull.points[left_hull.org(left_hull_outer_edge)])
        is_left_incorrect = utils.is_right_turn(left_reference_edge_origin, left_reference_edge_dest, right_hull.points[right_hull.org(right_hull_inner_edge)])
    return left_hull_outer_edge, right_hull_inner_edge

def connect_triangulations(combined_edges, handle_indices):
    ldi = handle_indices['left_inner']
    rdi = handle_indices['right_inner']
    ldo = handle_indices['left_outer']
    rdo = handle_indices['right_outer']
    base = combined_edges.connect(combined_edges.sym(ldi), rdi)
    if np.array_equal(combined_edges.points[combined_edges.org(ldi)], combined_edges.points[combined_edges.org(ldo)]):
        ldo = base
    if np.array_equal(combined_edges.points[combined_edges.org(rdi)], combined_edges.points[combined_edges.org(rdo)]):
        rdo = combined_edges.sym(base)
    combined_edges.set_extreme_edges(ldo, rdo)
    return base

def combine_triangulations(hull_left, hull_right, ldi, rdi):
    handle_indices = {
        'left_outer': hull_left.inner,
        'right_outer': hull_right.outer + hull_left.num_edges,
        'right_inner': rdi + hull_left.num_edges,
        'left_inner': ldi
    }
    combined_triangulation = hull_left.combine_data(hull_right)
    base_point_index = connect_triangulations(combined_triangulation, handle_indices)
    new_triangulation = weave_hulls(combined_triangulation, base_point_index)
    return new_triangulation

def weave_hulls(triangulation, base):
    base_origin = triangulation.points[triangulation.org(base)]
    base_dest = triangulation.points[triangulation.dest(base)]

    right_candidate = triangulation.rprev(base)
    right_candidate_dest = triangulation.points[triangulation.dest(right_candidate)]
    is_right_candidate_correct = utils.is_right_turn(base_origin, base_dest, right_candidate_dest)

    left_candidate = triangulation.oprev(base)
    left_candidate_dest = triangulation.points[triangulation.dest(left_candidate)]
    is_left_candidate_correct = utils.is_right_turn(base_origin, base_dest, left_candidate_dest)

    are_either_candidate_incorrect = True
    while are_either_candidate_incorrect:
        if not is_right_candidate_correct and not is_left_candidate_correct:
            break

        if is_right_candidate_correct:
            triangulation, right_candidate = generalized_candidate_finder(triangulation, right_candidate, base_origin, base_dest, left=False)

        if is_left_candidate_correct:
            triangulation, left_candidate = generalized_candidate_finder(triangulation, left_candidate, base_origin, base_dest, left=True)

        left_candidate_strong_valid = find_correct_edge(triangulation, left_candidate, right_candidate, is_left_candidate_correct)

        if not is_right_candidate_correct or left_candidate_strong_valid:
            base = triangulation.connect(left_candidate, triangulation.sym(base))
        else:
            base = triangulation.connect(triangulation.sym(base), triangulation.sym(right_candidate))

        base_origin = triangulation.points[triangulation.org(base)]
        base_dest = triangulation.points[triangulation.dest(base)]

        right_candidate = triangulation.rprev(base)
        right_candidate_dest = triangulation.points[triangulation.dest(right_candidate)]
        is_right_candidate_correct = utils.is_right_turn(base_origin, base_dest, right_candidate_dest)

        left_candidate = triangulation.oprev(base)
        left_candidate_dest = triangulation.points[triangulation.dest(left_candidate)]
        is_left_candidate_correct = utils.is_right_turn(base_origin, base_dest, left_candidate_dest)

        if not is_right_candidate_correct and not is_left_candidate_correct:
            are_either_candidate_incorrect = False

    return triangulation

def merge_triangulations(triangulations):
    if len(triangulations) == 1:
        return triangulations
    elif len(triangulations) == 2:
        ldi, rdi = lowest_common_tangent(triangulations[0], triangulations[1])
        d_triang = combine_triangulations(triangulations[0], triangulations[1], ldi, rdi)
        return d_triang
    else:
        new_triangulations = []
        number_of_triangulations = len(triangulations)
        for index in range(0, number_of_triangulations, 2):
            if (number_of_triangulations - index) == 1:
                new_triangulations.append(triangulations[index])
            else:
                new_triangulation = merge_triangulations(triangulations[index:index+2])
                new_triangulations.append(new_triangulation)
        return merge_triangulations(new_triangulations)

def partition_into_3s(points):
    number_of_points = len(points)
    if number_of_points < 4:
        raise ValueError("Point set size needs to be 4 or greater.")
    partitioned_list = []
    for index in range(0, number_of_points, 3):
        if index + 4 >= number_of_points:
            if number_of_points % 3 == 1:
                partitioned_list.append(points[index:index+2])
                partitioned_list.append(points[index+2:index+4])
                return partitioned_list
            elif number_of_points % 3 == 2:
                partitioned_list.append(points[index:index+2])
                return partitioned_list
        else:
            partitioned_list.append(points[index:index+3])
    return partitioned_list

# Not used. Could be used for an alternate implementation of Dwyer Divide and Conquer.
def partition_into_3s_cells(positions):
    q = []
    q.append(positions)
    new_q = []
    iteration_map = [-1]
    new_iteration_map = []
    loop_counter = 0

    while any([len(item) not in [2, 3] for item in q]):
        print(f'LC {loop_counter}')
        for index, item in enumerate(q):
            if len(item) not in [2, 3]:
                split_index = len(item)//2
                if loop_counter % 2 == 0:
                    split_1 = utils.lexigraphic_sort(item[:split_index])
                    split_2 = utils.lexigraphic_sort(item[split_index:])
                else:
                    split_1 = utils.lexigraphic_swap_sort(item[:split_index])
                    split_2 = utils.lexigraphic_swap_sort(item[split_index:])
                new_q.append(split_1)
                new_q.append(split_2)
                new_iteration_map.append(loop_counter)
                new_iteration_map.append(loop_counter)
            else:
                new_q.append(item)
                new_iteration_map.append(iteration_map[index])
        q = new_q
        iteration_map = new_iteration_map
        new_q = []
        new_iteration_map = []
        loop_counter += 1
    return q, iteration_map

def triangulate_dwyer_helper(triangulations):
    if len(triangulations)!=2:
        raise ValueError("len of triangulations must be 2")
    ldi, rdi = lowest_common_tangent(triangulations[0], triangulations[1])
    combined_triangulation = combine_triangulations(triangulations[0], triangulations[1], ldi, rdi)
    return combined_triangulation

def triangulate_dwyer(points, depth=0):
    if len(points) == 1:
        raise ValueError("Not accounted for case of 1 point")
    elif len(points) in [2, 3]:
        if depth % 2 == 0:
            sorted_points = utils.lexigraphic_sort(points)
        else:
            sorted_points = utils.lexigraphic_swap_sort(points)
        triangulation = edge_topology.TriangulationEdges(sorted_points)
        return triangulation
    else:
        if depth % 2 == 1:
            points_split_index = len(points)//2
            x_coords = [point[0] for point in points]
            median_x_coord = utils.quickSelect(x_coords, points_split_index)

            points_partition_1 = []
            points_partition_2 = []
            undecided = []
            for index, point in enumerate(points):
                if point[0] < median_x_coord:
                    points_partition_1.append(point)
                elif point[0] > median_x_coord:
                    points_partition_2.append(point)
                else:
                    undecided.append(point)

            y_coords = [undecided_point[1] for undecided_point in undecided]
            num_missing_from_partition_1 = points_split_index - len(points_partition_1)
            median_y_coord = utils.quickSelect(y_coords, num_missing_from_partition_1-1)

            for undecided_point in undecided:
                num_missing_from_partition_1 = points_split_index - len(points_partition_1)
                if num_missing_from_partition_1 > 0 and undecided_point[1] <= median_y_coord:
                    points_partition_1.append(undecided_point)
                else:
                    points_partition_2.append(undecided_point)

        else:
            points_split_index = len(points)//2
            y_coords = [point[1] for point in points]
            median_y_coord = utils.quickSelect(y_coords, points_split_index)
            points_partition_1 = []
            points_partition_2 = []
            undecided = []
            for index, point in enumerate(points):
                if point[1] < median_y_coord:
                    points_partition_1.append(point)
                elif point[1] > median_y_coord:
                    points_partition_2.append(point)
                else:
                    undecided.append(point)

            x_coords = [undecided_point[0] for undecided_point in undecided]
            num_missing_from_partition_1 = points_split_index - len(points_partition_1)
            median_x_coord = utils.quickSelect(x_coords, num_missing_from_partition_1-1)
            for undecided_point in undecided:
                num_missing_from_partition_1 = points_split_index - len(points_partition_1)
                if num_missing_from_partition_1 > 0 and undecided_point[0] <= median_x_coord:
                    points_partition_1.append(undecided_point)
                else:
                    points_partition_2.append(undecided_point)

        triangulation_part_1 = triangulate_dwyer(points_partition_1, depth+1)
        triangulation_part_2 = triangulate_dwyer(points_partition_2, depth+1)
        combined = triangulate_dwyer_helper([triangulation_part_1, triangulation_part_2])

        # Set handles for dwyer implementation
        left_e = combined.inner
        right_e = combined.outer

        combined_points = combined.points

        plo = combined_points[combined.edges[left_e].org]
        pld = combined_points[combined.edges[left_e].dest]

        pro = combined_points[combined.edges[right_e].org]
        prd = combined_points[combined.edges[right_e].dest]

        if depth % 2 == 0:
            while plo[0] != combined.get_min_x():
                left_e = combined.edges[combined.edges[left_e].onext].sym
                plo = combined_points[combined.edges[left_e].org]
                pld = combined_points[combined.edges[left_e].dest]

            while pro[0] != combined.get_max_x():
                right_e = combined.edges[combined.edges[right_e].sym].onext
                pro = combined_points[combined.edges[right_e].org]
                prd = combined_points[combined.edges[right_e].dest]
        else:
            while plo[1] != combined.get_min_y():
                left_e = combined.edges[combined.edges[left_e].onext].sym
                plo = combined_points[combined.edges[left_e].org]
                pld = combined_points[combined.edges[left_e].dest]

            while pro[1] != combined.get_max_y():
                right_e = combined.edges[combined.edges[right_e].sym].onext
                pro = combined_points[combined.edges[right_e].org]
                prd = combined_points[combined.edges[right_e].dest]

        combined.set_extreme_edges(left_e, right_e)
        return combined

def triangulate(positions):
    sorted_positions = utils.lexigraphic_sort(positions)
    split_positions = partition_into_3s(sorted_positions)
    groups = [edge_topology.TriangulationEdges(points_subset) for points_subset in split_positions]
    groups = merge_triangulations(groups)
    return groups

def plot_triangulation(triangulation):
    x_coords = []
    y_coords = []
    for edge in triangulation.edges:
        if not edge.deactivate:
            index_1 = edge.org
            index_2 = edge.dest
            coord_1, coord_2 = (triangulation.points[index_1], triangulation.points[index_2])
            x_coords.append([coord_1[0], coord_2[0]])
            y_coords.append([coord_1[1], coord_2[1]])
    for x_coord, y_coord in zip(x_coords, y_coords):
        plt.plot(x_coord, y_coord)
        plt.axis([0, 10, 0, 10])
