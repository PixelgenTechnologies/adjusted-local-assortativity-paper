"""
Higher order similarity measure
"""

import itertools
import numpy as np


def intersection_sets_list(set_1, list_of_sets):
    """
    Calculates the intersection of the given set and all sets in the list of sets.

    :param set_1: A set() containing the vertex indices that should be compared. 
    :param list_of_sets: A list of sets, containing at least one other set. These sets contain like set_1 vertex indices 
        that should be compared. E.g. list_of_sets = [set(), set(), ...]
    :return: The length of the intersection; the number of vertices shared by all sets as an integer.
    """

    # set sets to begin iteration 
    intersection_sets = set(set_1)

    # find intersection of all sets
    for set_2 in list_of_sets:
        intersection_sets = intersection_sets.intersection(set(set_2))

    # get number of elements in the intersection
    len_intersection = len(list(intersection_sets))

    return len_intersection



def calculate_higher_order_similarity(list_marker_scores):
    """
    Calculate higher order similarity based on the https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2373804/pdf/rsbl20060553.pdf
        paper.


    :param list_marker_scores: _description_
    """
    
    # number of markers to compare
    number_markers = len(list_marker_scores)

    # rescaling factor, is equivalent to T/T-1 from the paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2373804/pdf/rsbl20060553.pdf
    adjustement_by_number_markers = number_markers / (number_markers - 1)

    # initialize parameters 
    number_markers_detected = 0
    sets_of_markers_selected = []


    # get markers above the selected threshold
    for marker_scores in list_marker_scores:
        
        # set threshold for node selection 
        positive_marker_scores = [score for score in marker_scores if score > 0]
        if len(positive_marker_scores) > 0:
            threshold = np.mean(positive_marker_scores)
        else:
            threshold = 0

        # get positive markers as a discrete set
        marker_scores = [score for score in marker_scores if score > threshold]

        # add number of non zero markers to calculation factor
        number_markers_detected += sum(marker_scores)

        sets_of_markers_selected.append(marker_scores)


    # initalize S_T; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2373804/pdf/rsbl20060553.pdf
    inclusion_exclustion_term = number_markers_detected

    # set the length of the sets combinations, pairwise intersection, intersect three sets, four, five, ...
    # index_power + 1 give the signum of the inclusion exclusion principle 
    for index_power, number_of_combinations in enumerate(range(2, number_markers + 1)):

        # number of vertices in the current intersections
        inclusion_exclusion_sum = 0
        
        # get all intersection of two, three, .. sets, e.g \sum |A_1 \cap A_2 | + |A_2 \cap A_3| + |A_1 \cap A_3|
        for combination_markers in itertools.combinations(sets_of_markers_selected, number_of_combinations):
            
            # convert tuple to list
            marker_sets_combinations = list(combination_markers)
            
            # get intersection of the current length; e.g. |A_1 \cap A_2 \cap A_3 | 
            len_intersection = intersection_sets_list(marker_sets_combinations[0], marker_sets_combinations[1:])

            # sum uop the length of the intersection
            inclusion_exclusion_sum += len_intersection

        # add sum of intersections according to the  inclusion exclusion principle
        inclusion_exclustion_term += (-1) ** (index_power + 1) * inclusion_exclusion_sum

    # check not to divide by zero
    if number_markers_detected != 0:

        # calculate multiple-site similarity measure; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2373804/pdf/rsbl20060553.pdf
        higher_order_similarity_score = adjustement_by_number_markers * (1 - inclusion_exclustion_term / number_markers_detected)

    else:
        higher_order_similarity_score = 0


    return higher_order_similarity_score
