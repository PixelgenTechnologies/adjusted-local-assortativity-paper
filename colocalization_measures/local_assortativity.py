"""
Function to calculate the local assortativity for a given marker.

Based on the paper:     
https://doi.org/10.1073/pnas.1713019115

GitHub:
https://github.com/piratepeel/MultiscaleMixing/blob/master/localassort.py
"""

import scipy 
import igraph as ig
import numpy as np

import scipy.sparse as sparse


def calculateRWRrange(W, i, alphas, n, maxIter=1000):
    """
    Calculate the personalised TotalRank and personalised PageRank vectors.
    Parameters
    ----------
    W : array_like
        transition matrix (row normalised adjacency matrix)
    i : int
        index of the personalisation node
    alphas : array_like
        array of (1 - restart probabilties)
    n : int
        number of nodes in the network
    maxIter : int, optional
        maximum number of interations (default: 1000)
    cut_off: asd
    Returns
    -------
    pTotalRank : array_like
        personalised TotalRank (personalised PageRank with alpha integrated
        out)
    References
    ----------
    See [2]_ and [3]_ for further details.
    .. [2] Boldi, P. (2005). "TotalRank: Ranking without damping." In Special
        interest tracks and posters of the 14th international conference on
        World Wide Web (pp. 898-899).
    .. [3] Boldi, P., Santini, M., & Vigna, S. (2007). "A deeper investigation
        of PageRank as a function of the damping factor." In Dagstuhl Seminar
        Proceedings. Schloss Dagstuhl-Leibniz-Zentrum für Informatik.
    """

    alpha0 = np.max(alphas)
    WT = alpha0 * np.transpose(W)
    diff = 1
    it = 1

    # initialise PageRank vectors
    pPageRank = np.zeros(n)
    pPageRank[i] = 1
    pPageRank_old = pPageRank.copy()
    pTotalRank = pPageRank.copy()

    oneminusalpha0 = 1 - alpha0

    while diff > 1e-4:
        # calculate personalised PageRank via power iteration
        pPageRank = WT @ pPageRank
        pPageRank[i] += oneminusalpha0

        # calculate difference in pPageRank from previous iteration
        delta_pPageRank = pPageRank - pPageRank_old

        # Eq. [S23] Ref. [1]
        pTotalRank += (delta_pPageRank) / ((it + 1) * (alpha0 ** it))

        # calculate convergence criteria
        diff = np.sum((delta_pPageRank) ** 2) / n
        it += 1

        if it > maxIter:
            print(i, "max iterations exceeded")
            diff = 0

        pPageRank_old = pPageRank.copy()

   
    return pTotalRank


def localAssortF_numeric(graph, attribute, include_randomise=False):
    """
    Calculate local assortativity of A (undirected network) with respect to the values in attribute

    Parameters
    ----------
        A : array_like (e.g. from nx.to_scipy_sparse_array(G))
            adjacency matrix
        attribute : array_like
            array of numeric values 

    Returns
    -------
    loc_ass : array_like
        array of values representing the local assortativity of each node
    """
    permutations = 75

    adjacency = ig.Graph.get_adjacency_sparse(graph)
    
    # normalize attribute
    if np.std(attribute) != 0:
        attribute = (attribute - np.mean(attribute)) / np.std(attribute)
    else:
        attribute = np.array([0] * len(attribute))

    ## Construct transition matrix (row normalised adjacency matrix)
    # construct diagonal inverse degree matrix
    degree = np.array(adjacency.sum(1)).flatten()
    n_vertices = graph.vcount()
    dignal_matrix = scipy.sparse.diags(1./degree, 0, format='csc')

    # Construct transition matrix (row normalised adjacency matrix)
    W = dignal_matrix @ adjacency

    DM = sparse.diags(attribute, 0, format='csc')

    WY = np.array(DM.dot(W).dot(DM).sum(1)).flatten()

    # initialise the loc_ass array
    loc_ass = np.empty(n_vertices)
    pr = np.arange(0., 1., 0.1)

    per_pr = np.zeros((n_vertices, n_vertices))
    for index_vertex in range(n_vertices):
        # Calculate personalized totalrank
        pTotalRank = calculateRWRrange(W, index_vertex, pr, n_vertices)
        per_pr[index_vertex] = pTotalRank

        if include_randomise and permutations > 0:
            attribute_permutation = attribute.copy()
            permutation_testing = np.zeros(permutations)

            for index_permutation in range(permutations):
                np.random.shuffle(attribute_permutation)
                DMp = sparse.diags(attribute_permutation, 0, format='csc')
                WYp = np.array(DMp.dot(W).dot(DMp).sum(1)).flatten()
                permutation_testing[index_permutation] = (pTotalRank*WYp).sum()

            loc_ass_mean = np.mean(permutation_testing)

            # calculate local assortativity
            loc_ass[index_vertex] = (pTotalRank*WY).sum() - loc_ass_mean
        
        else:
            # calculate local assortativity
            loc_ass[index_vertex] = (pTotalRank*WY).sum()


    return loc_ass


def normalize_local_assortativity_scores(loc_ass):
    """
    Normalizes the assortativity score to have zero mean and unit standard deviation using the Z-score.
    This is done to make the scores comparable across samples.

    :param loc_ass: Local assorativity scores as a list or numpy.array(), e.g. [1, 1, 1, ...]
    :return: Returns ne normalized scores
    """

    # get weight of all positive in negative values to get mean 0 without moving the center
    sum_positiv_scores = sum([score for score in loc_ass if score > 0])
    sum_negative_scores = np.abs(sum([score for score in loc_ass if score < 0]))

    scores_adjusted_mean = []

    # adjust all scores to have center of mass (mean) 0 without moving the true value of 0
    for score in loc_ass:
        if score > 0:
            scores_adjusted_mean.append(score / sum_positiv_scores)
        elif score < 0:
            scores_adjusted_mean.append(score / sum_negative_scores)
        else:
            scores_adjusted_mean.append(0)

    # normalize assortativity using the zscore to have standard deviation 1
    normalized_loc_ass = scipy.stats.zscore(np.array(scores_adjusted_mean))

    # log scale
    normalized_loc_ass = [np.sign(score) * np.log(np.abs(score) + 1) for score in normalized_loc_ass]

    return normalized_loc_ass




def calculate_adjusted_local_assorativity(graph, attribute):
    """
    Calculates the local assortativity score for each vertex in the given graph based on the given attribute.

    :param graph: An iGraph Graph object, must contain vertex attribute "marker", e.g graph.vs[0]["marker"] = "CD45"
    :param attribute: String, Name of the marker for which the local assortativity should be computed. E.g. "CD45"

    :return: Normalized assortativity scores 
    """
    if type(attribute) == str: 

        # check if line graph or A-node projection
        if "marker" in graph.vs.attribute_names():

            # get a numerical list based on the count of the given marker
            numerical_attribute_list = [vertex["count"] if vertex["marker"] == attribute else 0 for vertex in graph.vs]

        else:
            numerical_attribute_list = [vertex["markers"][attribute] for vertex in graph.vs]

    else:
        numerical_attribute_list = attribute

    # get local assortativity based on the graph and given attribute
    loc_ass = localAssortF_numeric(graph, numerical_attribute_list, True)
    
    # normalize the given assortativity scores using the Z-score
    normalized_local_assorativity_scores = normalize_local_assortativity_scores(loc_ass)

    # round scores to 5 digits
    normalized_local_assorativity_scores = [np.round(num, 2) for num in normalized_local_assorativity_scores]

    # return the socres
    return normalized_local_assorativity_scores


