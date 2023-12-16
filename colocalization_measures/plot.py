"""
Plotting script for figures in the paper
"""

import matplotlib.pyplot as plt
import seaborn as sns
import igraph as ig
import numpy as np

import os 

from colocalization_measures.adjusted_local_assortativity import calculate_adjusted_local_assorativity


def local_assort_plot(graph, marker, marker_threshold=0):
    """_summary_

    :param graph: _description_
    :param marker: _description_
    """

    # get assortativity scores after normalization
    loc_rounded = calculate_adjusted_local_assorativity(graph, marker, marker_threshold)
    loc_rounded = [np.round(score, 2) for score in loc_rounded]

    # get color pallet with seaborn
    positive_scores = [score for score in loc_rounded if score > 0]
    negative_scores = [score for score in loc_rounded if score < 0]

    # create color dictionaries
    color_dict_pos = dict(zip(sorted(list(set(positive_scores))), sns.color_palette("Reds", n_colors=len(list(set(positive_scores)))))) 
    color_dict_neg = dict(zip(sorted(list(set(negative_scores))), sns.color_palette("Blues", n_colors=len(list(set(negative_scores)))))) 
    color_dict_zero = dict(zip([0], ["white"]))

    color_dict = color_dict_pos | color_dict_neg | color_dict_zero
                
    # color vertices
    if np.isnan(loc_rounded).any():
        colors_vertices = ["grey"] * graph.vcount()
    else:
        colors_vertices = [color_dict[num] for num in loc_rounded]

    # get mean of the assortativity scores
    mean_loc = np.mean(positive_scores)

    # define vertex sizes if score is above the mean
    size_vertices = [50 if loc > mean_loc else 15 for loc in loc_rounded]

    return colors_vertices, size_vertices


def marker_count_plot(graph, marker):
    """_summary_

    :param graph: _description_
    :param marker: _description_
    :return: _description_
    """

    # check wether the current graph is a line graph or an A-node-projection
    if "marker" in set(graph.vs.attribute_names()):

        # set vertex color
        colors_vertices = ["#FF9404" if vertex["marker"] else "#116178" for vertex in graph.vs]

        # define vertex sizes 
        size_vertices = [50 if vertex["marker"] == marker else 15  for vertex in graph.vs]

    else:
        # set vertex color
        counts = [np.round(vertex["markers"][marker] / sum(list(vertex["markers"].values())), 3) for vertex in graph.vs]
        counts_vertix = counts.copy()

        counts += list(range(-10,-1))
        counts.sort()

        colors = sns.dark_palette("#FF9404", reverse=False, n_colors=len(set(counts)))
                                
        color_dict = dict(zip(np.unique(counts), colors))

        colors_vertices = [color_dict[vertex_counts] if vertex_counts > 0 else "#116178" for vertex_counts in counts_vertix]

        # define vertex sizes 
        size_vertices = [50 if vertex["markers"][marker] > 0 else 15  for vertex in graph.vs]

    
    return colors_vertices, size_vertices



def plot_mutliple_markers(graph, list_markers, assortativity):
    """_summary_

    :param graph: _description_
    :param list_markers: _description_
    """

    # create figure
    fig = plt.figure(figsize=(20, 20))

    # number plots
    n_plots = len(list_markers)

    # set parameters for the multi plot
    number_grids = int(np.sqrt(n_plots)) + 1
    columns = number_grids
    rows = number_grids


    for index_place, marker in enumerate(list_markers):
        
        if assortativity:
            colors_vertices, size_vertices = local_assort_plot(graph, marker)
        else:
            colors_vertices, size_vertices = marker_count_plot(graph, marker)

        # plot graph and save as temporary file
        ig.plot(graph, vertex_size=size_vertices, vertex_color=colors_vertices, 
                layout=graph.layout_kamada_kawai(), bbox=(2000,2000)).save('temporary.png') 

        # create subfigure
        ax = fig.add_subplot(rows, columns, index_place + 1)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])

        # load png for matplot lib
        img = plt.imread('temporary.png')

        # plot the png in the panel
        plt.imshow(img)

        # remove tempory file
        os.remove('temporary.png')  

        # plot title
        if type(marker) == str:
            plt.title(marker)


    if assortativity:
        # Add a colorbar   
        sm = plt.cm.ScalarMappable(cmap="RdBu_r")
        sm.set_array([])

        v1 = np.array([0, 0.5, 1])

        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.88, 0.35, 0.04, 0.5])

        cbar = fig.colorbar(sm, ticks=v1, cax=cbar_ax)
        cbar.ax.set_yticklabels(["Disassortative", "Unifrom mixing", "Assortative"])

    # display the final plot
    plt.show()
