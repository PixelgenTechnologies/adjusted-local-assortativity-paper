import pandas as pd
import numpy as np

import os

import itertools
from scipy.stats import spearmanr

from .higher_order_similarity import calculate_higher_order_similarity

# remove performance warning from dataframes
import warnings
warnings.filterwarnings('ignore')


def create_colocalization_scores(marker_names: list, path_data: str, combinations=False):
    """
    Creates the coloclization scores based on the adjusted local assorativity dataframes saves as .csv files in the given path.    

    :param marker_names: The markers that should be analyzed. Should be a list of strings containing each marker. 
        E.g. ["CD50", "CD162", ...]
    :param path_data: The path to the saved dataframes containg all the adjustes local assorativity scores for each component.
    :return: Returns a pandas dataframe where each column contains the colocalization scores for the given markers, for each component.
    """

    # create new dataframe to keep track of similarity measures
    similarity_measure_df = pd.DataFrame()

    # iterate through all components saved in the given folder for the current sample
    for comp_index, local_assort_score_df in enumerate(sorted(os.listdir(path_data))):
        local_assort_score_df = pd.read_csv(path_data + "/" + local_assort_score_df)
        names, comp_values = [], []
        
        # create combinations from markers if not passed
        if not combinations:
            combinations = list(itertools.combinations(marker_names, 2))
  
        # iterate through all combinations of the current length 
        # this picks all combinations of dataframe columns
        for combination_df_columns in combinations:

                # append values to marker values list
                loc_assort_marker_0 = local_assort_score_df[combination_df_columns[0]].fillna(0).values
                loc_assort_marker_1 = local_assort_score_df[combination_df_columns[1]].fillna(0).values

                # name all combinations accordingly 
                combination_names = str(combination_df_columns[0]) + "_" + combination_df_columns[1]
                        
                # calculate pairwise similarity score
                if np.mean(loc_assort_marker_0) != 0 and np.mean(loc_assort_marker_1) != 0:
                    new_similarity_measure = spearmanr(loc_assort_marker_0, loc_assort_marker_1)[0]
                else:
                    new_similarity_measure = 0
                
                # append names and scores to list
                names.append(combination_names)
                comp_values.append(new_similarity_measure)

        # add higher order component values to list
        similarity_measure_df["CMP_" + str(comp_index)] = comp_values

    # defragment df
    similarity_measure_df = similarity_measure_df.copy()
            
    # reset index to the combination of names
    similarity_measure_df["Marker_combinations"] = names
    similarity_measure_df = similarity_measure_df.set_index("Marker_combinations")

    return similarity_measure_df


def create_higher_order_similarity_df(marker_names, path_data, order, preselected=False):
    """
    Calcalutes the higher order coloclaization for the given markers of any passed order.

    :param marker_names: The names of the markers for which the scores should be analyzed. E.g. ["CD50", "CD162", ...]
    :param path_data: The path to the saved dataframes containg all the adjustes local assorativity scores for each component.
    :param order: The order of marker combinations. Meaning the amout of markers that are compared at the same time. 
        E.g order=3 -> CD50_CD44_CD162, order=4 -> CD50_CD44_CD162_CD37, ...
    :param preselected: If the markers of the paper should be combined with the given set of markers or all combinations should be  
        computed.
    :return: Returns a pandas dataframe where each column contains the colocalization scores for the given markers, for each component.
    """

    # create new dataframe to keep track of similarity measures
    similarity_measure_df = pd.DataFrame()

    # iterate through the first 100 components in the current sample
    for comp_index, local_assort_score_df in enumerate(sorted(os.listdir(path_data))):
        
        # create list for results
        names, comp_values = [], []

        # load current assorativity dataframe
        local_assort_score_df = pd.read_csv(path_data + "/" + local_assort_score_df)

        if order == 3:
            # if preselected calculated the combinations of the given marker with CD50 and CD44 as in the paper
            if preselected:
                combinations = itertools.product(["CD50"], ["CD44"], marker_names)

                combinations = list(combinations)
                combinations.append(('CD50', 'mIgG1', 'mIgG2a'))
                combinations.append(('CD44', 'mIgG1', 'mIgG2a'))
                combinations.append(('mIgG2b', 'mIgG1', 'mIgG2a'))
                combinations.append(('CD50', 'CD162', 'CD37'))

            # otherwise calculates all combinations of order 3 for the given markers
            else:
                combinations = list(itertools.combinations(marker_names, 3))

        elif order == 4:
            # if preselected calculated the combinations of the given marker with CD50, CD44, CD54 as in the paper
            if preselected:
                combinations = itertools.product(["CD50"], ["CD44"], ["CD54"], marker_names)

                combinations = list(combinations)
                combinations.append(("CD50", 'mIgG2b', 'mIgG1', 'mIgG2a'))
                combinations.append(("CD44", 'mIgG2b', 'mIgG1', 'mIgG2a'))
                combinations.append(("CD43", 'mIgG2b', 'mIgG1', 'mIgG2a'))
                combinations.append(("CD44", 'CD50', 'CD162', 'CD37'))
                combinations.append(("CD43", 'CD50', 'CD162', 'CD37'))

            # otherwise calculates all combinations of order 4 for the given markers
            else:
                combinations = list(itertools.combinations(marker_names, 4))

        else:
             # in all other cases all combinations are calcualted
             combinations = list(itertools.combinations(marker_names, order))
        
        
        # iterate through all combinations of the current length 
        # this picks all combinations of dataframe columns
        for combination_df_columns in combinations:
                
                list_combination_columns = []
                combination_names = ""
                for column_name in combination_df_columns:

                    # append values to marker values list
                    values_column = local_assort_score_df[column_name].fillna(0).values
                    list_combination_columns.append(values_column)

                    # name all combinations accordingly 
                    combination_names += str(column_name) + "_"

                # remove last underscore
                combination_names = combination_names[:-1]
                        
                # calculate higher order similarity score
                new_similarity_measure = calculate_higher_order_similarity(list_combination_columns)
                
                # append names and scores to list
                names.append(combination_names)
                comp_values.append(new_similarity_measure)

        # add higher order component values to list
        similarity_measure_df["CMP_" + str(comp_index)] = comp_values
            
    # defragment df
    similarity_measure_df = similarity_measure_df.copy()

    # reset index to the combination of names
    similarity_measure_df["Marker_combinations"] = names
    similarity_measure_df = similarity_measure_df.set_index("Marker_combinations")

    return similarity_measure_df