import os
import pandas as pd
from pathlib import Path

# load functions from pixelator
from pixelator.pixeldataset import PixelDataset
from pixelator.graph import Graph

# load standard workflow single cell functions
from colocalization_measures.adjusted_local_assortativity import calculate_adjusted_local_assorativity


def load_pxl_data(path_name):
    """
    Load pixelator data from .pxl files and create an A-node projection for each file loaded. Save the mega-graphs and .pxl
        files into a list for later processing.

    Args:
        path_name (str): The path to the directory where all the .pxl are saved.

    Returns:
        List, List: List of .pxl data sets and a list containing the mega-graphs as iGraphs.Graphs
    """

    # new graph for all data
    edge_list_df_list = []
    pxl_data_list = []
    file_names = []

    for file in sorted(os.listdir(path_name)):
        
        # get name of current file in directory
        filename = os.fsdecode(file)

        # load marker edge list
        pxl_data = PixelDataset.from_file(path_name + "/" + filename)
        edge_list_df = pxl_data.edgelist   

        edge_list_df.columns = [col.lower() for col in edge_list_df.columns]

        edge_list_df_list.append(edge_list_df)
        pxl_data_list.append(pxl_data)
        file_names.append(filename)

        print(filename)

    sample_list = []

    for index_df, df in enumerate(edge_list_df_list):
        # build the graph
        graph = Graph.from_edgelist(edgelist=df,
                                    add_marker_counts=True,
                                    simplify=False,
                                    use_full_bipartite=False)
        mega_graph = graph._raw
        
        # add filename as a vertex attribute
        mega_graph.es["sample"] = [file_names[index_df]] * mega_graph.ecount()

        sample_list.append(mega_graph)

    return pxl_data_list, sample_list


def get_component_name_dict(pxl_data_list, sample_list, sample_index):
    """
    Creates a list of each component (graph) for the current sample (sample index in the sample list). Based on this
        a dictonary for A and B pixels to component name is created. This together with the set of pixel names is later
        used for nameing the component.

    Args:
        pxl_data_list (List[.pxl, .pxl, ...]): A list containing all loaded .pxl flies from the current sample.
        sample_list (List[iGraph.Graph, iGraph.Graph(), ...]]): A list containing all the mega-graphs for each sample as iGraph.Graphs()
        sample_index (int): The index of the currently analyzed sample. 

    Returns:
        List[Graphs(), ...], Dict, Dict, Set, Set: A list with all graphs of the current sample (1 graph = 1 cell). A and B UPIs dict
            translating UPI names to component names and two sets containing all UPIA and UPIB names respectivly. 
    """
    

    # subgraphs of the selected sample
    graph_component_sample = sample_list[sample_index].components().subgraphs()

    # load edge list from pixelator dataset
    df = pxl_data_list[sample_index].edgelist

    # create dictionary form all pixel names for each component
    dict_comp_A = dict(zip(df["upia"].values, df["component"].values))
    dict_comp_B = dict(zip(df["upib"].values, df["component"].values))

    # create a set for all UPIA and UPIB for fast assessments
    set_comp_A_pixels = set(df["upia"].values)
    set_comp_B_pixels = set(df["upib"].values)

    return graph_component_sample, dict_comp_A, dict_comp_B, set_comp_A_pixels, set_comp_B_pixels


def create_local_assorativity_dfs(directory_path, graph_component_sample, set_comp_A_pixels, dict_comp_A, dict_comp_B, marker_names):
    """
    Creates the local assorativity dfs for each component in the graph_component_sample list and saves the df as a .csv names after the
        component name. 

    Args:
        directory_path (str): Path where the .csv should be saved. Must end with ".../general_folder" and !not! with "/" as for each sample
            in the experiment a folder is added with "/S0", "/S1", ...
        graph_component_sample (List[Graph()]): List of all components (cells) in the current sample.
        set_comp_A_pixels (Set): Set of all UPIA names.
        dict_comp_A (Dict): Dict translating UPIAs to component names.
        dict_comp_B (Dict): Dict translating UPIBs to component names.
        marker_names (List[str, str, ...]): All markers for which the adjusted local assorativity should be computed. 
            E.g. ["CD20", "CD50" , ...]
    """


    # iterate through the first components in the current sample
    for component in graph_component_sample:

        # get the component name by vertex name
        if component.vs[0]["name"] in set_comp_A_pixels:
            comp_name = dict_comp_A[component.vs[0]["name"]]
        else:
            comp_name = dict_comp_B[component.vs[0]["name"]]

        # get file info
        component_file_name = directory_path + str(comp_name)
        component_file_path = Path(component_file_name)

        # check if file in directry 
        if not component_file_path.is_file():

            # create a df with the local assorativity scores for the current component
            local_assort_score_df = pd.DataFrame()
            for marker in marker_names:

                # calculate local assorativity scores for each marker in list
                loc_rounded = calculate_adjusted_local_assorativity(component, marker)
                local_assort_score_df[marker] = loc_rounded   

            # save df
            local_assort_score_df.to_csv(component_file_name)


if __name__ == "__main__":

    # definine all markers used
    marker_names = ['CD9', 'CD62P', 'CD27', 'CD36', 'mIgG2b', 'CD337', 'CD3E', 'CD274', 'CD55', 'CD161', 'CD14', 'CD48', 'CD82', 'CD64', 
                    'CD11c', 'CD54', 'CD11b', 'CD44', 'CD154', 'B2M', 'CD268', 'CD18', 'CD37', 'CD4', 'CD29', 'CD11a', 'CD47', 'CD7', 
                    'CD2', 'CD35', 'CD45', 'CD8', 'CD314', 'CD22', 'CD19', 'CD127', 'CD53', 'CD52', 'CD229', 'CD72', 'CD59', 'CD163', 'CD38', 
                    'CD25', 'CD41', 'CD150', 'CD278', 'ACTB', 'CD152', 'mIgG1',  'CD5', 'CD26', 'CD197', 'CD50', 'CD328', 'CD279', 'CD200',
                    'CD71', 'CD102', 'CD244', 'CD45RB', 'CD40', 'CD45RA', 'CD84', 'CD49D', 'CD162', 'CD1d', 'CD137', 'CD32', 'CD69', 'CD20', 'CD33', 
                    'CD158', 'HLA-ABC', 'mIgG2a', 'CD86', 'CD43', 'CD16', 'TCRb', 'HLA-DR']
    
    marker_names_retuximab = marker_names + ["Rituximab"]

    # select the current markers depending on the sample
    markers_selcted = marker_names_retuximab

    # set paths to load the data
    path_data_loaded = "...load data from path here..."
    
    # set directory where the component df should be saved
    path_data_saved = "...save data to path here..."

    # load pixel data and mega graphs into two lists
    pxl_data_list, sample_list = load_pxl_data(path_data_loaded)

    for sample_index in range(len(sample_list)):
        
        # seperate the samples into different folders
        path_data_saved += "/S" + str(sample_list)
        Path(path_data_loaded).mkdir(parents=True, exist_ok=True)

        # get a dictionary to name components accoringly in the dfs
        graph_component_sample, dict_comp_A, dict_comp_B, set_comp_A_pixels, set_comp_B_pixels = get_component_name_dict(pxl_data_list, sample_list, sample_index)

        # create and safe adjusted local assorativity dfs
        create_local_assorativity_dfs(path_data_saved, graph_component_sample, set_comp_A_pixels, dict_comp_A, dict_comp_B, markers_selcted)
