import numpy as np

class PPI_network(PPI, location_data):
    def __init__(self, PPI, location_data):
        self.PPI = PPI
        self.loc_data = location_data
        self.G = None


    def process_group_loc(self):
        #compress groups into keys and values as proteins and location
        return loc_dict

    def create_node_list(self):
        #group by ID, compress groups, generate ndoes
        return node_list

    def create_edge_list(self):
        #generate from PPI data frame
        return edge_list

    def starting_position(self):
        #choose starting position of proteins and set as pos
        return starting_pos

    def create_network(self):
        #use networkx to generate graph with add_x_from
        #save G into __init__
        return G

    def draw_edges(self):
        #draw edges with the positional arguments
        return 

    def draw_nodes(self):
        #draw nodes with positional arguments
        return

    def draw_network(self):
        #create plot, use draw_nodes and edges to draw network
        return network_plot

    def graph_history(self):
        #save generated graph temporarly
        return 

    def user_search(self):
        #create interactive terminal
        return search_request

    def network_response(self):
        #search_request response based on new or old request
        return new_plot
