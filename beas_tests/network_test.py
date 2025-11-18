import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
pd.set_option('display.max_colwidth', None) 

#def function for processing one group into {loc: (x,y), etc.}
	#create loc_dict
	#create tuple (loc, x, y) via zip
		#create entry in loc_dict with key=loc, value= (x,y)
	#return dict
def process_group_loc(group):
	loc_dict = {}
	for loc, x_val, y_val, z_val in zip(group['Location'], group['x'], group['y'], group['z']):
		loc_dict[loc] = (x_val, y_val, z_val)
	return loc_dict




if __name__ == '__main__':
	#Executable

	# SET-UP	
	PPI_df = pd.read_csv('test_PPI.csv')
	loc_df = pd.read_csv('test_loc_ppi.csv')
	grouped_loc_df = loc_df.groupby('ID')

	#processing each group from group.object into Protein {loc_dict}, dtype: object
	compressed_groups = grouped_loc_df.apply(process_group_loc)

	#create list = [(protein_node, location_dict),()]
	node_list = []
	for p_id, loc_dict in compressed_groups.items():
		node_list.append((p_id, loc_dict))

	# Initial positions for each node
	pos = {protein: next(iter(locs.values())) for protein, locs in compressed_groups.items()}

	
	#-------------------------------------------------------------
			# CHECK SET-UP

	#test PPI df
	print(PPI_df)
	#test loc_df
	print(loc_df)
	#object print
	print(grouped_loc_df)
	#summarized protein location coordinates
	print(compressed_groups)
	#protein list for networkx
	print(node_list)

	#-------------------------------------------------------------
			# CREATING NETWORK
	G = nx.Graph()
	G.add_nodes_from(node_list)
	edge_list = list(PPI_df.itertuples(index=False, name=None))  # [(p1, p2), ...]
	G.add_edges_from(edge_list)
	
	#-------------------------------------------------------------
			# DRAW NETWORK 
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	#draw edges using starting positions
	for p1, p2 in G.edges():	#edge between p1 and p2
		x = [pos[p1][0], pos[p2][0]] 	#first entry of tuple position pos
		y = [pos[p1][1], pos[p2][1]]
		z = [pos[p1][2], pos[p2][2]]
		ax.plot(x,y,z, color='black')

	#draw nodes with starting positions
	for node in G.nodes():
		x, y, z = pos[node]
		ax.scatter(x, y, z, s=100, color='skyblue')
		ax.text(x, y, z, node, size=10)

	#visualize plot
	#plt.draw()
	#plt.show()

	#----------------------------------------------------------
			#INTERACTIVE TERMINAL
	while True:
		inp = input("Enter protein and location (e.g., Q08379 loc2): ")
		protein, loc = inp.split(maxsplit=1)
		if protein in compressed_groups and loc in compressed_groups[protein]:
			pos[protein] = compressed_groups[protein][loc]  # update coordinates

			ax.clear()  # clear previous drawing
				    # redraw edges
			for u, v in G.edges():
				x = [pos[u][0], pos[v][0]]
				y = [pos[u][1], pos[v][1]]
				z = [pos[u][2], pos[v][2]]
				ax.plot(x, y, z, color='black')
			# redraw nodes
			for node in G.nodes():
				x, y, z = pos[node]
				ax.scatter(x, y, z, s=100, color='skyblue')
				ax.text(x, y, z, node, size=10)

			plt.draw()
			plt.pause(0.1)	

