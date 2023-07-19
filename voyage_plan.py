




# --------------------------------------------------
# --------------------------------------------------
#                       PACKAGES
# --------------------------------------------------
# --------------------------------------------------

import numpy as np
import pandas as pd
import geopandas  as gpd
from geopy.distance import geodesic as distance_GD 
import csv
from shapely.geometry import Point
import folium
import webbrowser


# --------------------------------------------------
# --------------------------------------------------
#                       Functions
# --------------------------------------------------
# --------------------------------------------------



#               MAP Related Functions
# --------------------------------------------------
# --------------------------------------------------


def find_nearest_point(Coords_Search,Coords_list):
    
    
    """
    Use to give coordinates to the function as inputs and the function returns the closest point on the available map
    
    Comments
    Requries the geopy package to compute the distances based on teh LAT, LON values
    
    Inputs
    Coords_search: array with [LAT, LON] values of the target to find
    Coords_list: n_nodes X 2 array with the [LAT, LON] values of all nodes in the grid
    
    
    Outputs
    closest_id: the id of the node that is closest to the coords given
    closest_coords: the coordinates of the closest node
    
    """
    
    node_search_tile = np.tile(Coords_Search,[len(Coords_list),1])
    
    
    
    temp_df = pd.DataFrame(data={'LAT_Search':node_search_tile[:,0],
                                 'LON_Search':node_search_tile[:,1],
                                 'LAT_List':Coords_list[:,0],
                                 'LON_List':Coords_list[:,1]}, index = range(len(Coords_list)))
    
    
    temp_df['Distance'] = temp_df.apply(lambda x: distance_GD((x['LAT_Search'], x['LON_Search']), (x['LAT_List'], x['LON_List'])).km, axis = 1)
    closest_id = temp_df['Distance'].values.argmin()
    closest_coords =  Coords_list[closest_id,:]
    
    
    return closest_id, closest_coords





def nd2nd_distance(target_Coords,Coords):
    
    """
    This Function is used to compute the distance between the node in the target_Coords and the nodes in Coords. 
        
    Inputs
    ----------
    target_Coords: n1 x 2 array where n1 is the number of target Coords. Columns are the LAT, LON values of each node
    Coords: n x 2 array where n is the number of Coords. 
    
    Comments
    -----------
    n1 = 1 or n1 = n, otherwise the function does not work
    
    Requries the geopy package to compute the distances based on teh LAT, LON values.
    
    
    Returns
    -------------
    D = [km] n x 1 array of distances
    
    """
    
    if target_Coords.shape[1]!=2 or Coords.shape[1]!=2:
        raise ValueError("The columns should be the Coords")
    
    if target_Coords.shape[0]!=1 and target_Coords.shape[0]!=Coords.shape[0]:
        raise ValueError("target_Coords should be one point or the same number of points as in Coords")
    
    if target_Coords.shape[0]==1:
        pt1=np.tile(target_Coords,target_Coords.shape[0])
    elif target_Coords.shape[0]==Coords.shape[0]:
        pt1=target_Coords
    
    # Create the dataframe with the coords
    temp_df = pd.DataFrame(data={'LAT_target':pt1[:,0],
                                 'LON_target':pt1[:,1],
                                 'LAT_List':Coords[:,0],
                                 'LON_List':Coords[:,1]}, index = range(len(Coords)))
    
    
    # Compute the distances
    temp_df['Distance'] = temp_df.apply(lambda x: distance_GD((x['LAT_target'], x['LON_target']), (x['LAT_List'], x['LON_List'])).km, axis = 1)
    D = temp_df['Distance'].values[:,None]
    
    return D




def get_row_and_column_lists(Connectivity):
    """
    Comments
    -----------
    This function generates the rows/columns when the connectivity list is available
    
    
    Parameters
    ----------
    Connectivity : List with length equal to the number of nodes. Each row is a node. Each element describes the connections between the nodes in the grid.
    
    Returns
    -------
    rows, columns: [n_nodes x 1] array containing the 
    
    """    
    
    # Find how many connections in Connectivity
    temp= [len(x) for x in Connectivity]
    total_connections = np.sum(temp)  
    
    # Initialize the rows and columns
    rows = np.zeros(total_connections,dtype=np.int64) 
    columns = np.zeros(total_connections,dtype=np.int64) 
    
    # Populates the rows/columns
    k=0
    for indx_current, current_connections in enumerate(Connectivity):
        
        no_connections_current = len(current_connections)
        temp = np.tile(indx_current, [no_connections_current,])
        
        rows[k:k+no_connections_current] = temp
        columns[k:k+no_connections_current] = current_connections
        
        k = k+no_connections_current
        
        
    return rows, columns





def get_points_within(Grid_Nodes_LAT, Grid_Nodes_LON, Target_Nodes_LAT, Target_Nodes_LON, buffer_distance):
    """
    Comments
    ----------
    This function finds all the points in Grid_Nodes that are within buffer_distance of Target_Nodes
    
    Parameters
    ----------
    Grid_Nodes_LAT and LON: list. contains the LAT/LON of each node in the grid.
    Target_Nodes_LAT and LON : list. contains the LAT/LON of each node that will be the center of the buffer
    buffer_distance : [METERS] the radius of the buffer
    
    Returns
    -------
    IDs : Indexes of the nodes in the grid that are within buffer_distance of the Target_Nodes_LAT
    
    """
    
    # Create the dataframes
    Grid_Nodes = pd.DataFrame(data={'LAT':Grid_Nodes_LAT, 'LON':Grid_Nodes_LON}, index=range(len(Grid_Nodes_LON)) )
    Target_Nodes = pd.DataFrame(data={'LAT':Target_Nodes_LAT, 'LON':Target_Nodes_LON}, index=range(len(Target_Nodes_LON)) )
    
    # Define the Coords of each node in a single column
    Grid_Nodes['Coords'] = list(zip(Grid_Nodes['LON'],Grid_Nodes['LAT'])) #the convention here is LON, LAT!!!!!
    Grid_Nodes['Coords'] = Grid_Nodes['Coords'].apply(Point)
    
    Target_Nodes['Coords'] = list(zip(Target_Nodes['LON'],Target_Nodes['LAT'])) #the convention here is LON, LAT!!!!!
    Target_Nodes['Coords'] = Target_Nodes['Coords'].apply(Point)
    
    # Create a geopandas dataframe
    gpd_Grid_Nodes   = gpd.GeoDataFrame(Grid_Nodes,   geometry='Coords', crs = 4326)
    gpd_Target_Nodes = gpd.GeoDataFrame(Target_Nodes, geometry='Coords', crs = 4326)
    
    # Convert the Coords to a different projection system
    gpd_Grid_Nodes['Coords_meters']   = gpd_Grid_Nodes['Coords'].to_crs(3857) #CHANGE THE PROJECTION HERE - The default is 4326 (WGS 84). run temp.crs to get projection info, 7789, 3035
    gpd_Target_Nodes['Coords_meters'] = gpd_Target_Nodes['Coords'].to_crs(3857)
    
    # Create the buffer around each node
    gpd_Target_Nodes['Buffer'] = gpd_Target_Nodes['Coords_meters'].buffer(distance = buffer_distance) #Create the region around each node
    
    # Get the nodes within par_buffer of the selected nodes
    a1 = gpd.GeoDataFrame(gpd_Target_Nodes['Buffer'],     geometry='Buffer')
    a2 = gpd.GeoDataFrame(gpd_Grid_Nodes['Coords_meters'],geometry='Coords_meters')
    
    ptWithin = gpd.tools.sjoin(left_df=a2,right_df=a1 ) #for each point index in the points, it stores the polygon index containing that point
    
    # Extract the final dataframe with the unique liste of indexes
    IDs = np.unique(ptWithin.index.values)
    
    return IDs
    
    








#                A-Star  Functions
# --------------------------------------------------
# --------------------------------------------------


# Define Design Values and calculations
# --------------------------------





# Define the Heuristic function
# --------------------------------

def get_heuristic_value(target_node_id,Coords):
    
    """
    This Function commutes the heuristic value between all nodes in Coords and node target_node_id (also in Coords) 
    
    Comments
    The definition of the heuristic value changes according the requiremetns of the problem
        
    Inputs
    target_node_id: the id of the target node
    Coords: n_nodes X 2 array with the [LAT, LON] values of all nodes in the grid
    
    Outputs
    heuristic_value: an array of all heuristic values from each node to the taget node
    
    """
    
    target_Coords = Coords[target_node_id][None,:] #Convert to 1 x 2 array to be compatible with the nd2nd_distance function
    heuristic_value = nd2nd_distance(target_Coords,Coords)
    
    return heuristic_value


# Define the g_n function
# --------------------------------

def get_gn_value(g_n, Coords, n_current, n_next):  
    """
    Parameters
    ----------
    g_n : list of g_n values for all nodes
    n_current : current node id in teh A* algorithm
    n_next : node connected to the n_current

    Returns
    -------
    g_value : value of the g function to reach from the start node the n_next node

    """
    g_value = g_n[n_current] + nd2nd_distance(Coords[n_current,:][None,:],Coords[n_next,:][None,:])
    
    return g_value



# Define the Cost function
# --------------------------------

def get_cost_value(h_n,g_n,open_list):
    """

    Parameters
    ----------
    h_n : heuristic function [n_nodes x 1]
    g_n : g_n function [n_nodes x 1]
    open_list : Current list of nodes that are considered for next step

    Returns
    -------
    f_n : Cost function value for each of the nodes in the open list [n_nodes_open x 1]
    n_current: the node with the smallest cost function

    """
    f_n = g_n[open_list] + h_n[open_list]
    n_current = open_list[f_n.argmin()]
    
    return f_n, n_current







#  A* Function
# --------------------------------

def A_star_Simple(Start_LAT, Start_LON, End_LAT, End_LON, Nodes, Connectivity):
    """

    Parameters
    ----------
    Start_LAT : float. Latitude of the start location
    Start_LON : float. Longitude of the start location
    End_LAT : float. Latitude of the end location.
    End_LON : float. Longitude of the end location
    Nodes : dataframe. columns:LAT, LON, Dimensions: [n x 2], n is the number of nodes in the grid
    Connectivity : List. len:[n]. each element contains the ids of the nodes that are connected to
    
    Returns
    -------
    reconst_path : list. sequence of nodes ids to reach the end location
    Nodes_Path : dataframe. columns:LAT, LON, Dimensions: [n1 x 2], n1 is the number of nodes in the path
    
    Credits
    ------------
    The A* algorithm is based on https://stackabuse.com/courses/graphs-in-python-theory-and-implementation/lessons/a-star-search-algorithm/
    
    """
    
    
    
    # --------------------------------------------------
    # --------------------------------------------------
    #                Initial Definitions
    # --------------------------------------------------
    # --------------------------------------------------
    
    # Define Start Coords
    # --------------------------------
    target_coords_start = np.asarray([Start_LAT, Start_LON])
    
    # Define End Coords
    # --------------------------------
    target_coords_end =   np.asarray([End_LAT, End_LON])
    
    
    # Get the Start and End nodes 
    # --------------------------------
    COORDS = np.vstack((Nodes['LAT'].values,Nodes['LON'].values)).T
    node_start_id, Coords_Start = find_nearest_point(target_coords_start,COORDS)
    node_end_id, Coords_End =     find_nearest_point(target_coords_end,COORDS)
    
    
    # --------------------------------------------------
    # --------------------------------------------------
    #                 Find Shortest Path
    # --------------------------------------------------
    # --------------------------------------------------
    
    
    # Initialize the Heuristic Function (Compute the heuristic distance of all nodes to the target)
    # --------------------------------
    h_n = get_heuristic_value(node_end_id,COORDS)
    
    
    # Initialize the G Function (Distance from the start to the current node - For initialization just give a large value)
    # --------------------------------
    g_n = np.ones(h_n.shape)*100000000
    
    
    
    # Initialize the lists
    # --------------------------------
    open_list = []
    closed_list =  []
    parent_list = [0]*h_n.shape[0]  #this list will keep track of the parent of each node
    
    # Start with the start node
    open_list.append(node_start_id)
    g_n[node_start_id] = 0
    parent_list[node_start_id] = node_start_id
    
    
    
    
    # A Star Loop
    # ---------------------------------------------------------------
    flag = 0
    counter = 1
    
    # # Loop until you find the end
    print('Initiate A* Search Algorithm')
    while len(open_list) > 0:
    # for ind in range(50):
        
        # Find the node in the open list with the smaller value of f = g + h (cost function)
        f_n, n_current = get_cost_value(h_n,g_n,open_list)
        
        if n_current == None:
            print('Path does not exist or none found!')
            break
        
        # If the selected node is the end node then construct the final path
        if n_current == node_end_id:
            print('Path Found!')
            
            reconst_path = []
            n_check = n_current
            while parent_list[n_check] != node_start_id:
                reconst_path.append(n_check)
                n_check = parent_list[n_check]
            reconst_path.append(n_check)
            reconst_path.append(node_start_id)
            
            
            # Print how many moves it took
            print('Moves taken: ' +str(counter))
            
            break
            # reconst_path.reverse()
            
            
        # Get the connected nodes to the selected node and check if they need to be added directly as children of this node
        current_neighbors = Connectivity[n_current]
        for n_next in current_neighbors:
            if n_next not in open_list and n_next not in closed_list:
                open_list.append(n_next)
                parent_list[n_next] = n_current
                g_n[n_next] = get_gn_value(g_n, COORDS, n_current, n_next)
                
            else:                                                                                 #if n_next is already in open or closed list( i.e. has been visited again), check if the current node n_current is a better candidate to be its parent.
                if g_n[n_next]>get_gn_value(g_n, COORDS, n_current, n_next):     #i.e. if the distance though node n is smaller compared to the previous path that leads to this node, then change the parent to n_current
                    g_n[n_next]=get_gn_value(g_n, COORDS, n_current, n_next)
                    parent_list[n_next] = n_current
                    
                    if n_next in closed_list:                                                      #If it was in the closed list and we found a better route, add it to the open list  
                        closed_list.remove(n_next)
                        open_list.append(n_next)
            
        
        
        # Remove n from the open list and add it to the closed list. 
        open_list.remove(n_current)
        closed_list.append(n_current)
        
        # Update Counter
        counter +=1
        
    
    
    # # Get the list of nodes for the path
    # # ---------------------------------------------------------------
    Nodes_Path_LAT = np.array(Nodes['LAT'].iloc[reconst_path])
    Nodes_Path_LON = np.array(Nodes['LON'].iloc[reconst_path])
    Nodes_Path = np.vstack((Nodes_Path_LAT,Nodes_Path_LON)).T
    
    
    # # Convert to Dataframe
    # # ---------------------------------------------------------------
    Nodes_Path = pd.DataFrame(Nodes_Path, columns=['LAT', 'LON'])
    
    
    
    
    return reconst_path, Nodes_Path




#         Compute Voyage Metrics
# --------------------------------


def get_Voyage_Metrics(Nodes_Path, dict_Parameters):
    """
    Parameters
    -------------
    Nodes_Path : n x 2 array of the LAT, LON of the nodes that define the shortest path
    dict_Parameters: dictionary that contains the following values V_design[km], Power[kW], SFC[kg/kWh], Fuel_Cost_ton[dollars/ton]
    
    Returns
    -------------
    Voyage_KPIs: a dictionary that contains: Distance_km[km], Time_h[h], Fuel_Cons_kg[kg], Fuel_Cost[dollars]
    
    """
    # Extract the values from the dictionary
    V_design = dict_Parameters['V_design']
    Power = dict_Parameters['Power']
    SFC = dict_Parameters['SFC']
    Fuel_Cost_ton = dict_Parameters['Fuel_Cost_ton']
    
    
    # Compute the distance traveled in this path
    d_Short_Path = 0.
    for indx in range(Nodes_Path.shape[0]-1):
        pt1 = Nodes_Path[indx,:][None,:]
        pt2 = Nodes_Path[indx+1,:][None,:]
        d_current = nd2nd_distance(pt1,pt2)[0][0]
        
        d_Short_Path = d_Short_Path + d_current
    
    # Get Time
    Voyage_time = d_Short_Path/V_design
    
    # Get Fuel Consumption
    Fuel_Cons_kg = Power*Voyage_time*SFC
    
    # Compute fuel cost
    FC = Fuel_Cons_kg * Fuel_Cost_ton/1000
    
    
    # Add everything to a dictionary    
    Voyage_KPIs = {'Distance_km':d_Short_Path,
                      'Time_h':Voyage_time,
                      'Fuel_Cons_kg':Fuel_Cons_kg,
                      "Fuel_Cost": FC,
                      }
    
    return Voyage_KPIs




#                Master Function
# --------------------------------------------------
# --------------------------------------------------

def main_voyage_planner(Start_LAT, Start_LON, End_LAT, End_LON, Nodes, Connectivity, dict_Parameters ):
    '''
    This function takes all the inputs and returs 3 list that contain all information necessary for the shortest path. 

    Parameters
    ----------
    Start_LAT : float, value of the Latitude for the start location
    Start_LON : float, value of the longitude for the start location
    End_LAT :  float, value of the latitude for the end location
    End_LON : float, value of the longitude for the end location
    dict_Parameters : a dictionary that contains the input values required for the computations of the shortest path and its KPIs
    Dir_MapGrid : directory of the csv files Map_Nodes and Map_Connectivity that are required for the A* algorithm

    Returns
    -------
    reconst_path : IDs of the nodes that denote the shortest path
    Nodes_Path : n x 2 array of the LAT, LON of the nodes that define the shortest path
    dict_KPIs : a dictionary that contains all the KPIs for this voyage
    '''
    
    
    # Run the A*
    # --------------------------------------------------
    reconst_path, Nodes_Path = A_star_Simple(Start_LAT, Start_LON, End_LAT, End_LON, Nodes, Connectivity)
    
    # Get the Voyage Metrics
    # --------------------------------------------------
    dict_Voyage_metrics = get_Voyage_Metrics(Nodes_Path.values, dict_Parameters)
    
    # Get the final voyage KPIs
    # -----------------------------------------------------
    dict_KPIs = dict_Voyage_metrics
        
    
    return reconst_path, Nodes_Path, dict_KPIs







# ------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------



if __name__ == "__main__":
    
       
    # Problem Inputs
    # --------------------------------------------------
    Start_LAT = 44.65 
    Start_LON = 37.85
    End_LAT = 45.65
    End_LON = 13.70
    
    # ME Design Values
    
    ME_Power = np.asarray([ 8000.,])        #[kW]
    ME_SFC =   np.asarray([ 170., ])/1000.  #[kg/kWh]
    V_design = np.asarray([ 10., ])         #[knots]
    par_Fuel_Cost_ton = np.array([860.])    #[USD/ton]
    
    
    dict_Parameters = {'Power':ME_Power[0],
                      'SFC':ME_SFC[0],
                      'V_design':V_design[0]*1.852,
                      "Fuel_Cost_ton": par_Fuel_Cost_ton[0],
                      }
    
    # Load the locations of the Nodes
    # --------------------------------------------------
    Nodes = pd.read_csv("Data\\Map_Nodes.csv",delimiter=',')
    
    # Load the connectivity list 
    # --------------------------------------------------
    with open("Data\\Map_Connectivity.csv", "r",newline='') as file:
        Connectivity = list(csv.reader(file, delimiter=","))
    file.close()
    
    # Initialy it reads them as strings so must be converted to ints
    for indx in range(len(Connectivity)):    # Skip the header row and convert first values to integers
        Connectivity[indx] = [int(i) for i in Connectivity[indx]]
    
    # Compute Route
    # --------------------------------------------------
    # Run A*
    reconst_path, Nodes_Path = A_star_Simple(Start_LAT, Start_LON, End_LAT, End_LON, Nodes, Connectivity)
    
    # Compute KPIs
    KPIs = get_Voyage_Metrics(Nodes_Path.values, dict_Parameters)


    # Plot Map
    # --------------------------------------------------
    folium_map = folium.Map(tiles="cartodbpositron", location = [43., 13.468777], zoom_start=3.5)
    
    folium.Marker(
      location=[Start_LAT, Start_LON],popup='Start',).add_to(folium_map)
    
    folium.Marker(
      location=[End_LAT,End_LON],popup='End',).add_to(folium_map)
    
    Coords = np.vstack((Nodes['LAT'].values,Nodes['LON'].values)).T
    for indx in range(len(reconst_path)-1):
        pt1 = Coords[reconst_path[indx],:][None,:][0]
        pt2 = Coords[reconst_path[indx+1],:][None,:][0]
        
        folium.PolyLine((pt1,pt2),
                        color= 'green',
                        weight=10,
                        opacity=0.8).add_to(folium_map)
        
    folium_map.save("world_map.html")
    webbrowser.open("world_map.html")



































