


# --------------------------------------------------
# --------------------------------------------------
#                       PACKAGES
# --------------------------------------------------
# -------------------------------------------------

import numpy as np
import streamlit as st
import pandas as pd
import csv
import folium
from streamlit_folium import st_folium
import voyage_plan

# cd C:\Users\ggFri\Desktop\Streamlit-VoyagePlan
# streamlit run main_streamlit_app.py


# --------------------------------------------------
# --------------------------------------------------
#                    FUNCTIONS
# --------------------------------------------------
# --------------------------------------------------







def head(title):
    
    temp = """ <h1 style='text-align: center; margin-bottom: -35px;'>""" + title+ """  </h1> """
    st.markdown(temp, unsafe_allow_html=True    )
    
    return


@st.cache_data
def load_map_data():
    
    # Load the Nodes
    Nodes = pd.read_csv("Data/Map_Nodes.csv",delimiter=',')

    # Load the connectivity list 
    fileName_Connectivity = "Data/Map_Connectivity.csv"
    with open(fileName_Connectivity, "r",newline='') as file:
        Connectivity = list(csv.reader(file, delimiter=","))
    file.close()
    
    for indx in range(len(Connectivity)):    # Skip the header row and convert first values to integers
        Connectivity[indx] = [int(i) for i in Connectivity[indx]]
        
    st.session_state.Nodes = Nodes
    st.session_state.Connectivity = Connectivity
    
    
    return Nodes, Connectivity



def display_map():
    
    folium_map = folium.Map(tiles="cartodbpositron", location = [43., 13.468777], zoom_start=3.5)
    
    
    #Plot start Location
    if st.session_state.plot_start_marker:
        folium.Marker(
          location=[st.session_state.Start_LAT, st.session_state.Start_LON],
          popup='Start',).add_to(folium_map)
    
    # Plot End Location
    if st.session_state.plot_end_marker:
        folium.Marker(
          location=[st.session_state.End_LAT, st.session_state.End_LON],
          popup='End',).add_to(folium_map)
    
    
    # Plot the computed route
    if st.session_state.route_computed:
        
        reconst_path = st.session_state.reconst_path
        Coords = st.session_state.Coords
        
        for indx in range(len(reconst_path)-1):
            pt1 = Coords[reconst_path[indx],:][None,:][0]
            pt2 = Coords[reconst_path[indx+1],:][None,:][0]
            
            folium.PolyLine((pt1,pt2),
                            color= 'green',
                            weight=10,
                            opacity=0.8).add_to(folium_map)


    
    # Plot the map
    st_map = st_folium(folium_map, width = 700, height = 450)
    

    
    return folium_map, st_map




def do_start_button(state):
    
    # Start by disabling the button
    st.session_state["button_disable_select_start"] = state
    
    # Check if compute route button should be activated
    if st.session_state.button_disable_select_start and  st.session_state.button_disable_select_end:
        st.session_state.button_disable_compute_route = False
        
    # Put the values into session_states
    st.session_state.Start_LAT = Start_LAT
    st.session_state.Start_LON = Start_LON
         
    
    st.session_state.plot_start_marker = True
    
    return






def do_end_button(state):
    
    # Start by disabling the button
    st.session_state["button_disable_select_end"] = state
    
    # Check if compute route button should be activated
    if st.session_state.button_disable_select_start and  st.session_state.button_disable_select_end:
        st.session_state.button_disable_compute_route = False
    
    # Put the values into session_states
    st.session_state.End_LAT = End_LAT
    st.session_state.End_LON = End_LON
         
    
    st.session_state.plot_end_marker = True
    
    
    
    return



def do_compute_route():
    
    # Load Map Values
    Nodes, Connectivity = load_map_data()
    
    # Extract Coords of the nodes
    Coords = np.vstack((Nodes['LAT'].values,Nodes['LON'].values)).T
    
    # Compute the Path
    reconst_path, Nodes_Path = voyage_plan.A_star_Simple(st.session_state.Start_LAT, st.session_state.Start_LON, 
                                                                st.session_state.End_LAT, st.session_state.End_LON,
                                                                Nodes, Connectivity)
    
    # Store Session States
    st.session_state.Coords= Coords
    st.session_state.Nodes_Path= Nodes_Path
    st.session_state.reconst_path= reconst_path
    st.session_state.Nodes= Nodes
    st.session_state.Connectivity= Connectivity
    
    st.session_state.route_computed = True

    return





# --------------------------------------------------
# --------------------------------------------------
#                    Main Body
# --------------------------------------------------
# --------------------------------------------------


# Create all state variables
# -----------------------------------------

# Create a state variable
if 'Power' not in st.session_state:
    st.session_state.Power = 8000.
if 'SFC' not in st.session_state:
    st.session_state.SFC = 170./1000.
if 'V_design' not in st.session_state:
    st.session_state.V_design = 14.
if 'Fuel_Cost_ton' not in st.session_state:
    st.session_state.Fuel_Cost_ton = 850.
if 'Start_LAT' not in st.session_state:
    st.session_state.Start_LAT = 44.65
if 'Start_LON' not in st.session_state:
    st.session_state.Start_LON = 37.85
if 'End_LAT' not in st.session_state:
    st.session_state.End_LAT = 45.65
if 'End_LON' not in st.session_state:
    st.session_state.End_LON = 13.70

if 'Distance_km' not in st.session_state:
    st.session_state.Distance_km = 0
if 'Time_h' not in st.session_state:
    st.session_state.Time_h = 0
if 'Fuel_Cons_kg' not in st.session_state:
    st.session_state.Fuel_Cons_kg = 0
if 'Fuel_Cost' not in st.session_state:
    st.session_state.Fuel_Cost = 0





if 'button_disable_select_start' not in st.session_state:
    st.session_state.button_disable_select_start = False
if 'button_disable_select_end' not in st.session_state:
    st.session_state.button_disable_select_end = False

if 'plot_start_marker' not in st.session_state:
    st.session_state.plot_start_marker = False
if 'plot_end_marker' not in st.session_state:
    st.session_state.plot_end_marker = False
if 'plot_route' not in st.session_state:
    st.session_state.plot_route = False


if 'route_computed' not in st.session_state:
    st.session_state.route_computed = False



if 'button_disable_compute_route' not in st.session_state:
    st.session_state.button_disable_compute_route = True




# Show the headtitle
# ---------------------------------------
app_title = 'Marine Vessel Voyage Planner'
head(app_title)


    
# Load the Map Data
# ---------------------------------------
load_map_data()





# Show the sidebar
# ---------------------------------------
st.sidebar.title("Navigation Parameters")

st.session_state.Power = st.sidebar.number_input("Main Engine Power [kW]", min_value=1., max_value=15000., value=st.session_state.Power)
st.session_state.SFOC = st.sidebar.number_input("Specific Fuel Oil Consumption [kg/kWh]", min_value=1./1000., max_value=250./1000., value=st.session_state.SFC)
st.session_state.V_design = st.sidebar.number_input("Sailing Speed [knots]", min_value=1., max_value=20., value=st.session_state.V_design)
st.session_state.Fuel_Cost_ton = st.sidebar.number_input("Fuel Cost [$/tonne]", min_value=1., max_value=3000., value=st.session_state.Fuel_Cost_ton)




# Select Start and end locations
# ---------------------------------------
st.subheader = ('Select Start and End Coordinates')
col1, col2 = st.columns(2)

with col1:
    # Provide Coords
    Start_LAT = st.number_input("Start Latitude", min_value=29., max_value=58., value = st.session_state.Start_LAT)
    Start_LON = st.number_input("Start Longitude", min_value=-12., max_value=45., value = st.session_state.Start_LON)

with col2:
    # Provide Coords
    End_LAT = st.number_input("End Latitude", min_value=29., max_value=58., value = st.session_state.End_LAT)
    End_LON = st.number_input("End Longitude", min_value=-12., max_value=45., value = st.session_state.End_LON)



# Display Map
# ---------------------------------------
st.subheader = ('Map Display')
f_map, st_map = display_map()


# Add select buttons
# ---------------------------------------
with col1:
    # Select Point
    Select_Start = st.button('Select Start Location', key='but_Start_Location', on_click=do_start_button, args=(True,), disabled = st.session_state.button_disable_select_start)

with col2:
    # Slect Point
    Select_End = st.button('Select End Location', key='but_End_Location', on_click=do_end_button, args=(True,), disabled = st.session_state.button_disable_select_end)



# Create a dictionary with the parameters from the side bar
# ---------------------------------------
dict_Parameters = {'Power':st.session_state.Power,
                   'SFC':st.session_state.SFC,
                   'V_design':st.session_state.V_design*1.852,
                   "Fuel_Cost_ton": st.session_state.Fuel_Cost_ton,
                  }

st.session_state.dict_Parameters=dict_Parameters





# Compute the route
# ---------------------------------------
Compute_Route = st.button('Compute Route', key='but_compute_route', on_click=do_compute_route, disabled = st.session_state.button_disable_compute_route)



if st.session_state.route_computed:
    # Compute the KPIs
    dict_Voyage_metrics = voyage_plan.get_Voyage_Metrics(st.session_state.Nodes_Path.values, st.session_state.dict_Parameters)
    
    st.session_state.Distance_km = dict_Voyage_metrics['Distance_km']
    st.session_state.Time_h = dict_Voyage_metrics['Time_h']
    st.session_state.Fuel_Cons_kg = dict_Voyage_metrics['Fuel_Cons_kg']
    st.session_state.Fuel_Cost = dict_Voyage_metrics['Fuel_Cost']
    


# Add KPI Metrics
# ---------------------------------------
st.subheader = ('KPIs')
KPI_col1, KPI_col2, KPI_col3, KPI_col4 = st.columns(4)
with KPI_col1:
    st.metric('Distance Traveled [km]', round(st.session_state.Distance_km, 2) )
with KPI_col2:
    st.metric('Voyage Duration [hours]', round(st.session_state.Time_h, 2))
with KPI_col3:
    st.metric('Fuel Consumed [ton]',  round(st.session_state.Fuel_Cons_kg/1000, 2))
with KPI_col4:
    st.metric('Fuel Cost [$]', round(st.session_state.Fuel_Cost, 2) )













