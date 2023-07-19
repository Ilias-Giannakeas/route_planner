# Marine Vessel Voyage Planer Using the A* Algorithm

## Introduction
This is the Marine Vessel Voyage Planer repository! This project aims to provide a comprehensive implementation of marine vessel navigation using the powerful A* algorithm. Whether you're a seafarer, a student, or just curious about maritime navigation, this codebase will be a valuable resource for you.

## Features
* The codebase leverages the A* algorithm, a widely-used pathfinding algorithm, to navigate marine vessels efficiently and effectively. 
* An interactive application has also been built to visualize vessel routes, providing a clear and intuitive representation of the navigational path.
* Land avoidance is achieved through a predefined grid of points. Users have the flexibility to provide custom grids, allowing them to define specific areas to avoid or include in the navigation path.
* Main voyage KPIs are provided for the journey. These include Distance, Time, Fuel Consumption and Fuel Cost.  

## Installation
To get started and plan your voyage, follow these steps:
* Clone the Repository: git clone  https://github.com/Ilias-Giannakeas/route_planner.git
* Install the required dependencies: 
```python 
pip install -r requirements.txt
```
    
## Use
In the script voyage_plan.py, an example is already included. 

The possible movement of the vessel is defined in the files Data\Nodes.csv and Data\Connectivity.csv. The first contains the coordinates (Lat, Lon) of a grid of points while the second describes how the grid points are connected with one another. Bot the grid and the connectivity between all points has been defined in such a way that major land masses are avoided.  

When specifying starting and end coordinates, the algorithm will search for the closest point in the available grid and pick that as the start/end.

Curently there is a limitation in the possible start and end location coordinates. The included grid contains mainly the mediterranean region and part of Europe. The latitudes are limited in [29., 58] while longitudes are constrained in [-12, 45].

Some reference values for main engine power, sailing speed, fuel consumption and fuel cost have been provided in the example to facilitate the computation of the KPIs. Make sure to change them according to your needs.


## Acknowledgements

 - For intuition behind the A* algorithm visit: [A* Search Algorithm](https://stackabuse.com/courses/graphs-in-python-theory-and-implementation/lessons/a-star-search-algorithm/) . The A* loop was based on this source. 



## Contact
If you have any questions, suggestions, or just want to say hi, feel free to reach out to us at iliasgiann@hotmail.com.

Happy sailing !
