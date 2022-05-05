import numpy as np 
from mpl_toolkits.basemap import Basemap
import csv
from calculations import ShipCalculations


class Grid(): 
    """
    Constructs the grid 
    """
    def __init__(self, uwind, vwind, latitude, longitude, waves):
        self.uwind = self.read_csv(uwind)
        self.vwind = self.read_csv(vwind) 
        self.lat = self.read_csv(latitude)
        self.long = self.read_csv(longitude)
        self.waves = self.read_csv(waves)
        self.grid = {}

        self.latMin=min(self.lat)[0] 
        self.latMax=max(self.lat)[0] 
        self.lonMin=min(self.long)[0] 
        self.lonMax=max(self.long)[-1] 

    def read_csv(self, filename): 
        """
        Reads from CSV files
        """
        newArr = []
        try:        
            with open(filename) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    tempArr = []
                    for number in row: 
                        if number == "--":
                            tempArr.append(float(0))
                        else: 
                            if isinstance(number, float):
                                tempArr.append(number)
                            else: 
                                tempArr.append(float(number))
                    newArr.append(tempArr)
            return newArr
        except: 
            print("Could not read from CSV file ", filename) 
            
    def construct_map(self):
        """
        Constructs the basemap projected map
        """
        latMin = self.latMin - 0.025
        latMax = self.latMax + 0.025
        lonMin = self.lonMin - 0.025
        lonMax = self.lonMax + 0.025
            
        m = Basemap(projection='merc',llcrnrlat=latMin,urcrnrlat=latMax,\
            llcrnrlon=lonMin,urcrnrlon=lonMax, area_thresh = 0.0001, lat_ts=20,resolution='h')
    
        m.drawcoastlines()
        m.fillcontinents(color='coral',lake_color='aqua')
        m.drawmapboundary(fill_color='aqua')

        return m

    def construct_forces_grid(self):
        """
        Adds wind and waves to nodes, so that it is easily accesable
        """
        for i in range(len(self.lat)):
            for j in range(len(self.long[i])): 
                self.grid[(self.lat[i][j],self.long[i][j])] = ForcesNode(self.uwind[i][j], self.vwind[i][j], self.waves[i][j])
               
    def get_node_values(self, x, y): 
        """
        Gets node values
        """
        for i in range(len(self.lat)):
            for j in range(len(self.long)): 
                if self.grid[i][j].x == x and self.grid[i][j].y == y:
                    return self.grid[i][j]
        return "Node could not be found"  

    def setup_meshgrid(self, m):
        """
        Sets up meshgrid for quiver plot, need to make long/lat into mapable projections
        """ 
        x_values = []
        y_values = []
        for i in range(len(self.lat)):
            temp_x = []
            temp_y = []
            for j in range(len(self.long[i])): 
                xpt,ypt = m(self.long[i][j],self.lat[i][j])
                temp_x.append(xpt)
                temp_y.append(ypt)
            x_values.append(temp_x)
            y_values.append(temp_y)
        return x_values,y_values
    
    def dijkstra(self, grid, start_node, end_node): 
        """
        Solution algorithm, takes in a map, start and end-node. 
        """
        x_step = abs(round(abs(self.long[0][1])- abs(self.long[0][0]), 2))
        y_step = abs(round(abs(self.lat[1][0])- abs(self.lat[0][0]),2))
        open = []
        closed = []

        #Adds startnode to open list
        startNode = GridNode(start_node, None)
        open.append(startNode)

        #While this is true, run
        while len(open) > 0:
            #Sorts the open list
            open.sort()
            currentNode = open.pop(0)
            closed.append(currentNode)
            
            #Checks if we reached the end
            if currentNode.pos == end_node:
                path = []
                resistance = []
                while currentNode.pos != start_node:
                    path.append(currentNode.pos)
                    resistance.append(currentNode.resistance)
                    currentNode = currentNode.parent
                #Returns both the path and resistance list
                return path[::-1], resistance

            (y,x) = currentNode.pos

            #Make a list of neighbouring nodes
            neighbours = self.get_neighbour_nodes((x,y),x_step,y_step)
            for next in neighbours:
                x,y = grid(next[0], next[1])
                isLand = grid.is_land(x,y)
                #Checks if current node is land
                if not isLand: 
                    neighbour = GridNode(next, currentNode)
                    if neighbour in closed:
                        continue
                    #Calculates cost of travel to next node
                    forces = self.calculate_resistance((y,x), next, end_node, x_step, y_step)
                    neighbour.h = forces[0]*0.15 +forces[1]*0.3+forces[2]*(forces[0]+forces[1])/2*0.55
                    neighbour.resistance = forces[0]+forces[1]
                    if (self.add_to_open(open, neighbour) == True): 
                        open.append(neighbour)

    def add_to_open(self, open, neighbour):
        for node in open:
            if (neighbour == node and neighbour.h >= node.h):
                return False 
        return True
                


    def get_neighbour_nodes(self, pos, x_step, y_step): 
        """
        Checks neighbouring nodes in all eight directions, which means
        vertical, horizontal and diagonal. If the current position does
        not satisfy the requirment to be within lon/lats min and max, it is
        not added to the neighbour list. 
        """
        neighbours = []
        x = pos[1]
        y = pos[0]
        if x-x_step >= self.lonMin:
            neighbours.append((round(x-x_step,2), y))
        if x-x_step >= self.lonMin and y-y_step >= self.latMin:
            neighbours.append((round(x-x_step,2), round(y-y_step,2)))
        if x-x_step >= self.lonMin and y+y_step <= self.latMax:
            neighbours.append((round(x-x_step,2), round(y+y_step,2)))
        if y+y_step <= self.latMax:
            neighbours.append((x, round(y+y_step,2)))
        if y-y_step >= self.latMin:
            neighbours.append((x, round(y-y_step,2)))
        if x+x_step<= self.lonMax:
            neighbours.append((round(x+x_step,2), y))
        if x+x_step <= self.lonMax and y+y_step <= self.latMax:
            neighbours.append((round(x+x_step,2), round(y+y_step,2)))
        if x+x_step <= self.lonMax and y-y_step >= self.latMin:
            neighbours.append((round(x+x_step,2), round(y-y_step,2)))
        return neighbours 

    def calculate_resistance(self, currentNode, neighbourNode, end_node, x_step, y_step): 
        """
        Calculates resistance
        """
        calculations = ShipCalculations()
        
        index_i = neighbourNode[1]
        index_j = neighbourNode[0]

        uwind = self.grid[(index_i,index_j)].uwind
        vwind = self.grid[(index_i,index_j)].vwind
        waveHeight = self.grid[(index_i,index_j)].waveHeight
        
        vec_boat = [neighbourNode[0]-currentNode[0], neighbourNode[1]-currentNode[1]]
        vec_wind = [uwind, vwind]

        wind_magnitude = calculations.get_vec_magnitude(vec_wind)
        angle = calculations.angle_between_vectors(vec_boat, vec_wind)
        #Checks if angle is between +- 45 deg or pi/4 radians
        if angle < np.pi/4 and angle > -np.pi/4: 
            wind_resistance = calculations.wind_resistance(calculations.rho_air,calculations.CD,calculations.AT, wind_magnitude)
        else:
            wind_resistance = 0

        #Wave resistance is not dependent on angle of the boat
        wave_resistance = calculations.STAwave_formulae(calculations.rho_water,calculations.g,waveHeight, calculations.breadth, calculations.waterlineLenght)
        
        #Gets distance to goal
        absolute_distance = self.get_distance_to_destination(neighbourNode, end_node)
        return [wind_resistance, wave_resistance, absolute_distance]


    def get_distance_to_destination(self, kord_1, kord_2):
        if kord_1[1]-kord_2[1] == 0: 
            return (kord_1[0]-kord_2[0])*87
        elif kord_1[0]-kord_2[0] == 0:
            return (kord_1[0]-kord_2[0])*111
        else:
            return np.sqrt((kord_1[0]-kord_2[0])**2*87+(kord_1[0]-kord_2[0])**2*111)
    

class GridNode(): 
    def __init__(self, pos, parent): 
        self.pos = pos
        self.parent = parent
        self.waveResistance = 0
        self.windResistance = 0 
        self.h = 0 
        self.resistance = 0
        
    # Compare nodes
    def __eq__(self, other):
        return self.pos == other.pos
    # Sort nodes
    def __lt__(self, other):
         return self.h < other.h


class ForcesNode(): 
    def __init__(self, uWind, vWind, waves):
        self.uwind = uWind
        self.vwind = vWind
        self.waveHeight = waves



