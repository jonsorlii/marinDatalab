from Hollenbach import hollenbach
import numpy as np 
import matplotlib.pyplot as plt
from grid import Grid
import datetime 

def main():
    #Ship dimensions
    L = 233.0 #[m]
    Lwl = 227.953 #[m]
    Los = 234.526 #[m]
    B = 32.2 #[m]
    TF = 11.0 #[m]
    TA = 11.0 #[m]
    CB = 0.6039
    S = 9203.9 #[m^2]
    Dp = 1.0 #[m]
    Nrud = 1.0
    NBrac = 0
    NBoss = 0
    NThr = 0
    Vvec = np.array([13.0, 23.0]) #[m/s]

    date = datetime.datetime.now()
    dateNow = str(date)[0:19]


    newGrid = Grid("u_wind.csv", "v_wind.csv", "latitude.csv", "longitude.csv", "wave_data.csv")
    m = newGrid.construct_map()
    newGrid.construct_forces_grid()
    solutionList, resistanceList = newGrid.dijkstra(m,(newGrid.long[0][58], newGrid.lat[1][0]),(newGrid.long[0][56],newGrid.lat[len(newGrid.lat) -1][0] ))
    x,y = newGrid.setup_meshgrid(m)
    m.quiver(x,y, newGrid.uwind, newGrid.vwind)
    print(solutionList)
    plt.figure()
    mapPlot = newGrid.construct_map()

    plt.figure()
    waveMap = newGrid.construct_map()
    x_zeros = np.zeros([len(newGrid.long[0]),len(newGrid.lat)], int)
    waveMap.quiver(x,y, x_zeros, newGrid.waves,color = 'r',headwidth = 2, headlength = 3, scale = None  )

    plt.figure()
    path_x = []
    path_y = []
    for coordinate in solutionList: 
        long = coordinate[0]
        lat = coordinate[1]
        x,y = m(long,lat)
        path_x.append(x)
        path_y.append(y)
    solutionMap = newGrid.construct_map()
    plt.plot(path_x,path_y, color = 'green', linewidth = 2)

    plt.figure()
    solutionWithWind = newGrid.construct_map()
    x,y = newGrid.setup_meshgrid(solutionWithWind)
    solutionWithWind.quiver(x,y, newGrid.uwind, newGrid.vwind)
    plt.plot(path_x,path_y, color = 'green', linewidth = 2)

    plt.figure()
    hollenbach_resistance = hollenbach(Vvec,L,Lwl,Los,B,TF,TA,CB,S,Dp,Nrud, NBrac, NBoss, NThr)
    plot_forces_vel1 = []
    plot_forces_vel2 = []
    for resistance in resistanceList:
        plot_forces_vel1.append(resistance + hollenbach_resistance[1][0]/1000) 
        plot_forces_vel2.append(resistance + hollenbach_resistance[1][1]/1000) 
        
    x_values = np.linspace(0,len(plot_forces_vel1),len(plot_forces_vel1))

    plt.subplot(1,2,1)
    plt.plot(x_values, plot_forces_vel1)
    plt.title("Resistance with V = " + str(Vvec[0]))

    plt.subplot(1,2,2)
    plt.plot(x_values, plot_forces_vel2)
    plt.title("Resistance with V = " + str(Vvec[1]))
   

    plt.figure()
    print(plot_forces_vel1, Vvec[0])
    velocity_1 = np.asarray(plot_forces_vel1)*float(Vvec[0])
    velocity_2 = np.asarray(plot_forces_vel2)*float(Vvec[1])

    plt.subplot(1,2,1)
    plt.plot(x_values, velocity_1)
    plt.title("Effective power P_E with V = " + str(Vvec[0]))

    plt.subplot(1,2,2)
    plt.plot(x_values, velocity_2)
    plt.title("Effective power P_E with V = " + str(Vvec[1]))
    
    plt.show()


main()