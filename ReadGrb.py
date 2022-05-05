from __future__ import print_function
import requests
import pygrib
import csv 
import datetime

from eccodes import *

date = datetime.datetime.now()
dateNow = str(date)[0:19]
newTime = str(int(dateNow[11:13])+1)+dateNow[13:19]

"""
Filename of Grb files
"""
filename = 'meps_weatherapi_west_norway.grb'
filename_2 = 'meps_wave_west_norway.grb'

"""
Urls to fetch from API
"""
url = 'https://api.met.no/weatherapi/gribfiles/1.1/?area=skagerrak&content=weather'
wave_url = 'https://api.met.no/weatherapi/gribfiles/1.1/?area=skagerrak&content=waves'

"""
Defining user-agent
"""
sitename = {
    'user-agent' : 'My user Agent 1.0',
    'From' : 'jonsorl@hotmial.com'
}

class API:
    """
    Class that calls the API and writes the weather data to a grib file.
    """
    def __init__(self, url, filename, sitename):
        self.url = url
        self.sitename = sitename
        self.filename = filename
        self.writeToGrib(self.filename, self.url)

    def writeToGrib(self, filename, url):
        response = requests.get(url, allow_redirects=True)
        open(filename, 'wb').write(response.content)


class writeToCSV():
    """
    Constructs the data into arrays and then writes it to a CSV file. 
    """
    def __init__(self, API):
        self.grbs = pygrib.open(API.filename)

    def get_u_wind(self):
        wind = self.grbs.select(name='10 metre U wind component')
        for i in range(len(wind)): 
            if str(wind[i].analDate)[0:10] >= dateNow[0:10] and str(wind[i].analDate)[11:19] > dateNow[11:19] and str(wind[i].analDate)[11:19] < newTime:
                return (wind[i].values)
        pass
    
    def get_v_wind(self): 
        wind = self.grbs.select(name='10 metre V wind component')
        for i in range(len(wind)): 
            if str(wind[i].analDate)[0:10] >= dateNow[0:10] and str(wind[i].analDate)[11:19] > dateNow[11:19] and str(wind[i].analDate)[11:19] < newTime:
                return (wind[i].values)
        pass
    
    def get_wave_data(self):
        wave_list = self.grbs.select(name = 'Experimental product')
        for i in range(len(wave_list)): 
            if str(wave_list[i].analDate)[0:10] >= dateNow[0:10] and str(wave_list[i].analDate)[11:19] > dateNow[11:19] and str(wave_list[i].analDate)[11:19] < newTime:
                return (wave_list[i].values) 
        pass
    def get_long_lat(self):
        grb = self.grbs.select(name='10 metre V wind component')[0]
        lat,long = grb.latlons()
        return lat,long
        
    def write_data_to_csv(self, data, filename):
        with open(filename, 'w') as f:
            writer = csv.writer(f, delimiter=',')
            for row in data: 
                if filename == "latitude.csv" or filename == "longitude.csv":
                    newRow = []
                    for number in row: 
                        newNum = round(float(number),2)
                        newRow.append(newNum)
                    writer.writerow(newRow)
                else:
                    writer.writerow(row) 

def main(): 
    """
    Get wind and waves
    """
    wind_reader = API(url, filename, sitename)
    wind_writer = writeToCSV(wind_reader)

    wind_writer.write_data_to_csv(wind_writer.get_u_wind(),'u_wind.csv')
    wind_writer.write_data_to_csv(wind_writer.get_v_wind(),'v_wind.csv')
    wind_writer.write_data_to_csv(wind_writer.get_long_lat()[0],'latitude.csv')
    wind_writer.write_data_to_csv(wind_writer.get_long_lat()[1],'longitude.csv')

    wave_reader = API(wave_url, filename_2, sitename)
    wave_writer = writeToCSV(wave_reader)
    wave_writer.write_data_to_csv(wave_writer.get_wave_data(), 'wave_data.csv')
main()
