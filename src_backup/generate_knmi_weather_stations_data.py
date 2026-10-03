# -*- coding: utf-8 -*-

import json
import pathlib
import numpy as np
import os
from osgeo.osr import SpatialReference, CoordinateTransformation




def generate_knmi_weather_stations_data(path, outputfilename = None):
    # Define the Rijksdriehoek projection system (EPSG 28992)
    epsg28992 = SpatialReference()
    epsg28992.ImportFromEPSG(28992)
     
    # correct the towgs84 # chob asks: ok, but where do these 7 parameters come from?
    epsg28992.SetTOWGS84(565.237,50.0087,465.658,-0.406857,0.350733,-1.87035,4.0812)
     
    # Define the wgs84 system (EPSG 4326)
    epsg4326 = SpatialReference()
    epsg4326.ImportFromEPSG(4326)
     
    rd2latlon = CoordinateTransformation(epsg28992, epsg4326)
    latlon2rd = CoordinateTransformation(epsg4326, epsg28992)
    
    
    F = open(path)
    all_lines = F.readlines()
    F.close()
    
    Nstations = len(all_lines)-1
    
    
    # option 1: use a structured array           
    weerstations_daggegevens_structured_array= np.empty(Nstations, dtype=[
        ('name', 'S20'),
        ('lon', 'f4'),
        ('lat', 'f4'),
        ('alt', 'f4'),
        ('code', 'S10'),
        ('X', 'f4'),
        ('Y', 'f4'),
        ('Z', 'f4')])
    
    # option 1: use a dict
    weerstations_daggegevens_dict = {}
    
    for i in range(0,Nstations):
        dataline_splitted = all_lines[i+1].rstrip().split(';')
        name = dataline_splitted[4]
        lon = float(dataline_splitted[1])
        lat = float(dataline_splitted[2])
        alt = float(dataline_splitted[3])
        STN = dataline_splitted[0]
        XYZ = latlon2rd.TransformPoint(lat,lon,alt)
        X = XYZ[0]
        Y = XYZ[1]
        Z = XYZ[2]
        
        
        # option 1
        weerstations_daggegevens_structured_array[i][0] = name
        weerstations_daggegevens_structured_array[i][1] = lon
        weerstations_daggegevens_structured_array[i][2] = lat
        weerstations_daggegevens_structured_array[i][3] = alt
        weerstations_daggegevens_structured_array[i][4] = STN
        weerstations_daggegevens_structured_array[i][5] = X
        weerstations_daggegevens_structured_array[i][6] = Y
        weerstations_daggegevens_structured_array[i][7] = Z
        
        
        #option 2
        weerstations_daggegevens_dict[name] = {}
        weerstations_daggegevens_dict[name]['name'] = name
        weerstations_daggegevens_dict[name]['lon'] = lon
        weerstations_daggegevens_dict[name]['lat'] = lat
        weerstations_daggegevens_dict[name]['alt'] = alt   
        weerstations_daggegevens_dict[name]['STN'] = STN
        weerstations_daggegevens_dict[name]['X'] = X
        weerstations_daggegevens_dict[name]['Y'] = Y  
        weerstations_daggegevens_dict[name]['Z'] = Z    
    
    weerstations_gegevens_json = json.dumps(weerstations_daggegevens_dict, indent = 4)
    
    
    if outputfilename is None:
        outputfilename = "test.json"
    else :
        outputfilename = outputfilename
    with open(outputfilename, "w") as outfile: 
        json.dump(weerstations_daggegevens_dict, outfile, indent = 4)
    
    
    
    #Just a technical note: do not forget to decode strings when retrieving 
    # strings from a structured array
    encoding = 'utf-8'
    weerstations_daggegevens_dict = {}
    for i in range(0,len(weerstations_daggegevens_structured_array)):
        name = weerstations_daggegevens_structured_array[i][0].decode(encoding)
        lon = weerstations_daggegevens_structured_array[i][1]
        lat = weerstations_daggegevens_structured_array[i][2]
        alt = weerstations_daggegevens_structured_array[i][3]
        STN = weerstations_daggegevens_structured_array[i][4].decode(encoding) # STN is the station number
        X = weerstations_daggegevens_structured_array[i][5]
        Y = weerstations_daggegevens_structured_array[i][5]
        Z = weerstations_daggegevens_structured_array[i][5] 
        
        
        weerstations_daggegevens_dict[name] = {}
        weerstations_daggegevens_dict[name]['name'] = name
        weerstations_daggegevens_dict[name]['lon'] = lon
        weerstations_daggegevens_dict[name]['lat'] = lat
        weerstations_daggegevens_dict[name]['alt'] = alt   
        weerstations_daggegevens_dict[name]['STN'] = STN
        weerstations_daggegevens_dict[name]['X'] = X
        weerstations_daggegevens_dict[name]['Y'] = Y  
        weerstations_daggegevens_dict[name]['Z'] = Z
        
    # with open("test.json", "w") as outfile: 
    #     json.dump(weerstations_daggegevens_dict, outfile, indent = 4)
    
    return weerstations_gegevens_json



if __name__ == "__main__":
    
    abs_path = os.path.dirname(os.path.abspath(__file__))
    # abs_path_splitted = abs_path.split('\\')
    # path_to_parent_folder_elements = abs_path_splitted[:-1]
    # path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
    # splitted=path_to_parent_folder.split(':')#trick  to  repair  path
    # path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
    # path_to_resources_folder = os.path.join(abs_path,'resources')
    curdir = os.getcwd()
    inputpath = os.path.join(curdir,"weerstations_daggegevens_input.txt")
    outputpath = os.path.join(curdir,"weerstations_daggegevens_input.json")
    weerstations_daggegevens_json = generate_knmi_weather_stations_data(inputpath, outputfilename = outputpath)

    inputpath = os.path.join(curdir,"weerstations_uurgegevens_input.txt")
    outputpath = os.path.join(curdir,"weerstations_uurgegevens_input.json")
    weerstations_daggegevens_json = generate_knmi_weather_stations_data(inputpath, outputfilename = outputpath)    




# # Opmerking: door stationsverplaatsingen en veranderingen in waarneemmethodieken zijn deze tijdreeksen van uurwaarden mogelijk inhomogeen! Dat betekent dat deze reeks van gemeten waarden niet geschikt is voor trendanalyse. Voor studies naar klimaatverandering verwijzen we naar de gehomogeniseerde dagreeksen <http://www.knmi.nl/nederland-nu/klimatologie/daggegevens> of de Centraal Nederland Temperatuur <http://www.knmi.nl/kennis-en-datacentrum/achtergrond/centraal-nederland-temperatuur-cnt>.
# # 
# # SOURCE: ROYAL NETHERLANDS METEOROLOGICAL INSTITUTE (KNMI)
# # Comment: These time series are inhomogeneous because of station relocations and changes in observation techniques. As a result these series are not suitable for trend analysis. For climate change studies we refer to the homogenized series of daily data <http://www.knmi.nl/nederland-nu/klimatologie/daggegevens> or the Central Netherlands Temperature <http://www.knmi.nl/kennis-en-datacentrum/achtergrond/centraal-nederland-temperatuur-cnt>.
# # 
# # STN         LON(east)   LAT(north)  ALT(m)      NAME
# # 209         4.518       52.465      0.00        IJmond      
# # 210         4.430       52.171      -0.20       Valkenburg Zh
# # 215         4.437       52.141      -1.10       Voorschoten 
# # 225         4.555       52.463      4.40        IJmuiden    
# # 235         4.781       52.928      1.20        De Kooy     
# # 240         4.790       52.318      -3.30       Schiphol    
# # 242         4.921       53.241      10.80       Vlieland    
# # 248         5.174       52.634      0.80        Wijdenes    
# # 249         4.979       52.644      -2.40       Berkhout    
# # 251         5.346       53.392      0.70        Hoorn Terschelling
# # 257         4.603       52.506      8.50        Wijk aan Zee
# # 258         5.401       52.649      7.30        Houtribdijk 
# # 260         5.180       52.100      1.90        De Bilt     
# # 265         5.274       52.130      13.90       Soesterberg 
# # 267         5.384       52.898      -1.30       Stavoren    
# # 269         5.520       52.458      -3.70       Lelystad    
# # 270         5.752       53.224      1.20        Leeuwarden  
# # 273         5.888       52.703      -3.30       Marknesse   
# # 275         5.873       52.056      48.20       Deelen      
# # 277         6.200       53.413      2.90        Lauwersoog  
# # 278         6.259       52.435      3.60        Heino       
# # 279         6.574       52.750      15.80       Hoogeveen   
# # 280         6.585       53.125      5.20        Eelde       
# # 283         6.657       52.069      29.10       Hupsel      
# # 285         6.399       53.575      0.00        Huibertgat  
# # 286         7.150       53.196      -0.20       Nieuw Beerta
# # 290         6.891       52.274      34.80       Twenthe     
# # 308         3.379       51.381      0.00        Cadzand     
# # 310         3.596       51.442      8.00        Vlissingen  
# # 311         3.672       51.379      0.00        Hoofdplaat  
# # 312         3.622       51.768      0.00        Oosterschelde
# # 313         3.242       51.505      0.00        Vlakte van De Raan
# # 315         3.998       51.447      0.00        Hansweert   
# # 316         3.694       51.657      0.00        Schaar      
# # 319         3.861       51.226      1.70        Westdorpe   
# # 323         3.884       51.527      1.40        Wilhelminadorp
# # 324         4.006       51.596      0.00        Stavenisse  
# # 330         4.122       51.992      11.90       Hoek van Holland
# # 331         4.193       51.480      0.00        Tholen      
# # 340         4.342       51.449      19.20       Woensdrecht 
# # 343         4.313       51.893      3.50        Rotterdam Geulhaven
# # 344         4.447       51.962      -4.30       Rotterdam   
# # 348         4.926       51.970      -0.70       Cabauw Mast 
# # 350         4.936       51.566      14.90       Gilze-Rijen 
# # 356         5.146       51.859      0.70        Herwijnen   
# # 370         5.377       51.451      22.60       Eindhoven   
# # 375         5.707       51.659      22.00       Volkel      
# # 377         5.763       51.198      30.00       Ell         
# # 380         5.762       50.906      114.30      Maastricht  
# # 391         6.197       51.498      19.50       Arcen      


