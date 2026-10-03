# -*- coding: utf-8 -*-

"""Main module."""
import json
import pathlib
import logging

import requests
import pandas as pd
import pytz
import tqdm
import dateutil
from datetime import datetime

from logging import getLogger
import numpy as np
import os
from datetime import datetime
from matplotlib.dates import date2num, num2date
from matplotlib.ticker import Formatter,FuncFormatter, NullLocator, FixedLocator
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure, show,savefig,clabel,contour,setp,gcf,getp,gca,close
# from utilities import orderOfMagnitude,generate_ticks,datestring2num,stdev_back_from_log
# from plotfunctions import simpleplot
from wellmetadata import WellMetaData

logger = getLogger(__name__)


BASE_URL = "https://waterwebservices.rijkswaterstaat.nl/"
BASE_URL = "https://publiek.broservices.nl"


url = BASE_URL + "/gm/gmw/v1/bro-ids"

url = "https://publiek.broservices.nl//gm//gmw//v1//bro-ids"

# url = 'https://www.w3schools.com/python/demopage.php'
# myobj = {'somekey': 'somevalue'}
myobj = {}

# x = requests.post(url, json = myobj)

broId = "GLD000000014660"
broId = "GMW000000019472"
broId = "GMW000000019472"
broId = "GMW000000026244"

	
url = f'https://publiek.broservices.nl//gm/gmw//v1//objects//{broId}'
# url = "https://publiek.broservices.nl//gm/gmw//v1//objects//GMW000000026244"
	

x = requests.get(url)

#print the response text (the content of the requested file):

# print(x.text)

#https://realpython.com/python-xml-parser/
from bs4 import BeautifulSoup
soup = BeautifulSoup(x.text)
print(soup.prettify())

# # from https://realpython.com/python-xml-parser/
# soup.find("stop", {"stop-color": "gold"})

# result1 = soup.find(text=lambda x: "pvc" in x).parent
result2 = soup.find_all("nitgcode", limit=1)

# print('result1', result1 )
print('result2', result2 )


abs_path = os.path.dirname(os.path.abspath(__file__))
abs_path_splitted = abs_path.split('\\')
path_to_parent_folder_elements = abs_path_splitted[:-1]
path_to_parent_folder = os.path.join(*path_to_parent_folder_elements)
splitted=path_to_parent_folder.split(':')#trick  to  repair  path
path_to_parent_folder = splitted[0]+':\\'+splitted[1]#trick  to  repair  path
path_to_resources_folder = os.path.join(path_to_parent_folder,'resources')



if path_to_parent_folder is not None:
    path_to_export_file =  os.path.join(path_to_resources_folder,broId+'.xml')
else:
    curdir = os.getcwd()
    path_to_export_file = curdir+'//'+broId+'.xml'

F=open(path_to_export_file,'w')
F.write(soup.prettify())
F.close()



