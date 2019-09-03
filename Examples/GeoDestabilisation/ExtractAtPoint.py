import yt
from yt import derived_field
import numpy as np
import time
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
yt.enable_parallelism()

#set up
input_files = '/mnt/extraspace/kclough/GeoInflation/GeoFieldplot_*'
output_filename = 'FieldsVsTime.txt'
grid_location =  np.array([0,0,0])
#grid_location =  ds[0].domain_right_edge/2
#variable = 'phi2'

from mpi4py import MPI
comm = MPI.COMM_WORLD

if(comm.rank==0) :
  datafile=open(output_filename,'a')
  datafile.write("# t_sim    scale_factor    phi1    phi2    K    \n")
  datafile.close()

ds = yt.load(input_files)
for i in ds:
  # find the center of the BH
  #value, location = i.find_min("chi")
  #center = [float(location[0]), float(location[1]), float(location[2])]
  #print ("New center ", center)
  time = float(i.current_time)
  my_point = i.point(grid_location)
  phi1 = float(my_point['phi1'])
  phi2 = float(my_point['phi2'])
  chi = float(my_point['chi'])
  K = float(my_point['K'])
  scale_factor = 1.0/np.sqrt(chi)

#  if (time > 100.0) :
#    break;

  print("scale factor ", scale_factor)
  print("chi ", phi2)

  if(comm.rank==0) :
    datafile=open(output_filename,'a')
    datafile.write("%f    %f    %f     %f   %f \n" % (time, scale_factor, phi1, phi2, K))
    datafile.close()
