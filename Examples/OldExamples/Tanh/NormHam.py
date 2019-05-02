import yt
import numpy as np
from yt import derived_field

#yt.enable_parallelism()

@derived_field(name="HamSquared", units="")
def _HamSquared(field, data) :
    value = data["Ham"] * data["Ham"]
    return value

filename = 'HamvsTime.txt'
loc = '/home/fmuia/GRChombo/Examples/StarobinskiPotentialRun/phi04+/Starobinski_phi04_a1_P_0000*'
ds = yt.load(loc)
		
#integral_data = []

volume = ( float(ds[0].domain_right_edge[0]) )**3.0

for i in ds :
    field = 'HamSquared'
    weight = 'cell_volume'
    ad = i.all_data()

    integral = ad.quantities.weighted_average_quantity(field, weight)
    #integral = ad.quantities.total_quantity(field)
    integral = integral * volume

    datafile=open(filename,'a')
    datafile.write("%f    %.10f \n" % (i.current_time, integral))
    datafile.close()
