import yt
import numpy as np

ds = yt.load("ScalarField_000000.3d.hdf5")
field = 'Ham'
weight = 'cell_volume'
ad = ds.all_data()

average_ham = ad.quantities.weighted_average_quantity(field, weight)
print("Average Ham is ", average_ham)

K = np.sqrt(1.5*np.abs(average_ham))

print("Value for K should be set to ", K)