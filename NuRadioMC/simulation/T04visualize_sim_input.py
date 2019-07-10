from __future__ import absolute_import, division, print_function
import numpy as np
from radiotools import helper as hp
from radiotools import plthelpers as php
from matplotlib import pyplot as plt
from NuRadioMC.utilities import units, plotting
import h5py
import argparse
import json
import time
import os

parser = argparse.ArgumentParser(description='Plot NuRadioMC event list input')
parser.add_argument('inputfilename', type=str,
                    help='path to NuRadioMC hdf5 simulation input')
# parser.add_argument('outputfilename', type=str,
#                     help='name of output file storing the electric field traces at detector positions')
args = parser.parse_args()

filename = os.path.splitext(os.path.basename(args.inputfilename))[0]
dirname = os.path.dirname(args.inputfilename)
plot_folder = os.path.join(dirname, 'plots', filename)
if(not os.path.exists(plot_folder)):
    os.makedirs(plot_folder)

fin = h5py.File(args.inputfilename, 'r')
V = None
if('xmax' in fin.attrs):
    dX = fin.attrs['xmax'] - fin.attrs['xmin']
    dY = fin.attrs['ymax'] - fin.attrs['ymin']
    dZ = fin.attrs['zmax'] - fin.attrs['zmin']
    V = dX * dY * dZ
elif('rmin' in fin.attrs):
    rmin = fin.attrs['rmin']
    rmax = fin.attrs['rmax']
    dZ = fin.attrs['zmax'] - fin.attrs['zmin']
    V = np.pi * (rmax**2 - rmin**2) * dZ

for flavor in np.unique(np.array(fin['flavors'])):
    print("{}: {} events".format(flavor, np.sum(np.array(fin['flavors']) == flavor)))
    

print("mean energy = {:.2g}eV +- {:.2g}eV".format(np.array(fin['energies']).mean(), np.array(fin['energies']).std()))

for ccnc in np.unique(np.array(fin['interaction_type'])):
    print("{}: {} events".format(ccnc, np.sum(np.array(fin['interaction_type']) == ccnc)))

fig, ax = plotting.plot_flavor_ratio(fin['flavors'], fin['interaction_type'])
fig.savefig(os.path.join(plot_folder, 'flavor.png'))

    
# plot vertex distribution
xx = np.array(fin['xx'])
yy = np.array(fin['yy'])
zz = np.array(fin['zz'])
rr = (xx**2 + yy**2)**0.5

fig, ax = php.get_histograms([rr, zz], xlabels=["r [m]", "z [m]"])
fig.savefig(os.path.join(plot_folder, 'vertex_distribution_r_z.png'))

plt.show()


fig, ax = plotting.plot_vertex_distribution(xx, yy, zz, weights=np.ones_like(xx), rmax=rmax, zmin=fin.attrs['zmin'],
                                           trigger_name="")
fig.savefig(os.path.join(plot_folder, 'vertex_distribution.png'), bbox='tight')


# plot incoming direction
zeniths = np.array(fin['zeniths'])
azimuths = np.array(fin['azimuths'])
fig, axs = php.get_histograms([zeniths / units.deg, azimuths / units.deg],
                              bins=[np.arange(0, 181, 10), np.arange(0, 361, 45)],
                              xlabels=['zenith [deg]', 'azimuth [deg]'],
                              stats=False)
fig.suptitle('neutrino direction')
plt.title('incoming direction')
plt.savefig(os.path.join(plot_folder, "simInputIncoming.pdf"))

# plot inelasticity
inelasticity = np.array(fin['inelasticity'])
fig, axs = php.get_histogram(inelasticity,
                             bins=np.logspace(np.log10(0.0001), np.log10(1.0), 50),
                             xlabel='inelasticity', figsize=(6, 6),
                             stats=True)
axs.semilogx(True)
plt.title('inelasticity')
plt.savefig(os.path.join(plot_folder, "simInputInelasticity.pdf"))
plt.show()
