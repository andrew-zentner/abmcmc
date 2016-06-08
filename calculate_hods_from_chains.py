#!/usr/bin/env python
""" This script is executable from the command-line and can be used to generate
the plot of our best-fit HODs.
"""
import argparse
import os
import numpy as np
from astropy.table import Table
from halotools.empirical_models import (AssembiasZheng07Cens,
        AssembiasZheng07Sats, TrivialPhaseSpace, NFWPhaseSpace, HodModelFactory)
from halotools.empirical_models import PrebuiltHodModelFactory
from halotools.sim_manager import CachedHaloCatalog
import h5py

# Determine the directory where the chain files are stored
parser = argparse.ArgumentParser()
parser.add_argument("chain_dirname", type=str,
    help="Absolute path to the directory containing the hdf5 files storing the chains")
parser.add_argument("-num_iterations", type=int, default=50,
    help="Number of times to repeatedly populate mocks to calculate the HOD."
    "Default is 50.")
parser.add_argument("-min_mass", type=int, default=11.5,
    help="Minimum mass to use in the logarithmically spaced mass bins used to calculate the HOD. "
    "Default is 11.5")
parser.add_argument("-max_mass", type=int, default=14.5,
    help="Minimum mass to use in the logarithmically spaced mass bins used to calculate the HOD. "
    "Default is 14.5")
parser.add_argument("-num_mass_bins", type=int, default=20,
    help="Number of logarithmically-spaced mass bins used to calculate the HOD. "
    "Default is 20")
parser.add_argument("-output_dirname", type=str, default='./',
    help="Directory where astropy tables storing the tabulated HODs will be stored. "
    "Default is the directory where the script is located.")

args = parser.parse_args()
chain_dirname = os.path.abspath(args.chain_dirname)
num_iterations = args.num_iterations
min_mass = args.min_mass
max_mass = args.max_mass
num_mass_bins = args.num_mass_bins
output_dirname = os.path.abspath(args.output_dirname)

mass_bins = np.logspace(min_mass, max_mass, num_mass_bins)
mass_bin_midpoints = (mass_bins[1:] + mass_bins[:-1])/2.
################################################################################
# Load the chains into memory and do some post-processing

wp20_ab_chain_fname = os.path.join(chain_dirname, 'wp20.0_ab_chain.hdf5')
wp20_std_chain_fname = os.path.join(chain_dirname, 'wp20.0_standard_chain.hdf5')

wp20_ab_chain = Table.read(wp20_ab_chain_fname, path='data')
wp20_std_chain = Table.read(wp20_std_chain_fname, path='data')


def correct_chi2_column_shape(chain):
    x = chain['chi2'].flatten()
    del chain['chi2']
    chain['chi2'] = x
    return chain
wp20_ab_chain = correct_chi2_column_shape(wp20_ab_chain)
wp20_std_chain = correct_chi2_column_shape(wp20_std_chain)

# pre-sort the chains by chi**2 for convenience
wp20_ab_chain.sort('chi2')
assert wp20_ab_chain['chi2'][0] == wp20_ab_chain['chi2'].min()
wp20_std_chain.sort('chi2')
assert wp20_std_chain['chi2'][0] == wp20_std_chain['chi2'].min()


def retrieve_random_param_dict(chain, delta_chi2_cutoff=1):
    idx_chi2_cuotff = np.searchsorted(chain['chi2'].data, chain['chi2'][0]+delta_chi2_cutoff)
    idx_random = np.random.randint(0, idx_chi2_cuotff)
    assert chain['chi2'][idx_random] <= chain['chi2'].min() + 1
    d = {}
    d['logMmin'] = chain['logMmin'][idx_random]
    d['logM0'] = chain['logM0'][idx_random]
    d['sigma_logM'] = chain['sigmalogM'][idx_random]
    d['logM1'] = chain['logM1'][idx_random]
    d['alpha'] = chain['alpha'][idx_random]
    try:
        d['mean_occupation_centrals_assembias_param1'] = chain['Acen'][idx_random]
        d['mean_occupation_satellites_assembias_param1'] = chain['Asat'][idx_random]
    except KeyError:
        pass
    return d


def update_param_dict(model, chain):
    new_param_dict = retrieve_random_param_dict(chain)
    for key, value in new_param_dict.iteritems():
        model.param_dict[key] = value

################################################################################
# build the two models

std_model = PrebuiltHodModelFactory('zheng07')

centrals_occupation = AssembiasZheng07Cens()
centrals_profile = TrivialPhaseSpace()
satellites_occupation = AssembiasZheng07Sats()
satellites_profile = NFWPhaseSpace()

satellites_occupation._suppress_repeated_param_warning = True

model_dict = ({'centrals_occupation': centrals_occupation,
    'centrals_profile': centrals_profile,
    'satellites_occupation': satellites_occupation,
    'satellites_profile': satellites_profile})
ab_model = HodModelFactory(**model_dict)

################################################################################
# Initially populate both models

halocat = CachedHaloCatalog(simname='bolplanck')
__  = halocat.halo_table['halo_x']

std_model.populate_mock(halocat)
ab_model.populate_mock(halocat)

################################################################################
np.seterr(divide='ignore', invalid='ignore')  # ignore divide by zero in e.g. DD/RR


def calculate_hod(halo_table, galaxy_table):
    halo_counts = np.histogram(halo_table['halo_mvir'].data, mass_bins)[0].astype('f4')

    central_mask = galaxy_table['gal_type'] == 'centrals'
    central_counts = np.histogram(galaxy_table['halo_mvir'][central_mask].data,
        mass_bins)[0].astype('f4')

    satellite_mask = galaxy_table['gal_type'] == 'satellites'
    satellite_counts = np.histogram(galaxy_table['halo_mvir'][satellite_mask].data,
        mass_bins)[0].astype('f4')

    return mass_bins, central_counts/halo_counts, satellite_counts/halo_counts


def calculate_mean_hods(model, chain):
    ncen = np.zeros((num_iterations, len(mass_bins)-1))
    nsat = np.zeros((num_iterations, len(mass_bins)-1))

    for i in range(num_iterations):
        update_param_dict(model, chain)
        model.mock.populate()
        __, ncen[i, :], nsat[i, :] = calculate_hod(model.mock.halo_table, model.mock.galaxy_table)

    mean_ncen = np.mean(ncen, axis=0)
    mean_nsat = np.mean(nsat, axis=0)

    stddev_ncen = np.std(ncen, axis=0)
    stddev_nsat = np.std(nsat, axis=0)

    return mean_ncen, stddev_ncen, mean_nsat, stddev_nsat


################################################################################

def write_result_to_disk(model, chain, output_fname):

    mean_ncen, stddev_ncen, mean_nsat, stddev_nsat = calculate_mean_hods(model, chain)

    hod_table = Table()
    hod_table['mean_ncen'] = mean_ncen
    hod_table['mean_nsat'] = mean_nsat
    hod_table['stddev_nsat'] = stddev_nsat
    hod_table['stddev_ncen'] = stddev_ncen
    hod_table['mass_bin_midpoints'] = mass_bin_midpoints

    try:
        os.remove(output_fname)
    except OSError:
        pass
    finally:
        hod_table.write(output_fname, path='data')

    f = h5py.File(output_fname)
    f.attrs.create('num_iterations', num_iterations)
    f.close()

################################################################################
################################################################################

output_fname_std = os.path.join(output_dirname, 'best_fit_hod_zheng_model.hdf5')
output_fname_ab = os.path.join(output_dirname, 'best_fit_hod_zentner_model.hdf5')

print("...working on calculation of HOD for standard case")
write_result_to_disk(std_model, wp20_std_chain, output_fname_std)
print("...working on calculation of HOD for assembly-biased case")
write_result_to_disk(ab_model, wp20_ab_chain, output_fname_ab)
