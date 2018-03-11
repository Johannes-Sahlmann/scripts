#!/usr/bin/env python
"""Script to download besancon model including kinematics using the astroquery interface.
Gaia DR2 parallax and proper motion sky-averaged uncertainties are derived using pygaia.
Random errors corresponding to Gaia DR2 expectations can be added to the parallaxes and distances of the
model output.

Authors
-------
    Johannes Sahlmann
    with inputs from Laura Watkins and Mark Fardal


Usage
-----

    python besancon_gaia_dr2.py youremail@email.com --options


Dependencies
------------

    astroquery
    pygaia
    astropy

"""

from __future__ import print_function
import sys
import os
import numpy as np
import argparse

from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
from astropy import units as u

from astroquery.besancon import Besancon

from pygaia.photometry import transformations
from pygaia.errors.astrometric import properMotionErrorSkyAvg, parallaxErrorSkyAvg


def vminusi_to_gaiacolor(v_minus_i):
    """Return Gaia BP-RP color given V-I colour
    formula is from Jordi+10

    :param v_minus_i: float
    :return: float
    """
    return -0.0660 + 1.2061 * v_minus_i - 0.0614 * v_minus_i**2 + 0.0041 * v_minus_i**3


def vminusi_to_gminusbp(v_minus_i):
    """
    Return Gaia G-BP color given V-I colour
    formula is from Jordi+10, Table 3

    :param v_minus_i:
    :return:
    """

    return 0.0387 - 0.4191 * v_minus_i - 0.0736 * v_minus_i ** 2 + 0.0040 * v_minus_i ** 3


def vminusi_to_gminusrp(v_minus_i):
    """
    Return Gaia G-BP color given V-I colour
    formula is from Jordi+10, Table 3

    :param v_minus_i:
    :return:
    """

    return -0.0274 + 0.7870 * v_minus_i - 0.1350 * v_minus_i ** 2 + 0.0082 * v_minus_i ** 3


def get_gaia_bp_and_rp(v, v_minus_i):

    gminusv = transformations.gminvFromVmini(v_minus_i)
    g = gminusv + v

    gminusbp = vminusi_to_gminusbp(v_minus_i)
    gminusrp = vminusi_to_gminusrp(v_minus_i)

    bp = g - gminusbp
    rp = g - gminusrp

    return bp, rp


def main(argv):
    parser = argparse.ArgumentParser(
        description="This Script downloads a besancon model including include_kinematics using the astroquery interface. Gaia DR2 parallax and proper motion sky-averaged uncertainties are derived using pygaia. Random errors corresponding to Gaia DR2 expectations can be added to the parallaxes and distances of the model output.")
    parser.add_argument('email', help='Your valid email address.')
    parser.add_argument('--ra_deg', type=float, default=10., help='RA in degrees')
    parser.add_argument('--dec_deg', type=float, default=0., help='Dec in degrees')
    parser.add_argument('--field_size_squaredegree', type=float, default=0.05, help='Field size in squaredegrees.')
    parser.add_argument('--data_dir', type=str, default='', help='Path to directory for table saving')
    parser.add_argument('--overwrite', type=bool, default=False,
                        help='Overwrite download from besancon server.')
    parser.add_argument('--add_random_gaia_errors', type=bool, default=False,
                        help='Add random errors to parallax and PM assuming Gaia DR2 performances.')
    parser.add_argument('--random_seed', type=int, default=None,
                        help='numpy random seed for error simulation. Allows for repeatability.')
    # parser.add_argument('--argument_dict', type=, default=None,
    #                     help='Dictionary of arguments passed to Besancon.query. See Besancon.query documentation for allowed values and format.')
    parser.add_argument('--vmag_upper_limit', type=float, default=None,
                        help='V mag upper limit passed to Besancon.query. Default is 18')


    args = parser.parse_args(argv)

    email = args.email

    # central pointing and field size
    ra_deg = args.ra_deg
    dec_deg = args.dec_deg
    field_size_squaredegree = args.field_size_squaredegree

    data_dir = args.data_dir
    overwrite = args.overwrite
    random_seed = args.random_seed
    add_random_gaia_errors = args.add_random_gaia_errors
    # argument_dict = args.argument_dict
    # vmag_limits = args.vmag_limits
    vmag_upper_limit = args.vmag_upper_limit

    galactic_latitude_deg=None
    galactic_longitude_deg=None

    if (galactic_latitude_deg is None) and (galactic_longitude_deg is None):
        equatorial_coordinates = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
        galactic_longitude_deg = equatorial_coordinates.galactic.l.value
        galactic_latitude_deg = equatorial_coordinates.galactic.b.value

    include_kinematics = True
    pm_in_lb = False # flag to obtain proper motions in galactic coordinates instead of equatorial

    # include_kinematics options
    kwd = {}
    if (include_kinematics):
        kwd['cinem'] = 1
        kwd[ 'sc'] = [[0, 0, 0]] * 15 # changed default *9 to *15 to allow for kinematics to be downloaded
        # note: if you leave klee at zero you get no kinematics
        if (pm_in_lb):
            kwd['klee'] = 2  #l/b
        else:
            kwd['klee'] = 1  #ra/dec

    # merge two dictionaries of arguments
    # if argument_dict is not None:
    #     kwd = {**kwd, **argument_dict}
    # if vmag_limits is not None:
    #     kwd['mag_limits'] = {'V': vmag_limits}

    # execute the query and save into astropy table on disk
    table_file = os.path.join(data_dir, 'besancon_table_{:3.2f}_{:3.2f}_{:3.2f}.csv'.format(galactic_longitude_deg, galactic_latitude_deg, field_size_squaredegree))

    if (not os.path.isfile(table_file)) | (overwrite):
        if vmag_upper_limit is not None:
            # print(vmag_upper_limit)
            mag_limits = {'V': (10, vmag_upper_limit)}
            besancon_model = Besancon.query(glon=galactic_longitude_deg, glat=galactic_latitude_deg, email=email, retrieve_file=False, area=field_size_squaredegree, mag_limits=mag_limits, **kwd)
            # print('Faintest V magnitude = {}'.format(np.max(besancon_model['V'])))
        else:
            besancon_model = Besancon.query(glon=galactic_longitude_deg, glat=galactic_latitude_deg, email=email, retrieve_file=False, area=field_size_squaredegree, **kwd)
        besancon_model.write(table_file, format='ascii.fixed_width', delimiter=',', bookend=False)
    else:
        besancon_model = Table.read(table_file, format='ascii.basic', delimiter=',')

    dr2_table_file = table_file.replace('besancon_', 'besancon_gaia_dr2_')

    # compute Gaia magnitude and color
    besancon_model['G-V'] = transformations.gminvFromVmini(besancon_model['V-I'])
    besancon_model['Gmag'] = besancon_model['G-V'] + besancon_model['V']
    besancon_model['BP-RP'] = vminusi_to_gaiacolor(besancon_model['V-I'])
    besancon_model['BP'], besancon_model['RP'] = get_gaia_bp_and_rp(besancon_model['V'], besancon_model['V-I'])
    # print(besancon_model['BP']-besancon_model['RP'] - besancon_model['BP-RP'])

    # compute expected DR2 errors, parallaxErrorSkyAvg is in microarcsec
    dr2_offset = -(60. - 22.) / 12.
    besancon_model['properMotionErrorSkyAvg_dr2_ra_maspyr'], besancon_model['properMotionErrorSkyAvg_dr2_dec_maspyr'] = np.array(properMotionErrorSkyAvg(besancon_model['Gmag'].data, besancon_model['V-I'].data, extension=dr2_offset))/1000.
    besancon_model['parallaxErrorSkyAvg_dr2_mas'] = np.array(parallaxErrorSkyAvg(besancon_model['Gmag'], besancon_model['V-I'], extension=dr2_offset))/1000.

    besancon_model['parallax_mas'] = 1./besancon_model['Dist'] # 'Dist' is in kpc

    # coordinates are constant and equal to the central pointing in the standard output
    galactic_coords = SkyCoord('galactic', l=besancon_model['l'], b=besancon_model['b'], unit='deg')
    besancon_model['ra_deg'] = galactic_coords.icrs.ra.value
    besancon_model['dec_deg'] = galactic_coords.icrs.dec.value

    if pm_in_lb is False:
        # mux,muy are in arcsec/century
        besancon_model['pm_ra*_maspyr'] = besancon_model['mux'] * 1000/100.
        besancon_model['pm_dec_maspyr'] = besancon_model['muy'] * 1000/100.

    if add_random_gaia_errors:
        # add random errors to parallax and PM assuming Gaia DR2 performances
        if random_seed is not None:
            np.random.seed(random_seed)
        besancon_model['parallax_mas'] += (np.random.normal(0., 1, len(besancon_model)) * besancon_model['parallaxErrorSkyAvg_dr2_mas'])
        if random_seed is not None:
            np.random.seed(random_seed+1)
        besancon_model['pm_ra*_maspyr'] += (np.random.normal(0., 1, len(besancon_model)) * besancon_model['properMotionErrorSkyAvg_dr2_ra_maspyr'])
        if random_seed is not None:
            np.random.seed(random_seed+2)
        besancon_model['pm_dec_maspyr'] += (np.random.normal(0., 1, len(besancon_model)) * besancon_model['properMotionErrorSkyAvg_dr2_dec_maspyr'])

    besancon_model.meta={}
    besancon_model.write(dr2_table_file, format='ascii.fixed_width', delimiter=',', bookend=False)
    print('Besancon model and Gaia DR2 performance simulation written to {}'.format(dr2_table_file))

if __name__ == '__main__':
    main(sys.argv[1:])
    sys.exit(0)
