"""Script to download a selection of Gaia sources from the Gaia archive that are also present in a
crossmatch table, here 2MASS.

Authors
-------
    Johannes Sahlmann

"""
from collections import OrderedDict
import os

from astroquery.gaia import Gaia
import astropy.units as u
from astropy.table import Table

overwrite = True
# overwrite = True

field_radius = 10*u.arcmin

field_number = 1
ra_deg = 240.634
dec_deg = -11.012


data = OrderedDict()
data['gaia'] = OrderedDict()
data['tmass'] = OrderedDict()
data['xmatch'] = OrderedDict()

data['gaia']['query'] = """SELECT * FROM gaiadr2.gaia_source AS gaia
                        WHERE 1=CONTAINS(POINT('ICRS',gaia.ra,gaia.dec), CIRCLE('ICRS',{}, {}, {}))
                        """.format(ra_deg, dec_deg, field_radius.to(u.deg).value)

data['tmass']['query'] = """SELECT * FROM gaiadr1.tmass_original_valid AS tmass
                        WHERE 1=CONTAINS(POINT('ICRS',tmass.ra,tmass.dec), CIRCLE('ICRS',{}, {}, {}))
                        """.format(ra_deg, dec_deg, field_radius.to(u.deg).value)

data['xmatch']['query'] = """SELECT * from
            (SELECT gaia.*
            FROM gaiadr2.gaia_source AS gaia
            WHERE 1=CONTAINS(POINT('ICRS',gaia.ra,gaia.dec), CIRCLE('ICRS',{}, {}, {})))
            AS field
            INNER JOIN gaiadr2.tmass_best_neighbour AS xmatch
                ON field.source_id = xmatch.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass
                ON tmass.tmass_oid = xmatch.tmass_oid
        """.format(ra_deg, dec_deg, field_radius.to(u.deg).value)


for key in data.keys():
    output_file = os.path.join(os.environ['HOME'], '{}_query_field{:02d}_{}.vot'.format(key, field_number, field_radius).replace(' ',''))


    if (not os.path.isfile(output_file)) or (overwrite):
        job = Gaia.launch_job_async(data[key]['query'], dump_to_file=True, output_file=output_file)
        table = job.get_results()
    else:
        table = Table.read(output_file)
    print('Retrieved {} sources for catalog {}'.format(len(table), key))
    # table.pprint()
