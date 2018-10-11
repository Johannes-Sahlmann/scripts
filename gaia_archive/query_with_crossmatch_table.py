"""Script to download a selection of Gaia sources from the Gaia archive that are also present in a
crossmatch table, here 2MASS.

Authors
-------
    Johannes Sahlmann

"""
import os

from astroquery.gaia import Gaia
import astropy.units as u
from astropy.table import Table

field_radius = 5*u.arcmin

field_number = 1
ra_deg = 240.634
dec_deg = -11.012

output_file = os.path.join(os.environ['HOME'], 'gaia_query_field{:02d}.vot'.format(field_number))

query = """SELECT * from
            (SELECT gaia.*
            FROM gaiadr2.gaia_source AS gaia
            WHERE 1=CONTAINS(POINT('ICRS',gaia.ra,gaia.dec), CIRCLE('ICRS',{}, {}, {})))
            AS field
            INNER JOIN gaiadr1.tmass_best_neighbour AS xmatch
                ON field.source_id = xmatch.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass
                ON tmass.tmass_oid = xmatch.tmass_oid
        """.format(ra_deg, dec_deg, field_radius.to(u.deg).value)


if not os.path.isfile(output_file):
    job = Gaia.launch_job_async(query, dump_to_file=True, output_file=output_file)
    table = job.get_results()
else:
    table = Table.read(output_file)
print('Retrieved {} sources'.format(len(table)))
table.pprint()
