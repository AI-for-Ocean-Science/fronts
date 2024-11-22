""" Module to generate a prototype for the SST-SSH U-Net with
LLC4320 data """


import os
import numpy as np

import pandas
import xarray

#from ulmo.llc import extract 
#from ulmo.llc import uniform
#from ulmo import io as ulmo_io
#from ulmo.analysis import evaluate as ulmo_evaluate
#from ulmo.nflows import nn 
#from ulmo.preproc import plotting as pp_plotting

from fronts.llc import table
from fronts import io as fronts_io

from IPython import embed

local_out_path = os.path.join(os.getenv('OS_OGCM'), 'LLC', 'Fronts')

def generate_super_table(debug=False, resol=0.5, plot=False,
               minmax_lat=(-90.,57.), field_size=(64,64),
               local:bool=True):
    """ Get the show started by sampling uniformly
    in space and and time

    Args:
        tbl_file (str): _description_
        debug (bool, optional): _description_. Defaults to True.
        resol (float, optional): 
            Typical separation of images in deg
        plot (bool, optional): Plot the spatial distribution?
        max_lat (float, optional): Restrict on latitude
    """

    if debug:
        tbl_file = os.path.join(local_out_path, 'LLC4320_SST144_SSS40_super_test.parquet')
    else:
        tbl_file = os.path.join(local_out_path, 'LLC4320_SST144_SSS40_super.parquet')

    # Begin 
    llc_table = table.uniform_coords(resol=resol, minmax_lat=minmax_lat,
                               field_size=field_size)#, outfile=tbl_file)
    # Plot
    if plot:
        table.plot_extraction(llc_table, s=1, resol=resol)

    # Temporal sampling
    if debug:
        # Extract 6 days across the full range;  ends of months
        dti = pandas.date_range('2011-09-13', periods=6, freq='2M')
    else:
        # Extract 24 days across the full range;  ends of months; every 2 weeks
        dti = pandas.date_range('2011-09-13', periods=24, freq='2W')
    llc_table = table.add_days(llc_table, dti, outfile=tbl_file)

    print(f"Wrote: {tbl_file} with {len(llc_table)} unique cutouts.")
    print("All done with init")

    if debug:
        embed(header='63 of llc4320_sst_ssh_proto.py')


def main(flg:str):
    flg= int(flg)

    # Generate the LLC Table
    if flg == 1:
        # Debug
        generate_super_table(debug=True, plot=True)

# Command line execution
if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1:

        flg = 1 # Generate super table
        #flg += 2 ** 1  # 2 -- Extract
        #flg += 2 ** 2  # 4 -- Evaluate (with noise)
        #flg += 2 ** 3  # 8 -- Evaluate (without noise)
    else:
        flg = sys.argv[1]

    main(flg)