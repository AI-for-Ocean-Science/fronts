""" Module to generate a prototype for the SST-SSH U-Net with
LLC4320 data """


import os
import numpy as np

import pandas
import xarray
import json
import h5py

#from ulmo.llc import extract 
#from ulmo.llc import uniform
#from ulmo import io as ulmo_io
#from ulmo.analysis import evaluate as ulmo_evaluate
#from ulmo.nflows import nn 
#from ulmo.preproc import plotting as pp_plotting

from fronts.llc import table
from fronts.llc import extract
from fronts.preproc import process
from fronts.tables import catalog
from fronts import io as fronts_io

from IPython import embed

local_out_path = os.path.join(os.getenv('OS_OGCM'), 'LLC', 'Fronts')
super_tbl_file = os.path.join(local_out_path, 'LLC4320_SST144_SSS40_super_test.parquet')
super_preproc_file = os.path.join(local_out_path, 'LLC4320_SST144_SSS40_super.h5')

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
        tbl_file = super_tbl_file

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

    # Add days
    llc_table = table.add_days(llc_table, dti, outfile=tbl_file, to_s3=False)

    # Add UIDs
    _ = table.add_uid(llc_table, outfile=tbl_file, to_s3=False)

    print(f"Wrote: {tbl_file} with {len(llc_table)} unique cutouts.")
    print("All done with init")

    if debug:
        embed(header='71 of llc4320_sst_ssh_proto.py')

def preproc_super(extract_file:str, debug:bool=False):

    outfile = super_preproc_file 

    # Load the table
    llc_table = pandas.read_parquet(super_tbl_file)

    # Debug?
    if debug:
        llc_table = llc_table.iloc[:100].copy()
        outfile = os.path.join(local_out_path, 'LLC4320_SST144_SSS40_super_test.h5')

    # Extract dict
    # Load JSON
    #with open(preproc_file, 'r') as infile:
    #    preproc_dict = json.load(infile)


    #extract_dict = {'fields': ['SSS'],
    #extract_dict = {'fields': ['Divb2'],
    extract_dict = {'fields': ['SST','SSS','Divb2'],
             'field_size': 64,
             'pdicts': 
                 {
                     'SST': 
                        {
                        'fixed_km': 144.,
                        'field_size': 64,
                        "quality_thresh": 2,
                        "nrepeat": 1,
                        "downscale": False,
                        "inpaint": False,
                        "median": False,
                        "only_inpaint": False
                        }
                    ,
                     'SSS': 
                        {
                        'fixed_km': 144.,
                        'field_size': 64,
                        'smooth_km': 40., 
                        "quality_thresh": 2,
                        "nrepeat": 1,
                        "downscale": False,
                        "inpaint": False,
                        "de_mean": False,
                        "median": False,
                        "only_inpaint": False
                        }
                    ,
                     'Divb2': 
                        {
                        'fixed_km': 144.,
                        'field_size': 64,
                        'dx': 144./64,
                        }
                    ,
                 }
             }

    # Prep LLC Table
    llc_table = process.prep_table_for_preproc(
        llc_table, extract_file, field_size=
        (extract_dict['field_size'], extract_dict['field_size']))


    # Open HDF5 file
    f = h5py.File(outfile, 'w')
    for field in extract_dict['fields']:
        # Preprocess
        llc_table, success, pp_fields, meta = extract.preproc_field(
            llc_table, field, extract_dict['pdicts'][field],
            fixed_km=extract_dict['pdicts'][field]['fixed_km'],
            n_cores=10, dlocal=True,
            test_failures=False,
            test_process=False)

        # Write data
        pp_fields = np.array(pp_fields).astype(np.float32)
        f.create_dataset(field, data=pp_fields)

        # Add meta
        for key in meta.keys():
            llc_table[field+key] = meta[key]

        # Deal with failures in Table
        if np.any(~success):
            print(f"Failed to preprocess some {field} fields")
            fail = np.where(~success)[0]
            llc_table.loc[fail,'pp_type'] = -999
        # Good
        good = success & (llc_table['pp_type'] != -999)
        llc_table.loc[np.where(good)[0],'pp_type'] = 0
        #
        del pp_fields

    # Close
    f.close()    

    # Write table
    assert catalog.vet_main_table(llc_table)
    if not debug:
        fronts_io.write_main_table(llc_table, tbl_file, to_s3=False)
    else:
        tbl_file = os.path.join(local_out_path, 'LLC4320_SST144_SSS40_super_test.parquet')
        fronts_io.write_main_table(llc_table, tbl_file, to_s3=False)
        embed(header='preproc_super 118')


def main(flg:str):
    flg= int(flg)

    # Generate the LLC Table
    if flg == 1:
        generate_super_table()#debug=True, plot=True)

    # Generate the Super Preproc File
    if flg == 2:
        #preproc_super('dummy_file.json', debug=True)
        preproc_super('llc4320_sst144_sss40_proto.json')#, debug=True)

# Command line execution
if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1:

        flg = 1 # Generate super table
        flg = 2 # Preproc super table

        #flg += 2 ** 1  # 2 -- Extract
        #flg += 2 ** 2  # 4 -- Evaluate (with noise)
        #flg += 2 ** 3  # 8 -- Evaluate (without noise)
    else:
        flg = sys.argv[1]

    main(flg)