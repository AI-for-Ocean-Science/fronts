""" Routines to extract cutouts from LLC data """

import numpy as np

import pandas

from functools import partial
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

from fronts.llc import io as llc_io
from fronts.preproc import process
from fronts.tables import catalog

from IPython import embed

def preproc_field(llc_table:pandas.DataFrame, 
                  field:str,
                  pdict:str,
                  field_size=(64,64), 
                  fixed_km=None,
                  n_cores=10,
                  dlocal:bool=True,
                  override_RAM=False,
                  test_failures:bool=False,
                  test_process:bool=False,
                  debug=False):
    """Main routine to extract and pre-process LLC data for later SST analysis
    The llc_table is modified in place (and also returned).

    Args:
        llc_table (pandas.DataFrame): cutout table
        field (str): Field to extract
        preproc_dict (dict): Preprocessing steps. 
        field_size (tuple, optional): Defines cutout shape. Defaults to (64,64).
        fixed_km (float, optional): Require cutout to be this size in km
        n_cores (int, optional): Number of cores for parallel processing. Defaults to 10.
        valid_fraction (float, optional): [description]. Defaults to 1..
        dlocal (bool, optional): Data files are local? Defaults to False.
        override_RAM (bool, optional): Over-ride RAM warning?

    Raises:
        IOError: [description]

    Returns:
        pandas.DataFrame: Modified in place table

    """
    # Load coords?
    if fixed_km is not None:
        coords_ds = llc_io.load_coords()
        R_earth = 6371. # km
        circum = 2 * np.pi* R_earth
        km_deg = circum / 360.
    
    # Setup for parallel
    map_fn = partial(process.preproc_image, pdict=pdict)

    # Setup for dates
    uni_date = np.unique(llc_table.datetime)
    if len(llc_table) > 1000000 and not override_RAM:
        raise IOError("You are likely to exceed the RAM.  Deal")

    # Init
    pp_fields, meta, img_idx, all_sub = [], [], [], []

    # Loop

    for udate in uni_date:
        # Parse filename
        filename = llc_io.grab_llc_datafile(udate, local=dlocal)

        # 
        ds = llc_io.load_llc_ds(filename, local=dlocal)

        # Field
        if field == 'SST':
            data = ds.Theta.values
        elif field == 'SSS':
            data = ds.Salt.values
        else:
            raise IOError(f"Not ready for this field {field}")
        # Parse 
        gd_date = llc_table.datetime == udate
        sub_idx = np.where(gd_date)[0]
        all_sub += sub_idx.tolist()  # These really should be the indices of the Table
        coord_tbl = llc_table[gd_date]

        # Add to table
        llc_table.loc[gd_date, 'filename'] = filename

        # Load up the cutouts
        fields, smooth_pixs = [], []
        for r, c in zip(coord_tbl.row, coord_tbl.col):
            if fixed_km is None:
                dr = field_size[0]
                dc = field_size[1]
            else:
                dlat_km = (coords_ds.lat.data[r+1,c]-coords_ds.lat.data[r,c]) * km_deg
                dr = int(np.round(fixed_km / dlat_km))
                dc = dr
            # Deal with smoothing
            if 'smooth_km' in pdict.keys():
                smooth_pix = int(np.round(pdict['smooth_km'] / dlat_km))
                pad = 2*smooth_pix
                #
                use_r = r - pad
                dr += 2*pad
                use_c = c - pad
                dc += 2*pad
                smooth_pixs.append(smooth_pix)
            else:
                use_r, use_c = r, c
            # Off the image?
            if (r+dr >= data.shape[0]) or (c+dc > data.shape[1]) or (
                use_r < 0) or (use_c < 0):
                fields.append(None)
            else:
                fields.append(data[use_r:use_r+dr, use_c:use_c+dc])
        print("Cutouts loaded for {}".format(filename))

        # Prep items
        if 'smooth_km' in pdict.keys():
            myzip = zip(fields, sub_idx, smooth_pixs)
        else:
            myzip = zip(fields, sub_idx)
        items = [item for item in myzip]

        # Test processing
        if test_process:
            idx = 75
            img, tmeta = process.preproc_field(fields[idx], None,
                                  smooth_pix=smooth_pixs[idx], 
                                  **pdict)
            # 
            from matplotlib import pyplot as plt
            fig = plt.figure(figsize=(8,8))
            plt.clf()
            plt.imshow(img, origin='lower')
            plt.show()
            embed(header='extract.py/preproc_field 130')

        # Multi-process time
        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            chunksize = len(items) // n_cores if len(items) // n_cores > 0 else 1
            answers = list(tqdm(executor.map(map_fn, items,
                                             chunksize=chunksize), total=len(items)))

        # Debuggin
        if test_failures:
            answers[50] = [None, answers[50][1], None]

        # Deal with failures
        #answers = [f for f in answers if f is not None]
        cur_img_idx = [item[1] for item in answers]

        # Find missing items and replace with -1.

        # Slurp
        pp_fields += [item[0] for item in answers]
        img_idx += cur_img_idx
        meta += [item[2] for item in answers]

        del answers, fields, items

    # Fuss with indices
    ex_idx = np.array(all_sub)
    img_idx = np.array(img_idx)
    ppf_idx = catalog.match_ids(img_idx, ex_idx, require_in_match=True)

    #llc_table = pp_utils.write_pp_fields(
    #    pp_fields, meta, llc_table, 
    #    ex_idx, ppf_idx, 
    #    valid_fraction, s3_file, local_file,
    #    kin_meta=kin_meta, debug=debug, write_cutouts=write_cutouts)


    # Clean up time

    # Find the bad ones (if any)
    bad_idx = [int(item[1]) for item in zip(pp_fields, img_idx) if item[0] is None]
    good_idx = np.array([int(item[1]) for item in zip(pp_fields, img_idx) if item[0] is not None])

    success = np.ones(len(pp_fields), dtype=bool)

    # Replace with -1 images
    if len(bad_idx) > 0:
        bad_img = -1*np.ones(field_size)
        for ii in bad_idx:
            pp_fields[ii] = bad_img.copy()
            success[ii] = False
    pp_fields = np.array(pp_fields)[ppf_idx]

    # Meta time
    good_meta = pandas.DataFrame([item for item in meta if item is not None])
    final_meta = pandas.DataFrame()
    for key in good_meta.keys():
        final_meta[key] = np.zeros(len(ppf_idx))
        #
        final_meta.loc[good_idx, key] = good_meta[key].values

    # Return
    return llc_table, success, pp_fields, final_meta