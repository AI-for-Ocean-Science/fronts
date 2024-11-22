""" Routines to extract cutouts from LLC data """

import pandas

def preproc_for_analysis(llc_table:pandas.DataFrame, 
                         local_file:str,
                         preproc_root:str,
                         fields:list,
                         field_size=(64,64), 
                         fixed_km=None,
                         n_cores=10,
                         valid_fraction=1., 
                         dlocal=False,
                         write_cutouts:bool=True,
                         override_RAM=False,
                         s3_file=None, debug=False):
    """Main routine to extract and pre-process LLC data for later SST analysis
    The llc_table is modified in place (and also returned).

    Args:
        llc_table (pandas.DataFrame): cutout table
        local_file (str): path to PreProc file
        preproc_root (str, optional): Preprocessing steps. Defaults to 'llc_std'.
        fields (list): Fields to extract
        field_size (tuple, optional): Defines cutout shape. Defaults to (64,64).
        fixed_km (float, optional): Require cutout to be this size in km
        n_cores (int, optional): Number of cores for parallel processing. Defaults to 10.
        valid_fraction (float, optional): [description]. Defaults to 1..
        calculate_kin (bool, optional): Perform frontogenesis calculations?
        extract_kin (bool, optional): Extract kinematic cutouts too!
        kin_stat_dict (dict, optional): dict for guiding FS stats
        dlocal (bool, optional): Data files are local? Defaults to False.
        override_RAM (bool, optional): Over-ride RAM warning?
        s3_file (str, optional): s3 URL for file to write. Defaults to None.
        write_cutouts (bool, optional): 
            Write the cutouts to disk?

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
    
    # Preprocess options
    pdict = pp_io.load_options(preproc_root)

    # Setup for parallel
    map_fn = partial(pp_utils.preproc_image, pdict=pdict)

    # Kinematics
    if calculate_kin:
        if kin_stat_dict is None:
            raise IOError("You must provide kin_stat_dict with calculate_kin")
        # Prep
        if 'calc_FS' in kin_stat_dict.keys() and kin_stat_dict['calc_FS']:
            map_kin = partial(kinematics.cutout_kin, 
                         kin_stats=kin_stat_dict,
                         extract_kin=extract_kin,
                         field_size=field_size[0])

    # Setup for dates
    uni_date = np.unique(llc_table.datetime)
    if len(llc_table) > 1000000 and not override_RAM:
        raise IOError("You are likely to exceed the RAM.  Deal")

    # Init
    pp_fields, meta, img_idx, all_sub = [], [], [], []
    if calculate_kin:
        kin_meta = []
    else:
        kin_meta = None
    if extract_kin:  # Cutouts of kinematic information
        Fs_fields, divb_fields = [], []

    # Prep LLC Table
    llc_table = pp_utils.prep_table_for_preproc(
        llc_table, preproc_root, field_size=field_size)
    # Loop
    #if debug:
    #    uni_date = uni_date[0:1]

    for udate in uni_date:
        # Parse filename
        filename = llc_io.grab_llc_datafile(udate, local=dlocal)

        # Allow for s3
        ds = llc_io.load_llc_ds(filename, local=dlocal)
        sst = ds.Theta.values
        # Parse 
        gd_date = llc_table.datetime == udate
        sub_idx = np.where(gd_date)[0]
        all_sub += sub_idx.tolist()  # These really should be the indices of the Table
        coord_tbl = llc_table[gd_date]

        # Add to table
        llc_table.loc[gd_date, 'filename'] = filename

        # Load up the cutouts
        fields, rs, cs, drs = [], [], [], []
        for r, c in zip(coord_tbl.row, coord_tbl.col):
            if fixed_km is None:
                dr = field_size[0]
                dc = field_size[1]
            else:
                dlat_km = (coords_ds.lat.data[r+1,c]-coords_ds.lat.data[r,c]) * km_deg
                dr = int(np.round(fixed_km / dlat_km))
                dc = dr
                # Save for kinematics
                drs.append(dr)
                rs.append(r)
                cs.append(c)
            #
            if (r+dr >= sst.shape[0]) or (c+dc > sst.shape[1]):
                fields.append(None)
            else:
                fields.append(sst[r:r+dr, c:c+dc])
        print("Cutouts loaded for {}".format(filename))

        # Multi-process time
        # 
        items = [item for item in zip(fields,sub_idx)]

        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            chunksize = len(items) // n_cores if len(items) // n_cores > 0 else 1
            answers = list(tqdm(executor.map(map_fn, items,
                                             chunksize=chunksize), total=len(items)))

        # Deal with failures
        answers = [f for f in answers if f is not None]
        cur_img_idx = [item[1] for item in answers]

        # Slurp
        pp_fields += [item[0] for item in answers]
        img_idx += cur_img_idx
        meta += [item[2] for item in answers]

        del answers, fields, items

    # Fuss with indices
    ex_idx = np.array(all_sub)
    ppf_idx = []
    ppf_idx = catalog.match_ids(np.array(img_idx), ex_idx)

    # Write
    llc_table = pp_utils.write_pp_fields(
        pp_fields, meta, llc_table, 
        ex_idx, ppf_idx, 
        valid_fraction, s3_file, local_file,
        kin_meta=kin_meta, debug=debug, write_cutouts=write_cutouts)

    # Write kin?
    if extract_kin:
        # F_s
        Fs_local_file = local_file.replace('.h5', '_Fs.h5')
        pp_utils.write_extra_fields(Fs_fields, llc_table, Fs_local_file)
        # divb
        divb_local_file = local_file.replace('.h5', '_divb.h5')
        pp_utils.write_extra_fields(divb_fields, llc_table, divb_local_file)
    
    # Clean up
    del pp_fields

    # Upload to s3? 
    if s3_file is not None:
        ulmo_io.upload_file_to_s3(local_file, s3_file)
        print("Wrote: {}".format(s3_file))
        # Delete local?

    # Return
    return llc_table 