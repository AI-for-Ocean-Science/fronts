""" Script to grab/plot data from LLC """

from IPython import embed

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Grab/plot LLC data')
    parser.add_argument("table_file", type=str, help="File+path to Preproc file")
    parser.add_argument("UID", type=int, help="Evaluation table filename (.csv)")
    parser.add_argument("--fields", type=str, default="SST,SSS,Divb2", help="Evaluation table filename (.csv)")
    parser.add_argument("--fig_file", type=str, help="Filename for figure")
    parser.add_argument("--data_file", type=str, help="Filename for figure")
    parser.add_argument("--field_size", type=int, default=64, help="Filename for figure")
    parser.add_argument("--fixed_km", type=float, default=144., help="Filename for figure")
    parser.add_argument("-s","--show", default=False, action="store_true", help="Show pre-processed image?")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):
    """ Run
    """
    import numpy as np
    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec
    import seaborn as sns

    import pandas

    from fronts.llc import io as llc_io
    from fronts.llc import extract as llc_extract

    # Load the table
    print("Loading the table..")
    df = pandas.read_parquet(pargs.table_file)
    df_idx = np.where(df.UID == pargs.UID)[0][0]
    idf = df.iloc[df_idx]

    # Fields
    fields = pargs.fields.split(",")
    nfields = len(fields)

    # Grab the data
    udate = idf.datetime
    filename = llc_io.grab_llc_datafile(udate, local=True)
    ds = llc_io.load_llc_ds(filename, local=True)

    images = []
    for field in fields:
        full = llc_extract.field_from_ds(ds, field)
        # Cut out
        cutout = full[idf.row:idf.row+pargs.field_size, 
                      idf.col:idf.col+pargs.field_size]
        images.append(cutout)
        
    # Show?
    if pargs.show or pargs.fig_file is not None:
        vmnx=(-5,5)
        # Plot
        fig = plt.figure(figsize=(10, 4))
        #pal, cm = plotting.load_palette()
        plt.clf()
        gs = gridspec.GridSpec(1, nfields)

        # Original
        ax1 = plt.subplot(gs[0])
        sns.heatmap(field[0, ...], ax=ax1, xticklabels=[], yticklabels=[], cmap=cm,
                    vmin=vmnx[0], vmax=vmnx[1])

        # Reconstructed
        ax2 = plt.subplot(gs[1])
        sns.heatmap(recons[0, 0, ...], ax=ax2, xticklabels=[], yticklabels=[], cmap=cm,
                    vmin=vmnx[0], vmax=vmnx[1])
        if pargs.fig_file is not None:
            plt.savefig(pargs.fig_file, dpi=300)
        if pargs.show:
            plt.show()

