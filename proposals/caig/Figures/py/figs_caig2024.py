# imports
from importlib import reload

import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec


import pandas
import seaborn as sns
import xarray as xr

from skimage import filters as sf

# ulmo -- on f_s branch
from ulmo.llc import kinematics
from ulmo.llc import io as llc_io
from ulmo import io as ulmo_io
from ulmo.plotting import plotting 

from IPython import embed

def fig_fronts(outfile:str='MNIST.png'):

    _, cm = plotting.load_palette()

    # Load the data
    llc_tfile = os.path.join(os.getenv('OS_OGCM'), 'LLC', 'Ulmo', 'Tables', 'LLC_uniform144_r0.5_nonoise.parquet')
    llc_u = ulmo_io.load_main_table(llc_tfile)


    # Grab Brazil-Malvanis
    # Assume 1.5km
    lat_BM = -42.  # 42S
    lon_BM = -50.  # 50W
    d = np.sqrt((llc_u.lat-lat_BM)**2 + (llc_u.lon-lon_BM)**2)

    # Pick one
    srt = np.argsort(d)
    idx = srt[50]
    # 
    cutout = llc_u.iloc[idx].copy()


    cutout.filename = 'LLC4320_2011-09-30T00_00_00.nc'
    ds = xr.open_dataset(os.path.join(os.getenv('OS_OGCM'), 'LLC', 'data', 
                                      'ThetaUVWSaltEta', cutout.filename))

    # Grab T, S
    cutout.field_size = 128
    xoffset = 100
    yoffset = 350
    #cutout.field_size = 500
    #xoffset = 0
    #yoffset = 0
    row = cutout.row + yoffset
    col = cutout.col + xoffset


    Theta = ds.Theta[row:row+cutout.field_size, 
                col:col+cutout.field_size].values
    Salt = ds.Salt[row:row+cutout.field_size, 
                col:col+cutout.field_size].values

    # Coords
    R_earth = 6371. # km
    circum = 2 * np.pi* R_earth
    km_deg = circum / 360.
    coords_ds = llc_io.load_coords()
    lat0, lon0 = coords_ds.lat.data[row,col], coords_ds.lon.data[row,col]
    erow, ecol = row+cutout.field_size, col+cutout.field_size
    late, lone = coords_ds.lat.data[erow,ecol], coords_ds.lon.data[erow,ecol]


    # Front Intensity
    divb_2 = kinematics.calc_gradb(Theta, Salt)

    cbar_kws = {} 
    cbar_kws['pad']=0. 
    cbar_kws['fraction']=0.040

    # Figure
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(2,2)

    all_ax = []
    # Theta
    ax_T = plt.subplot(gs[0])
    cbar_kws['label'] = 'SST (deg C)'
    img = plt.(np.flipud(Theta), xticklabels=[], 
                     #vmin=vmnx[0], vmax=vmnx[1], 
                     ax=ax_T,
                     yticklabels=[], cmap=cm, cbar=True, 
                     cbar_kws=cbar_kws)
    all_ax.append(ax_T)

    # Front intensity
    ax_b = plt.subplot(gs[3])
    all_ax.append(ax_b)
    img = np.log10(divb_2) 
    vmnx=(-14., -11.3)
    cbar_kws['label'] = r'$\log |\nabla b|^2$' 
    xticks = np.linspace(0, cutout.field_size, 6)
    yticks = np.linspace(0, cutout.field_size, 6)
    _ = sns.heatmap(np.flipud(img), 
                    #xticklabels=[], yticklabels=[], 
                    vmin=vmnx[0], vmax=vmnx[1], 
                    ax=ax_b,
                    cmap='Greys', cbar=True, 
                    cbar_kws=cbar_kws)

    # Label
    ax_b.set_xlabel('Longitude')
    ax_b.set_ylabel('Latitude')

    # Set ticks
    xlim = ax_b.get_xlim()
    ylim = ax_b.get_ylim()
    
    #embed(header='100 of figs_caig2024.py')
    


    # Set Aspect ratio
    for ax in all_ax:
        ax.set_aspect('equal')

    plt.tight_layout()#pad=0.0, h_pad=0.0, w_pad=0.3)
    plt.savefig(outfile, dpi=300)
    print(f"Saved: {outfile}")

# Command line execution
if __name__ == '__main__':
    fig_fronts()