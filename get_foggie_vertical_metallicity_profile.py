#!/usr/bin/env python3

"""
    Title :      get_foggie_vertical_metallicity_profile
    Notes :      Plot metallicity profile of a given list of filenames (FOGGIE snapshots)
    Output :     Pandas dataframe
    Author :     Ayan Acharyya
    Started :    29-04-2026
    Examples :   run get_foggie_vertical_metallicity_profile.py --system ayan_ssd --halo 8508 --run nref11c_nref9f --output RD0027 --upto_kpc 20
                 run get_foggie_vertical_metallicity_profile.py --system ayan_ssd --halo 2392 --run 2162_39 --output DD1360 --upto_kpc 20
"""
from craft_header import *
from craft_utils import *
setup_plot_style()
from foggie_header import *
from make_3D_FRB_electron_density import get_AM_vector
from get_foggie_metallicity_profile import get_halo_coords, my_foggie_load, get_projection_frb
start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------
def make_hdu(data, quant, name_suffix, args):
    '''
    Makes a HDU, with given data, and some header keywords, to be saved as aprt of an HDULIST later
    Returns HDU
    '''
    hdu = fits.ImageHDU(data=data.astype(np.float32))
    header = hdu.header
    header['EXTNAME'] = f'{quant_dict[quant][1]} {name_suffix}'
    header['FIELD'] = quant_dict[quant][1]
    header['BUNIT'] = quant_dict[quant][8]
    header['WIDTH'] = (args.plot_width, 'kpc')
    
    for ind in range(2):
        header[f'CRPIX{ind+1}'] = np.shape(data)[ind] / 2. + 0.5
        header[f'CRVAL{ind+1}'] = 0.
        header[f'CDELT{ind+1}'] = (args.kpc_per_pix_proj, 'kpc')

    header['REDSHIFT'] = args.current_redshift
    header['TIME'] = (args.current_time, 'Gyr')

    return hdu

# --------------------------------------------------------------------------------------------------------
def write_projected_frb(quant_arr, args):
    '''
    Write projected edge on and face on views as FRB on to fits file
    Return fits hdulist
    '''
    # ----------reading in snapshot-----------
    halos_df_name = Path(args.code_path) / f'halo_infos/00{args.halo}/{args.run}/halo_c_v'
    halo_center = get_halo_coords(halos_df_name, args.output)
    snap_name = Path(args.foggie_dir) / f'halo_{int(args.halo):06d}' / args.run / args.output / args.output        
    
    ds, _ = my_foggie_load(snap_name, halo_center)
    args.current_redshift = ds.current_redshift
    args.current_time = ds.current_time.in_units('Gyr').tolist()

    # ----------determining extent-----------------------
    if args.upto_kpc is not None:
        if args.docomoving:
            args.galrad = args.upto_kpc / (1 + args.current_redshift) / 0.695  # fit within a fixed comoving kpc h^-1, 0.695 is Hubble constant
        else:
            args.galrad = args.upto_kpc  # fit within a fixed physical kpc
    else:
        args.re = get_re_from_coldgas(args) if args.use_gasre else get_re_from_stars(ds, args)
        args.galrad = args.re * args.upto_re  # kpc

    # -------extract the required box------------
    box_center = ds.halo_center_kpc
    box_width = args.galrad * 2  # in kpc
    box_width_kpc = ds.arr(box_width, 'kpc')
    box = ds.r[box_center[0] - box_width_kpc / 2.: box_center[0] + box_width_kpc / 2., box_center[1] - box_width_kpc / 2.: box_center[1] + box_width_kpc / 2., box_center[2] - box_width_kpc / 2.: box_center[2] + box_width_kpc / 2.]

    # ----------------setting up FITS file-------------------
    norm_L = get_AM_vector(ds, radius_kpc=3., use_particles=False)
    img_hdu_list = []
    primary_hdu = fits.PrimaryHDU()
    img_hdu_list.append(primary_hdu)

    args.plot_width = box_width/1.44
    args.ncell_buff = 800 # this buffer size is purely for visualising and saving the projected map
    args.kpc_per_pix_proj = args.plot_width / args.ncell_buff

    for quant in quant_arr:
        # ----------------making projection plots-------------------
        data_faceon, data_edgeon = get_projection_frb(box, quant_dict[quant][0], args.plot_width, norm_L, args, 
                                                quant_label=quant_dict[quant][0], 
                                                unit='Msun/pc**2' if quant == 'density' else quant_dict[quant][2], 
                                                return_frb=True,
                                                do_plot=False,
                                                ncell_buff = args.ncell_buff,
                                                )
    
        # --------making the FITS ImageHDU for 2D projected FRB: face on---------------
        hdu_faceon = make_hdu(data_faceon, quant, 'FACE ON PROJ', args)
        img_hdu_list.append(hdu_faceon)

        # --------making the FITS ImageHDU for 2D projected FRB: face on---------------
        hdu_edgeon = make_hdu(data_edgeon, quant, 'EDGE ON PROJ', args)
        img_hdu_list.append(hdu_edgeon)

    # -----------------saving the projections as fits file------------
    hdul = fits.HDUList(img_hdu_list)
    hdul.writeto(fitsname, overwrite=True)
    print(f'Successfully saved {len(hdul)-1} extensions to {fitsname}"')

    return hdul

# --------------------------------------------------------------------------------------------------------
def plot_projections_from_frb(data_proj, args, ax=None, crpix=None, cdelt=None, crval=0, quant_label='density', clim=None,  cmap='viridis', takelog=True, clabel=None, alpha=1, hide_cbar=False):
    '''
    Plots the input face on and edge on FRBs onto the same figure
    Returns figure handle
    '''
    if takelog:
        data_proj = np.log10(data_proj)

    if crpix is None:
        crpix = np.shape(data_proj)[0]/2.
    
    # ------plotting onto a matplotlib figure--------------
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        right, top, bottom, left, wspace = 0.85, 0.98, 0.1, 0.12, 0.0
        fig.subplots_adjust(right=right, top=top, bottom=bottom, left=left, wspace=wspace)
    
    fontsize = args.fontsize

    # ----------plotting axes----------------
    p = ax.imshow(data_proj, origin='lower', cmap=cmap, vmin=clim[0] if clim is not None else None, vmax=clim[1] if clim is not None else None, alpha=alpha)

    # ---------------making colorbar------------------------
    if not hide_cbar:
        fig = ax.figure
        cbar = fig.colorbar(p, ax=ax, orientation='vertical', fraction=0.047, pad=0.01)
        cbar.ax.tick_params(labelsize=fontsize, width=2.5, length=5)
        cbar.set_label(clabel, fontsize=fontsize) # this is to have more control on the display of the color label

    # ---------------prepping axes------------------------
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    if cdelt is None: ax.set_xticklabels(['%.1F' % item for item in ax.get_xticks()], fontsize=fontsize)
    else: ax.set_xticklabels(['%.1F' % ((item - crpix) * cdelt + crval) for item in ax.get_xticks()], fontsize=fontsize)
    ax.set_xlabel('Offset (kpc)', fontsize=fontsize)
    if cdelt is None: ax.set_yticklabels(['%.1F' % item for item in ax.get_yticks()], fontsize=fontsize)
    else: ax.set_yticklabels(['%.1F' % ((item - crpix) * cdelt + crval) for item in ax.get_yticks()], fontsize=fontsize)
    ax.set_ylabel('Offset (kpc)', fontsize=fontsize)
  
    # ---------------making annotations------------------------
    ax.text(0.97, 0.95, 'z = %.2F' % args.current_redshift, c='white', ha='right', va='top', transform=ax.transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))

    # ---------------saving fig------------------------
    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try: mplcyberpunk.make_lines_glow()
        except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    return ax

# --------------------------------------------------------------------------------------------------------
def bin_cones(data, center, opening_angle_deg, bin_edges):
    '''
    Computes the binned mean value in conical bins around the center of the image, within a given angle (in degrees)
    Returns the bin distances, binned mean values, conical mask and binned_map
    '''
    # Create coordinate grids
    yc, xc = center
    y, x = np.indices(data.shape)
    r = np.sqrt((x - xc)**2 + (y - yc)**2)
    # theta is 0 at +x axis, goes from -pi to pi
    theta = np.arctan2(y - yc, x - xc)
    
    # Define the Conical Mask (Vertical outflows)
    # Top cone is centered at pi/2 (90 deg), Bottom at -pi/2 (-90 deg)
    half_angle = np.radians(opening_angle_deg / 2)
    mask_top = np.abs(theta - np.pi/2) < half_angle
    mask_bottom = np.abs(theta + np.pi/2) < half_angle
    cone_mask = mask_top | mask_bottom
    
    # 3. Extract data within the cones
    r_in_cone = r[cone_mask]
    metal_in_cone = data[cone_mask]
    
    # 4. Perform Binning
    bin_means, edges, _ = binned_statistic(r_in_cone, metal_in_cone, statistic='mean', bins=bin_edges)
    bin_centers = (edges[:-1] + edges[1:]) / 2
    
    # Create a 2D map of the binned metallicity for plotting
    binned_map = np.full_like(data, np.nan)
    # Digitize the entire radial grid
    bin_indices = np.digitize(r, edges)
    
    for i in range(1, len(edges)):
        # Apply the mean value to all pixels in that radial bin within the cone
        bin_mask = (bin_indices == i) & cone_mask
        binned_map[bin_mask] = bin_means[i-1]
        
    return bin_centers, bin_means, cone_mask, binned_map

# --------------------------------------------------------------------------------------------------------
def overplot_bin_mask(xc, yc, angle, bins, ax):
    '''
    Overplot binned mask on data
    Returns axis handle
    '''
    # Draw cone boundaries
    for a in [90-angle/2, 90+angle/2, -90-angle/2, -90+angle/2]:
        ax.plot([xc, xc + 150*np.cos(np.radians(a))], [yc, yc + 150*np.sin(np.radians(a))], color='white', lw=1, ls='--')
    
    # Draw radial arcs
    for edge in bins:
        circle = plt.Circle((xc, yc), edge, color='white', fill=False, alpha=0.3)
        ax.add_artist(circle)

    return ax

# --------------------------------------------------------------------------------------------------------
def overplot_binned_data(binned_map, ax, args, clim=None, clabel='Mean Metallicity', cmap='viridis'):
    '''
    Overplot binned data on data
    Returns axis handle
    '''
    p = ax.imshow(binned_map, origin='lower', vmin=clim[0] if clim is not None else None, vmax=clim[1] if clim is not None else None, cmap=cmap)

    fig = ax.figure
    cbar = fig.colorbar(p, ax=ax, orientation='vertical', fraction=0.047, pad=0.01)
    cbar.ax.tick_params(labelsize=args.fontsize, width=2.5, length=5)
    cbar.set_label(clabel, fontsize=args.fontsize) # this is to have more control on the display of the color label
    
    return ax

# --------------------------------------------------------------------------------------------------------
def plot_bin_profile(centers, means, ax, args, ylabel='Mean Metallicity'):
    '''
    Plot radial binned profile
    Returns axis handle
    '''
    ax.scatter(centers, means, color='black', lw=2)
    ax = annotate_axes(ax, 'Radius (kpc)', ylabel, args=args)

    return ax

# ----------global variables------------
quant_dict = {'density':['density', 'Gas density', 'Msun/pc**3', -2.5, 2.5, 'cornflowerblue', density_color_map, True, 'Msun/pc**2', r'Gas density [M$_\odot$ pc$^{-2}$]'], 
              'metal':['metallicity', 'Gas metallicity', r'Zsun', -1.7, 0.7, 'cornflowerblue', metal_color_map, True, r'Zsun', r'Z/Z$_\odot$'],
              } # for each quantity: [yt field, label in plots, units, lower limit in log, upper limit in log, color for scatter plot, colormap, whether to take log, units for projection plot, units to display in projection plot]

# ---------------------main code-----------------------------
if __name__ == '__main__':
    args = parse_args()  # default simulation to work upon when comand line args not provided
    if not args.keep: plt.close('all')
    args.fontfactor = 1.
    if not 'pleiades' in args.system: args.output_dir = args.output_dir.replace('CRAFT', 'sharing')

    quant_arr = ['metal', 'density']

    # -----------determining directories----------------
    args.fig_dir = args.output_dir + 'metallicity_plots/'
    Path(args.fig_dir).mkdir(parents=True, exist_ok=True)
    
    args.data_dir = Path(args.output_dir) / 'metallicity_data'
    args.data_dir.mkdir(exist_ok=True, parents=True)
    
    args.upto_text = '_upto%.1Fckpchinv' % args.upto_kpc if args.docomoving and args.upto_kpc is not None else '_upto%.1Fkpc' % args.upto_kpc if args.upto_kpc is not None else f'_upto{args.re:.1f}re'
    fitsname = args.data_dir / f'{args.halo}_{args.run}_{args.output}_projected_met_map{args.upto_text}.fits'

    # ---------------getting the metallicity profile----------------
    if not os.path.exists(fitsname) or args.clobber:
        print(f']\n{fitsname} does not exist. Creating afresh..')
        hdul = write_projected_frb(quant_arr, args)
    else:
        print(f'\nReading from existing file {fitsname}')
        hdul = fits.open(fitsname)

    # --------preparing data------------
    quant = 'metal'
    hdu = hdul[f'{quant_dict[quant][1].upper()} EDGE ON PROJ']
    header = hdu.header
    crpix, crval, cdelt = header['CRPIX1'], header['CRVAL1'], header['CDELT1']
    data = hdu.data

    args.current_redshift = 1#header['REDSHIFT']
    
    # --------preparing cone bins------------
    bin_upto_kpc = args.upto_kpc / 4
    bin_width_kpc = 0.5
    angle = 45 # total opening angle

    bin_upto_pix = bin_upto_kpc / cdelt
    nbins = int(bin_upto_kpc / bin_width_kpc)
    rad_bins = np.linspace(0, bin_upto_pix, nbins) 
    
    centers, means, cone_mask, binned_map = bin_cones(data, (crpix, crpix), angle, rad_bins)
    if quant_dict[quant][7]:
        binned_map = np.log10(binned_map)
        means = np.log10(means)
        quant_dict[quant][9] = 'Log' + quant_dict[quant][9]

    # ----------------making projection plots-------------------
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))#, constrained_layout=True)
    fig.subplots_adjust(left=0.07, bottom=0.1, right=0.98, top=0.98, wspace=0.4)
    
    axes[0] = plot_projections_from_frb(data, args, ax=axes[0],
                                            crpix=crpix, cdelt=cdelt, crval=crval, 
                                            quant_label=quant_dict[quant][0], 
                                            clim=[quant_dict[quant][3], quant_dict[quant][4]] if quant_dict[quant][3] is not None else None,  
                                            cmap=quant_dict[quant][6], 
                                            takelog=quant_dict[quant][7], 
                                            clabel=quant_dict[quant][9],
                                            alpha=1,
                                            hide_cbar=False,
                                            )

    axes[1] = plot_projections_from_frb(data, args, ax=axes[1],
                                            crpix=crpix, cdelt=cdelt, crval=crval, 
                                            quant_label=quant_dict[quant][0], 
                                            clim=[quant_dict[quant][3], quant_dict[quant][4]] if quant_dict[quant][3] is not None else None,  
                                            cmap=quant_dict[quant][6], 
                                            takelog=quant_dict[quant][7], 
                                            clabel=quant_dict[quant][9],
                                            alpha=0.5,
                                            hide_cbar=True,
                                            )

    # ----------------overplotting bins-------------------
    axes[0] = overplot_bin_mask(crpix, crpix, angle, rad_bins, axes[0])
    axes[1] = overplot_binned_data(binned_map, axes[1], args, clim=[quant_dict[quant][3], quant_dict[quant][4]] if quant_dict[quant][3] is not None else None, clabel=quant_dict[quant][9], cmap=quant_dict[quant][6])
    
    # ----------------making profile plots-------------------
    centers_kpc = centers * cdelt + crval
    axes[2] = plot_bin_profile(centers_kpc, means, axes[2], args, ylabel=quant_dict[quant][9])

    # ------------saving figure--------------
    figname = f'{args.output}_{args.halo}_binned_{quant}{args.upto_text}.png'
    save_fig(fig, Path(args.fig_dir), figname, args)
    
    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
