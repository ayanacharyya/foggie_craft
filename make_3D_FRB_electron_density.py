#!/usr/bin/env python3

"""

    Title :      make_3D_FRB_electron_density
    Notes :      Make a 3D Fixed Resolution Buffer (FRB) for gas density and electron number density and save as a multi-extension fits file, and optionally plot along a chosen line of sight
    Output :     3D data cube as fits file, and optionally png figures
    Author :     Ayan Acharyya
    Started :    Aug 2024
    Examples :   run make_3D_FRB_electron_density.py --system ayan_pleiades --halo 8508 --res 1 --upto_kpc 50 --do_all_sims --use_cen_smoothed
                 run make_3D_FRB_electron_density.py --system ayan_hd --halo 4123 --res 1 --upto_kpc 10 --output RD0038 --clobber --plot_3d --use_cen_smoothed
                 run make_3D_FRB_electron_density.py --system ayan_hd --halo 8508 --res 1 --upto_kpc 200 --output RD0030,RD0042 --clobber --use_cen_smoothed
                 run make_3D_FRB_electron_density.py --system ayan_local --halo 8508 --res 0.5 --upto_kpc 100 --output RD0027 --clobber --use_cen_smoothed --do_only_plot
"""
from foggie_header import *
from yt.visualization.fits_image import FITSImageData
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams['axes.linewidth'] = 1

start_time = datetime.now()

# -----------------------------------------------------------------------------
def get_sfr_df(args):
    '''
    Reads in the dataframe that contains SFR of each snapshot of a given halo
    Returns pandas dataframe
    '''
    sfr_filename = args.code_path + 'halo_infos/00' + args.halo + '/' + args.run + '/sfr'
    if os.path.exists(sfr_filename):
        print('Reading SFR history from', sfr_filename)
        sfr_df = pd.read_table(sfr_filename, names=('output', 'redshift', 'sfr'), comment='#', delim_whitespace=True)
    else:
        print('Did not find', sfr_filename, ', therefore will not include SFR')
        sfr_df = pd.DataFrame()

    return sfr_df

# -----------------------------------------------------------------------------
def get_AM_vector(ds, radius_kpc=15., use_particles='young_stars'):
    '''
    Computes the orientation vector of angular momentum of the disk, in the given dataset, considering young star particles
    Based on foggie_load()
    Returns the unit vector as a list
    '''
    start_time = datetime.now()

    print('Starting to derive angular momentum vector. This can take a while..')
    sphere = ds.sphere(ds.halo_center_kpc, (radius_kpc, 'kpc'))
    if use_particles == 'gas' or use_particles == False:
        L = sphere.quantities.angular_momentum_vector(use_gas=True, use_particles=False)
    else:
        L = sphere.quantities.angular_momentum_vector(use_gas=False, use_particles=True, particle_type=use_particles)
    print('Completed deriving angular momentum vector, in %s'% timedelta(seconds=(datetime.now() - start_time).seconds))
    norm_L = L / np.sqrt((L ** 2).sum())
    norm_L = np.array(norm_L.value)

    return norm_L

# --------------------------------------------------------------------------
def plot_3d_frb(data, ax, args, label=None, unit=None, clim=None,  cmap='viridis'):
    '''
    Function to make a 3D plot given a 3D numpy array in a given axis
    Returns axis handle
    '''
    z, x, y = data.nonzero()
    ax.scatter3D(z, x, y, c=np.log10(data), cmap=cmap, alpha=0.3)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    return ax

# --------------------------------------------------------------------------
def plot_proj_frb(data, ax, args, label='', unit='', clim=None,  cmap='viridis', hidex=False, hidey=False):
    '''
    Function to make a 2D projection plot (along one line of sight) given a 3D numpy array, in a given axis
    Returns axis handle
    '''
    data_proj = np.sum(data, axis=2)
    p = ax.imshow(np.log10(data_proj), cmap=cmap)

    # -----------making the axis labels etc--------------
    if hidex:
        ax.set_xticklabels(['' % item for item in ax.get_xticks()])
        ax.set_xlabel('')
    else:
        ax.set_xticklabels(['%.1F' % ((item - central_pixel) * args.kpc_per_pix) for item in ax.get_xticks()], fontsize=args.fontsize)
        ax.set_xlabel('Offset (kpc)', fontsize=args.fontsize)

    if hidey:
        ax.set_yticklabels(['' % item for item in ax.get_yticks()])
        ax.set_ylabel('')
    else:
        ax.set_yticklabels(['%.1F' % ((item - central_pixel) * args.kpc_per_pix) for item in ax.get_yticks()], fontsize=args.fontsize)
        ax.set_ylabel('Offset (kpc)', fontsize=args.fontsize)

    # ---------making the colorbar axis once, that will correspond to all projections--------------
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('top', size='5%', pad=0.05)
    cbar = fig.colorbar(p, orientation='horizontal', cax=cax)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.tick_params(labelsize=args.fontsize, width=2.5, length=5)
    cbar.set_label(f'log(LoS summed {label} ({unit}))', fontsize=args.fontsize/1.)
    cbar.set_label(f'log(LoS summed {label} ({unit}))', fontsize=args.fontsize/1.)

    # ---------------making annotations------------------------
    ax.text(0.97, 0.95, 'z = %.2F' % args.current_redshift, c='white', ha='right', va='top', transform=ax.transAxes, fontsize=args.fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))
    ax.text(0.97, 0.85, 't = %.1F Gyr' % args.current_time, c='white', ha='right', va='top', transform=ax.transAxes, fontsize=args.fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))

    return ax

# --------------------------------------------------------------------------
def plot_projection_diskrel(box, field, box_width, norm_L, args, quant_label='density', unit='', clim=None,  cmap='viridis', takelog=True, clabel=None, annotate_labels=None, dpi=300, do_plot=True, return_frb=False, ncell_buff=800, field_type='gas'):
    '''
    Function to make a 2D projection plot along edge-on and face-on views given a dataset
    Borrowed a little from foggie_load()
    Returns figure handle
    '''
    if not do_plot: return_frb = True # if not plotting then has to return FRB objects
    
    x = np.array([1., 0., 0.])  # take a random vector
    x -= x.dot(norm_L) * norm_L  # make it orthogonal to L
    x /= np.linalg.norm(x)  # normalize it
    y = np.cross(norm_L, x)  # cross product with L

    field = (field_type, field)
    fontsize = args.fontsize

    # ---------------making face on and edge on projections------------------------
    p_faceon = yt.OffAxisProjectionPlot(box.ds, box.ds.arr(norm_L), field, data_source=box, width=(box_width, 'kpc'), weight_field=None, center=box.ds.halo_center_kpc, north_vector=box.ds.arr(x))
    p_edgeon = yt.OffAxisProjectionPlot(box.ds, box.ds.arr(x), field, data_source=box, width=(box_width, 'kpc'), weight_field=None, center=box.ds.halo_center_kpc, north_vector=box.ds.arr(norm_L))

    # ------------set projection parameters------------
    p_faceon.set_buff_size((ncell_buff, ncell_buff))
    p_edgeon.set_buff_size((ncell_buff, ncell_buff))

    p_faceon.set_unit(field, unit)
    p_edgeon.set_unit(field, unit)

    # ----------------make the FRB if required---------------------
    if return_frb:
        data_faceon = p_faceon.frb[field].v
        data_edgeon = p_edgeon.frb[field].v
    
    # ----------------make the plot if required---------------------
    if do_plot:
        # ---------------setting up units, colormaps, etc------------------------
        p_faceon.set_log(field, takelog)
        if clim is not None: p_faceon.set_zlim(field, zmin=10**clim[0] if takelog else clim[0], zmax=10**clim[1] if takelog else clim[1])
        p_faceon.set_cmap(field, cmap)

        p_edgeon.set_log(field, takelog)
        if clim is not None: p_edgeon.set_zlim(field, zmin=10**clim[0] if takelog else clim[0], zmax=10**clim[1] if takelog else clim[1])
        p_edgeon.set_cmap(field, cmap)

        # ------plotting onto a matplotlib figure--------------
        fig, axes = plt.subplots(1, 2, figsize=(8, 4))

        p_faceon.plots[field].axes = axes[0]
        p_faceon._setup_plots()
        p_edgeon.plots[field].axes = axes[1]
        p_edgeon._setup_plots()

        right, top, bottom, left, wspace = 0.85, 0.98, 0.1, 0.12, 0.0
        fig.subplots_adjust(right=right, top=top, bottom=bottom, left=left, wspace=wspace)

        # ---------------making colorbar------------------------
        offset, width = 0.075, 0.02
        cax = fig.add_axes([right, bottom + offset, width, top - bottom - 2 * offset])
        cbar = fig.colorbar(p_edgeon.plots[field].cb.mappable, orientation='vertical', cax=cax)
        cbar.ax.tick_params(labelsize=fontsize, width=2.5, length=5)
        if clabel is None: cbar.set_label(p_edgeon.plots[field].cax.get_ylabel(), fontsize=fontsize)
        else: cbar.set_label(clabel, fontsize=fontsize) # this is to have more control on the display of the color label

        # ---------------prepping axes------------------------
        for index in range(len(axes)):
            ax = axes[index]
            ax.xaxis.set_major_locator(plt.MaxNLocator(5))
            ax.yaxis.set_major_locator(plt.MaxNLocator(5))
            ax.set_xticklabels(['%.1F' % item for item in ax.get_xticks()], fontsize=fontsize)
            ax.set_xlabel('Offset (kpc)', fontsize=fontsize)
            if index == 0:
                ax.set_yticklabels(['%.1F' % item for item in ax.get_yticks()], fontsize=fontsize)
                ax.set_ylabel('Offset (kpc)', fontsize=fontsize)
            else:
                ax.set_yticklabels(['' % item for item in ax.get_yticks()])
                ax.set_ylabel('')
            
            if 'massrad' in args:
                ax.add_patch(plt.Circle((0, 0), args.massrad, color='r', fill=False, lw=1)) # for over-plotting radius within which mass was computed

        # ---------------making annotations------------------------
        if annotate_labels is None:
            axes[0].text(0.97, 0.95, 'z = %.2F' % args.current_redshift, c='white', ha='right', va='top', transform=axes[0].transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))
            axes[0].text(0.97, 0.85, 't = %.1F Gyr' % args.current_time, c='white', ha='right', va='top', transform=axes[0].transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))
        else:
            for index, label in enumerate(annotate_labels):
                axes[index % 2].text(0.97, 0.95 - (index // 2) * 0.1, label, c='white', ha='right', va='top', transform=axes[index % 2].transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))
    
        axes[0].text(0.98, 0.02, 'Face on', c='white', ha='right', va='bottom', transform=axes[0].transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))
        axes[1].text(0.98, 0.02, 'Edge on', c='white', ha='right', va='bottom', transform=axes[1].transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))

        # ---------------saving fig------------------------
        outfile_rootname = '%s_%s_diskrel_%s%s.pdf' % (args.output, args.halo, quant_label, args.upto_text)
        if args.do_all_sims: outfile_rootname = 'z=*_' + outfile_rootname[len(args.output) + 1:]
        if 'el' in quant_label.lower(): fig_dir = args.fig_dir + 'electron_density/'
        else: fig_dir = args.fig_dir + 'gas_density/'
        figname = fig_dir + outfile_rootname.replace('*', '%.5F' % (args.current_redshift))

        if args.fortalk:
            mplcyberpunk.add_glow_effects()
            try: mplcyberpunk.make_lines_glow()
            except: pass
            try: mplcyberpunk.make_scatter_glow()
            except: pass

        plt.savefig(figname, dpi=dpi, transparent=args.fortalk)
        myprint('Saved figure ' + figname, args)
        if not ('pleiades' in args.system or args.hide_plot): plt.show()

    if return_frb and do_plot:
        return fig, data_faceon, data_edgeon
    elif return_frb:
        return data_faceon, data_edgeon
    else:
        return fig

# ---------------global dictionary--------------------------
quant_dict = {'density':['density', 'Gas density', 'Msun/pc**3', -2.5, 2.5, 'cornflowerblue', density_color_map, True, 'Msun/pc**2', r'Gas density [M$_\odot$ pc$^{-2}$]'], 
              'el_density':['El_number_density', 'Electron density', 'cm**-3', 0, 220, 'cornflowerblue', 'viridis', False, 'pc*cm**-3', r'DM [pc cm$^{-3}$]'],
              } # for each quantity: [yt field, label in plots, units, lower limit in log, upper limit in log, color for scatter plot, colormap, whether to take log, units for projection plot, units to display in projection plot]

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')

    quant_arr = ['el_density', 'density']

    # ------------reading SFR and mstar df-----------------
    sfr_df  = get_sfr_df(args)
    smooth_over_snap = 20
    sfr_df[f'sfr_smooth{smooth_over_snap}'] = sfr_df['sfr'].rolling(window=smooth_over_snap, center=True).mean()

    # --------domain decomposition; for mpi parallelisation-------------
    if args.do_all_sims: list_of_sims = get_all_sims_for_this_halo(args) # all snapshots of this particular halo
    else: list_of_sims = list(itertools.product([args.halo], args.output_arr))
    total_snaps = len(list_of_sims)

    comm = MPI.COMM_WORLD
    ncores = comm.size
    rank = comm.rank
    print_master('Total number of MPI ranks = ' + str(ncores) + '. Starting at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.now()), args)
    comm.Barrier() # wait till all cores reached here and then resume

    split_at_cpu = total_snaps - ncores * int(total_snaps/ncores)
    nper_cpu1 = int(total_snaps / ncores)
    nper_cpu2 = nper_cpu1 + 1
    if rank < split_at_cpu:
        core_start = rank * nper_cpu2
        core_end = (rank+1) * nper_cpu2 - 1
    else:
        core_start = split_at_cpu * nper_cpu2 + (rank - split_at_cpu) * nper_cpu1
        core_end = split_at_cpu * nper_cpu2 + (rank - split_at_cpu + 1) * nper_cpu1 - 1

    # -------------loop over snapshots-----------------
    print_mpi('Operating on snapshots ' + str(core_start + 1) + ' to ' + str(core_end + 1) + ', i.e., ' + str(core_end - core_start + 1) + ' out of ' + str(total_snaps) + ' snapshots', args)

    for index in range(core_start + args.start_index, core_end + 1):
        start_time_this_snapshot = datetime.now()
        this_sim = list_of_sims[index]
        args.output = this_sim[1]
        print_mpi('Doing snapshot ' + this_sim[1] + ' of halo ' + this_sim[0] + ' which is ' + str(index + 1 - core_start) + ' out of the total ' + str(core_end - core_start + 1) + ' snapshots...', args)

        # -------loading in snapshot-------------------
        halos_df_name = args.code_path + 'halo_infos/00' + this_sim[0] + '/' + args.run + '/'
        halos_df_name += 'halo_cen_smoothed' if args.use_cen_smoothed else 'halo_c_v'
        ds, refine_box = load_sim(args, region='refine_box', do_filter_particles=True, disk_relative=False, halo_c_v_name=halos_df_name)

        norm_L = get_AM_vector(ds) #computing disk orientation #
        #norm_L = np.array([-0.51443095, -0.62174905, -0.59058354]) # this is just for faster testing; this is for halo 8505 snap RD0027

        # --------assigning additional keyword args-------------
        args.upto_text = '_upto%.1Fckpchinv' % args.upto_kpc if args.docomoving else '_upto%.1Fkpc' % args.upto_kpc
        args.current_redshift = ds.current_redshift
        args.current_time = ds.current_time.in_units('Gyr').tolist()
        args.res = args.res_arr[0]
        args.res_text = f'_res{args.res:.1f}kpc'
        if args.docomoving: args.res = args.res / (1 + args.current_redshift) / 0.695  # converting from comoving kcp h^-1 to physical kpc
        args.fontsize = 15

        # -----------determining SFR amd stellar mass--------------------
        try: sfr = sfr_df[sfr_df['output'] == args.output][f'sfr_smooth{smooth_over_snap}'].values[0]
        except: sfr = -99

        if args.do_only_plot: args.diskrad, log_mstar = np.nan, np.nan
        else: args.diskrad, log_mstar = get_stellar_mass(args, refine_box=refine_box)
        #else: args.diskrad, log_mstar = 5, 10.46  # this is just for faster testing; this is for halo 8505 snap RD0027

        # --------determining corresponding text suffixes and figname-------------
        args.fig_dir = args.output_dir + 'plots/'
        Path(args.fig_dir).mkdir(parents=True, exist_ok=True)

        args.fits_dir = args.output_dir + 'data/'
        Path(args.fits_dir).mkdir(parents=True, exist_ok=True)

        outfile_rootname = '%s_%s_FRB_%s%s%s.png' % (args.output, args.halo, quant_dict[quant_arr[0]][0], args.upto_text, args.res_text)
        if args.do_all_sims: outfile_rootname = 'z=*_' + outfile_rootname[len(args.output) + 1:]
        figname = args.fig_dir + outfile_rootname.replace('*', '%.5F' % (args.current_redshift))
        fitsname = args.fits_dir + outfile_rootname.replace('*', '%.5F' % (args.current_redshift)).replace('.png', '.fits')

        if not os.path.exists(fitsname) or args.clobber:
            try:
                # ------tailoring the simulation box for individual snapshot analysis--------
                if args.upto_kpc is not None:
                    if args.docomoving: args.galrad = args.upto_kpc / (1 + args.current_redshift) / 0.695  # fit within a fixed comoving kpc h^-1, 0.695 is Hubble constant
                    else: args.galrad = args.upto_kpc  # fit within a fixed physical kpc
                else:
                    args.re = get_re_from_coldgas(args) if args.use_gasre else get_re_from_stars(ds, args)
                    args.galrad = args.re * args.upto_re  # kpc

                # extract the required box
                box_width = args.galrad * 2  # in kpc
                box_width_kpc = ds.arr(box_width, 'kpc')
                args.ncells = int(box_width / args.res)
                box_center = ds.halo_center_kpc
                box = ds.r[box_center[0] - box_width_kpc / 2.: box_center[0] + box_width_kpc / 2., box_center[1] - box_width_kpc / 2.: box_center[1] + box_width_kpc / 2., box_center[2] - box_width_kpc / 2.: box_center[2] + box_width_kpc / 2., ]

                central_pixel = args.ncells / 2
                args.kpc_per_pix = 2 * args.galrad / args.ncells

                # -------setting up fig--------------
                if args.plot_3d or args.plot_proj:
                    fig = plt.figure(figsize=(2 + 4 * len(quant_arr), 5))
                    fig.subplots_adjust(top=0.88, bottom=0.12, left=0.07, right=0.92, wspace=0.4 if args.plot_3d else 0.02, hspace=0.)

                # -------making and plotting the 3D FRBs--------------
                all_data = ds.arbitrary_grid(left_edge=box.left_edge, right_edge=box.right_edge, dims=[args.ncells, args.ncells, args.ncells])
                img_hdu_list = []

                plot_width = box_width/5
                ncell_buff = 800 # this buffer size is purely for visualising and saving the projected map
                kpc_per_pix_proj = plot_width / ncell_buff

                for index, quant in enumerate(quant_arr):
                    myprint(f'Making and plotting FRB for {quant} which is {index+1} out of {len(quant_arr)} quantities..', args)

                    # --------making the projection plots------------
                    fig_diskrel, data_faceon, data_edgeon = plot_projection_diskrel(box, quant_dict[quant][0], plot_width, norm_L, args, 
                                                            quant_label=quant_dict[quant][0], 
                                                            unit=quant_dict[quant][8], 
                                                            clim=[quant_dict[quant][3], quant_dict[quant][4]] if quant_dict[quant][3] is not None else None,  
                                                            cmap=quant_dict[quant][6], 
                                                            takelog=quant_dict[quant][7], 
                                                            clabel=quant_dict[quant][9],
                                                            annotate_labels = [rf'SFR = {sfr:.1f} M$_\odot$/yr', rf'$\log$ (M$_*$/M$_\odot$) = {log_mstar:.1f}'],
                                                            return_frb=True,
                                                            do_plot=True,
                                                            ncell_buff = ncell_buff,
                                                            )
                    # --------making the 3D FRB------------
                    if not args.do_only_plot:
                        FRB = all_data[('gas', quant_dict[quant][0])].in_units(quant_dict[quant][2]).astype(np.float32)

                        # --------making the FITS ImageHDU for 3D FRB---------------
                        img_hdu = FITSImageData(FRB, ('gas', quant_dict[quant][1]))
                        header = img_hdu[0].header     
                        
                        for ind in range(3):
                            header[f'CDELT{ind+1}'] = args.kpc_per_pix
                            header[f'CUNIT{ind+1}'] = 'kpc'
                            header[f'NORMAL_UNIT_VECTOR{ind+1}'] = norm_L[ind]

                        header[f'SFR'] = 'NaN' if np.isnan(sfr) else sfr
                        header[f'SFRUNIT'] = 'Msun/yr'

                        header[f'LOG_MSTAR'] = 'NaN' if np.isnan(log_mstar) else log_mstar
                        header[f'MSTARUNIT'] = 'Msun'

                        header[f'REDSHIFT'] = args.current_redshift
                        
                        img_hdu_list.append(img_hdu[0])
                        
                        # --------making the FITS ImageHDU for 2D projected FRB: face on---------------                    
                        hdu_faceon = FITSImageData(data_faceon, f'{quant_dict[quant][1]} FACE ON PROJ')
                        header = hdu_faceon[0].header                        
        
                        header['WIDTH'] = (plot_width, 'kpc')
                        for ind in range(2):
                            header[f'CDELT{ind+1}'] = kpc_per_pix_proj
                            header[f'CUNIT{ind+1}'] = 'kpc'
                        
                        img_hdu_list.append(hdu_faceon[0])

                        # --------making the FITS ImageHDU for 2D projected FRB: face on---------------
                        hdu_edgeon = FITSImageData(data_edgeon, f'{quant_dict[quant][1]} EDGE ON PROJ')
                        header = hdu_edgeon[0].header       
                        
                        header['WIDTH'] = (plot_width, 'kpc')
                        for ind in range(2):
                            header[f'CDELT{ind+1}'] = kpc_per_pix_proj
                            header[f'CUNIT{ind+1}'] = 'kpc'
                        
                        img_hdu_list.append(hdu_edgeon[0])

                        # ------making the plots-----------
                        if args.plot_3d or args.plot_proj:
                            ax = fig.add_subplot(1, len(quant_arr), index + 1, projection='3d' if args.plot_3d else None)
                            if args.plot_3d: ax = plot_3d_frb(FRB, ax, args, label=quant_dict[quant][1], unit=quant_dict[quant][2], clim=[quant_dict[quant][3], quant_dict[quant][4]], cmap=quant_dict[quant][6])
                            elif args.plot_proj: ax = plot_proj_frb(FRB, ax, args, label=quant_dict[quant][1], unit=quant_dict[quant][2], clim=[quant_dict[quant][3], quant_dict[quant][4]], cmap=quant_dict[quant][6], hidey=index > 0)

                # ------saving fits file------------------
                if not args.do_only_plot:
                    primary = fits.PrimaryHDU(header=img_hdu_list[0].header, data=img_hdu_list[0].data)
                    combined_img_hdu = fits.HDUList([primary] + [fits.ImageHDU(h.data, header=h.header) for h in img_hdu_list[1:]])

                    combined_img_hdu.writeto(fitsname, overwrite=args.clobber)
                    myprint('Saved fits file as ' + fitsname, args)

                    # ------saving fig------------------
                    if args.plot_3d or args.plot_proj:
                        fig.savefig(figname)
                        myprint('Saved plot as ' + figname, args)

                if not ('pleiades' in args.system or args.hide_plot): plt.show(block=False)
                print_mpi('This snapshots completed in %s' % timedelta(seconds=(datetime.now() - start_time_this_snapshot).seconds), args)
            
            except Exception as e:
                print_mpi('Skipping ' + this_sim[1] + ' because ' + str(e), args)
                continue
                
        else:
            print('Skipping snapshot %s as %s already exists. Use --clobber to remake figure.' %(args.output, fitsname))
            continue

    # -----------------------------------------------------------------------------------
    if ncores > 1: print_master('Parallely: time taken for ' + str(total_snaps) + ' snapshots with ' + str(ncores) + ' cores was %s' % timedelta(seconds=(datetime.now() - start_time).seconds), args)
    else: print_master('Serially: time taken for ' + str(total_snaps) + ' snapshots with ' + str(ncores) + ' core was %s' % timedelta(seconds=(datetime.now() - start_time).seconds), args)
