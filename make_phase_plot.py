#!/usr/bin/env python3

"""

    Title :      make_phase_plot.py
    Notes :      Make a phase plot (T vs rho) color-coded by electron density using datashader for a given list of FOGGIE snapshots
                 Run this script in astroconda env
    Output :     Plot saved as png file
    Author :     Ayan Acharyya
    Started :    Nov 2026
    Examples :   run make_phase_plot.py --system ayan_pleiades --halo 8508 --upto_kpc 100 --do_all_sims --use_cen_smoothed
                 run make_phase_plot.py --system ayan_hd --halo 4123 --upto_kpc 100 --output RD0038,RD0042 --use_cen_smoothed --xcol density --ycol temp --colorcol metal
                 run make_phase_plot.py --system ayan_local --halo 8508 --upto_kpc 100 --output RD0027 --use_cen_smoothed
                 run make_phase_plot.py --system ayan_local --do_all_halos --upto_kpc 100 --use_cen_smoothed --weight edensity
                 run make_phase_plot.py --system ayan_local --halo 2878 --output RD0037 --upto_kpc 100 --use_cen_smoothed --weight edensity
"""
from foggie_header import *
from craft_utils import annotate_axes, save_fig
from yt.visualization.fits_image import FITSImageData
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpi4py import MPI

plt.rcParams['axes.linewidth'] = 1

start_time = datetime.now()

# -------------------------------------------------------------------------------
def get_df_from_ds(ds, args):
    '''
    Function to make a pandas dataframe from the yt dataset based on the given field list and color category
    The dataframe can be written to file for faster access in future
    This function is somewhat based on foggie.utils.prep_dataframe.prep_dataframe()
    :return: dataframe
    '''
    print('\nCreating datframe..')
    all_fields = [args.xcol, args.ycol, args.colorcol]
    if args.weight is not None:
        all_fields += [args.weight]
    df = pd.DataFrame()
    for index,field in enumerate(all_fields):
        print('Doing property: ' + field + ', which is ' + str(index + 1) + ' of the ' + str(len(all_fields)) + ' fields..')
        arr = ds[field_dict[field]].in_units(unit_dict[field]).ndarray_view()
        df[field] = arr

    return df

# --------------------------------------------------------------------------------
def make_datashader_plot_mpl(df, args):
    '''
    Function to make data shader plot of y_field vs x_field, colored in bins of color_field
    This function is based on foggie.render.shade_maps.render_image()
    :return figure
    '''
    print(f'Making datashader plot..')
    # --------to make the main datashader plot--------------------------
    fig, ax = plt.subplots(1, figsize=(6, 5), constrained_layout=True)
    color_key = [to_hex(item) for item in args.colormap.colors]
    if args.weight is not None and args.weight in df.columns:
        #if islog_dict[args.ycol]: df[args.ycol] = 10 ** df[args.ycol] # converting back to linear because they were already in log
        if islog_dict[args.weight]: df[args.weight] = 10 ** df[args.weight] # converting back to linear because they were already in log
        #df[args.ycol] = df[args.ycol] * df[args.weight] * len(df[args.weight]) / np.sum(df[args.weight])
        #if islog_dict[args.ycol]: df[args.ycol] = np.log10(df[args.ycol]) # converting back to log because the plot needs log

        artist = dsshow(df, dsh.Point(args.xcol, args.ycol), dsh.sum(args.weight), 
                        norm='log', cmap='rainbow', aspect = 'auto', ax=ax,
                        x_range=(args.xmin, args.xmax), y_range=(args.ymin, args.ymax)
                    )
    else:
        artist = dsshow(df, dsh.Point(args.xcol, args.ycol), dsh.mean(args.colorcol), 
                        norm='linear', cmap=color_key, aspect = 'auto', ax=ax,
                        x_range=(args.xmin, args.xmax), y_range=(args.ymin, args.ymax), vmin=args.cmin, vmax=args.cmax 
                    )
    # ------to make the axes-------------
    xlabel = f'{label_dict[args.xcol]} [{unit_dict[args.xcol]}]'
    ylabel = f'{label_dict[args.ycol]} [{unit_dict[args.ycol]}]'
    clabel = f'{label_dict[args.colorcol]} [{unit_dict[args.colorcol]}]'
    if islog_dict[args.xcol]: xlabel = 'Log ' + xlabel
    if islog_dict[args.ycol]: ylabel = 'Log ' + ylabel
    if islog_dict[args.colorcol]: clabel = 'Log ' + clabel

    ax = annotate_axes(ax, xlabel, ylabel, args=args)
    logT = 4.
    ax.axhline(logT, lw=1, ls='--', c='k', label=f'T=10^{logT:.0f} K')
    ax.legend(loc='lower left', fontsize = args.fontsize / args.fontfactor)

    # ------to make the colorbar axis-------------
    if args.weight is None or args.weight not in df.columns:
        cax_xpos, cax_ypos, cax_width, cax_height = 0.6, 0.9, 0.3, 0.035
        cax = fig.add_axes([cax_xpos, cax_ypos, cax_width, cax_height])
        plt.colorbar(artist, cax=cax, orientation='horizontal')

        cax.set_xticklabels(['%.0F' % index for index in cax.get_xticks()], fontsize=args.fontsize / args.fontfactor)
        fig.text(cax_xpos + cax_width / 2, cax_ypos + cax_height + 0.005, clabel, ha='center', va='bottom', fontsize=args.fontsize / args.fontfactor)

    # ---------to annotate and save the figure----------------------
    colorby_text = f'_colorby_{args.colorcolname}' if args.weight is None or args.weight not in df.columns else ''
    weightby_text = f'_weightby_{args.colorcolname}' if args.weight is not None and args.weight in df.columns else ''
    figname = Path(args.output_dir) / 'plots' / f'{args.output}_{args.halo}_phaseplot_{args.ycolname}_vs_{args.xcolname}{colorby_text}{weightby_text}{args.upto_text}.png'
    fig.savefig(figname)
    print(f'Saved plot as {figname}')    
    plt.show(block=False)

    return fig

# ----------to initialise dicts-----------------------
field_dict = {'rad': ('gas', 'radius_corrected'), 'metal': ('gas', 'metallicity'), 'vrad': ('gas', 'radial_velocity_corrected'), 'density':('gas', 'density'), 'temp':('gas', 'temperature'), 'edensity':('gas', 'El_number_density')}
unit_dict = {'rad': 'kpc', 'metal': r'Zsun', 'vrad': 'km/s', 'density':'g/cm**3', 'temp':'K', 'edensity':'cm**-3'}
label_dict = {'rad': 'Radius (kpc)', 'metal': 'Log Metallicity (Zsun)', 'vrad': 'Radial velocity (km/s)', 'density':'Gas density', 'temp':'Temperature', 'edensity':'Electron number density'}
bounds_dict = {'metal': (-2, 1), 'vrad': (-400, 400), 'density': (1e-30, 1e-22), 'temp': (1e1, 1e8), 'edensity': (1e-6, 1e-1)}
colormap_dict = {'vrad': outflow_inflow_discrete_cmap, 'temp': temperature_discrete_cmap, 'density': density_discrete_cmap, 'metal': metal_discrete_cmap, 'edensity': edensity_discrete_cmap}
islog_dict = {'rad': False, 'metal': True, 'vrad': False, 'density': True, 'temp': True, 'edensity': True}

# -----main code-----------------
if __name__ == '__main__':
    '''
    args = parse_args()
    if not args.keep: plt.close('all')

    # --------getting all snapshot IDs-------------

    total_snaps = len(list_of_sims)
    '''
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.system == "ayan_pleiades": args.code_dir = '/nobackupp19/aachary2/ayan_codes/foggie/foggie/'
    else: args.code_dir = '/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/'
    smooth_over_snap = 20
    
    # ----------getting list of snapshots-----------
    if args.do_all_halos:
        halos = ['8508', '5036', '5016', '4123', '2392', '2878']
    else:
        halos = args.halo_arr


    # ---------looping over halos-----------
    for index, thishalo in enumerate(halos):
        print_master(f'Starting halo {thishalo}..', args)
        args.halo = thishalo
        # --------preparing to read in this halo-----------------
        if args.do_all_sims:
            list_of_sims = get_all_sims_for_this_halo(args) # all snapshots of this particular halo
        elif args.output == 'RD0042': # i.e. the default, no output has been specified
            z_min, z_max, delta_z = 0, 2, 0.05
            redshift_list = np.arange(z_min, z_max + delta_z, delta_z)

            df_halo = pd.read_csv(args.code_dir + f'halo_infos/00{args.halo}/nref11c_nref9f/halo_cen_smoothed', sep=r'\s*\|\s*', engine='python')
            df_halo = df_halo.dropna(axis=1, how='all')[['snap', 'redshift']]
            output_list = []
            for redshift in redshift_list:
                idx = (df_halo['redshift'] - redshift).abs().idxmin()
                output_list.append(df_halo.loc[idx, 'snap'])
            output_list.sort()
            list_of_sims = list(itertools.product([args.halo], output_list))    
        else:
            list_of_sims = list(itertools.product([args.halo], args.output_arr))    

        # --------domain decomposition; for mpi parallelisation-------------
        total_snaps = len(list_of_sims)

        comm = MPI.COMM_WORLD
        ncores = comm.size
        rank = comm.rank
        print(f"RAW DEBUG: Process system PID {os.getpid()} reports Rank {rank} out of {ncores} total cores.")
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
            print('Doing snapshot ' + this_sim[1] + ' of halo ' + this_sim[0] + ' which is ' + str(index + 1) + ' out of the total ' + str(total_snaps) + ' snapshots...')

            # --------determining corresponding text suffixes and figname-------------
            args.fig_dir = args.output_dir + 'plots/'
            Path(args.fig_dir).mkdir(parents=True, exist_ok=True)

            args.fits_dir = args.output_dir + 'data/'
            Path(args.fits_dir).mkdir(parents=True, exist_ok=True)

            args.upto_text = '_upto%.1Fckpchinv' % args.upto_kpc if args.docomoving else '_upto%.1Fkpc' % args.upto_kpc

            # ------tailoring the simulation box for individual snapshot analysis--------
            if args.upto_kpc is not None:
                if args.docomoving: args.galrad = args.upto_kpc / (1 + args.current_redshift) / 0.695  # fit within a fixed comoving kpc h^-1, 0.695 is Hubble constant
                else: args.galrad = args.upto_kpc  # fit within a fixed physical kpc
            else:
                args.re = get_re_from_coldgas(args) if args.use_gasre else get_re_from_stars(ds, args)
                args.galrad = args.re * args.upto_re  # kpc

            # ----------to initialise axes limits and colormaps-----------------------
            args.xcolname = 'log_' + args.xcol if islog_dict[args.xcol] else args.xcol
            args.ycolname = 'log_' + args.ycol if islog_dict[args.ycol] else args.ycol
            args.colorcolname = 'log_' + args.colorcol if islog_dict[args.colorcol] else args.colorcol
            args.outfilename = Path(args.output_dir) / 'data' / f'{args.output}_{args.halo}_list_{args.xcol}_{args.ycol}_{args.colorcol}{args.upto_text}.fits'

            bounds_dict.update(rad=(0, np.ceil(args.galrad * (1 + args.current_redshift)) if args.docomoving else args.galrad)) # if comoving is used, put plot limits to maximum radius, i.e. physical radius at z = 0
            args.xmin = np.log10(bounds_dict[args.xcol][0]) if islog_dict[args.xcol] else bounds_dict[args.xcol][0]
            args.xmax = np.log10(bounds_dict[args.xcol][1]) if islog_dict[args.xcol] else bounds_dict[args.xcol][1]
            args.ymin = np.log10(bounds_dict[args.ycol][0]) if islog_dict[args.ycol] else bounds_dict[args.ycol][0]
            args.ymax = np.log10(bounds_dict[args.ycol][1]) if islog_dict[args.ycol] else bounds_dict[args.ycol][1]
            args.cmin = np.log10(bounds_dict[args.colorcol][0]) if islog_dict[args.colorcol] else bounds_dict[args.colorcol][0]
            args.cmax = np.log10(bounds_dict[args.colorcol][1]) if islog_dict[args.colorcol] else bounds_dict[args.colorcol][1]
            args.colormap = colormap_dict[args.colorcol]
            
            # -----------loading/building the dataframe-----------------
            if not os.path.exists(args.outfilename) or args.clobber:
                # -------loading in snapshot-------------------
                halos_df_name = args.code_path + 'halo_infos/00' + this_sim[0] + '/' + args.run + '/'
                halos_df_name += 'halo_cen_smoothed' if args.use_cen_smoothed else 'halo_c_v'
                try:
                    ds, refine_box = load_sim(args, region='refine_box', do_filter_particles=True, disk_relative=False, halo_c_v_name=halos_df_name)
                except FileNotFoundError:
                    print(f'This snapshot {args.halo}:{args.output} was not found, so skipping it..')
                    continue

                # --------assigning additional keyword args-------------
                args.current_redshift = ds.current_redshift
                args.current_time = ds.current_time.in_units('Gyr').tolist()

                # extract the required box
                box_width = args.galrad * 2  # in kpc
                box_width_kpc = ds.arr(box_width, 'kpc')
                box_center = ds.halo_center_kpc
                box = ds.r[box_center[0] - box_width_kpc / 2.: box_center[0] + box_width_kpc / 2., 
                        box_center[1] - box_width_kpc / 2.: box_center[1] + box_width_kpc / 2., 
                        box_center[2] - box_width_kpc / 2.: box_center[2] + box_width_kpc / 2., ]

                # ----------to make dataframe from dataset-----------------------
                df = get_df_from_ds(box, args)
                Table.from_pandas(df).write(args.outfilename, format='fits', overwrite=True)
                #df.to_parquet(str(Path(args.outfilename).with_suffix('.parquet')))
                print(f'Saved dataframe in {args.outfilename}')
            else:
                print(f'Reading from existing {args.outfilename}')
                df = Table.read(args.outfilename, format='fits', memmap=True).to_pandas()
                #df.read_parquet(str(Path(args.outfilename).with_suffix('.parquet')))

            # -----------prepping dataframe--------------
            print(f'Prepping dataframe..')
            for col in df.columns:
                if islog_dict[col]:
                    df[col] = np.log10(df[col])

            # ----------to make datashader plot-----------------------
            fig = make_datashader_plot_mpl(df, args)

            # --------------------
            if not ('pleiades' in args.system or args.hide_plot): plt.show(block=False)
            print('This snapshots completed in %s' % timedelta(seconds=(datetime.now() - start_time_this_snapshot).seconds))

    # -----------------------------------------------------------------------------------
    print('\nCompleted ' + str(total_snaps) + ' snapshots in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
