#!/usr/bin/env python3

"""
    Title :      get_foggie_metallicity_profile
    Notes :      Plot metallicity profile of a given list of filenames (FOGGIE snapshots)
    Output :     Pandas dataframe
    Author :     Ayan Acharyya
    Started :    29-04-2026
    Examples :   run get_foggie_metallicity_profile.py --system ayan_ssd --halo 8508 --run nref11c_nref9f --output RD0027 --upto_kpc 20
                 run get_foggie_metallicity_profile.py --system ayan_ssd --halo 2392 --run 2162_39 --output DD1360 --upto_kpc 20
"""
from craft_header import *
from craft_utils import *
setup_plot_style()
from foggie_header import *
from make_3D_FRB_electron_density import plot_projection_diskrel, get_AM_vector
start_time = datetime.now()

# ------------------------------------------------------------------
def get_halo_coords(halo_c_v_filename, snap):
    '''
    Reads in a given halo_c_v file, and finds the central coordinate of a given snapshot
    Returns center as tuple
    '''
    df = pd.read_csv(halo_c_v_filename, sep='|', skipinitialspace=True)

    # We drop columns that are all NaN and strip the whitespace from headers
    df = df.dropna(axis=1, how='all')
    df.columns = df.columns.str.strip()
    df['name'] = df['name'].str.strip()

    row = df[df['name'] == snap]

    if not row.empty:
        # Return as a tuple (x_c, y_c, z_c)
        return (row.iloc[0]['x_c'], row.iloc[0]['y_c'], row.iloc[0]['z_c'])
    else:
        print(f"Output {snap} not found in {halo_c_v_filename}")
        return None

# ------------------------------------------------------------------
def my_foggie_load(snap, halo_center, refine_width_chkpc=200.):
    '''
    Reads in a given FOGGIE output
    Returns yt dataset and a refine region that is <refine_width_chkpc> h^-1 comoving kpc on a side
    '''
    print (f'Opening snapshot {snap}')
    # ------------load dataset and estimate domain size--------------
    ds = yt.load(snap)
    refine_box = ds.all_data()

    proper_box_size = get_proper_box_size(ds) # get the proper size of the computational domain (NOT the refine region) 
    ds.omega_baryon = ds.parameters['CosmologyOmegaMatterNow']-ds.parameters['CosmologyOmegaDarkMatterNow']

    # ----------estimate center------------
    if halo_center[0] < 1: # the coordinates are in code units
        ds.halo_center_code = ds.arr(np.array(halo_center), 'code_length')
        ds.halo_center_kpc  = ds.arr(np.array(halo_center)*proper_box_size, 'kpc')
    else: # the coordinates must be in kpc units
        ds.halo_center_kpc = ds.arr(np.array(halo_center), 'kpc')
        ds.halo_center_code = ds.arr(np.array(halo_center)/proper_box_size, 'code_length')
    
    # --------estimate bulk velocity----------------
    sp = ds.sphere(ds.halo_center_kpc, (3., 'kpc'))
    bulk_vel = sp.quantities.bulk_velocity(use_gas=False,use_particles=True,particle_type='all').to('km/s')
    ds.halo_velocity_kms = bulk_vel
  
    # --------compute additional fields--------------
    ds.add_field(('gas','vx_corrected'), function=vx_corrected, units='km/s', take_log=False, sampling_type='cell')
    ds.add_field(('gas', 'vy_corrected'), function=vy_corrected, units='km/s', take_log=False, sampling_type='cell')
    ds.add_field(('gas', 'vz_corrected'), function=vz_corrected, units='km/s', take_log=False, sampling_type='cell')
    ds.add_field(('gas', 'vel_mag_corrected'), function=vel_mag_corrected, units='km/s', take_log=False, sampling_type='cell')
    ds.add_field(('gas', 'radius_corrected'), function=radius_corrected, units='kpc', take_log=False, force_override=True, sampling_type='cell')
    ds.add_field(('gas', 'cell_mass_msun'), function=cell_mass_msun, units='Msun', take_log=True, force_override=True, sampling_type='cell')

    # # filter particles into star and dm
    # filter_particles(refine_box, filter_particle_types = ['stars', 'young_stars', 'dm'])

    # ds.add_field(('stars', 'radius_corrected'), function=radius_corrected_stars, units='kpc', take_log=False, force_override=True, sampling_type='particle')
    # ds.add_field(('young_stars', 'radius_corrected'), function=radius_corrected_young_stars, units='kpc', take_log=False, force_override=True, sampling_type='particle')
    # ds.add_field(('dm', 'radius_corrected'), function=radius_corrected_dm, units='kpc', take_log=False, force_override=True, sampling_type='particle')

    # -------------------get refine region-----------------
    zsnap = ds.get_parameter('CosmologyCurrentRedshift') 
    ds.refine_width_kpc = refine_width_chkpc / (1 + zsnap) / 0.695
    ds.refine_width_code = ds.arr(ds.refine_width_kpc/ proper_box_size, 'code_length')

    refine_box = ds.r[ds.halo_center_code[0] - ds.refine_width_code/2 : ds.halo_center_code[0] + ds.refine_width_code/2, 
                      ds.halo_center_code[1] - ds.refine_width_code/2 : ds.halo_center_code[1] + ds.refine_width_code/2, 
                      ds.halo_center_code[2] - ds.refine_width_code/2 : ds.halo_center_code[2] + ds.refine_width_code/2]

    return ds, refine_box

# ----------------------------------------------------------------
def make_df_from_box(box, args):
    '''
    Function to convert a given dataset to pandas dataframe
    :return: modified dataframe
    '''
    rho_cut = get_density_cut(args.current_time)  # based on Cassi's CGM-ISM density cut-off
    box = box.cut_region(['obj["gas", "density"] > %.1E' % rho_cut])
    print('Imposing a density criteria to get ISM above density', rho_cut, 'g/cm^3')

    df = pd.DataFrame()
    fields = ['radius', 'metal', 'density'] # only the relevant properties

    for index, field in enumerate(fields):
        print(f'Doing property: {field} ({(index + 1)}/{len(fields)})')
        df[field] = box[field_dict[field]].in_units(unit_dict[field]).ndarray_view()

    df = df.rename(columns={'radius': 'radius_kpc'})
    df = df.dropna()
    df = df.sort_values(by='radius_kpc')

    return df

# ----------global variables------------
quant_dict = {'density':['density', 'Gas density', 'Msun/pc**3', -2.5, 2.5, 'cornflowerblue', density_color_map], 
              'el_density':['El_number_density', 'Electron density', 'cm**-3', -6, -1, 'cornflowerblue', e_color_map],
              'metal':['metallicity', 'Gas metallicity', r'Zsun', -1.7, 0.7, 'cornflowerblue', metal_color_map],
              } # for each quantity: [yt field, label in plots, units, lower limit in log, upper limit in log, color for scatter plot, colormap]

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()  # default simulation to work upon when comand line args not provided
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2
    if not 'pleiades' in args.system: args.output_dir = args.output_dir.replace('CRAFT', 'sharing')

    # -----------determining directories----------------
    args.fig_dir = args.output_dir + 'metallicity_plots/'
    Path(args.fig_dir).mkdir(parents=True, exist_ok=True)
    
    args.text_dir = Path(args.output_dir) / 'metallicity_txtfiles'
    args.text_dir.mkdir(exist_ok=True, parents=True)
    
    args.upto_text = '_upto%.1Fckpchinv' % args.upto_kpc if args.docomoving and args.upto_kpc is not None else '_upto%.1Fkpc' % args.upto_kpc if args.upto_kpc is not None else f'_upto{args.re:.1f}re'
    df_snap_filename = args.text_dir / f'met_profile_{args.halo}_{args.run}_{args.output}{args.upto_text}.txt'

    # ---------------getting the metallicity profile----------------
    if not os.path.exists(df_snap_filename) or args.clobber:
        print(f']\n{df_snap_filename} does not exist. Creating afresh..')

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

        # ----------------making projection plots-------------------
        norm_L = get_AM_vector(ds, radius_kpc=3., use_particles=False)
        quant_arr = ['density', 'metal']
        for quant in quant_arr:
            fig_diskrel = plot_projection_diskrel(box, quant_dict[quant][0], box_width, norm_L, args, quant_label=quant_dict[quant][0], unit='Msun/pc**2' if quant == 'density' else quant_dict[quant][2], clim=[quant_dict[quant][3], quant_dict[quant][4]],  cmap=quant_dict[quant][6])

        # ---------get profile-------------
        df = make_df_from_box(box, args)

        df.to_csv(df_snap_filename, sep='\t', index=None)
        print(f'Saved profile to file {df_snap_filename}')
    else:
        print(f'\nReading from existing file {df_snap_filename}')
        df = pd.read_table(df_snap_filename, delim_whitespace=True, comment='#')

    # --------preparing profile------------
    quant = 'metal'
    weight_by = None # choose between 'density' and None
    met_profile_upto = 10. # kpc

    if weight_by is not None:
        weight_sum = np.sum(df[weight_by])
        colname = f'weighted_{quant}'
        df[colname] = df[quant] * df[weight_by] / weight_sum
        weightby_text = f'_wtby_{weight_by}'
    else:
        colname = quant
        weightby_text = ''
    
    df[f'log_{colname}'] = np.log10(df[colname])
    df = df[df['radius_kpc'] < met_profile_upto]
    
    # ---------------smoothing the profile----------------------
    ncells = len(df) // 100
    if ncells %2 == 0: ncells += 1
    df[f'smoothed_log_{colname}'] = df[f'log_{colname}'].rolling(window=ncells, center=True, min_periods=1).mean()

    # ----------------making profile plots-------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(df['radius_kpc'], df[f'smoothed_log_{colname}'], c=quant_dict[quant][5], lw=1)
    ax = annotate_axes(ax, 'Radius (kpc)', f'Log {quant_dict[quant][1]} ({quant_dict[quant][2]})', args=args)
    save_fig(fig, Path(args.fig_dir), f'met_profile_{args.halo}_{args.run}_{args.output}_upto{met_profile_upto:.1f}kpc{weightby_text}.png')

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
