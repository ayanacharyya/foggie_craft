#!/usr/bin/env python3

"""

    Title :      get_mass_sfr
    Notes :      Derive stellar masses and SFRs for a list of FOGGIE snapshots
    Output :     Pandas dataframe
    Author :     Ayan Acharyya
    Started :    Jan 2026
    Examples :   run get_mass_sfr.py --system ayan_pleiades --do_all_halos
                 run get_mass_sfr.py --system ayan_pleiades --halo 2878
                 run get_mass_sfr.py --system ayan_local --plot_sfh --fontsize 8
"""
from foggie_header import *
from craft_utils import annotate_axes, save_fig
start_time = datetime.now()

# -----------------------------------------------------------------------------
def get_sfr_df(args):
    '''
    Reads in the dataframe that contains SFR of each snapshot of a given halo
    Returns pandas dataframe
    '''
    sfr_filename = args.code_path + 'halo_infos/00' + args.halo + '/' + args.run + '/sfr'
    if os.path.exists(sfr_filename):
        print_master(f'Reading SFR history from {sfr_filename}', args)
        sfr_df = pd.read_table(sfr_filename, names=('output', 'redshift', 'sfr'), comment='#', delim_whitespace=True)
    else:
        print_master(f'Did not find {sfr_filename}, therefore will not include SFR', args)
        sfr_df = pd.DataFrame()

    return sfr_df

# -----------------------------------------------------------------------------
def get_mass_profile(args):
    '''
    Function to get the mass profile for a stars and gas, for a given output
    '''
    mass_filename = args.code_path + 'halo_infos/00' + args.halo + '/' + args.run + '/masses_z-less-2.hdf5'

    if os.path.exists(mass_filename):
        print_mpi(f'Reading in {mass_filename}', args)
        alldata = pd.read_hdf(mass_filename, key='all_data')
        thisdata = alldata[alldata['snapshot'] == args.output]

        if len(thisdata) == 0: # snapshot not found in masses less than z=2 file, so try greater than z=2 file
            mass_filename = mass_filename.replace('less', 'gtr')
            print_mpi(f'Could not find snapshot in previous file, now reading in {mass_filename}', args)
            alldata = pd.read_hdf(mass_filename, key='all_data')
            thisdata = alldata[alldata['snapshot'] == args.output]

            if len(thisdata) == 0: # snapshot still not found in file
                print_mpi('Snapshot not found in either file. Returning bogus mass', args)
                return np.nan

        mass_profile = thisdata[thisdata['radius'] <= args.diskrad]
        if len(mass_profile) == 0: # the smallest shell available in the mass profile is larger than the necessary radius within which we need the stellar mass
            print_mpi('Smallest shell avialable in mass profile is too large compared to args.diskrad. Returning nearest shell', args)
            return thisdata.head(1)
    else:
        print_mpi(f'File not found: {mass_filename}', args)
        return np.nan

    return mass_profile

# -----------------------------------------------------------------------------
def get_masses_and_re(args, get_re_using='gas_HI_mass'):
    '''
    Function to determine the stellar, gas, and halo (total) mass, for a given snapshot, which is defined as the mass contained within args.diskrad, which can either be a fixed absolute size in kpc OR = args.upto_re*Re
    Also determines the effective radius of stellar disk, based on the stellar mass profile
    Returns halo mass, gas mass, stellar mass, and stellar half-light radius
    '''

    mass_profile = get_mass_profile(args)

    if type(mass_profile) == pd.DataFrame:    
        mass_profile = mass_profile.sort_values('radius')

        mhalo = mass_profile['total_mass'].iloc[-1] # radius is in kpc, mass in Msun
        mgas = mass_profile['gas_mass'].iloc[-1] # radius is in kpc, mass in Msun
        mstar = mass_profile['stars_mass'].iloc[-1] # radius is in kpc, mass in Msun

        total_mass = mass_profile[get_re_using].iloc[-1]
        mask = mass_profile[get_re_using] <= total_mass/2
        if len(mass_profile[mask]) > 0:
            half_mass_radius = mass_profile[mask]['radius'].iloc[-1]
        else:
            print_mpi('Smallest shell avialable in mass profile is larger than half-mass. So returning the smallest shell as half-mass radius', args)
            half_mass_radius = mass_profile['radius'].iloc[0]
    else:
        mhalo, mgas, mstar, half_mass_radius = np.nan, np.nan, np.nan, np.nan

    return mhalo, mgas, mstar, half_mass_radius

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.system == "ayan_pleiades": args.code_dir = '/nobackupp19/aachary2/ayan_codes/foggie/foggie/'
    else: args.code_dir = '/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/'
    smooth_over_snap = 20
    
   # ---------initialising output dataframe-------------
    output_dfname = args.output_dir + 'data/lsm_sfr_masses_upto_disk_dummy.csv'
    columns = ['halo', 'snap', 'redshift', 'sfr', f'sfr_{int(5 * smooth_over_snap)}Myr', 'disk_rad', 'log_star_mass_from_snap', 'log_star_mass_from_profile', 'log_gas_mass_from_profile', 'half_mass_rad']

    # ----------getting list of snapshots-----------
    if args.do_all_halos:
        halos = ['8508', '5036', '5016', '4123', '2392', '2878']
    else:
        halos = args.halo_arr
    z_min, z_max, delta_z = 0, 2, 0.05
    redshift_list = np.arange(z_min, z_max + delta_z, delta_z)

    # ------------------reading existing file if any-------------------
    if os.path.exists(output_dfname):
        print_master(f'Reading existing df from {output_dfname}', args)
        df_existing = pd.read_csv(output_dfname, comment='#')
        existing_combos = set(zip(df_existing['halo'].astype(str), df_existing['snap'].astype(str)))
    else:
        existing_combos = set()

    # -------setup SFH figure-------
    if args.plot_sfh:
        print_master(f'Will plot SFH, so will not compute masses; to compute masses run without the --plot_sfh option', args)
        fig, axes = plt.subplots(len(halos), 1, figsize=(10, 7), sharex=True, constrained_layout=True)
        fig2, axes2 = plt.subplots(len(halos), 1, figsize=(10, 7), sharex=True, constrained_layout=True)
        smooth_over_nsnap_arr = [1, 20]
        col_arr = ['purple', 'orange', 'darkolivegreen', 'cornflowerblue', 'salmon']
    
    # ---------looping over halos-----------
    for index, thishalo in enumerate(halos):
        print_master(f'Starting halo {thishalo}..', args)
        args.halo = thishalo
        sfr_df  = get_sfr_df(args) # reading SFR df
        sfr_df[f'sfr_smooth{smooth_over_snap}'] = sfr_df['sfr'].rolling(window=smooth_over_snap, center=True).mean()

        # ---------make SFH plots-------------
        if args.plot_sfh:
            sfr_df['time'] = cosmo.age(sfr_df['redshift']).value # in Gyr

            # ------smoothing the sfh----------
            for index2, smooth_over_snap in enumerate(smooth_over_nsnap_arr):
                sfr_df[f'sfr_smooth{smooth_over_snap}'] = sfr_df['sfr'].rolling(window=smooth_over_snap, center=True).mean()
                
                # ---------plotting smoothed sfh---------
                axes[index].plot(sfr_df['time'], sfr_df[f'sfr_smooth{smooth_over_snap}'], c=col_arr[index2], lw=0.5 + 0.1*index2, label=f'{int(smooth_over_snap * 5)} Myr scale')
                
                # -------annotating plots-----------
                axes[index] = annotate_axes(axes[index], r'Time [Gyr]', r'SFR [M$_\odot$/yr]', fontsize=args.fontsize, label=f'{args.halo}: {halo_dict[args.halo]}', hide_xaxis=index < len(halos) - 1)
            
                # -----plotting SFR at our redshifts----------
                sfr_list = [sfr_df.loc[(sfr_df['redshift'] - this_redshift).abs().idxmin(), f'sfr_smooth{smooth_over_snap}'] for this_redshift in redshift_list]
                axes2[index].plot(redshift_list, np.log10(sfr_list), 'o-', c=col_arr[index2], lw=0.5, label=f'{int(smooth_over_snap * 5)} Myr scale')

                # -------annotating plots-----------
                axes2[index] = annotate_axes(axes2[index], r'Redshift', r'$\log$ SFR' +'\n'+ r'(M$_\odot$ yr$^{-1}$)', fontsize=args.fontsize, label=f'{args.halo}: {halo_dict[args.halo]}', hide_xaxis=index < len(halos) - 1)

            # -----making legends---------
            if index == 3: axes[index].legend(loc='upper right', fontsize = args.fontsize)
            if index == 5: axes2[index].legend(loc='lower right', fontsize = args.fontsize, ncol=2)

            # ------plotting redshifts of interest---------
            for this_redshift in redshift_list:
                this_time = cosmo.age(this_redshift).value
                axes[index].axvline(this_time, c='k', ls='dotted', lw=0.5)
                if index == 0:
                    axes[index].text(this_time, axes[index].get_ylim()[1], f'z={this_redshift:.2f}', c='k', va='top', ha='right', rotation=90, fontsize=args.fontsize)

        # --------determining snapshots-----------
        else:
            # --------preparing to read in this halo-----------------
            df = pd.read_csv(args.code_dir + f'halo_infos/00{args.halo}/nref11c_nref9f/halo_cen_smoothed', sep=r'\s*\|\s*', engine='python')
            df = df.dropna(axis=1, how='all')[['snap', 'redshift']]
            output_list = []
            for redshift in redshift_list:
                idx = (df['redshift'] - redshift).abs().idxmin()
                output_list.append(df.loc[idx, 'snap'])
            
            outputs_existing = [snap for snap in output_list if (str(thishalo), str(snap)) in existing_combos]
            outputs_todo = list(set(output_list) - set(outputs_existing))

            # --------domain decomposition; for mpi parallelisation-------------
            total_snaps = len(outputs_todo)

            comm = MPI.COMM_WORLD
            ncores = comm.size
            rank = comm.rank
            print_master(f'For halo {args.halo}: need to do only {len(outputs_todo)} of original {len(output_list)} outputs because {len(outputs_existing)} already exist', args)
            if len(outputs_todo) == 0:
                print_master(f'Therefore continuing to next halo..', args)
                continue
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
            #print_mpi('Operating on snapshots ' + str(core_start + 1) + ' to ' + str(core_end + 1) + ', i.e., ' + str(core_end - core_start + 1) + ' out of ' + str(total_snaps) + ' snapshots', args)
            print_mpi(f'Operating on snapshots {thishalo:} {outputs_todo[core_start : core_end + 1]} ' + str(core_start + 1) + ' to ' + str(core_end + 1) + ', i.e., ' + str(core_end - core_start + 1) + ' out of ' + str(total_snaps) + ' snapshots', args) ##
            for index in range(core_start + args.start_index, core_end + 1):
                start_time_this_snapshot = datetime.now()
                args.output = outputs_todo[index]
                args.halo = thishalo
                print_mpi('Doing snapshot ' + args.output + ' of halo ' + args.halo + ' which is ' + str(index + 1 - core_start) + ' out of the total ' + str(core_end - core_start + 1) + ' snapshots...', args)

                # -----------determining SFR and redshift--------------------
                if args.output in sfr_df['output'].values:
                    sfr = sfr_df[sfr_df['output'] == args.output][f'sfr'].values[0]
                    sfr_smooth = sfr_df[sfr_df['output'] == args.output][f'sfr_smooth{smooth_over_snap}'].values[0]
                    args.current_redshift = sfr_df[sfr_df['output'] == args.output]['redshift'].values[0]
                    df_row = pd.DataFrame([[args.halo, args.output, args.current_redshift, sfr, sfr_smooth]], columns=['halo', 'snap', 'z', 'sfr', 'sfr_smooth'])
                    '''
                    try:
                        # ------determining extent for computing mass--------
                        args.diskrad, log_mstar_from_snap = get_stellar_mass(args)
                        if np.isnan(log_mstar_from_snap):
                            raise ValueError

                        # ------determining stellar mass--------                
                        _, mgas, mstar, half_mass_radius = get_masses_and_re(args, get_re_using='gas_HI_mass')
                        log_mstar_from_profile = np.log10(mstar)
                        log_mgas_from_profile = np.log10(mgas)

                        # ------appending to dataframe----------
                        df_row = pd.DataFrame([[args.halo, args.output, args.current_redshift, sfr, sfr_smooth, args.diskrad, log_mstar_from_snap, log_mstar_from_profile, log_mgas_from_profile, half_mass_radius]], columns=columns)
                    except Exception as e:
                        print_mpi(f'Snapshot {args.halo}:{args.output} failed due to {e}, therefore skipping, and not adding to dataframe', args)
                        continue
                    '''
                else:
                    print_mpi(f'Snapshot {args.halo}:{args.output} is not in sfr_df, therefore skipping, and not adding to dataframe', args)
                    continue

                # -----------saving dataframe---------------
                file_exists = os.path.exists(output_dfname)
                df_row.to_csv(output_dfname, index=False, mode='a' if file_exists else 'w', header=not file_exists)
                print_mpi(f'Completed snapshot {args.halo}:{args.output} in {timedelta(seconds=(datetime.now() - start_time_this_snapshot).seconds)}', args)
    
    # --------saving final figure---------
    if args.plot_sfh:
        save_fig(fig, Path(args.output_dir) / 'plots', 'all_halos_sfh.png', args)
        save_fig(fig2, Path(args.output_dir) / 'plots', 'all_halos_sfr_to_pick.png', args)

        print_master('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds), args)
    else:
        print_master(f'Saved dataframe as {output_dfname}', args)

        if ncores > 1: print_master('Parallely: time taken for ' + str(total_snaps) + ' snapshots with ' + str(ncores) + ' cores was %s' % timedelta(seconds=(datetime.now() - start_time).seconds), args)
        else: print_master('Serially: time taken for ' + str(total_snaps) + ' snapshots with ' + str(ncores) + ' core was %s' % timedelta(seconds=(datetime.now() - start_time).seconds), args)
