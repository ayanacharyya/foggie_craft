#!/usr/bin/env python3

"""

    Title :      get_hotgas_sfr
    Notes :      Derive hot gas masses and SFRs for a list of FOGGIE snapshots
    Output :     Pandas dataframe
    Author :     Ayan Acharyya
    Started :    Nov 2026
    Examples :   run get_hotgas_sfr.py --system ayan_pleiades --do_all_halos
                 run get_hotgas_sfr.py --system ayan_pleiades --halo 2878
                 run get_hotgas_sfr.py --system ayan_local --do_all_halos
"""
from foggie_header import *
from craft_utils import annotate_axes, save_fig
from mpi4py import MPI

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
def get_cgm_gas_masses(args):
    '''
    Computes the gas mass and hot (T > 1e4 K) gas mass within the refine box for a given FOGGIE snapshot
    Returns the two masses (in linear scale)
    '''
    halos_df_name = args.code_path + 'halo_infos/00' + args.halo + '/' + args.run + '/'
    halos_df_name += 'halo_cen_smoothed' if args.use_cen_smoothed else 'halo_c_v'
    ds, refine_box = load_sim(args, region='refine_box', do_filter_particles=True, disk_relative=False, halo_c_v_name=halos_df_name)

    gas_masses = refine_box[('gas', 'mass')]
    gas_mass = np.sum(gas_masses).in_units('Msun')
    gas_temps = refine_box[('gas', 'temperature')]
    hot_mask = gas_temps > ds.quan(1e4, 'K')
    hot_gas_mass = np.sum(gas_masses[hot_mask]).in_units('Msun')
    del gas_masses, gas_temps, hot_mask
    '''
    gas_mass = refine_box.quantities.total_quantity(('gas', 'mass')).in_units('Msun')
    hot_gas = refine_box.cut_region(["obj[('gas', 'temperature')].in_units('K') > 1e4"])
    hot_gas_mass = hot_gas.quantities.total_quantity(('gas', 'mass')).in_units('Msun')
    '''
    return gas_mass, hot_gas_mass

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.system == "ayan_pleiades": args.code_dir = '/nobackupp19/aachary2/ayan_codes/foggie/foggie/'
    else: args.code_dir = '/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/'
    smooth_over_snap = 20
    
   # ---------initialising output dataframe-------------
    output_dfname = args.output_dir + 'data/hgm_sfr_upto_cgm.csv'
    columns = ['halo', 'snap', 'redshift', 'sfr', f'sfr_{int(5 * smooth_over_snap)}Myr', 'log_gas_mass', 'log_hot_gas_mass']

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
        df_existing['time'] = cosmo.age(df_existing['redshift']).value # in Gyr
        existing_combos = set(zip(df_existing['halo'].astype(str), df_existing['snap'].astype(str)))
    else:
        existing_combos = set()

    # -------setup SFH figure-------
    if args.plot_sfh and os.path.exists(output_dfname):
        print_master(f'Will plot SFH and hot gas mass, so will not compute masses; to compute masses run without the --plot_sfh option', args)
        fig, axes = plt.subplots(len(halos), 1, figsize=(10, 1.5 + 1 * len(halos)), sharex=True)
        fig.subplots_adjust(top=0.98, bottom=0.2 if len(halos) == 1 else 0.07, left=0.07, right=0.93, wspace=0., hspace=0.05)
        axes = np.atleast_1d(axes)
        smooth_over_nsnap_arr = [1, 20]
        col_arr = ['purple', 'darkolivegreen', 'orange', 'cornflowerblue', 'salmon']
    
    # ---------looping over halos-----------
    for index, thishalo in enumerate(halos):
        print_master(f'Starting halo {thishalo}..', args)
        args.halo = thishalo

        # ---------make SFH plots-------------
        if args.plot_sfh:
            df = df_existing[df_existing['halo'].astype(str) == thishalo]
            ax = axes[index]
            # ---------plotting smoothed sfh---------
            ax.plot(df['time'], df['sfr'], c=col_arr[0], lw=0.5 , label=f'5 Myr scale')
            ax.plot(df['time'], df[f'sfr_{int(5 * smooth_over_snap)}Myr'], c=col_arr[1], lw=1.5, label=f'{int(smooth_over_snap * 5)} Myr scale')
            ax.scatter(df['time'], df['sfr'], c=col_arr[0], lw=0.5, s=50)

            # -------annotating plots-----------
            ax = annotate_axes(ax, r'Time [Gyr]', '', args=args, label=f'{args.halo}: {halo_dict[args.halo]}', hide_xaxis=index < len(halos) - 1)
            #ax.set_ylim(0.1, 300)
            ax.set_ylim(1, 30)
            ax.set_yscale('log')

            # -----making legends---------
            if index == 0:
                ax.legend(loc='upper right', fontsize = args.fontsize)
                fig.text(0.01, 0.5, r'SFR [M$_\odot$/yr]', ha='left', va='center', fontsize=args.fontsize, rotation=90)
            if index == len(halos) - 1:
                ax.set_xlim(1, 14.)

            # -----plotting gas mass with SFH----------
            ax2 = ax.twinx()
            ax2.plot(df['time'], df['log_hot_gas_mass'], c=col_arr[2], lw=0.5)
            ax2.scatter(df['time'], df['log_hot_gas_mass'], c=col_arr[2], lw=0.5, s=50)

            # -------annotating plots-----------
            ax2.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, prune='both'))
            ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
            ax2.tick_params(axis='y', which='major', labelsize=args.fontsize, labelright=True)
            ax2.set_ylim(8, 11)
            if index == 0:
                fig.text(0.99, 0.5, r'$\log$ Gas mass [M$_*$/M$_\odot$] at T > 10$^4$ K', ha='right', va='center', fontsize=args.fontsize, rotation=90)

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
            outputs_todo.sort()

            sfr_df  = get_sfr_df(args) # reading SFR df
            sfr_df[f'sfr_smooth{smooth_over_snap}'] = sfr_df['sfr'].rolling(window=smooth_over_snap, center=True).mean()

            # --------domain decomposition; for mpi parallelisation-------------
            total_snaps = len(outputs_todo)

            comm = MPI.COMM_WORLD
            ncores = comm.size
            rank = comm.rank
            print(f"RAW DEBUG: Process system PID {os.getpid()} reports Rank {rank} out of {ncores} total cores.")
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
            print_mpi(f'halo {thishalo}: outputs_todo = {outputs_todo}', args) ##
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
                    print_mpi(f'{args.halo}, {args.output}, {args.current_redshift}, {sfr}, {sfr_smooth}', args) ##
                    
                    try:
                        # ------determining gas mass--------                
                        gas_mass, hot_gas_mass = get_cgm_gas_masses(args)
                        log_gas_mass = np.log10(gas_mass)
                        log_hot_gas_mass = np.log10(hot_gas_mass)

                        # ------appending to dataframe----------
                        df_row = pd.DataFrame([[args.halo, args.output, args.current_redshift, sfr, sfr_smooth, log_gas_mass, log_hot_gas_mass]], columns=columns)
                    except Exception as e:
                        print_mpi(f'Snapshot {args.halo}:{args.output} failed due to {e}, therefore skipping, and not adding to dataframe', args)
                        continue
                    
                else:
                    print_mpi(f'Snapshot {args.halo}:{args.output} is not in sfr_df, therefore skipping, and not adding to dataframe', args)
                    continue

                # -----------saving dataframe---------------
                file_exists = os.path.exists(output_dfname)
                df_row.to_csv(output_dfname, index=False, mode='a' if file_exists else 'w', header=not file_exists)
                print_mpi(f'Completed snapshot {args.halo}:{args.output} in {timedelta(seconds=(datetime.now() - start_time_this_snapshot).seconds)}', args)
    
    # --------saving final figure---------
    if args.plot_sfh:
        save_fig(fig, Path(args.output_dir) / 'plots', 'all_halos_sfh_hgm.png', args)
        print_master('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds), args)
    else:
        print_master(f'Saved dataframe as {output_dfname}', args)

        if ncores > 1: print_master('Parallely: time taken for ' + str(total_snaps) + ' snapshots with ' + str(ncores) + ' cores was %s' % timedelta(seconds=(datetime.now() - start_time).seconds), args)
        else: print_master('Serially: time taken for ' + str(total_snaps) + ' snapshots with ' + str(ncores) + ' core was %s' % timedelta(seconds=(datetime.now() - start_time).seconds), args)
