#
#	Script for estimating FRB host galaxy DM from simulated electron density cubes
#
#								AB, August 2024
#   Modified by AA in Apr 2026
#   
#   Examples of how to run (from within ipython):   
#   run dmplot.py --mode lsmzsfr --rangekpc 200 --reskpc 0.5 --z_range 0,6 --inc 0,90 --lsm 9.5,10.5 --lsfr 0,0.5 --resfile_prefix all_lsm
#   run dmplot.py --mode lsmzsfr --lsm 9.5,10.5 --lsfr 0,0.5 --resfile_prefix all_lsm
#   run dmplot.py --mode indi --lsm 9.5,10.0 --lsfr=-1,0
#   run dmplot.py --mode indi --lsm 9.5,10.0 --lsfr=-1,0 --multi_panel
#   run dmplot.py --mode halo --halo 5036 --lsm 9.5,10.0 --lsfr=-1,0
#   run dmplot.py --mode lsmzsfr --lsm 9.5,10.5
#   run dmplot.py --mode lsmzsfr --lsm all --multi_panel
#   run dmplot.py --mode lsmzsfr --lsm all
#   run dmplot.py --mode lsmzsfr --lsm all --inc 0,30
#   run dmplot.py --mode lsmzsfr --lsm all --inc 80,90

#	--------------------------	Import modules	---------------------------
from craft_utils import *
import plotfns as pfns

start_time = datetime.now()

# -----------------------------------------------------------------------------
def print_instructions ():

	#	Print instructions to terminal
	
	print("\n            You probably need some assistance here!\n")
	print("\n Arguments are     --- <mode> <zmin/zmax> <args.lsm_range[0]/args.lsm_range[1]> <args.lsfr_range[0]/args.lsfr_range[1]> <rangkpc> <reskpc> <args.inc_range[0]/args.inc_range[1]> <halo> <figname> <args.resfile_prefix>\n")
	print(" Supported Modes are --- indi     (Plot individual files)")
	print("                     --- halo     (Track a particular halo)")
	print("                     --- lsmzsfr  (Combine all wihin the lsm and z range)")
	
	
	print("\n            Now let's try again!\n")
	
	return(0)

# -----------------------------------------------------------------------------
def create_tiled_layout(nrows, ncols, fig_size=4):
    # Create the multi-multi-panel figure layout, to be used by pltdm_ind_imf
    fig = plt.figure(figsize=(ncols * fig_size * 1.5, nrows * fig_size))
    
    master_gs = fig.add_gridspec(nrows, ncols, wspace=0.3, hspace=0.3) # hspace/wspace control the gaps BETWEEN the triplets
    
    all_axes = []

    for r in range(nrows):
        for c in range(ncols):
            inner_gs = master_gs[r, c].subgridspec(1, 3, width_ratios=[40, 40, 2], wspace=0.1) # Gap inside the triplet
            
            ax1 = fig.add_subplot(inner_gs[0, 0])
            ax2 = fig.add_subplot(inner_gs[0, 1])
            ax3 = fig.add_subplot(inner_gs[0, 2])
            
            all_axes.append((ax1, ax2, ax3))
            
    return fig, all_axes

# -----------------------------------------------------------------------------
def execute_mode_indi(df_snap, args):
    '''
    Function to execute mode indi
    Returns nothing
    '''
    print("\nPloting individual snaps...\n")

    if args.multi_panel:
        nrows, ncols = get_grid_size(len(df_snap))
        fig, axes = plt.subplots(nrows, ncols, figsize=(10, 8))
        fig.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.98, wspace=0.01, hspace=0.01)

    for i, snap in df_snap.iterrows():
        thisfile = args.los_dir / f'{snap["snap"]}_{snap["halo"]}_FRB_El_number_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_150.npy'
        dm_arr	= np.load(thisfile)
        this_df = pd.DataFrame(dm_arr, columns=['inc', 'impf', 'distmaj', 'losdm'])
        print(f'{snap["snap"]}_{snap["halo"]}: Total number of LoS = {len(this_df)}')
        
        print("Plotting DMs within inclination ",args.inc_range[0], args.inc_range[1])
        this_df = this_df[this_df['inc'].between(args.inc_range[0], args.inc_range[1])]

        outfile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{snap["halo"]}_{snap["snap"]}'
        pars, epars, ax	= pfns.pltdm_ind_imf_1d(this_df, snap['log_star_mass'], snap['sfr'], args.lsm_range, outfile + '_1d', 3.0, hide=args.hide, bin_col='impf', data_col='losdm', given_ax=axes[i // ncols][i % ncols] if args.multi_panel else None)

        if args.multi_panel:
            if i // ncols < nrows - 1:
                ax.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                ax.set_xlabel('')
            if i % ncols > 0:
                ax.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                ax.set_ylabel('')

        # --------------initialise dataframe------------------
        lgsm	= np.log10(10.0 ** snap["log_gas_mass"] + 10.0 ** snap["log_star_mass"])
        df_out = pd.DataFrame({'lsm_bin': pd.Interval(args.lsm_range[0], args.lsm_range[1]),
                            'lsfr_bin': pd.Interval(args.lsfr_range[0], args.lsfr_range[1]),
                            'inc_bin': pd.Interval(args.inc_range[0], args.inc_range[1]),
                            'medlsm': snap["log_star_mass"],
                            'medsfr': snap["sfr"],
                            'medlgm': snap["log_gas_mass"],
                            'medlgsm': lgsm,
                            'ledssfr9': np.log10(snap["sfr"]) - snap["log_star_mass"] + 9,
                            'medgsfr9': np.log10(snap["sfr"]) - snap["log_gas_mass"] + 9,
                            'medgssfr9': np.log10(snap["sfr"]) - lgsm + 9,
                            'r0': pars[0],
                            'er0': epars[0],
                            'D0': pars[1],
                            'eD0': epars[1],
                            }, index=[0])
        
        outfile = f'{args.resfile_prefix}_indiv_allinc.txt'
        df_out.to_csv(outfile, mode='a', sep='\t', header=not os.path.exists(outfile), index=None)

    if args.multi_panel: save_fig(fig, args.fig_dir, f'{args.mode}_inc_{args.inc_range[0]}_{args.inc_range[1]}_multipanel_1d.pdf', args)

    return

# -----------------------------------------------------------------------------
def execute_mode_halo(df_snap, args):
    '''
    Function to execute mode halo
    Returns nothing
    '''
    print(f"\nTracking halo {args.halo} ...\n")
    df_snap = df_snap[df_snap["halo"].astype(str) == args.halo]

    if args.multi_panel:
        nrows, ncols = get_grid_size(len(df_snap))
        fig_1d, axes_1d = plt.subplots(nrows, ncols, figsize=(10,8))
        fig_1d.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.98, wspace=0.01, hspace=0.01)

        fig, axes = create_tiled_layout(nrows, ncols, fig_size=8)
        fig.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.98, wspace=0.01, hspace=0.01)

    for i, snap in df_snap.iterrows():
        thisfile = args.los_dir / f'{snap["snap"]}_{snap["halo"]}_FRB_El_number_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_150.npy'
        dm_arr	= np.load(thisfile)
        this_df = pd.DataFrame(dm_arr, columns=['inc', 'impf', 'distmaj', 'losdm'])
        this_df = this_df[this_df['inc'].between(args.inc_range[0], args.inc_range[1])]

        print(f'{snap["snap"]}_{snap["halo"]}: Total number of LoS = {len(this_df)}')
        print("Plotting DMs within inclination ",args.inc_range[0], args.inc_range[1])

        outfile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{snap["halo"]}_{snap["snap"]}'
        
        ax_list = pfns.pltdm_ind_imf(this_df, snap['log_star_mass'], snap['sfr'], args.inc_range, snap['redshift'], outfile, 3.0, hide=False, bin_col1='impf', bin_col2='distmaj', data_col='losdm', given_ax=axes[i] if args.multi_panel else None)
        pars, epars, ax_1d	= pfns.pltdm_ind_imf_1d(this_df, snap['log_star_mass'], snap['sfr'], args.lsm_range, outfile + '_1d', 3.0, hide=args.hide, bin_col='impf', data_col='losdm', given_ax=axes_1d[i // ncols][i % ncols] if args.multi_panel else None)

        if args.multi_panel:
            if i // ncols < nrows - 1:
                for ax in ax_list:
                    ax.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                    ax.set_xlabel('')
                ax_1d.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                ax_1d.set_xlabel('')
            if i % ncols > 0:
                for ax in ax_list:
                    ax.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                    ax.set_ylabel('')
                ax_1d.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                ax_1d.set_ylabel('')

    if args.multi_panel:
        save_fig(fig, args.fig_dir, f'{args.mode}_inc_{args.inc_range[0]}_{args.inc_range[1]}_{args.halo}_multipanel_imf_heatmap.pdf', args)
        save_fig(fig_1d, args.fig_dir, f'{args.mode}_inc_{args.inc_range[0]}_{args.inc_range[1]}_{args.halo}_multipanel_1d.pdf', args)

    return

# -----------------------------------------------------------------------------
def execute_mode_lsmzsfr(df_snap, args, given_ax=None):
    '''
    Function to execute mode lsmzsfr
    Returns nothing
    '''
    print(f"\nCombining all within lsm range {args.lsm_range} and z range {args.z_range} \n")

    combined_df = pd.DataFrame()
    
    for i, snap in df_snap.iterrows():
        thisfile = args.los_dir / f'{snap["snap"]}_{snap["halo"]}_FRB_El_number_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_150.npy'
        dm_arr	= np.load(thisfile)
        this_df = pd.DataFrame(dm_arr, columns=['inc', 'impf', 'distmaj', 'losdm'])
        this_df = this_df[this_df['inc'].between(args.inc_range[0], args.inc_range[1])]
        combined_df = pd.concat([combined_df, this_df], ignore_index=True)
        
    print(f"Total number of LoS = {len(combined_df)}")	
    print("Plotting DMs within inclination ",args.inc_range[0], args.inc_range[1])

    median_lsm = np.nanmedian(df_snap["log_star_mass"])
    median_lgm = np.nanmedian(df_snap["log_gas_mass"])
    median_sfr = np.nanmedian(df_snap["sfr"])
    lgsm	= np.log10(10.0 ** df_snap["log_gas_mass"] + 10.0 ** df_snap["log_star_mass"])
    
    outfile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{args.mode}_lsm_{args.lsm_range[0]}_{args.lsm_range[1]}_lsfr_{args.lsfr_range[0]}_{args.lsfr_range[1]}_1d'
    pars, epars, ax	= pfns.pltdm_ind_imf_1d(combined_df, median_lsm, median_sfr, args.lsm_range, outfile, 2.6, hide=args.hide, bin_col='impf', data_col='losdm', given_ax=given_ax)

    # --------------initialise dataframe------------------
    df_out = pd.DataFrame({'lsm_bin': pd.Interval(args.lsm_range[0], args.lsm_range[1]),
                           'lsfr_bin': pd.Interval(args.lsfr_range[0], args.lsfr_range[1]),
                           'inc_bin': pd.Interval(args.inc_range[0], args.inc_range[1]),
                           'medlsm': median_lsm,
                           'medsfr': median_sfr,
                           'medlgm': median_lgm,
                           'medlgsm': np.nanmedian(lgsm),
                           'ledssfr9': np.nanmedian(np.log10(df_snap["sfr"]) - df_snap["log_star_mass"]) + 9,
                           'medgsfr9': np.nanmedian(np.log10(df_snap["sfr"]) - df_snap["log_gas_mass"]) + 9,
                           'medgssfr9': np.nanmedian(np.log10(df_snap["sfr"]) - lgsm) + 9,
                           'r0': pars[0],
                           'er0': epars[0],
                           'D0': pars[1],
                           'eD0': epars[1],
                           }, index=[0])
    
    outfile = f'{args.resfile_prefix}_allinc.txt'
    df_out.to_csv(outfile, mode='a', sep='\t', header=not os.path.exists(outfile), index=None)

    return ax

# -----main code-----------------
if __name__ == '__main__':
    #	--------------------------	Read inputs	-------------------------------
    args = parse_args()
    if not args.keep: plt.close('all')
    out_dir = Path(f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}')
    out_dir.mkdir(exist_ok=True, parents=True)

    # ------------looping over inclination bins----------
    for index2, this_inc_bin in enumerate(args.inc_bins):
        print(f'\n\nRunning ({index2 + 1}/{len(args.inc_bins)}) for inclination bin {this_inc_bin}..\n')
        args.inc_range = this_inc_bin

        # ------------setup multi-panel figure if needed-------
        if args.multi_panel and args.mode == 'lsmzsfr':
            nrows, ncols = get_grid_size(len(args.lsm_bins))
            fig, axes = plt.subplots(nrows, ncols, figsize=(10, 8))
            fig.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.98, wspace=0.01, hspace=0.01)

        # ------------looping over stellar mass bins----------
        for index, this_lsm_bin in enumerate(args.lsm_bins):
            print(f'\t\nRunning ({index + 1}/{len(args.lsm_bins)}) for stellar mass bin {this_lsm_bin}..\n')
            args.lsm_range = this_lsm_bin
            #	-------------------------	Initialize	-----------------------------------
            filespecs =	args.data_dir / "lsm_sfr.txt"
            df_snap = pd.read_csv(filespecs, sep=r'\s+', engine='python')
            df_snap = df_snap.rename(columns={f'{df_snap.columns[0]}':f'{df_snap.columns[0][1:]}'})
            df_snap['log_sfr'] = np.log10(df_snap['sfr'])
            df_snap = df_snap[(df_snap['redshift'].between(args.z_range[0], args.z_range[1])) & 
                            (df_snap['log_star_mass'].between(args.lsm_range[0], args.lsm_range[1])) & 
                            (df_snap['log_sfr'].between(args.lsfr_range[0], args.lsfr_range[1]))].reset_index(drop=True)
            print (f'\tFound {len(df_snap)} snapshots, within log mass range {args.lsm_range}, log sfr range {args.lsfr_range} and redshift range {args.z_range}')
            if (len(df_snap) < 1): sys.exit('\tExiting because no snapshot found... ')
            
            #	-------------------------	Execute tasks	-------------------------------
            if (args.mode=='indi'):
                execute_mode_indi(df_snap, args)
            elif (args.mode=='halo'):
                execute_mode_halo(df_snap, args)            
            elif (args.mode=='lsmzsfr'):                
                # ---------------make the plots-----------------
                ax = execute_mode_lsmzsfr(df_snap, args, given_ax=axes[index // ncols][index % ncols] if args.multi_panel else None)

                if args.multi_panel:
                    if index // ncols < nrows - 1:
                        ax.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                        ax.set_xlabel('')
                    if index % ncols > 0:
                        ax.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                        ax.set_ylabel('')
            else:
                print("\n\tHmm...What mode is that again...?\n")

        if args.mode == 'lsmzsfr' and args.multi_panel:
            save_fig(fig, args.fig_dir, f'{args.mode}_inc_{args.inc_range[0]}_{args.inc_range[1]}_multipanel_1d.pdf', args)
    
    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))











































