'''
	Script for estimating FRB host galaxy DM from simulated electron density cubes

								AB, August 2024
  Modified by AA in Apr 2026
  
  Examples of how to run (from within ipython):   
  run dmplot.py --mode lsmzsfr --rangekpc 200 --reskpc 0.5 --z_range 0,6 --inc 0,90 --lsm 9.5,10.5 --lsfr 0,0.5 --resfile_prefix all_lsm
  run dmplot.py --mode lsmzsfr --lsm 9.5,10.5 --lsfr 0,0.5 --resfile_prefix all_lsm
  run dmplot.py --mode indi --lsm 9.5,10.0 --lsfr=-1,0
  run dmplot.py --mode indi --lsm 9.5,10.0 --lsfr=-1,0 --multi_panel
  run dmplot.py --mode halo --halo 5036 --lsm 9.5,10.0 --lsfr=-1,0
  run dmplot.py --mode halo --halo 5036 --z_range 0,2 --fontsize 15
  run dmplot.py --mode plot_halo --halo 5036 --z_range 0,2 --fontsize 15
  run dmplot.py --mode indi --lsm 10.5,11.0 --z_range 0,2 --fontsize 15
  run dmplot.py --mode plot_indi --lsm 10.5,11.0 --lsfr=-0.3,0.4,1.1,1.8 --z_range 0,2 --fontsize 15
  run dmplot.py --mode lsmzsfr --lsm 9.5,10.5
  run dmplot.py --mode lsmzsfr --lsm all --multi_panel
  run dmplot.py --mode lsmzsfr --lsm all
  run dmplot.py --mode lsmzsfr --lsm all --inc 0,30
  run dmplot.py --mode lsmzsfr --lsm all --inc 80,90
  run dmplot.py --mode lsmzsfr --lsm 8.5,9.0,9.5,10.0,10.5,11.0,11.5 --lsfr=-3.0,0.0,0.5,1.0,1.5,2.0 --resfile_prefix binby_lsm_lsfr --multi_panel --fontsize 6
'''

#	--------------------------	Import modules	---------------------------
from craft_utils import *
setup_plot_style()
import plotfns as pfns
from plot_sfms import read_snap_list

start_time = datetime.now()

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
    param_outfile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_indiv_allinc.txt'
    if os.path.exists(param_outfile) and args.clobber:
        os.remove(param_outfile)
        print(f'Removed existing {param_outfile} because --clobber was used.')

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

        outfile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{snap["halo"]}_{snap["snap"]}'
        
        #multifit_par_filename = f'{args.fig_dir}/{Path(args.resfile_prefix).stem}_z_{args.z_range[0]}_{args.z_range[1]}_DM0_r0_vs_lsm_inc_{args.inc_range[0]}_{args.inc_range[1]}_multifit_params.txt'
        multifit_par_filename = f'{args.fig_dir}/{Path(args.resfile_prefix).stem}_z_0.0_1.0_DM0_r0_vs_lsm_inc_0.0_90.0_multifit_params.txt'
        if not os.path.exists(multifit_par_filename): multifit_par_filename = None
        #try:
        pars, epars, ax	= pfns.pltdm_ind_imf_1d(this_df, snap['log_star_mass'], snap['sfr'], args.lsm_range, outfile + '_1d', 2.6, hide=args.hide, bin_col='impf', data_col='losdm', given_ax=axes[i // ncols][i % ncols] if args.multi_panel else None, fortalk=args.fortalk, multifit_par_filename=multifit_par_filename)
        '''
        except:
            print(f'Failing DM profile fit for {snap["snap"]}')
            continue
        '''
        if args.multi_panel:
            if i // ncols < nrows - 1:
                ax.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                ax.set_xlabel('')
            if i % ncols > 0:
                ax.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                ax.set_ylabel('')

        # --------------initialise dataframe------------------
        lgsm	= np.log10(10.0 ** snap["log_gas_mass"] + 10.0 ** snap["log_star_mass"])
        df_out = pd.DataFrame({'inc_bin': pd.Interval(args.inc_range[0], args.inc_range[1]),
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
        
        df_out.to_csv(param_outfile, mode='a', sep='\t', header=not os.path.exists(param_outfile), index=None)

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

        outfile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{snap["halo"]}_{snap["snap"]}'
        
        ax_list = pfns.pltdm_ind_imf(this_df, snap['log_star_mass'], snap['sfr'], args.inc_range, snap['redshift'], outfile, 3.0, hide=False, bin_col1='impf', bin_col2='distmaj', data_col='losdm', given_ax=axes[i] if args.multi_panel else None, fortalk=args.fortalk)
        pars, epars, ax_1d	= pfns.pltdm_ind_imf_1d(this_df, snap['log_star_mass'], snap['sfr'], args.lsm_range, outfile + '_1d', 3.0, hide=args.hide, bin_col='impf', data_col='losdm', given_ax=axes_1d[i // ncols][i % ncols] if args.multi_panel else None, fortalk=args.fortalk)

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
    
    outfile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{args.mode}_lsm_{args.lsm_range[0]}_{args.lsm_range[1]}_lsfr_{args.lsfr_range[0]}_{args.lsfr_range[1]}_1d'
    pars, epars, ax	= pfns.pltdm_ind_imf_1d(combined_df, median_lsm, median_sfr, args.lsm_range, outfile, 2.6, hide=args.hide, bin_col='impf', data_col='losdm', given_ax=given_ax, nobj=len(df_snap), lsfr_lims=args.lsfr_range if len(args.lsfr_bins) > 1 else None, fortalk=args.fortalk)

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
                           'ngal': len(df_snap),
                           }, index=[0])
    
    outfile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_allinc.txt'
    df_out.to_csv(outfile, mode='a', sep='\t', header=not os.path.exists(outfile), index=None)

    return ax

# ------------------------------------------------------------------------------------------------
def plot_dm_impfac_halo_combined(df_snap, args, cmap='viridis'):
    '''
    Plot DM vs Impact factor in a single panel, for a given halo
    Saves plot
    Returns axis handle
    '''
    # ---------setup figure---------------
    fig, ax = plt.subplots(1, figsize=(8, 5))
    fig.subplots_adjust(left=0.12, bottom=0.12, right=0.99, top=0.98)

    norm = mplcolors.Normalize(vmin=df_snap['redshift'].min(), vmax=df_snap['redshift'].max())
    sm = mpl_cm.ScalarMappable(cmap=plt.get_cmap(cmap), norm=norm)

    # -----------loop through mass bins--------------------
    for index, snap in df_snap.iterrows(): 
        infile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{snap["halo"]}_{snap["snap"]}_1d.npy'
        col = sm.to_rgba(snap['redshift'])

        data_arr = np.load(infile) # data_arr is of the format [impx, dmavg, dmlower, dmhier]
        ax.errorbar(data_arr[0], data_arr[1], yerr=[data_arr[2], data_arr[3]], c=col, fmt='o-', lw=2, markersize=15, capsize=4, alpha=0.5)

    ax.set_xscale("log")
    ax.set_xticks(impbinegs[1:],impbinegs[1:])

    if not args.set_ylin:
        ax.set_yscale("log")
        ax.set_yticks(dm_ticks, dm_ticks)
        ax.set_ylim(ymin=0.7)
    
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label('Redshift', fontsize=args.fontsize)
    cbar.ax.tick_params(labelsize=args.fontsize)

    ax = annotate_axes(ax, "Impact factor (kpc)", "DM (pc cm$^{-3}$)", args=args, clabel='Redshift', set_ticks=False)
    ax.text(0.95, 0.95, f'Halo {args.halo}', c='k', fontsize=args.fontsize, ha='right', va='top', transform=ax.transAxes)

    save_fig(fig, args.fig_dir, f'DM_vs_impfact_halo_{args.halo}_inc{args.inc_range[0]}-{args.inc_range[1]}.pdf', args)
    plt.show(block=False)

    return ax

# ------------------------------------------------------------------------------------------------
def plot_dm_impfac_indi_combined(df_snap, args, cmap='viridis', colorcol='redshift'):
    '''
    Plot DM vs Impact factor in a single panel, for a list of stellar mass and sfr ranges
    Saves plot
    Returns axis handle
    '''
    # ---------setup figure---------------
    fig, ax = plt.subplots(1, figsize=(8, 5))
    fig.subplots_adjust(left=0.12, bottom=0.12, right=0.99, top=0.98)

    norm = mplcolors.Normalize(vmin=df_snap[colorcol].min(), vmax=df_snap[colorcol].max())
    sm = mpl_cm.ScalarMappable(cmap=plt.get_cmap(cmap), norm=norm)

    # -----------loop through mass bins--------------------
    for index, snap in df_snap.iterrows(): 
        infile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{snap["halo"]}_{snap["snap"]}_1d.npy'
        col = sm.to_rgba(snap[colorcol])

        data_arr = np.load(infile) # data_arr is of the format [impx, dmavg, dmlower, dmhier]
        ax.errorbar(data_arr[0], data_arr[1], yerr=[data_arr[2], data_arr[3]], c=col, fmt='o-', lw=2, markersize=15, capsize=4, alpha=0.5)

    ax.set_xscale("log")
    ax.set_xticks(impbinegs[1:],impbinegs[1:])

    if not args.set_ylin:
        ax.set_yscale("log")
        ax.set_yticks(dm_ticks, dm_ticks)
        ax.set_ylim(ymin=0.7)
    
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(colorcol, fontsize=args.fontsize)
    cbar.ax.tick_params(labelsize=args.fontsize)

    ax = annotate_axes(ax, "Impact factor (kpc)", "DM (pc cm$^{-3}$)", args=args, clabel=colorcol, set_ticks=False)
    ax.text(0.95, 0.95, rf'{args.lsm_range[0]} < $\log$(M/M$_\odot$) < {args.lsm_range[1]}', c='k', fontsize=args.fontsize, ha='right', va='top', transform=ax.transAxes)
    ax.text(0.95, 0.85, rf'{args.lsfr_range[0]} < $\log$(SFR/M$_\odot$ yr$^{-1}$) < {args.lsfr_range[1]} [{len(df_snap)}]', c='k', fontsize=args.fontsize, ha='right', va='top', transform=ax.transAxes)

    save_fig(fig, args.fig_dir, f'DM_vs_impfact_indi_lsm_bin_{args.lsm_range[0]}_{args.lsm_range[1]}_lsfr_bin_{args.lsfr_range[0]}_{args.lsfr_range[1]}_zrange_{args.z_range[0]}_{args.z_range[1]}_inc{args.inc_range[0]}-{args.inc_range[1]}.pdf', args)
    plt.show(block=False)

    return ax

# -----main code-----------------
if __name__ == '__main__':
    #	--------------------------	Read inputs	-------------------------------
    args = parse_args()
    if not args.keep: plt.close('all')
    out_dir = Path(f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{args.inc_range[0]}_{args.inc_range[1]}')
    out_dir.mkdir(exist_ok=True, parents=True)

    # ------------looping over inclination bins----------
    for index3, this_inc_bin in enumerate(args.inc_bins):
        print(f'\n\nRunning ({index3 + 1}/{len(args.inc_bins)}) for inclination bin {this_inc_bin}..\n')
        args.inc_range = this_inc_bin

        # ------------setup multi-panel figure if needed-------
        if args.multi_panel and args.mode == 'lsmzsfr':
            nrows, ncols = get_grid_size(len(args.lsm_bins) * len(args.lsfr_bins))
            fig, axes = plt.subplots(nrows, ncols, figsize=(8, 8))
            fig.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.98, wspace=0.01, hspace=0.01)

        # ------------looping over SFR bins----------
        for index2, this_lsfr_bin in enumerate(args.lsfr_bins):
            print(f'\n\tRunning ({index2 + 1}/{len(args.lsfr_bins)}) for SFR bin {this_lsfr_bin}..\n')
            args.lsfr_range = this_lsfr_bin

            # ------------looping over stellar mass bins----------
            for index, this_lsm_bin in enumerate(args.lsm_bins):
                print(f'\n\t\tRunning ({index + 1}/{len(args.lsm_bins)}) for stellar mass bin {this_lsm_bin}..\n')
                args.lsm_range = this_lsm_bin

                if args.multi_panel:
                    nrow = index if len(args.lsfr_bins) > 1 else index // ncols
                    ncol = index2 if len(args.lsfr_bins) > 1 else index % ncols
                #	-------------------------	Initialize	-----------------------------------
                df_snap = read_snap_list(args)
                if (len(df_snap) < 1):
                    print('\t\tNo snapshot found. Continuing to next loop iteration... ')
                    if args.multi_panel:
                        axes[nrow][ncol].remove()
                    continue
                
                #	-------------------------	Execute tasks	-------------------------------
                if (args.mode=='indi'):
                    execute_mode_indi(df_snap, args)
                
                elif (args.mode=='halo') :
                    execute_mode_halo(df_snap, args)     
                    
                elif (args.mode=='lsmzsfr'):                
                    # ---------------make the plots-----------------
                    ax = execute_mode_lsmzsfr(df_snap, args, given_ax=axes[nrow][ncol] if args.multi_panel else None)

                    if args.multi_panel:
                        if nrow < nrows - 1:
                            ax.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                            ax.set_xlabel('')
                        if ncol > 0:
                            ax.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                            ax.set_ylabel('')
                elif args.mode != 'plot_halo':
                    print("\n\tHmm...What mode is that again...?\n")

                if args.mode == 'halo' or args.mode == 'plot_halo':
                    df_snap = df_snap[df_snap["halo"].astype(str) == args.halo]
                    ax = plot_dm_impfac_halo_combined(df_snap, args)

                elif args.mode == 'indi' or args.mode == 'plot_indi':
                    #ax = plot_dm_impfac_indi_combined(df_snap, args, colorcol='redshift')
                    #ax = plot_dm_impfac_indi_combined(df_snap, args, colorcol='log_star_mass')
                    ax = plot_dm_impfac_indi_combined(df_snap, args, colorcol='log_sfr')


        if args.mode == 'lsmzsfr' and args.multi_panel:
            save_fig(fig, args.fig_dir, f'{Path(args.resfile_prefix).stem}_z_{args.z_range[0]}_{args.z_range[1]}_{args.mode}_inc_{args.inc_range[0]}_{args.inc_range[1]}_multipanel_1d.pdf', args)
        
    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))











































