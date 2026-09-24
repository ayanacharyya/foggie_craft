'''
	Script for estimating FRB host galaxy DM from simulated electron density cubes

								AB, August 2024
  Modified by AA in Apr 2026
  
  Examples of how to run (from within ipython):   
  run dmplot.py --mode lsmzsfr --rangekpc 200 --reskpc 0.5 --z_range 0,6 --inc 0,90 --lsm 9.5,10.5 --lsfr 0,0.5 --resfile_prefix all_lsm
  run dmplot.py --mode lsmzsfr --lsm 9.5,10.5 --lsfr 0,0.5 --resfile_prefix all_lsm

  run dmplot.py --mode indi --lsm 9.5,10.0 --lsfr=-1,0
  run dmplot.py --mode indi --lsm 9.5,10.0 --lsfr=-1,0 --multi_panel
  run dmplot.py --mode indi --z_range 0,2 --inc 0,90 --hide --clobber
  run dmplot.py --mode plot_indi --lsm 10.5,11.0 --lsfr=-0.3,0.4,1.1,1.8 --z_range 0,2 --fontsize 15

  run dmplot.py --mode halo --halo 5036 --lsm 9.5,10.0 --lsfr=-1,0
  run dmplot.py --mode halo --halo 5036 --z_range 0,2 --fontsize 15
  run dmplot.py --mode plot_halo --halo 5036 --z_range 0,2 --fontsize 15

  run dmplot.py --mode proj --halo 5036 --z_range 0,0.5
  run dmplot.py --mode proj --z_range 0,2 --inc 80,90 --plot_all --hide --clobber
  run dmplot.py --mode plot_2d_radius_ratio --z_range 0,2 --fontsize 12
  run dmplot.py --mode plot_2d_param_comp --z_range 0,2 --fontsize 12
  run dmplot.py --mode plot_1d_param_comp --z_range 0,2 --fontsize 12
  run dmplot.py --mode plot_1d_2d_comp --z_range 0,2 --fontsize 12

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
    param_outfile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_indiv_allinc.csv'
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
        
        multifit_par_filename = f'{args.fig_dir}/{Path(args.resfile_prefix).stem}_z_{args.z_range[0]}_{args.z_range[1]}_DM0_r0_vs_lsm_inc_{args.inc_range[0]}_{args.inc_range[1]}_multifit_params.txt'
        if not os.path.exists(multifit_par_filename): multifit_par_filename = None
        
        pars, epars, ax	= pfns.pltdm_ind_imf_1d(this_df, snap['log_star_mass'], snap['sfr'], args.lsm_range, outfile + '_1d', 2.6, hide=args.hide, bin_col='impf', data_col='losdm', given_ax=axes[i // ncols][i % ncols] if args.multi_panel else None, fortalk=args.fortalk, multifit_par_filename=multifit_par_filename, redshift=snap['redshift'])

        if args.multi_panel:
            if i // ncols < nrows - 1:
                ax.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                ax.set_xlabel('')
            if i % ncols > 0:
                ax.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                ax.set_ylabel('')

        # --------------initialise dataframe------------------        
        fit_dict = {'inc_bin': f'{args.inc_range[0]}_{args.inc_range[1]}',
                    'r0': pars[0],
                    'er0': epars[0],
                    'D0': pars[1],
                    'eD0': epars[1],
                    }
        
        combined_row = {**snap.to_dict(), **fit_dict}
        df_out = pd.DataFrame(combined_row, index=[0])

        df_out.to_csv(param_outfile, index=None, mode='a', header=not os.path.exists(param_outfile))
        if len(df_snap) > 10 and not args.multi_panel: plt.close('all')

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
        this_df['distmin'] = np.sqrt(this_df['impf'] ** 2 - this_df['distmaj'] ** 2)
        this_df = this_df[this_df['inc'].between(args.inc_range[0], args.inc_range[1])]

        print(f'{snap["snap"]}_{snap["halo"]}: Total number of LoS = {len(this_df)}')
        print("Plotting DMs within inclination ",args.inc_range[0], args.inc_range[1])

        outfile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{snap["halo"]}_{snap["snap"]}'
        
        ax_list = pfns.pltdm_ind_imf_2d(this_df, snap['log_star_mass'], snap['sfr'], args.inc_range, snap['redshift'], outfile, 3.0, hide=False, bin_col1='distmin', bin_col2='distmaj', data_col='losdm', given_ax=axes[i] if args.multi_panel else None, fortalk=args.fortalk)
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

# -----------------------------------------------------------------------------
def execute_mode_projection(df_snap, args):
    '''
    Function to execute mode projection
    Returns dataframe with fitted parameters
    '''
    outdir = Path(f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{args.inc_range[0]}_{args.inc_range[1]}')

    if not args.plot_all:
        df_snap = df_snap[df_snap["halo"].astype(str) == args.halo]
        print (f'\t\tFound {len(df_snap)} snapshots, for halo {args.halo}')

        df_snap = df_snap[df_snap["snap"].astype(str) == args.output] ##
        print (f'\t\tFound {len(df_snap)} snapshots, for output {args.output}')

    compiled_rows = [] # to store the fitted parameters later in a separate file

    for i, snap in df_snap.iterrows():
        thisfile = args.los_dir / f'{snap["snap"]}_{snap["halo"]}_FRB_El_number_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_150.npy'
        dm_arr	= np.load(thisfile)
        this_df = pd.DataFrame(dm_arr, columns=['inc', 'impf', 'distmaj', 'losdm'])
        this_df['distmin'] = np.sqrt(this_df['impf'] ** 2 - this_df['distmaj'] ** 2)
        this_df = this_df[this_df['inc'].between(args.inc_range[0], args.inc_range[1])]

        print(f'{snap["snap"]}_{snap["halo"]}: Total number of LoS = {len(this_df)}')
        print("Plotting DMs within inclination ",args.inc_range[0], args.inc_range[1])

        outfile = str(outdir / f'{snap["halo"]}_{snap["snap"]}')
        
        data_col, bin_col1, bin_col2 = 'losdm', 'distmin', 'distmaj'
        popt, perr, rx0, e_rx0, ry0, e_ry0 = pfns.pltdm_ind_imf_2d(this_df, snap['log_star_mass'], snap['sfr'], args.inc_range, snap['redshift'], outfile, 3.2, hide=args.hide, bin_col1=bin_col1, bin_col2=bin_col2, data_col=data_col, given_ax=axes[i] if args.multi_panel else None, fortalk=args.fortalk)

        fit_dict = {f'{bin_col1}0': popt[0], f'e_{bin_col1}0': perr[0],
                    f'{bin_col2}0': popt[1], f'e_{bin_col2}0': perr[1],
                    f'{data_col}0': popt[2], f'e_{data_col}0': perr[2],
                    f'{bin_col1}0_indep': rx0, f'e_{bin_col1}0_indep': e_rx0,
                    f'{bin_col2}0_indep': ry0, f'e_{bin_col2}0_indep': e_ry0,
                    }
        
        combined_row = {**snap.to_dict(), **fit_dict}
        compiled_rows.append(combined_row)
    
    df_results = pd.DataFrame(compiled_rows)
    output_df = outdir / 'projection_fit.csv'
    
    if not os.path.exists(output_df) or args.clobber:
        df_results.to_csv(output_df, index=None)
    else:
        df_results.to_csv(output_df, index=None, mode='a', header=False)

    return df_results

# -----------------------------------------------------------------------------
def execute_mode_plot_2d_fit_param_comparison(args):
    '''
    Function to execute mode plot_param_comp, which reads in the 2D fitted parameter file and plots rx0 vs ry0, or D0 vs r0
    Saves the plots
    Returns nothing
    '''
    # -----------setting up figure-----------------
    inc_ranges = [[0, 30], [60, 90], [80, 90]]
    #inc_ranges = [[60, 90]]

    col1, col2, figlabel = 'distmin0', 'distmaj0', 'rmin_rmaj' # this is for the scenario where the minor and major r0 have been fitted together as a 2D profile
    #col1, col2, figlabel = 'distmin0_indep', 'distmaj0_indep', 'rmin_rmaj_indep' # this is for the scenario where the minor and major r0 have been fitted independently (individually)    
    rlimits = [0.5, 300]

    # col1, col2, figlabel = 'losdm0', 'r_sq', 'rsq_dm' # 
    # dlimits = [7e0, 2e4]
    rlimits = [1e-1, 5e2]

    # --------------dictionaries for labels and colors-------------
    fixed_color = 'cornflowerblue'
    label_dict = {'distmaj0': r'r$_{0,maj}$ (kpc)',
                    'distmin0': r'r$_{0,min}$ (kpc)',
                    'distmaj0_indep': r'r$_{0,maj}$ (kpc)',
                    'distmin0_indep': r'r$_{0,min}$ (kpc)',
                    'losdm0': r'D$_{0,2D}$ pc cm$^{-3}$',
                    'r_sq': r'$\sqrt{r_{0,maj}^2 + r_{0,min}^2}$ (kpc)',
                    'log_star_mass': r'$\log{(M_*/M_\odot)}$',
                    'log_sfr': r'$\log{SFR/M_\odot yr^{-1}}$',
                    'log_ssfr': r'$\log{sSFR/yr^{-1}}$',
                    'redshift': 'Redshift',
                    }

    #colorby_col_arr = [None]
    colorby_col_arr = [None, 'log_star_mass', 'log_sfr', 'log_ssfr', 'redshift']

    # ---------looping over color cols---------------
    for colorby_col in colorby_col_arr:
        print(f'\n\nDoing color by {colorby_col}..')

        # -----------setting up figure-----------------
        fig, axes = plt.subplots(1, len(inc_ranges), figsize=(4 * len(inc_ranges), 4), layout='constrained', sharey=True)
        axes = np.atleast_1d(axes)

        # ---------looping over inc ranges----------
        for index, inc_range in enumerate(inc_ranges):
            outdir = Path(f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{inc_range[0]:.1f}_{inc_range[1]:.1f}')
            input_filename = args.data_dir / outdir / 'projection_fit.csv'
            print(f'For inc_range {inc_range} ({index + 1}/{len(inc_ranges)}): Reading projection fit parameters from {input_filename}..')
            df = pd.read_csv(input_filename)
            
            # ----computing r sqaured---------
            quant = ((unp.uarray(df['distmin0'], df['e_distmin0']) ** 2 + unp.uarray(df['distmaj0'], df['e_distmaj0']) ** 2) ** 0.5)
            df['r_sq'] = unp.nominal_values(quant)
            df['e_r_sq'] = unp.std_devs(quant)

            # -------computing sSFR---
            df['log_ssfr'] = df['log_sfr'] - df['log_star_mass']

            # ---------doing the plot------
            axes[index].errorbar(df[col1], df[col2], xerr=df[f'e_{col1}'], yerr=df[f'e_{col2}'], fmt='none', color=fixed_color, lw=0.5, alpha=0.8, capsize=2, zorder=-10)
            im = axes[index].scatter(df[col1], df[col2], c=fixed_color if colorby_col is None else df[colorby_col], s=10, lw=0.5, ec='k', alpha=0.8)

            # -------axes limits-----------
            axes[index].set_xscale('log')
            axes[index].set_yscale('log')

            if not ('losdm' in col1 or 'losdm' in col2):
                axes[index].plot([axes[index].get_xlim()[0], axes[index].get_xlim()[1]], [axes[index].get_xlim()[0], axes[index].get_xlim()[1]], c='k', ls='dashed', lw=1)
                
                mad = median_abs_deviation(df[col2] - df[col1])
                axes[index].text(0.05, 0.8, f'MAD={mad:.2f}', color='k', ha='left', va='top', transform=axes[index].transAxes, fontsize=args.fontsize)
            
                axes[index].set_xlim(rlimits[0], rlimits[1])
                axes[index].set_ylim(rlimits[0], rlimits[1])
            else:
                axes[index].set_xlim(dlimits[0], dlimits[1])
                axes[index].set_ylim(rlimits[0], rlimits[1])
            
            axes[index] = annotate_axes(axes[index], label_dict[col1], label_dict[col2], args=args, xloc=0.65 if 'losdm' in col1 or 'losdm' in col2 else 0.05, label=rf'{inc_range[0]}$^\circ$ $< i <$ {inc_range[1]}$^\circ$', hide_xaxis=False, hide_yaxis=index, bbox=False, set_ticks=False)

        # --------color axis-----------
        if colorby_col is not None:
            cbar = fig.colorbar(
                im, 
                ax=axes,          # Pass the entire array/list of axes here
                location='top',   # Forces it above the subplots
                orientation='horizontal', # Ensures the colorbar orientation is horizontal
                shrink=1.,       # Optional: scales width (1.0 = 100% width of the axes grid)
                pad=0.02,          # Optional: spacing between colorbar and subplots top edge
                aspect = 70,       # higher value for thinner colorbar
            )
            cbar.set_label(label_dict[colorby_col], labelpad=10, fontsize=args.fontsize)

        # ------------saving the figure----------------------
        figname = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_{len(inc_ranges)}_inc_ranges_2D_fit_{figlabel}.png'
        if colorby_col is not None:
            figname = figname.replace('.png', f'_colby_{colorby_col}.png')

        save_fig(fig, args.plot_dir, figname, args=args)
             
    return

# -----------------------------------------------------------------------------
def execute_mode_plot_1d_fit_param_comparison(args):
    '''
    Function to execute mode plot_param_comp, which reads in the 2D fitted parameter file and plots rx0 vs ry0, or D0 vs r0
    Saves the plots
    Returns nothing
    '''
    # -----------setting up figure-----------------
    #inc_ranges = [[0, 30], [60, 90], [80, 90]]
    inc_ranges = [[0, 90]]
    col1, col2, figlabel = 'D0', 'r0', 'r_dm' # 
    dlimits = [7e0, 1e3]
    rlimits = [5e-1, 5e2]

    # --------------dictionaries for labels and colors-------------
    fixed_color = 'cornflowerblue'
    label_dict = {'log_star_mass': r'$\log{(M_*/M_\odot)}$',
                    'log_sfr': r'$\log{SFR/M_\odot yr^{-1}}$',
                    'log_ssfr': r'$\log{sSFR/yr^{-1}}$',
                    'redshift': 'Redshift',
                    'r0': r'r$_{0}$ (kpc)',
                    'D0': r'D$_{0}$ pc cm$^{-3}$',                    
                    }

    colorby_col_arr = [None]
    #colorby_col_arr = [None, 'log_star_mass', 'log_sfr', 'log_ssfr', 'redshift']

    # ---------looping over color cols---------------
    for colorby_col in colorby_col_arr:
        print(f'\n\nDoing color by {colorby_col}..')

        # -----------setting up figure-----------------
        fig, axes = plt.subplots(1, len(inc_ranges), figsize=(4 * len(inc_ranges), 3.6), layout='constrained', sharey=True)
        axes = np.atleast_1d(axes)

        # ---------reading in input file with fitted parameters for 1D fit (with all inc)-----------
        input_filename_1dfit = args.data_dir / f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_indiv_allinc.csv'
        df_1d_allinc = pd.read_csv(input_filename_1dfit)
        df_1d_allinc = df_1d_allinc.rename(columns={'er0':'e_r0', 'eD0':'e_D0'})

        # ---------looping over inc ranges----------
        for index, inc_range in enumerate(inc_ranges):
            print(f'For inc_range {inc_range} ({index + 1}/{len(inc_ranges)})..')
            df = df_1d_allinc[df_1d_allinc['inc_bin'].astype(str) == f'{inc_range[0]:.1f}_{inc_range[1]:.1f}'].drop(columns='inc_bin')

            # -------computing sSFR---
            df['log_ssfr'] = df['log_sfr'] - df['log_star_mass']

            # ---------doing the plot------
            axes[index].errorbar(df[col1], df[col2], xerr=df[f'e_{col1}'], yerr=df[f'e_{col2}'], fmt='none', color=fixed_color, lw=0.5, alpha=0.8, capsize=2, zorder=-10)
            im = axes[index].scatter(df[col1], df[col2], c=fixed_color if colorby_col is None else df[colorby_col], s=10, lw=0.5, ec='k', alpha=0.8)

            # -------axes limits-----------
            axes[index].set_xscale('log')
            axes[index].set_yscale('log')

            axes[index].set_xlim(dlimits[0], dlimits[1])
            axes[index].set_ylim(rlimits[0], rlimits[1])
            
            axes[index] = annotate_axes(axes[index], label_dict[col1], label_dict[col2], args=args, xloc=0.65 if 'losdm' in col1 or 'losdm' in col2 else 0.65, label=rf'{inc_range[0]}$^\circ$ $< i <$ {inc_range[1]}$^\circ$', hide_xaxis=False, hide_yaxis=index, bbox=False, set_ticks=False)

        # --------color axis-----------
        if colorby_col is not None:
            cbar = fig.colorbar(
                im, 
                ax=axes,          # Pass the entire array/list of axes here
                location='top',   # Forces it above the subplots
                orientation='horizontal', # Ensures the colorbar orientation is horizontal
                shrink=1.,       # Optional: scales width (1.0 = 100% width of the axes grid)
                pad=0.02,          # Optional: spacing between colorbar and subplots top edge
                aspect = 70,       # higher value for thinner colorbar
            )
            cbar.set_label(label_dict[colorby_col], labelpad=10, fontsize=args.fontsize)

        # ------------saving the figure----------------------
        figname = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_{len(inc_ranges)}_inc_ranges_1D_fit_{figlabel}.png'
        if colorby_col is not None:
            figname = figname.replace('.png', f'_colby_{colorby_col}.png')

        save_fig(fig, args.plot_dir, figname, args=args)
             
    return

# -----------------------------------------------------------------------------
def execute_mode_plot_1d_2d_comparison(args):
    '''
    Function to execute mode plot_1d_2d_comp, which reads in the 1D and 2D fitted parameter files and plots rx0 vs ry0, or D0 vs r0
    Saves the plots
    Returns nothing
    '''
    # -----------setting up figure-----------------
    inc_ranges = [[0, 30], [60, 90], [80, 90]]
    #inc_ranges = [[60, 90]]
    #col1, col2, figlabel = 'r_sq', 'r0', 'r2d_r1d' # this is for the scenario where the minor and major r0 have been fitted together as a 2D profile
    #limits = [0.5, 300]

    col1, col2, figlabel = 'losdm0', 'D0', 'dm2d_dm1d' # 
    limits = [7e0, 2e4]

    # --------------dictionaries for labels and colors-------------
    fixed_color = 'cornflowerblue'
    label_dict = {'distmaj0': r'r$_{0,maj}$ (kpc)',
                    'distmin0': r'r$_{0,min}$ (kpc)',
                    'distmaj0_indep': r'r$_{0,maj}$ (kpc)',
                    'distmin0_indep': r'r$_{0,min}$ (kpc)',
                    'losdm0': r'D$_{0, 2D}$ pc cm$^{-3}$',
                    'r_sq': r'$\sqrt{r_{0,maj}^2 + r_{0,min}^2}$ (kpc)',
                    'log_star_mass': r'$\log{(M_*/M_\odot)}$',
                    'log_sfr': r'$\log{SFR/M_\odot yr^{-1}}$',
                    'log_ssfr': r'$\log{sSFR/yr^{-1}}$',
                    'redshift': 'Redshift',
                    'r0': r'2 $\times$ r$_{0}$ (kpc)',
                    'D0': r'D$_{0}$ pc cm$^{-3}$',
                    }

    colorby_col_arr = [None]
    #colorby_col_arr = [None, 'log_star_mass', 'log_sfr', 'log_ssfr', 'redshift']

    # ---------looping over color cols---------------
    for colorby_col in colorby_col_arr:
        print(f'\n\nDoing color by {colorby_col}..')

        # -----------setting up figure-----------------
        fig, axes = plt.subplots(1, len(inc_ranges), figsize=(4 * len(inc_ranges), 4), layout='constrained', sharey=True)
        axes = np.atleast_1d(axes)

        # ---------reading in input file with fitted parameters for 1D fit (with all inc)-----------
        input_filename_1dfit = args.data_dir / f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_indiv_allinc.csv'
        df_1d_allinc = pd.read_csv(input_filename_1dfit)
        df_1d_allinc = df_1d_allinc.rename(columns={'er0':'e_r0', 'eD0':'e_D0'})

        # ---------looping over inc ranges----------
        for index, inc_range in enumerate(inc_ranges):
            outdir = Path(f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{inc_range[0]:.1f}_{inc_range[1]:.1f}')
            input_filename_2dfit = args.data_dir / outdir / 'projection_fit.csv'
            print(f'For inc_range {inc_range} ({index + 1}/{len(inc_ranges)}): Reading projection fit parameters from {input_filename_2dfit}..')
            df_2d = pd.read_csv(input_filename_2dfit)

            # ----------merging with 1D fit param dataframe------------------
            df_1d = df_1d_allinc[df_1d_allinc['inc_bin'].astype(str) == f'{inc_range[0]:.1f}_{inc_range[1]:.1f}'].drop(columns='inc_bin')
            df_1d['r0'] = 2 * df_1d['r0'] # factor of 2 before comparing with 2D r0
            df_1d['e_r0'] = 2 * df_1d['e_r0']
            merge_on = ['halo', 'snap']
            common_cols = df_1d.columns.intersection(df_2d.columns)
            cols_to_drop = [col for col in common_cols if col not in merge_on]
            df = pd.merge(df_1d, df_2d.drop(columns=cols_to_drop), on=merge_on, how='inner')

            # ----computing r sqaured---------
            quant = ((unp.uarray(df['distmin0'], df['e_distmin0']) ** 2 + unp.uarray(df['distmaj0'], df['e_distmaj0']) ** 2) ** 0.5)
            df['r_sq'] = unp.nominal_values(quant)
            df['e_r_sq'] = unp.std_devs(quant)

            # -------computing sSFR---
            df['log_ssfr'] = df['log_sfr'] - df['log_star_mass']

            # ---------doing the plot------
            axes[index].errorbar(df[col1], df[col2], xerr=df[f'e_{col1}'], yerr=df[f'e_{col2}'], fmt='none', color=fixed_color, lw=0.5, alpha=0.8, capsize=2, zorder=-10)
            im = axes[index].scatter(df[col1], df[col2], c=fixed_color if colorby_col is None else df[colorby_col], s=10, lw=0.5, ec='k', alpha=0.8)

            # -------axes limits-----------
            axes[index].set_xscale('log')
            axes[index].set_yscale('log')

            axes[index].plot([axes[index].get_xlim()[0], axes[index].get_xlim()[1]], [axes[index].get_xlim()[0], axes[index].get_xlim()[1]], c='k', ls='dashed', lw=1)
            
            mad = median_abs_deviation(df[col2] - df[col1])
            axes[index].text(0.05, 0.8, f'MAD={mad:.2f}', color='k', ha='left', va='top', transform=axes[index].transAxes, fontsize=args.fontsize)
        
            axes[index].set_xlim(limits[0], limits[1])
            axes[index].set_ylim(limits[0], limits[1])
            
            axes[index] = annotate_axes(axes[index], label_dict[col1], label_dict[col2], args=args, xloc=0.05, label=rf'{inc_range[0]}$^\circ$ $< i <$ {inc_range[1]}$^\circ$', hide_xaxis=False, hide_yaxis=index, bbox=False, set_ticks=False)

        # --------color axis-----------
        if colorby_col is not None:
            cbar = fig.colorbar(
                im, 
                ax=axes,          # Pass the entire array/list of axes here
                location='top',   # Forces it above the subplots
                orientation='horizontal', # Ensures the colorbar orientation is horizontal
                shrink=1.,       # Optional: scales width (1.0 = 100% width of the axes grid)
                pad=0.02,          # Optional: spacing between colorbar and subplots top edge
                aspect = 70,       # higher value for thinner colorbar
            )
            cbar.set_label(label_dict[colorby_col], labelpad=10, fontsize=args.fontsize)

        # ------------saving the figure----------------------
        figname = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_{len(inc_ranges)}_inc_ranges_2D_vs_1D_fit_{figlabel}.png'
        if colorby_col is not None:
            figname = figname.replace('.png', f'_colby_{colorby_col}.png')

        save_fig(fig, args.plot_dir, figname, args=args)
             
    return

# -----------------------------------------------------------------------------
def execute_mode_plot_2d_fit_radius_ratio(args):
    '''
    Function to execute mode plot_2d_fit_radius_ratio, which reads in the 2D fitted parameter file and plots rx0 / ry0 vs global properties
    Saves the plots
    Returns nothing
    '''
    # -----------setting up figure-----------------
    inc_ranges = [[0, 30], [60, 90], [80, 90]]
    #inc_ranges = [[60, 90]]
    
    #col1, col2, figlabel = 'distmin0', 'distmaj0', 'rmin_rmaj_ratio' # this is for the scenario where the minor and major r0 have been fitted together as a 2D profile
    col1, col2, figlabel = 'distmin0_indep', 'distmaj0_indep', 'rmin_rmaj_ratio_indep' # this is for the scenario where the minor and major r0 have been fitted independently (individually)    

    # --------------dictionaries for labels and colors-------------
    fixed_color = 'cornflowerblue'
    label_dict = {'distmaj0': r'r$_{0,maj}$ (kpc)',
                    'distmin0': r'r$_{0,min}$ (kpc)',
                    'distmaj0_indep': r'r$_{0,maj}$ (kpc)',
                    'distmin0_indep': r'r$_{0,min}$ (kpc)',
                    'losdm0': r'D$_{0,2D}$ pc cm$^{-3}$',
                    'r_sq': r'$\sqrt{r_{0,maj}^2 + r_{0,min}^2}$ (kpc)',
                    'log_star_mass': r'$\log{(M_*/M_\odot)}$',
                    'log_sfr': r'$\log{SFR/M_\odot yr^{-1}}$',
                    'log_ssfr': r'$\log{sSFR/yr^{-1}}$',
                    'redshift': 'Redshift',
                    'log_r_ratio': r'$\log$ (r$_{0,maj}$ / r$_{0,min}$)'
                    }

    colorby_col_arr = ['log_sfr'] #[None, 'log_sfr', 'log_ssfr', 'redshift']
    xcol = 'log_star_mass'

    # ---------looping over color cols---------------
    for colorby_col in colorby_col_arr:
        print(f'\n\nDoing color by {colorby_col}..')

        # -----------setting up figure-----------------
        fig, axes = plt.subplots(1, len(inc_ranges), figsize=(3.6 * len(inc_ranges), 4), layout='constrained', sharey=True)
        axes = np.atleast_1d(axes)

        # ---------looping over inc ranges----------
        for index, inc_range in enumerate(inc_ranges):
            outdir = Path(f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{inc_range[0]:.1f}_{inc_range[1]:.1f}')
            input_filename = args.data_dir / outdir / 'projection_fit.csv'
            print(f'For inc_range {inc_range} ({index + 1}/{len(inc_ranges)}): Reading projection fit parameters from {input_filename}..')
            df = pd.read_csv(input_filename)

            # -------computing sSFR---
            df['log_ssfr'] = df['log_sfr'] - df['log_star_mass']

            # ----computing r_min / r_maj ratio--------
            quant = unp.log10(unp.uarray(df[col1], df[f'e_{col1}']) / unp.uarray(df[col2], df[f'e_{col2}']))
            df['log_r_ratio'] = unp.nominal_values(quant)
            df['e_log_r_ratio'] = unp.std_devs(quant)

            # ---------doing the plot------
            axes[index].errorbar(df[xcol], df['log_r_ratio'], yerr=df['e_log_r_ratio'], fmt='none', color=fixed_color, lw=0.5, alpha=0.8, capsize=2, zorder=-10)
            im = axes[index].scatter(df[xcol], df['log_r_ratio'], c=fixed_color if colorby_col is None else df[colorby_col], s=10, lw=0.5, ec='k', alpha=0.8)

            # -------axes limits-----------
            axes[index].axhline(0, c='k', ls='--')
            axes[index].set_ylim(-1.1, 1.1)
            
            axes[index] = annotate_axes(axes[index], label_dict[xcol], label_dict['log_r_ratio'], args=args, xloc=0.05, label=rf'{inc_range[0]}$^\circ$ $< i <$ {inc_range[1]}$^\circ$', hide_xaxis=False, hide_yaxis=index, bbox=False, set_ticks=False)

        # --------color axis-----------
        if colorby_col is not None:
            cbar = fig.colorbar(
                im, 
                ax=axes,          # Pass the entire array/list of axes here
                location='top',   # Forces it above the subplots
                orientation='horizontal', # Ensures the colorbar orientation is horizontal
                shrink=1.,       # Optional: scales width (1.0 = 100% width of the axes grid)
                pad=0.02,          # Optional: spacing between colorbar and subplots top edge
                aspect = 70,       # higher value for thinner colorbar
            )
            cbar.set_label(label_dict[colorby_col], labelpad=10, fontsize=args.fontsize)

        # ------------saving the figure----------------------
        figname = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_{len(inc_ranges)}_inc_ranges_2D_fit_{figlabel}.png'
        if colorby_col is not None:
            figname = figname.replace('.png', f'_colby_{colorby_col}.png')

        save_fig(fig, args.plot_dir, figname, args=args)
             
    return

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
    label_dict = {'log_star_mass': r'$\log{(M_*/M_\odot)}$',
                    'log_sfr': r'$\log{SFR/M_\odot yr^{-1}}$',
                    'log_ssfr': r'$\log{sSFR/yr^{-1}}$',
                    'redshift': 'Redshift',
                    }

    # ---------setup figure---------------
    fig, ax = plt.subplots(1, figsize=(7, 5))
    fig.subplots_adjust(left=0.12, bottom=0.12, right=0.93, top=0.98)

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

    colorlabel = label_dict[colorcol]
    cbar = fig.colorbar(sm, ax=ax, pad=0.01)
    cbar.set_label(colorlabel, fontsize=args.fontsize)
    cbar.ax.tick_params(labelsize=args.fontsize)

    ax = annotate_axes(ax, "Impact factor (kpc)", "DM (pc cm$^{-3}$)", args=args, clabel=colorlabel, set_ticks=False)
    ax.text(0.95, 0.95, rf'{args.z_range[0]} $\leq z <$ {args.z_range[1]}', c='k', fontsize=args.fontsize, ha='right', va='top', transform=ax.transAxes)
    #ax.text(0.95, 0.95, rf'{args.lsm_range[0]} < $\log$(M/M$_\odot$) < {args.lsm_range[1]}', c='k', fontsize=args.fontsize, ha='right', va='top', transform=ax.transAxes)
    #ax.text(0.95, 0.85, rf'{args.lsfr_range[0]} < $\log$(SFR/M$_\odot$ yr$^{-1}$) < {args.lsfr_range[1]} [{len(df_snap)}]', c='k', fontsize=args.fontsize, ha='right', va='top', transform=ax.transAxes)

    save_fig(fig, args.fig_dir, f'DM_vs_impfact_indi_lsm_bin_{args.lsm_range[0]}_{args.lsm_range[1]}_lsfr_bin_{args.lsfr_range[0]}_{args.lsfr_range[1]}_zrange_{args.z_range[0]}_{args.z_range[1]}_inc{args.inc_range[0]}-{args.inc_range[1]}_colorby_{colorcol}.pdf', args)
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
                    
                elif (args.mode=='proj' or args.mode=='projection') :
                    df_fit = execute_mode_projection(df_snap, args)     

                elif (args.mode=='plot_2d_param_comp') : # to compare fitted parameters (r0, D0, etc) from 2D fitting
                    execute_mode_plot_2d_fit_param_comparison(args)     

                elif (args.mode=='plot_1d_param_comp') : # to compare fitted parameters (r0, D0, etc) from 1D fitting
                    execute_mode_plot_1d_fit_param_comparison(args)     

                elif (args.mode=='plot_1d_2d_comp') : # to compare fitted parameters (r0, D0, etc) across 1D and 2D fitting
                    execute_mode_plot_1d_2d_comparison(args)

                elif (args.mode=='plot_2d_radius_ratio') : # to compare r_min/r_max from 2D fitting
                    execute_mode_plot_2d_fit_radius_ratio(args)

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

                if args.mode == 'plot_indi':
                    #ax = plot_dm_impfac_indi_combined(df_snap, args, colorcol='redshift')
                    ax = plot_dm_impfac_indi_combined(df_snap, args, colorcol='log_star_mass')
                    ax = plot_dm_impfac_indi_combined(df_snap, args, colorcol='log_sfr')
                    ax = plot_dm_impfac_indi_combined(df_snap, args, colorcol='log_ssfr')


        if args.mode == 'lsmzsfr' and args.multi_panel:
            save_fig(fig, args.fig_dir, f'{Path(args.resfile_prefix).stem}_z_{args.z_range[0]}_{args.z_range[1]}_{args.mode}_inc_{args.inc_range[0]}_{args.inc_range[1]}_multipanel_1d.pdf', args)
        
    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))











































