#!/usr/bin/env python3
"""
    Title :      plots_for_frb_paper
    Notes :      Make various plots to be used for the FRB paper
    Output :     Plots as PDF
    Author :     Ayan Acharyya
    Started :    31-03-26
    Examples :   run plots_for_frb_paper.py --plot_dm_lsm --lsm 9.5,9.75
                 run plots_for_frb_paper.py --plot_dm_lsm --lsm 10.75,11,11.25,11.5 --inc 0,30,80,90
                 run plots_for_frb_paper.py --plot_dm_lsm --lsm 10.75,11,11.25,11.5 --inc 0,30,80,90 --multi_panel
                 run plots_for_frb_paper.py --plot_radprof
                 run plots_for_frb_paper.py --plot_dm_fit --fit_robust
                 run plots_for_frb_paper.py --plot_dm_fit --regres linear
                 run plots_for_frb_paper.py --plot_dm_fit --mode indi
                 run plots_for_frb_paper.py --plot_dm_fit --fit_robust --resfile_prefix binby_lsm_lsfr
                 run plots_for_frb_paper.py --plot_dm_all_lsm --cmap tab10 --set_ylin
                 run plots_for_frb_paper.py --plot_dm_all_lsm
                 run plots_for_frb_paper.py --make_latex_table --resfile_prefix binby_lsm_lsfr
                 run plots_for_frb_paper.py --make_latex_table --resfile_prefix all_lsm
                 run plots_for_frb_paper.py --plot_foggie_snaps --system ayan_local --halo 2878 --upto_kpc 100 --reskpc 0.5

"""
from craft_header import *
from craft_utils import *
setup_plot_style()
import plotfns as pfns

start_time = datetime.now()

# ------------------------------------------------------------------------------------------------
def read_dataframe_txt(filename, interval_cols=['lsm_bin', 'lsfr_bin', 'inc_bin']):
    '''
    Function to read txt file as pandas dataframe and properly parse intervals
    Returns dataframe
    '''
    df = pd.read_csv(filename, sep='\t')
    if len(interval_cols) > 1: df = df.drop_duplicates(subset=interval_cols, keep='last')

    for col in interval_cols:
        temp_df = df[col].str.strip('()[]').str.split(',', expand=True).astype(float)
        df[col] = temp_df.apply(lambda x: pd.Interval(x[0], x[1], closed='right'), axis=1)
    
    return df

# ------------------------------------------------------------------------------------------------
def read_dataframe_csv(filename, interval_cols=['lsm_bin', 'lsfr_bin', 'inc_bin']):
    '''
    Function to read csv file as pandas dataframe and properly parse intervals
    Returns dataframe
    '''
    df = pd.read_csv(filename, comment='#')

    for col in interval_cols:
        temp_df = df[col].str.split('_', expand=True).astype(float)
        df[col] = temp_df.apply(lambda x: pd.Interval(x[0], x[1], closed='right'), axis=1)

    df = df.rename(columns={'e_r0':'er0', 
                            'e_D0':'eD0',
                            'log_star_mass': 'medlsm',
                            'sfr': 'medsfr',
                            'sfr_100Myr': 'medsfr_100Myr',
                            'log_gas_mass': 'medlgsm',
                            })    
    return df

# ------------------------------------------------------------------------------------------------
def plot_dm_impfac_one_lsm_bin(df_dmpars, args, given_ax=None):
    '''
    Plot DM vs Impact factor in a single panel, for a given stellar mass and inclination range
    Saves plot
    Returns axis handle
    '''
    # ---------setup figure---------------
    if given_ax is None:
        fig, ax = plt.subplots(1, figsize=(5, 4), layout='constrained')
    else:
        ax = given_ax
    
    face_col_arr    = ['b', 'lightblue', 'r']
    mark_arr        = ['o', 'x', 's']
    fill_arr        = ['full', 'full', 'none']

    # --------loop over inclination bins--------------
    for index, this_inc_bin in enumerate(args.inc_bins):
        print(f'\t\nRunning ({index + 1}/{len(args.inc_bins)}) for inclination bin {this_inc_bin}..\n')
        args.inc_range = this_inc_bin

        # -----------get required bin--------------------
        df_dmpars_sub = df_dmpars[(df_dmpars['inc_bin'] == pd.Interval(args.inc_range[0], args.inc_range[1])) &
                                  (df_dmpars['lsm_bin'] == pd.Interval(args.lsm_range[0], args.lsm_range[1])) & 
                                  (df_dmpars['lsfr_bin'] == pd.Interval(args.lsfr_range[0], args.lsfr_range[1]))]
        if len(df_dmpars_sub) == 0:
            print(f'\tNo entries found for inclination range {args.inc_range}, continuing..')
            continue
        
        dmpars = df_dmpars_sub.iloc[0]
        infile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{args.inc_range[0]}_{args.inc_range[1]}/lsmzsfr_lsm_{dmpars["lsm_bin"].left}_{dmpars["lsm_bin"].right}_lsfr_{dmpars["lsfr_bin"].left}_{dmpars["lsfr_bin"].right}_1d.npy'

        data_arr = np.load(infile)

        ax.errorbar(data_arr[0], data_arr[1], yerr=[data_arr[2], data_arr[3]], fmt=mark_arr[index], fillstyle=fill_arr[index], mfc=face_col_arr[index], mec=face_col_arr[index], ecolor=face_col_arr[index], lw=1, markersize=6, capsize=4)
        
        dm_arr = schechter(data_arr[0], dmpars['r0'], dmpars['D0'])
        ax.plot(data_arr[0], dm_arr, color=face_col_arr[index], lw=1, ls='dashed')

        ax.text(x=impbinegs[-1], y=320 / (index+1), s=f"$D_0$ = {dmpars['D0']:.1f}", c=face_col_arr[index], fontsize=args.fontsize / args.fontfactor, ha='right', va='top')
        ax.text(x=impbinegs[-1], y=200 / (index+1), s=f"$r_0$ = {dmpars['r0']:.1f}", c=face_col_arr[index], fontsize=args.fontsize / args.fontfactor, ha='right', va='top')

        if len(args.inc_bins) > 1:
            ax.text(x=0.4 * impbinegs[1], y=0.8 + index * 0.3, s=f'{args.inc_range[0]}' + r' < $i$ < ' + f'{args.inc_range[1]}', fontsize=args.fontsize / args.fontfactor, color=face_col_arr[index])

    # ------------plot based on expected DM from scaling relation-----------
    # if len(args.inc_bins) < 2:
    #     dm_expected_arr = 10 ** logradialexp3(data_arr[0], 10.0 ** (0.61 -0.53 * (dmpars['medlsm'] - 10)), 10.0 ** (2.15 + 0.24 * (dmpars['medlsm'] - 10)))
    #     ax.plot(data_arr[0], dm_expected_arr, color='r', lw=1, ls='dotted')

    # -------annotating plot----------------
    ax.set_xscale("log")
    ax.set_xticks(impbinegs[1:],impbinegs[1:])
    ax.set_xlim([0.25 * impbinegs[1], 1.5 * impbinegs[-1]])

    if not args.set_ylin:
        ax.set_yscale("log")
        ax.set_yticks(dm_ticks, dm_ticks)
        ax.set_ylim([0.5, 3 * maxdmcol])
    
    ax = annotate_axes(ax, "Impact factor (kpc)", "DM (pc cm$^{-3}$)", args=args, set_ticks=False)

    ax.text(x=0.4*impbinegs[1], y=300, s=rf'{args.z_range[0]} $< z <$ {args.z_range[1]}', fontsize=args.fontsize / args.fontfactor)
    #ax.text(x=0.4*impbinegs[1], y=300, s="%.2f < log ($M_* / M_{\odot}$) < %.2f"%(dmpars['lsm_bin'].left, dmpars['lsm_bin'].right), fontsize=args.fontsize / args.fontfactor)
    #ax.text(x=0.4*impbinegs[1], y=1.6, s="log ($M_* / M_{\odot}$) = %.2f"% dmpars['medlsm'], fontsize=args.fontsize / args.fontfactor)
    #ax.text(x=0.4*impbinegs[1], y=0.8, s="SFR = %.2f $M_{\odot} yr^{-1}$"% dmpars['medsfr'], fontsize=args.fontsize / args.fontfactor)	
    #ax.text(x=1.0*impbinegs[-4], y=150, s="$D_0$ = %d $\pm$ %d"%(dmpars['D0'], dmpars['eD0']), fontsize=args.fontsize / args.fontfactor)
    #ax.text(x=1.0*impbinegs[-4], y=75, s="$r_0$ = %.1f $\pm$ %.1f"%(dmpars['r0'], dmpars['er0']), fontsize=args.fontsize / args.fontfactor)

    if given_ax is None:
        save_fig(fig, args.fig_dir, f'DM_vs_impfact_inc{",".join(np.array(args.inc_bins).flatten().astype(str))}_lsm_{args.lsm_range[0]}_{args.lsm_range[1]}_lsfr_{args.lsfr_range[0]}_{args.lsfr_range[1]}_zrange_{args.z_range[0]}_{args.z_range[1]}.pdf', args)
        plt.show(block=False)

    return ax

# ------------------------------------------------------------------------------------------------
def plot_dm_impfac_all_lsm_bin(df_dmpars, args, cmap='tab10'):
    '''
    Plot DM vs Impact factor in a single panel, for a list of stellar mass ranges and a given inclination range
    Saves plot
    Returns axis handle
    '''
    # ---------setup figure---------------
    fig, ax = plt.subplots(1, figsize=(8, 5), layout='constrained')
    color_list = plt.get_cmap(cmap)(range(12))

    # -----------get required bin--------------------
    df_dmpars = df_dmpars[(df_dmpars['lsfr_bin'] == pd.Interval(args.lsfr_range[0], args.lsfr_range[1]))].reset_index(drop=True)

    # -----------loop through mass bins--------------------
    for index, dmpars in df_dmpars.iterrows(): 
        infile = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_inc_{args.inc_range[0]}_{args.inc_range[1]}/lsmzsfr_lsm_{dmpars["lsm_bin"].left}_{dmpars["lsm_bin"].right}_lsfr_{dmpars["lsfr_bin"].left}_{dmpars["lsfr_bin"].right}_1d.npy'

        data_arr = np.load(infile)
        col = color_list[index]

        #ax.plot(data_arr[0], data_arr[1], color=col, lw=1)
        dm_arr = schechter(data_arr[0], dmpars['r0'], dmpars['D0'])

        ax.plot(data_arr[0], dm_arr, color=col, lw=2, ls=lslist[index % len(lslist)], label=f'{dmpars["lsm_bin"].left:.1f}-{dmpars["lsm_bin"].right:.1f}')
        ax.fill_between(data_arr[0], data_arr[1] - data_arr[2], data_arr[1] + data_arr[3], color=col,alpha=0.1)
		
        if args.set_ylin: ax.text(0.98, 0.98 - index * 0.04, f'{dmpars["lsm_bin"].left:.2f}' + r' $\leq \log$ M/M$_\odot \leq$ ' + f'{dmpars["lsm_bin"].right:.2f}', c=col, ha='right', va='top', transform=ax.transAxes, fontsize=args.fontsize / args.fontfactor)
        else: ax.text(0.02, 0.02 + index * 0.04, f'{dmpars["lsm_bin"].left:.2f}' + r' $\leq \log$ M/M$_\odot \leq$ ' + f'{dmpars["lsm_bin"].right:.2f}', c=col, ha='left', va='bottom', transform=ax.transAxes, fontsize=args.fontsize / args.fontfactor)

    #ax.legend(ncol=2, loc='upper right', fontsize=args.fontsize / args.fontfactor)
    ax.set_xscale("log")
    ax.set_xticks(impbinegs[1:],impbinegs[1:])

    if not args.set_ylin:
        ax.set_yscale("log")
        ax.set_yticks(dm_ticks, dm_ticks)
    
    ax = annotate_axes(ax, "Impact factor (kpc)", "DM (pc cm$^{-3}$)", args=args, set_ticks=False)
        
    save_fig(fig, args.fig_dir, f'DM_vs_impfact_all_lsm_inc{args.inc_range[0]}-{args.inc_range[1]}_zrange_{args.z_range[0]}_{args.z_range[1]}.pdf', args)
    plt.show(block=False)

    return fig

# ------------------------------------------------------------------------------------------------
def plot_dm_fit(df_dmpars, args):
    '''
    Plot DM0 and r0 vs stellar mass in a single panel
    Saves plot
    Returns axis handle
    '''
    outfilename = f'{args.fig_dir}/{Path(args.resfile_prefix).stem}_z_{args.z_range[0]}_{args.z_range[1]}_DM0_r0_vs_lsm_inc_{args.inc_range[0]}_{args.inc_range[1]}'
    df_dmpars['log_ssfr'] = np.log10((10 ** df_dmpars['medlsm']) / df_dmpars['medsfr']) + 20
    df_dmpars['medlgsm_offset'] = df_dmpars['medlgsm'] - 10
    df_dmpars['medlsm_offset'] = df_dmpars['medlsm'] - 10
    df_dmpars['log_medsfr'] = np.log10(df_dmpars['medsfr'])
    df_dmpars['log_medsfr_100Myr'] = np.log10(df_dmpars['medsfr_100Myr'])
    
    res = pfns.plt_dmpars(df_dmpars, outfilename, 3.0, xcol='medlsm_offset', y1col='D0', y2col='r0', x2col='log_medsfr', fit_robust=args.fit_robust, regres=args.regres, fortalk=args.fortalk,scale_fit_thresh=10)
    res = pfns.plt_dmpars(df_dmpars, outfilename, 3.0, xcol='log_medsfr_100Myr', y1col='D0', y2col='r0', x2col='medlsm_offset', fit_robust=args.fit_robust, regres=args.regres, fortalk=args.fortalk,scale_fit_thresh=2)
    res = pfns.plt_dmpars(df_dmpars, outfilename, 3.0, xcol='log_medsfr', y1col='D0', y2col='r0', x2col='medlsm_offset', fit_robust=args.fit_robust, regres=args.regres, fortalk=args.fortalk,scale_fit_thresh=2)
    #res = pfns.plt_dmpars(df_dmpars, outfilename, 3.0, xcol='log_medsfr', y1col='D0', y2col='r0', x2col='medlgsm_offset', fit_robust=args.fit_robust, fortalk=args.fortalk)
    #res = pfns.plt_dmpars(df_dmpars, outfilename, 3.0, xcol='log_ssfr', y1col='D0', y2col='r0', x2col='log_medsfr', fit_robust=args.fit_robust, fortalk=args.fortalk)

    return res

# ------------------------------------------------------------------------------------------------
def make_latex_table(df_dmpars, args, columns=['lsm_bin', 'ngal', 'medlsm', 'medsfr', 'D0', 'r0']):
    '''
    Convert the input dataframe into a latex table
    Saves latex table
    Returns latex dataframe
    '''    
    colnames_dict = {'lsm_bin':r'\makecell{$\log(M_*/M_\odot$)\\range}', 
                     'lsfr_bin':r'\makecell{$\log SFR (M_\odot/yr$)\\range}', 
                     'ngal':r'N$_{\rm snapshot}$', 
                     'medlsm':r'\makecell{Median\\$\log(M_*/M_\odot$)}', 
                     'medsfr':r'\makecell{Median SFR\\($M_\odot\: yr^{-1}$)}', 
                     'r0':r'\makecell{$r_0$\\(kpc)}', 
                     'D0':r'\makecell{$D_0$\\($pc\: cm^{-3}$)}',
                     }

    columns_to_publish = columns
    if 'lsfr_bin' in df_dmpars and len(pd.unique(df_dmpars['lsfr_bin'])) > 1:
        try: columns_to_publish.insert(columns_to_publish.index('lsm_bin') +1, 'lsfr_bin') # if 'lsm_bin' column exists, then insert 'lsfr_bin' column immediately next to it
        except: columns_to_publish += ['lsfr_bin'] # otherwise append 'lsfr_bin' column at the end of the table

    columns_with_err = ['r0', 'D0']
    columns_with_interval = [item for item in columns_to_publish if item.endswith('_bin')]

    df_latex = df_dmpars[np.hstack([columns_to_publish, ['e' + item for item in columns_with_err]])]
    df_mread = df_latex.copy()

    for col in columns_with_interval:
        df_latex[col] = df_latex[col].apply(lambda x: f"{x.left:.2f} -- {x.right:.2f}")
    
    for col in columns_with_err:
        if 'D0' in col: df_latex[col] = df_latex.apply(lambda x: f"{x[col]:.0f} $\pm$ {x['e' + col]:.0f}", axis=1) # 0 floating point precision for D0
        else: df_latex[col] = df_latex.apply(lambda x: f"{x[col]:.1f} $\pm$ {x['e' + col]:.1f}", axis=1)
        df_latex.drop(columns=['e' + col], inplace=True)
    
    for col in (set(df_latex.columns) - set(np.hstack([columns_with_err, columns_with_interval, ['ngal']]))):
        df_latex[col] = df_latex[col].map('{:.2f}'.format)

    df_latex = df_latex.rename(columns=colnames_dict)

    outfilename = f'{args.fig_dir}/{Path(args.resfile_prefix).stem}_z_{args.z_range[0]}_{args.z_range[1]}_table_DM0_r0_vs_lsm_inc_{args.inc_range[0]}_{args.inc_range[1]}.tex'

    df_mread.to_csv(outfilename.replace('.tex', '.txt'), index=None, sep='\t')
    df_latex.to_latex(outfilename, index=False, escape=False, column_format='l' * 1 + 'c' * (len(df_latex.columns) - 1))
    
    # -----------to insert lines between SFR groups------------
    insert_line_in_file('\\toprule\n', 1, outfilename) # to insert an additioal \toprule
    if 'lsfr_bin' in columns_to_publish:
        pos = 4
        for lsfr_bin in pd.unique(df_dmpars['lsfr_bin'])[:-1]:
            df_sub = df_dmpars[df_dmpars['lsfr_bin']==lsfr_bin]
            pos = pos + len(df_sub) + 1
            insert_line_in_file('\\midrule\n', pos, outfilename) # to insert an additioal \midrule

    print(f'Saved latex table as {outfilename} and as .txt')
    print(df_latex)

    return df_latex

# ------------------------------------------------------------------------------------------------
def plot_multipanel_foggie(args):
    '''
    Make a 3-row multi-panel plot for a given FOGGIE halo with gas projection, electron density projection, and electron density radial profile, for a series of redshifts
    Saves the plot
    Returns figure handle
    '''
    # --------setup plot parameters----------
    redshift_arr = [0.05, 0.10, 0.15, 0.20, 0.25]
    inc_ranges	=	np.array([[0,10], [80,90]])
    colist	= ['cornflowerblue','salmon', 'k']
    shlist	= ['aqua', 'coral','grey']

    quant_dict = {'density':['density', 'Gas density', 'Msun/pc**3', -1.5, 3.5, 'cornflowerblue', 'cividis', True, 'Msun/pc**2', r'Projected gas density / M$_\odot$ pc$^{-2}$'], 
                'el_density':['El_number_density', 'Electron density', 'cm**-3', 0, 220, 'cornflowerblue', 'viridis', False, 'pc*cm**-3', r'Projected electron density / pc cm$^{-3}$']
                } # for each quantity: [yt field, label in plots, units, lower limit in log, upper limit in log, color for scatter plot, colormap, whether to take log, units for projection plot, units to display in projection plot]
    quant_arr = ['density', 'el_density']

    # ------plotting onto a matplotlib figure--------------
    ncols = len(redshift_arr)
    fig = plt.figure(figsize=(1.6 * ncols, 7.))

    # ------- 1. Outer & Sub GridSpec Setup -------
    outer_gs = fig.add_gridspec(
        nrows=3, ncols=1, 
        height_ratios=[1.1, 1.1, 1.0], 
        hspace=0.2, 
        left=0.07, right=0.98, top=0.92, bottom=0.07
    )

    # Top GridSpec has 4 rows: [Cbar 0, Plot Row 0, Cbar 1, Plot Row 1]
    top_gs1 = outer_gs[0].subgridspec(
        nrows=2, ncols=ncols, 
        height_ratios=[0.05, 1.0], 
        wspace=0.0, hspace=0.0
    )

    top_gs2 = outer_gs[1].subgridspec(
        nrows=2, ncols=ncols, 
        height_ratios=[0.05, 1.0], 
        wspace=0.0, hspace=0.0
    )

    bot_gs = outer_gs[2].subgridspec(1, ncols, wspace=0.0)

    # Axes mapping to preserve axes[row, col] indexing
    axes = np.empty((3, ncols), dtype=object)
    for col in range(ncols):
        axes[0, col] = fig.add_subplot(top_gs1[1, col]) # Plot Row 0
        axes[1, col] = fig.add_subplot(top_gs2[1, col]) # Plot Row 1
        axes[2, col] = fig.add_subplot(bot_gs[0, col]) # Plot Row 2 (Radial Profiles)

    # Colorbar axes spanning all columns
    cax0 = fig.add_subplot(top_gs1[0, :])
    cax1 = fig.add_subplot(top_gs2[0, :])

    # -----------determine snapshot list from redshift list---------
    args.code_dir = '/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/'
    df = pd.read_csv(args.code_dir + f'halo_infos/00{args.halo}/nref11c_nref9f/halo_cen_smoothed', sep=r'\s*\|\s*', engine='python')
    df = df.dropna(axis=1, how='all')[['snap', 'redshift']]
    output_arr = []
    for redshift in redshift_arr:
        idx = (df['redshift'] - redshift).abs().idxmin()
        output_arr.append(df.loc[idx, 'snap'])

    # -------looping over redshift snapshots--------
    for index, args.output in enumerate(output_arr):
        print(f'\nDoing output {args.output} which is {index + 1} of {len(output_arr)}..')

        # -----------read in FRB data-------------
        fitsname = Path(args.fits_dir) / f'{args.output}_{args.halo}_FRB_{quant_dict["el_density"][0]}{args.upto_text}{args.res_text}.fits'
        print(f'Trying to read {fitsname}..')
        hdul = fits.open(fitsname)
        sfr = hdul[0].header['SFR']
        log_mstar = hdul[0].header['LOG_MSTAR']
        if sfr == 'NaN': sfr = np.nan

        # ----------plot gas and electron density projections-------
        for index2, quant in enumerate(quant_arr):
            ax = axes[index2, index]
            clim = [quant_dict[quant][3], quant_dict[quant][4]]
            cmap = quant_dict[quant][6]

            # --------read in the quantity------
            data = hdul[f'{quant_dict[quant][1]} FACE ON PROJ'].data
            crpix = hdul[f'{quant_dict[quant][1]} FACE ON PROJ'].header['CRPIX1']
            cdelt = hdul[f'{quant_dict[quant][1]} FACE ON PROJ'].header['CDELT1']
            crval = hdul[f'{quant_dict[quant][1]} FACE ON PROJ'].header['CRVAL1']
            if quant_dict[quant][7]: data = np.log10(data)

            # ---------plot projection-----------
            p = ax.imshow(data, cmap=cmap, vmin=clim[0], vmax=clim[1])

            # ---------------prepping axes------------------------
            ax.xaxis.set_major_locator(plt.MaxNLocator(5))
            ax.yaxis.set_major_locator(plt.MaxNLocator(5))

            ax = annotate_axes(ax, 'Offset (kpc)', 'Offset (kpc)', args=args, xloc=0.05, 
                               #label=rf'$\log$(M/M$_\odot$) = {log_mstar:.2f}' if index2 else rf'$\log$ SFR = {sfr:.2f}',# rf'$\log$(SFR/M$_\odot$ yr$^{-1}$) = {sfr:.2f}', 
                               hide_xaxis=index2 == 0, hide_yaxis=index, bbox=True, set_ticks=False, 
                               p=p, hide_cbar=True)

            if index2 > 0:
                ax.set_xticklabels(['%.1F' % ((item - crpix) * cdelt + crval) for item in ax.get_xticks()], fontsize=args.fontsize)

            if index == 0:
                ax.set_yticklabels(['%.1F' % ((item - crpix) * cdelt + crval) for item in ax.get_yticks()], fontsize=args.fontsize)

            # ---------------making colorbar------------------------
            if index == len(output_arr) - 1:
                cax = [cax0, cax1][index2]
                cbar = fig.colorbar(p, cax=cax, orientation='horizontal')
                cax.xaxis.set_ticks_position('top')
                cax.xaxis.set_label_position('top')
                cbar.ax.tick_params(labelsize=args.fontsize, width=1.2, length=3)
                clabel = fr'Log ({quant_dict[quant][9]})' if quant_dict[quant][7] else quant_dict[quant][9]
                cbar.set_label(clabel, fontsize=args.fontsize)

        # ----------read in radial electron density profile-------
        fitsname = fitsname.stem
        if 'El' in fitsname: profsubdir = 'electron_density/'
        else: profsubdir = 'gas_density/'
        profdir = radialdir + profsubdir
        Path(profdir).mkdir(exist_ok=True, parents=True)
        profile_pkl_filename = profdir + fitsname + '_radprof.pkl'

        print("\nPlotting radial electron density profiles...\n")            
        with open(profile_pkl_filename, 'rb') as file_obj:
            cube = pkl.load(file_obj) # load the pickle file

        # ----------plot gas and electron density projections-------
        ax = axes[2, index]
        
        for i in range(0,len(inc_ranges)):
            inclim	= inc_ranges[i]
            relinds	= np.where((cube.inclination - np.deg2rad(inclim[0]))*(cube.inclination - np.deg2rad(inclim[1])) <= 0.0)	
            relne	= cube.neincrad[relinds]
            binned_ne	= np.zeros((len(radbins)-1, 6), dtype=np.float32)

            for k in range (0,len(radbins)-1):
                rel2inds		= np.where((cube.radkpc - radbins[k])*(cube.radkpc - radbins[k+1]) <= 0.0)
                binned_ne[k,0]	= (radbins[k]+radbins[k+1])/2.0
                binned_ne[k,1:6]= np.percentile(relne[:,rel2inds], (16, 25, 50, 75, 84))

            ax.fill_between(binned_ne[:,0], binned_ne[:,1], binned_ne[:,5], color=shlist[i],alpha=0.2)
            ax.plot(binned_ne[:,0], binned_ne[:,3], c=colist[i], marker='s', markersize=6, label=str(inclim[0])+"$^{\circ}$ < $i$ < "+str(inclim[1])+"$^{\circ}$")

        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylim(1e-5, 1e0)
        ax = annotate_axes(ax, 'Radius (kpc)', '$n_e$ (cm$^{-3}$)', args=args, xloc=0.6, 
                               label=f'z ={redshift_arr[index]:.2f}', 
                               hide_xaxis=False, hide_yaxis=index, bbox=False, set_ticks=False)

    # ---------------saving fig------------------------
    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try: mplcyberpunk.make_lines_glow()
        except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    save_fig(fig, Path(args.fig_dir), f'{args.halo}_multi_z_snapshots.pdf', args)
    plt.show(block=False)
    return fig

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')

    if args.plot_foggie_snaps:
        # -----------determining directories----------------
        args.fig_dir = root_dir + 'plots/'
        Path(args.fig_dir).mkdir(parents=True, exist_ok=True)
        
        args.fits_dir = root_dir + 'data/'
        Path(args.fits_dir).mkdir(parents=True, exist_ok=True)

        args.res_text = f'_res{args.reskpc:.1f}kpc'
        args.upto_text = '_upto%.1Fckpchinv' % args.upto_kpc if args.docomoving and args.upto_kpc is not None else '_upto%.1Fkpc' % args.upto_kpc if args.upto_kpc is not None else f'_upto{args.re:.1f}re'

        fig = plot_multipanel_foggie(args)

    else:    
        # -------------determining directories------------------
        if args.mode == 'indi': catalog_name = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_indiv_allinc.csv'
        else: catalog_name = f'{args.resfile_prefix}_z_{args.z_range[0]}_{args.z_range[1]}_allinc.csv'
        #df_dmpars = read_dataframe_txt(catalog_name, interval_cols=['inc_bin'] if args.mode == 'indi' else ['lsm_bin', 'lsfr_bin', 'inc_bin'])
        df_dmpars = read_dataframe_csv(catalog_name, interval_cols=['inc_bin'] if args.mode == 'indi' else ['lsm_bin', 'lsfr_bin', 'inc_bin'])

        # -------------calling plotting functions------------------
        if args.plot_dm_lsm:

            # ------------setup multi-panel figure if needed-------
            if args.multi_panel:
                nrows, ncols = 1, len(args.lsm_bins)
                fig, axes = plt.subplots(nrows, ncols, figsize=(3.0*np.array([len(args.lsm_bins),1])))
                axes = np.atleast_2d(axes)
                fig.subplots_adjust(left=0.1, bottom=0.15, right=0.98, top=0.98, wspace=0.01, hspace=0.01)

            # --------loop over log stellar mass bins------------
            for index, this_lsm_bin in enumerate(args.lsm_bins):
                print(f'\nRunning ({index + 1}/{len(args.lsm_bins)}) for stellar mass bin {this_lsm_bin}..\n')
                args.lsm_range = this_lsm_bin
                ax = plot_dm_impfac_one_lsm_bin(df_dmpars, args, given_ax=axes[index // ncols][index % ncols] if args.multi_panel else None)
        
                if args.multi_panel:
                    if index // ncols < nrows - 1:
                        ax.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                        ax.set_xlabel('')
                    if index % ncols > 0:
                        ax.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                        ax.set_ylabel('')
            if args.multi_panel:
                save_fig(fig, args.fig_dir, f'DM_vs_impfact_all_inc_all_lsm_multipanel_zrange_{args.z_range[0]}_{args.z_range[1]}.pdf', args)

        if args.plot_dm_all_lsm:
            df_dmpars = df_dmpars[df_dmpars['inc_bin'] == pd.Interval(args.inc_range[0], args.inc_range[1])].reset_index(drop=True) # choosing the correct inclination bin from the dataframe
            ax = plot_dm_impfac_all_lsm_bin(df_dmpars, args, cmap=args.cmap)

        if args.plot_dm_fit:
            df_dmpars = df_dmpars[df_dmpars['inc_bin'] == pd.Interval(args.inc_range[0], args.inc_range[1])].reset_index(drop=True) # choosing the correct inclination bin from the dataframe
            ax = plot_dm_fit(df_dmpars, args)
        
        if args.make_latex_table:
            df_dmpars = df_dmpars[df_dmpars['inc_bin'] == pd.Interval(args.inc_range[0], args.inc_range[1])].reset_index(drop=True) # choosing the correct inclination bin from the dataframe
            ax = make_latex_table(df_dmpars, args)

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
