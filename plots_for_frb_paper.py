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
                 run plots_for_frb_paper.py --plot_dm_fit --fit_lsm_range 8.5,11.0 --fit_robust
                 run plots_for_frb_paper.py --plot_dm_fit --fit_lsm_range 8.5,11.0
                 run plots_for_frb_paper.py --plot_dm_all_lsm --cmap tab10 --set_ylin
                 run plots_for_frb_paper.py --plot_dm_all_lsm
                 run plots_for_frb_paper.py --make_latex_table --resfile_prefix binby_lsm_lsfr
                 run plots_for_frb_paper.py --make_latex_table --resfile_prefix all_lsm
"""
from craft_header import *
from craft_utils import *
setup_plot_style()
import plotfns as pfns

start_time = datetime.now()

# ------------------------------------------------------------------------------------------------
def read_dataframe(filename, interval_cols=['lsm_bin', 'lsfr_bin', 'inc_bin']):
    '''
    Function to read txt file as pandas dataframe and properly parse intervals
    Returns dataframe
    '''
    df = pd.read_csv(filename, sep='\t')
    df = df.drop_duplicates(subset=interval_cols, keep='last')

    for col in interval_cols:
        temp_df = df[col].str.strip('()[]').str.split(',', expand=True).astype(float)
        df[col] = temp_df.apply(lambda x: pd.Interval(x[0], x[1], closed='right'), axis=1)
    
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
        infile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}/lsmzsfr_lsm_{dmpars["lsm_bin"].left}_{dmpars["lsm_bin"].right}_lsfr_{dmpars["lsfr_bin"].left}_{dmpars["lsfr_bin"].right}_1d.npy'

        data_arr = np.load(infile)

        ax.errorbar(data_arr[0], data_arr[1], yerr=[data_arr[2], data_arr[3]], fmt=mark_arr[index], fillstyle=fill_arr[index], mfc=face_col_arr[index], mec=face_col_arr[index], ecolor=face_col_arr[index], lw=1, markersize=6, capsize=4)
        
        dm_arr = 10 ** logradialexp3(data_arr[0], dmpars['r0'], dmpars['D0'])
        ax.plot(data_arr[0], dm_arr, color=face_col_arr[index], lw=1, ls='dashed')
            
        if len(args.inc_bins) > 1:
            ax.text(x=0.4 * impbinegs[1], y=0.8 + index * 0.3, s=f'{args.inc_range[0]}' + r' < $i$ < ' + f'{args.inc_range[1]}', fontsize=args.fontsize / args.fontfactor, color=face_col_arr[index])

    # ------------plot based on expected DM from scaling relation-----------
    if len(args.inc_bins) < 2:
        dm_expected_arr = 10 ** logradialexp3(data_arr[0], 10.0 ** (0.61 -0.53 * (dmpars['medlsm'] - 10)), 10.0 ** (2.15 + 0.24 * (dmpars['medlsm'] - 10)))
        ax.plot(data_arr[0], dm_expected_arr, color='r', lw=1, ls='dotted')

    # -------annotating plot----------------
    ax.set_xscale("log")
    ax.set_xticks(impbinegs[1:],impbinegs[1:])
    ax.set_xlim([0.25 * impbinegs[1], 1.5 * impbinegs[-1]])

    if not args.set_ylin:
        ax.set_yscale("log")
        ax.set_yticks(dm_ticks, dm_ticks)
        ax.set_ylim([0.5, 3 * maxdmcol])
    
    ax = annotate_axes(ax, "Impact factor (kpc)", "DM (pc cm$^{-3}$)", args=args, set_ticks=False)

    ax.text(x=0.4*impbinegs[1], y=300, s="%.2f < log ($M_* / M_{\odot}$) < %.2f"%(dmpars['lsm_bin'].left, dmpars['lsm_bin'].right), fontsize=args.fontsize / args.fontfactor)
    #ax.text(x=0.4*impbinegs[1], y=1.6, s="log ($M_* / M_{\odot}$) = %.2f"% dmpars['medlsm'], fontsize=args.fontsize / args.fontfactor)
    #ax.text(x=0.4*impbinegs[1], y=0.8, s="SFR = %.2f $M_{\odot} yr^{-1}$"% dmpars['medsfr'], fontsize=args.fontsize / args.fontfactor)	
    #ax.text(x=1.0*impbinegs[-4], y=150, s="$D_0$ = %d $\pm$ %d"%(dmpars['D0'], dmpars['eD0']), fontsize=args.fontsize / args.fontfactor)
    #ax.text(x=1.0*impbinegs[-4], y=75, s="$r_0$ = %.1f $\pm$ %.1f"%(dmpars['r0'], dmpars['er0']), fontsize=args.fontsize / args.fontfactor)

    if given_ax is None:
        save_fig(fig, args.fig_dir, f'DM_vs_impfact_inc{",".join(np.array(args.inc_bins).flatten().astype(str))}_lsm_{args.lsm_range[0]}_{args.lsm_range[1]}_lsfr_{args.lsfr_range[0]}_{args.lsfr_range[1]}.pdf', args)
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
        infile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}/lsmzsfr_lsm_{dmpars["lsm_bin"].left}_{dmpars["lsm_bin"].right}_lsfr_{dmpars["lsfr_bin"].left}_{dmpars["lsfr_bin"].right}_1d.npy'

        data_arr = np.load(infile)
        col = color_list[index]

        #ax.plot(data_arr[0], data_arr[1], color=col, lw=1)
        dm_arr = 10 ** logradialexp3(data_arr[0], dmpars['r0'], dmpars['D0'])
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
        
    save_fig(fig, args.fig_dir, f'DM_vs_impfact_all_lsm_inc{args.inc_range[0]}-{args.inc_range[1]}.pdf', args)
    plt.show(block=False)

    return fig

# ------------------------------------------------------------------------------------------------
def plot_dm_fit(df_dmpars, args):
    '''
    Plot DM0 and r0 vs stellar mass in a single panel
    Saves plot
    Returns axis handle
    '''
    outfilename = f'{args.fig_dir}/{Path(args.resfile_prefix).stem}_DM0_r0_vs_lsm_inc_{args.inc_range[0]}_{args.inc_range[1]}'
    res = pfns.plt_dmpars(df_dmpars, outfilename, 3.0, xcol='medlsm', y1col='D0', y2col='r0', x2col='medsfr', fit_robust=args.fit_robust)

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

    outfilename = f'{args.fig_dir}/{Path(args.resfile_prefix).stem}_table_DM0_r0_vs_lsm_inc_{args.inc_range[0]}_{args.inc_range[1]}.tex'

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

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')

    # -------------determining directories------------------
    catalog_name = f'{args.resfile_prefix}_allinc.txt'
    df_dmpars = read_dataframe(catalog_name)

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
            save_fig(fig, args.fig_dir, f'DM_vs_impfact_all_inc_all_lsm_multipanel.pdf', args)

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
