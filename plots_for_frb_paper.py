#!/usr/bin/env python3
"""
    Title :      plots_for_frb_paper
    Notes :      Make various plots to be used for the FRB paper
    Output :     Plots as PDF
    Author :     Ayan Acharyya
    Started :    31-03-26
    Examples :   run plots_for_frb_paper.py --plot_dm_lsm --lsm 9.5,9.75
                 run plots_for_frb_paper.py --plot_radprof
                 run plots_for_frb_paper.py --plot_dm_fit --fit_lsm_range 8.5,11.0
                 run plots_for_frb_paper.py --plot_dm_all_lsm --cmap tab10 --set_ylin
                 run plots_for_frb_paper.py --plot_dm_all_lsm
"""
from craft_header import *
from craft_utils import *
import plotfns as pfns

start_time = datetime.now()

# ------------------------------------------------------------------------------------------------
def read_dataframe(filename, interval_cols=['lsm_bin', 'lsfr_bin']):
    '''
    Function to read txt file as pandas dataframe and properly parse intervals
    Returns dataframe
    '''
    df = pd.read_csv(filename, sep='\t')
    df = df.drop_duplicates(subset=['lsm_bin', 'lsfr_bin'], keep='last')

    for col in interval_cols:
        temp_df = df[col].str.strip('()[]').str.split(',', expand=True).astype(float)
        df[col] = temp_df.apply(lambda x: pd.Interval(x[0], x[1], closed='right'), axis=1)
    
    return df

# ------------------------------------------------------------------------------------------------
def plot_dm_impfac_one_lsm_bin(df_dmpars, args):
    '''
    Plot DM vs Impact factor in a single panel, for a given stellar mass and inclination range
    Saves plot
    Returns figure handle
    '''
    # ---------setup figure---------------
    fig, ax = plt.subplots(1, figsize=(5, 4), layout='constrained')

    # -----------get required bin--------------------
    dmpars = df_dmpars[(df_dmpars['lsm_bin'] == pd.Interval(args.lsm_range[0], args.lsm_range[1])) & 
                       (df_dmpars['lsfr_bin'] == pd.Interval(args.lsfr_range[0], args.lsfr_range[1]))].iloc[0]
    infile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}/lsmzsfr_lsm_{dmpars["lsm_bin"].left}_{dmpars["lsm_bin"].right}_lsfr_{dmpars["lsfr_bin"].left}_{dmpars["lsfr_bin"].right}_1d.npy'

    data_arr = np.load(infile)

    ax.errorbar(data_arr[0], data_arr[1], fmt='bo', yerr=[data_arr[2], data_arr[3]], lw=1, ls='dashed')
    dm_arr = 10 ** logradialexp3(data_arr[0], dmpars['r0'], dmpars['D0'])
    ax.plot(data_arr[0], dm_arr, color='r', lw=1, ls='dotted')
    
    ax.set_xscale("log")
    ax.set_xticks(impbinegs[1:],impbinegs[1:])

    if not args.set_ylin:
        ax.set_yscale("log")
        ax.set_yticks(dm_ticks, dm_ticks)
    
    ax = annotate_axes(ax, "Impact factor (kpc)", "DM (pc cm$^{-3}$)", args=args, set_ticks=False)
        
    save_fig(fig, args.fig_dir, f'DM_vs_impfact_inc{args.inc_range[0]}-{args.inc_range[1]}_lsm_{args.lsm_range[0]}_{args.lsm_range[1]}_lsfr_{args.lsfr_range[0]}_{args.lsfr_range[1]}.pdf', args)
    plt.show(block=False)

    return fig

# ------------------------------------------------------------------------------------------------
def plot_dm_impfac_all_lsm_bin(df_dmpars, args, cmap='tab10'):
    '''
    Plot DM vs Impact factor in a single panel, for a list of stellar mass ranges and a given inclination range
    Saves plot
    Returns figure handle
    '''
    # ---------setup figure---------------
    fig, ax = plt.subplots(1, figsize=(12, 7), layout='constrained')
    color_list = plt.get_cmap(cmap)(range(12))

    # -----------loop through mass bins--------------------
    for index, dmpars in df_dmpars.iterrows(): 
        infile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}/lsmzsfr_lsm_{dmpars["lsm_bin"].left}_{dmpars["lsm_bin"].right}_lsfr_{dmpars["lsfr_bin"].left}_{dmpars["lsfr_bin"].right}_1d.npy'

        data_arr = np.load(infile)
        col = color_list[index]

        #ax.plot(data_arr[0], data_arr[1], color=col, lw=1)
        dm_arr = 10 ** logradialexp3(data_arr[0], dmpars['r0'], dmpars['D0'])
        ax.plot(data_arr[0], dm_arr, color=col, lw=2, ls=lslist[index % len(lslist)], label=f'{dmpars["lsm_bin"].left:.1f}-{dmpars["lsm_bin"].right:.1f}')
        ax.fill_between(data_arr[0], data_arr[1] - data_arr[2], data_arr[1] + data_arr[3], color=col,alpha=0.1)
		
        if args.set_ylin: ax.text(0.98, 0.98 - index * 0.03, f'{dmpars["lsm_bin"].left:.1f}' + r' $\leq \log$ M/M$_\odot \leq$ ' + f'{dmpars["lsm_bin"].right:.1f}', c=col, ha='right', va='top', transform=ax.transAxes, fontsize=args.fontsize / args.fontfactor)
        else: ax.text(0.02, 0.02 + index * 0.03, f'{dmpars["lsm_bin"].left:.1f}' + r' $\leq \log$ M/M$_\odot \leq$ ' + f'{dmpars["lsm_bin"].right:.1f}', c=col, ha='left', va='bottom', transform=ax.transAxes, fontsize=args.fontsize / args.fontfactor)

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
    Returns figure handle
    '''
    outfilename = f'{args.fig_dir}/DM0_r0_vs_lsm_inc_{args.inc_range[0]}_{args.inc_range[1]}'
    res = pfns.plt_dmpars(df_dmpars, argsfit_lsm_range, outfilename, 3.0, xcol='medlsm', y1col='D0', y2col='r0', x2col='medsfr')

    return res

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fig_dir = args.plot_dir / 'plots_for_paper'
    args.fig_dir.mkdir(exist_ok=True, parents=True)

    # -------------determining directories------------------
    catalog_name = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}.txt'
    df_dmpars = read_dataframe(catalog_name)

    # -------------calling plotting functions------------------
    if args.plot_dm_lsm:
        for index, this_lsm_bin in enumerate(args.lsm_bins):
            print(f'\nRunning ({index + 1}/{len(args.lsm_bins)}) for stellar mass bin {this_lsm_bin}..\n')
            args.lsm_range = this_lsm_bin
            fig = plot_dm_impfac_one_lsm_bin(df_dmpars, args)
    if args.plot_dm_all_lsm:
        fig = plot_dm_impfac_all_lsm_bin(df_dmpars, args, cmap=args.cmap)
    if args.plot_dm_fit:
        fig = plot_dm_fit(df_dmpars, args)
    
    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
