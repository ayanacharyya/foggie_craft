#!/usr/bin/env python3
"""
    Title :      plot_sfms.py
    Notes :      Make stellar mass vs SFR of FOGGIE snapshots
    Output :     Plots as PDF
    Author :     Ayan Acharyya
    Started :    10-05-26
    Examples :   run plot_sfms.py --sample high_mass
                 run plot_sfms.py --sample high_mass --z_range 0,2
"""
from craft_header import *
from craft_utils import *
setup_plot_style()

start_time = datetime.now()

# ------------------------------------------------------------------------------------------------
def read_snap_list(args, filename="lsm_sfr_masses_upto_disk.txt"):
    '''
    Reads in the list of FOGGIE simulation snapshots from <filename>
    Returns pandas dataframe
    '''
    filespecs =	args.data_dir / filename
    df = pd.read_csv(filespecs, sep=r'\s+', engine='python', comment='#')

    if 'sfr_100Myr' in df.columns:
        print('\nReplacing sfr column with sfr_100My column')
        df = df.drop('sfr', axis=1)
        df = df.rename(columns={'sfr_100Myr':'sfr'})
    
    df['log_sfr'] = np.log10(df['sfr'])

    df = df.rename(columns={'log_star_mass_from_snap': 'log_star_mass', 'log_gas_mass_from_profile': 'log_gas_mass'})
    df = df[(df['redshift'].between(args.z_range[0], args.z_range[1])) & 
                    (df['log_star_mass'].between(args.lsm_range[0], args.lsm_range[1])) & 
                    (df['log_sfr'].between(args.lsfr_range[0], args.lsfr_range[1]))].reset_index(drop=True)
    print (f'\t\tFound {len(df)} snapshots, within log mass range {args.lsm_range}, log sfr range {args.lsfr_range} and redshift range {args.z_range}')

    return df

# ------------------------------------------------------------------------------------------------
def plot_sfms(df, args, xcol='log_star_mass', ycol='sfr', colorcol=None):
    '''
    Plot mass vs SFR in a single panel
    Saves plot
    Returns axis handle
    '''
    label_dict = {'log_star_mass': r'$\log$ [M$_*$/M$_{\odot}$]', 
                  'sfr': r'SFR [M$_{\odot}$ yr$^{-1}$]', 
                  'log_sfr': r'$\log$ [SFR/M$_{\odot}$ yr$^{-1}$]', 
                  'redshift': 'Redshift'}
    
    # --------------setting up figure---------
    fig, ax = plt.subplots(1, figsize=(6, 5))
    fig.subplots_adjust(left=0.12, bottom=0.12, right=0.98 if colorcol is None else 0.88, top=0.98)

    marker_arr = ['o', 's', '^', 'P', 'd', '*']
    color_arr = ['cornflowerblue', 'salmon', 'sienna', 'goldenrod', 'darkgreen', 'teal']
    df = df.sort_values(by=xcol)

    popt = np.polyfit(df[xcol] - 10, df[ycol], 1)
    print(popt)
    ax.plot(df[xcol], np.poly1d(popt)(df[xcol].values - 10), c='k')

    # -------------plotting-------------
    for index, this_halo in enumerate(pd.unique(df['halo'])):
        df_sub = df[df['halo'] == this_halo]
        color = df_sub[colorcol] if colorcol is not None else 'cornflowerblue'
        p = ax.scatter(df_sub[xcol], df_sub[ycol], c=color, s=70, lw=1, ec='k', marker=marker_arr[index], alpha=0.8)
        ax.plot(df_sub[xcol], df_sub[ycol], lw=0.5, c=color_arr[index])

    # ----------annotating and saving-----------
    ax = annotate_axes(ax, label_dict[xcol], label_dict[ycol], args=args, label='', clabel=label_dict[colorcol] if colorcol is not None else '', hide_cbar=colorcol is None, p=p, cticks_integer=True)
    ax.set_ylim(-0.4, 2)

    figname = f'{xcol}_vs_{ycol}_lsm_range_{args.lsm_range[0]}_{args.lsm_range[1]}_lsfr_range_{args.lsfr_range[0]}_{args.lsfr_range[1]}_zrange_{args.z_range[0]}_{args.z_range[1]}.png'
    save_fig(fig, args.plot_dir, figname, args=args, dpi=300)

    return ax

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')
    args.lsm_range = args.lsm_bins[0]

    # -----------------------read in snapshot list----------------------------------
    df_snap = read_snap_list(args)
    if (len(df_snap) < 1):
        print('\t\tNo snapshot found. Continuing to next loop iteration... ')

    # -----------make SFMS plot--------------------
    ax = plot_sfms(df_snap, args, xcol='log_star_mass', ycol='log_sfr', colorcol='redshift')

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
