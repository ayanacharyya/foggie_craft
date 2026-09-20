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
from compute_host_dm import read_obs_catalog

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_SFMS_Popesso23(log_mass_min, log_mass_max, redshift, nbins=40):
    '''
    Computes an empirical SFMS based on Popesso+23 (https://arxiv.org/abs/2203.10487) Eq 10, for given redshift
    Then returns a two-part log SFR array based on an input minimum and maximum log mass; the two part are based on the lower limit of this empirical relation and any extrapolation below it
    Returns two tuples: (log_mass1, log_SFR1) and (log_mass2, log_SFR2)
    '''
    a0, a1, b0, b1, b2 = ufloat(0.2,0.02), ufloat(-0.034, 0.002), ufloat(-26.134,0.015), ufloat(4.722, 0.012), ufloat(-0.1925, 0.0011)  # Table 2, Eq 10
    log_mass_low_lim = 8.7 # lower limit of mass they fitted up to

    age_at_z = cosmo.age(redshift).value # Gyr

    log_mass1 = np.linspace(log_mass_min, log_mass_low_lim, nbins//2)
    log_SFR1 = (a1 * age_at_z + b1) * log_mass1 + b2 * (log_mass1) ** 2 + b0 + a0 * age_at_z

    log_mass2 = np.linspace(log_mass_low_lim, log_mass_max, nbins//2)
    log_SFR2 = (a1 * age_at_z + b1) * log_mass2 + b2 * (log_mass2) ** 2 + b0 + a0 * age_at_z

    return (log_mass1, log_SFR1), (log_mass2, log_SFR2)

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Popesso23(ax, redshift, color='cornflowerblue'):
    '''
    Computes an empirical SFMS based on Popesso+23 (https://arxiv.org/abs/2203.10487) Eq 10, for given redshift
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''
    (log_mass1, log_SFR1), (log_mass2, log_SFR2) = get_SFMS_Popesso23(ax.get_xlim()[0], ax.get_xlim()[1], redshift)

    #ax.plot(log_mass1, unp.nominal_values(log_SFR1), ls='dashed', c=color, lw=2)
    #ax.fill_between(log_mass1, unp.nominal_values(log_SFR1) - unp.std_devs(log_SFR1)/2, unp.nominal_values(log_SFR1) + unp.std_devs(log_SFR1)/2, alpha=0.3, facecolor=color)

    ax.plot(log_mass2, unp.nominal_values(log_SFR2), ls='dashed', c=color, lw=2, label=f'Popesso+23: z = {redshift}', zorder=-10)
    #ax.fill_between(log_mass2, unp.nominal_values(log_SFR2) - unp.std_devs(log_SFR2)/2, unp.nominal_values(log_SFR2) + unp.std_devs(log_SFR2)/2, alpha=0.3, facecolor=color)

    return ax

# ------------------------------------------------------------------------------------------------
def read_snap_list(args, filename="lsm_sfr_masses_upto_disk.csv", filepath=None):
    '''
    Reads in the list of FOGGIE simulation snapshots from <filename>
    Returns pandas dataframe
    '''
    if filepath is None:
        filepath =	args.data_dir / filename
    df = pd.read_csv(filepath, comment='#')
    df = df.drop_duplicates(subset=['halo', 'snap'], keep='last')

    # if 'sfr_100Myr' in df.columns:
    #     print('\nReplacing sfr column with sfr_100My column')
    #     df = df.drop('sfr', axis=1)
    #     df = df.rename(columns={'sfr_100Myr':'sfr'})
    
    df['log_sfr'] = np.log10(df['sfr'])

    df = df.rename(columns={'log_star_mass_from_snap': 'log_star_mass', 'log_gas_mass_from_profile': 'log_gas_mass'})
    df = df[(df['redshift'].between(args.z_range[0], args.z_range[1], inclusive='left')) & 
                    (df['log_star_mass'].between(args.lsm_range[0], args.lsm_range[1], inclusive='left')) & 
                    (df['log_sfr'].between(args.lsfr_range[0], args.lsfr_range[1], inclusive='left'))].reset_index(drop=True)
    df['log_ssfr'] = df['log_sfr'] - df['log_star_mass']

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
    fig, ax = plt.subplots(1, figsize=(3.6, 2.6))
    fig.subplots_adjust(left=0.12, bottom=0.13, right=0.98 if colorcol is None else 0.8, top=0.98)

    marker_arr = ['o', 's', '^', 'P', 'd', '*']
    color_arr = ['cornflowerblue', 'salmon', 'sienna', 'goldenrod', 'darkgreen', 'teal']
    df = df.sort_values(by=xcol)

    # popt = np.polyfit(df[xcol] - 10, df[ycol], 1)
    # print(popt)
    # ax.plot(df[xcol], np.poly1d(popt)(df[xcol].values - 10), c='k')

    # -------------plotting-------------
    for index, this_halo in enumerate(pd.unique(df['halo'])):
        df_sub = df[df['halo'] == this_halo]
        color = df_sub[colorcol] if colorcol is not None else 'cornflowerblue'
        p = ax.scatter(df_sub[xcol], df_sub[ycol], c=color, s=30, lw=1, ec='k', marker=marker_arr[index], alpha=0.8)
        #ax.plot(df_sub[xcol], df_sub[ycol], lw=0.5, c=color_arr[index])

    # ----------annotating and saving-----------
    ax.text(x=9.7, y=2.4, s=rf'{args.z_range[0]} $\leq z \less$ {args.z_range[1]}', ha='left', va='top', fontsize=args.fontsize / args.fontfactor)
    ax.text(x=9.7, y=2.2, s=f'#snapshots = {len(df)}', ha='left', va='top', fontsize=args.fontsize / args.fontfactor)
    ax = annotate_axes(ax, label_dict[xcol], label_dict[ycol], args=args, label='', clabel=label_dict[colorcol] if colorcol is not None else '', hide_cbar=colorcol is None, p=p, cticks_integer=True)
    ax.set_ylim(-0.6, 2.6)
    ax.set_xlim(9.6, 11.6)
    ax.set_xticks([10, 10.5, 11.0, 11.5])
    ax.set_yticks([0, 1, 2])

    # ------plotting literature MS-----------
    z = np.mean(args.z_range)
    ax = plot_SFMS_Popesso23(ax, z, color='cornflowerblue')

    # ------plotting observed FRBs-----------------
    df_obs, args.input_cat = read_obs_catalog(args.input_cat, args)
    df_obs = df_obs[(df_obs['z'].between(args.z_range[0], args.z_range[1], inclusive='left'))]
    q = ax.scatter(df_obs['log_mass'], df_obs['sfr_med'], s=70, lw=0.1, ec='k', marker='*', alpha=0.8, facecolor='goldenrod')

    # --------saving figure----------
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
