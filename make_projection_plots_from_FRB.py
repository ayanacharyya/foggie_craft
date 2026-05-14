#!/usr/bin/env python3
"""
    Title :      make_projection_plots_from_FRB
    Notes :      Make various plots to be used for the FRB paper
    Output :     Plots as PDF
    Author :     Ayan Acharyya
    Started :    31-03-26
    Examples :   run make_projection_plots_from_FRB.py --system ayan_ssd --halo 2392 --run 2162_39 --output DD1360 --upto_kpc 50 --reskpc 0.5
                 run make_projection_plots_from_FRB.py --system ayan_pleiades --foggie_dir /nobackupp19/aachary2/LowZRuns/ --halo 2392 --run 2162_39 --output DD1360 --upto_kpc 50 --reskpc 0.5
                 run make_projection_plots_from_FRB.py --system ayan_local --halo 8508 --output RD0027 --upto_kpc 100 --reskpc 0.5
"""
from craft_header import *
from craft_utils import *
setup_plot_style()
import plotfns as pfns

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------
def plot_projections_from_frb(data_faceon, data_edgeon, args, crpix=None, cdelt=None, crval=0, quant_label='density', clim=None,  cmap='viridis', takelog=True, clabel=None, annotate_labels=None, dpi=300):
    '''
    Plots the input face on and edge on FRBs onto the same figure
    Returns figure handle
    '''
    if takelog:
        data_faceon = np.log10(data_faceon)
        data_edgeon = np.log10(data_edgeon)

    if crpix is None:
        crpix = np.shape(data_edgeon)[0]/2.
    
    # ------plotting onto a matplotlib figure--------------
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    right, top, bottom, left, wspace = 0.85, 0.98, 0.1, 0.12, 0.0
    fig.subplots_adjust(right=right, top=top, bottom=bottom, left=left, wspace=wspace)
    fontsize = args.fontsize

    # ----------plotting axes----------------
    p = axes[0].imshow(data_faceon, cmap=cmap, vmin=clim[0] if clim is not None else None, vmax=clim[1] if clim is not None else None)
    p = axes[1].imshow(data_edgeon, cmap=cmap, vmin=clim[0] if clim is not None else None, vmax=clim[1] if clim is not None else None)

    # ---------------making colorbar------------------------
    offset, width = 0.075, 0.02
    cax = fig.add_axes([right, bottom + offset, width, top - bottom - 2 * offset])
    cbar = fig.colorbar(p, orientation='vertical', cax=cax)
    cbar.ax.tick_params(labelsize=fontsize, width=2.5, length=5)
    cbar.set_label(clabel, fontsize=fontsize) # this is to have more control on the display of the color label

    # ---------------prepping axes------------------------
    for index in range(len(axes)):
        ax = axes[index]
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        if cdelt is None: ax.set_xticklabels(['%.1F' % item for item in ax.get_xticks()], fontsize=fontsize)
        else: ax.set_xticklabels(['%.1F' % ((item - crpix) * cdelt + crval) for item in ax.get_xticks()], fontsize=fontsize)
        ax.set_xlabel('Offset (kpc)', fontsize=fontsize)
        if index == 0:
            if cdelt is None: ax.set_yticklabels(['%.1F' % item for item in ax.get_yticks()], fontsize=fontsize)
            else: ax.set_yticklabels(['%.1F' % ((item - crpix) * cdelt + crval) for item in ax.get_yticks()], fontsize=fontsize)
            ax.set_ylabel('Offset (kpc)', fontsize=fontsize)
        else:
            ax.set_yticklabels(['' % item for item in ax.get_yticks()])
            ax.set_ylabel('')
        
    # ---------------making annotations------------------------
    if annotate_labels is None:
        axes[0].text(0.97, 0.95, 'z = %.2F' % args.current_redshift, c='white', ha='right', va='top', transform=axes[0].transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))
        axes[0].text(0.97, 0.85, 't = %.1F Gyr' % args.current_time, c='white', ha='right', va='top', transform=axes[0].transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))
    else:
        for index, label in enumerate(annotate_labels):
            axes[index % 2].text(0.97, 0.95 - (index // 2) * 0.1, label, c='white', ha='right', va='top', transform=axes[index % 2].transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))

    axes[0].text(0.98, 0.02, 'Face on', c='white', ha='right', va='bottom', transform=axes[0].transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))
    axes[1].text(0.98, 0.02, 'Edge on', c='white', ha='right', va='bottom', transform=axes[1].transAxes, fontsize=fontsize, bbox=dict(facecolor='k', alpha=0.3, edgecolor='k'))

    # ---------------saving fig------------------------
    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try: mplcyberpunk.make_lines_glow()
        except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    figname = f'{args.output}_{args.halo}_diskrel_{quant_label}{args.upto_text}.png'
    if 'el' in quant_label.lower(): fig_dir = args.fig_dir + 'electron_density/'
    else: fig_dir = args.fig_dir + 'gas_density/'
    save_fig(fig, Path(fig_dir), figname, args, dpi=dpi)

    return fig

# ---------------main code------------------------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.fontsize == 8: args.fontsize = 15

    quant_dict = {'density':['density', 'Gas density', 'Msun/pc**3', -2.5, 2.5, 'cornflowerblue', 'cividis', True, 'Msun/pc**2', r'Gas density [M$_\odot$ pc$^{-2}$]'], 
                'el_density':['El_number_density', 'Electron density', 'cm**-3', 0, 220, 'cornflowerblue', 'viridis', False, 'pc*cm**-3', r'DM [pc cm$^{-3}$]']
                } # for each quantity: [yt field, label in plots, units, lower limit in log, upper limit in log, color for scatter plot, colormap, whether to take log, units for projection plot, units to display in projection plot]
    quant_arr = ['el_density']#, 'density']

    # -----------determining directories----------------
    args.fig_dir = args.output_dir + 'plots/'
    Path(args.fig_dir).mkdir(parents=True, exist_ok=True)
    
    args.fits_dir = args.output_dir + 'data/'
    Path(args.fits_dir).mkdir(parents=True, exist_ok=True)

    args.res_text = f'_res{args.reskpc:.1f}kpc'
    args.upto_text = '_upto%.1Fckpchinv' % args.upto_kpc if args.docomoving and args.upto_kpc is not None else '_upto%.1Fkpc' % args.upto_kpc if args.upto_kpc is not None else f'_upto{args.re:.1f}re'

    fitsname = Path(args.fits_dir) / f'{args.output}_{args.halo}_FRB_{quant_dict[quant_arr[0]][0]}{args.upto_text}{args.res_text}.fits'

    # ---------reading in fits file----------------------
    print(f'Trying to read {fitsname}..')
    hdul = fits.open(fitsname)
    sfr = hdul[0].header['SFR']
    log_mstar = hdul[0].header['LOG_MSTAR']

    # ---------looping over quants----------------------
    for index, quant in enumerate(quant_arr):
        print(f'Making and plotting FRB for {quant} which is {index+1} out of {len(quant_arr)} quantities..')
        data_faceon = hdul[f'{quant_dict[quant][1]} FACE ON PROJ'].data
        data_edgeon = hdul[f'{quant_dict[quant][1]} EDGE ON PROJ'].data

        cdelt = hdul[f'{quant_dict[quant][1]} FACE ON PROJ'].header['CDELT1']
        crpix = hdul[f'{quant_dict[quant][1]} FACE ON PROJ'].header['CRPIX1']
        crval = hdul[f'{quant_dict[quant][1]} FACE ON PROJ'].header['CRVAL1']

        # -------------make the plots--------------
        fig_diskrel = plot_projections_from_frb(data_faceon, data_edgeon, args,
                                                crpix=crpix, cdelt=cdelt, crval=crval, 
                                                quant_label=quant_dict[quant][0], 
                                                clim=[quant_dict[quant][3], quant_dict[quant][4]] if quant_dict[quant][3] is not None else None,  
                                                cmap=quant_dict[quant][6], 
                                                takelog=quant_dict[quant][7], 
                                                clabel=quant_dict[quant][9],
                                                annotate_labels = [rf'SFR = {sfr:.1f} M$_\odot$/yr', rf'$\log$ (M$_*$/M$_\odot$) = {log_mstar:.1f}'],
                                                )


    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))

