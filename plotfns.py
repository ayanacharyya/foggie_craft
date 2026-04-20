#
#	Plotting functions
#
#								AB, August 2024
#
#	Function list
#
#
#	--------------------------	Import modules	---------------------------
from craft_utils import *
setup_plot_style()

#	----------------------------------------------------------------------------------------------------------
def plot_nerad(cubes, inc_ranges, title, outfile, fig_size, hide=False, subtitle='', given_ax=None, fortalk=False):
    #	Plot Radial profile of ne within the inclination range

    if given_ax is None:
        fig 	= plt.figure(figsize=(1.2 * fig_size, fig_size))
        ax	 	= fig.add_axes([0.17,0.15,0.82,0.84])
    else:
        ax = given_ax
    
    # ------------loop over inclination ranges-----------------
    for i, inc_range in enumerate(inc_ranges):
        binned_ne	= np.zeros((len(radbins)-1, 6), dtype=np.float32)		
        binrelnes	= []

        for cube in cubes:
            relinds	= np.where((cube.inclination - np.deg2rad(inc_range[0]))*(cube.inclination - np.deg2rad(inc_range[1])) <= 0.0)	
            relne	= cube.neincrad[relinds]
            binrelnes.append(relne)				

        binrelne	= np.concatenate(binrelnes, axis=0)

        #	The ugly binning in radius
        for k in range (0,len(radbins)-1):
            rel2inds		= np.where((cubes[0].radkpc - radbins[k]) * (cubes[0].radkpc - radbins[k+1]) <= 0.0)
            binned_ne[k,0]	= (radbins[k]+radbins[k+1])/2.0
            binned_ne[k,1:6]= np.percentile(binrelne[:,rel2inds], (16, 25, 50, 75, 84))

        ax.fill_between(binned_ne[:,0], binned_ne[:,1], binned_ne[:,5], color=shlist[i],alpha=0.1)
        #ax.plot(binned_ne[:,0], binned_ne[:,1], colist[i]+'--', lw=0.5)
        #ax.plot(binned_ne[:,0], binned_ne[:,5], colist[i]+'--', lw=0.5)
        ax.plot(binned_ne[:,0], binned_ne[:,3], clist[i]+marklist[i], markersize=6, fillstyle="none", \
        label="%d$^{\circ}-$%d$^{\circ}$"%(inc_range[0],inc_range[1]))

    ax.legend(loc="lower left",frameon=False)
    #ax.set_xlim([-1,205])
    #ax.set_ylim([3.5e-6,1.5])

    ax.text(x=0.3, y=0.8, s=title, ha='left', va='bottom', transform=ax.transAxes)
    ax.text(x=0.3, y=0.9, s=subtitle, ha='left', va='bottom', transform=ax.transAxes)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Radius (kpc)')
    ax.set_ylabel('$n_e$ (cm$^{-3}$)')

    if given_ax is None:
        figname = Path(outfile)
        save_fig(fig, figname.parent, figname.name, fortalk=fortalk)
        if hide: plt.close()
        else: plt.show(block=False)

    return ax

#	----------------------------------------------------------------------------------------------------------	
def pltdm_ind_imf(df, lsm, sfr, inc_range, redshift, outfilename, fig_size, hide=False, bin_col1='impf', bin_col2='distmaj', data_col='losdm', given_ax=None, fortalk=False):
	#	Plot LoSDM vs impact factor and dmaj_proj for a given inclination range 

    percentiles = [16, 25, 50, 75, 84]
    maps = {}
    
    for p in percentiles:
        ret = binned_statistic_2d(df[bin_col1], df[bin_col2], df[data_col], statistic=lambda x, p_val=p: np.nanpercentile(x, p_val) if len(x) > 0 else np.nan, bins=[impbinegs, impbinegs])
        maps[f'p{p}'] = 0.5 * ret.statistic

    dmavg	= maps['p50']	
    dmdiff	= maps['p84'] - maps['p16']

    # ------------plot the profile----------
    if given_ax is None:
        fig 	= plt.figure(figsize=(1.2 * fig_size, fig_size))
        ax1	 	= fig.add_axes([0.07,0.15,0.40,0.83])
        ax2	 	= fig.add_axes([0.49,0.15,0.40,0.83])
        ax3	 	= fig.add_axes([0.90,0.15,0.02,0.83])
    else:
        ax1, ax2, ax3 = given_ax
    
    medim	= ax1.imshow(dmavg, origin='lower', interpolation='none', aspect='auto', cmap="Blues", vmin=0, vmax=maxdmcol)

    ax1.text(x=0, y=len(impbinegs) - 2, s="Median")
    ax1.text(x=0, y=len(impbinegs) - 3, s="z = %.1f"%redshift)
    ax1.text(x=0, y=len(impbinegs) - 3.8, s="log ($M_* / M_{\odot}$) = %.2f"%lsm)
    ax1.text(x=0, y=len(impbinegs) - 4.6, s="SFR = %.2f $M_{\odot} yr^{-1}$"%sfr)
    ax1.set_xticks(np.arange(0,len(impbinegs) - 1, 1) - 0.5,impbinegs[:-1])
    ax1.set_yticks(np.arange(0,len(impbinegs) - 1, 1) - 0.5,impbinegs[:-1])
    ax1.set_xlabel("Impact factor (kpc)")
    ax1.set_ylabel("Projected distance from major axis (kpc)")

    ax2.imshow(dmdiff, origin='lower', interpolation='none', aspect='auto', cmap="Blues", vmin=0, vmax=maxdmcol)
    ax2.text(x=0, y=len(impbinegs) - 1-1, s="Scatter")
    ax2.text(x=0, y=len(impbinegs) - 1-2, s="%.1f$^{\circ}$ < i < %.1f$^{\circ}$"%(inc_range[0],inc_range[1]))
    ax2.set_xticks(np.arange(0,len(impbinegs) - 1, 1) - 0.5,impbinegs[:-1])
    ax2.set_yticks(np.arange(0,len(impbinegs) - 1, 1) - 0.5,impbinegs[:-1])
    ax2.set_yticklabels([])
    ax2.set_xlabel("Impact factor (kpc)")

    fig.colorbar(medim, cax=ax3, label="DM (pc cm$^{-3}$)")	

    if given_ax is None:
        figname = Path(outfilename + ".png")
        save_fig(fig, figname.parent, figname.name, fortalk=fortalk)
        if hide: plt.close()
        else: plt.show(block=False)

    return ax1, ax2, ax3

#	----------------------------------------------------------------------------------------------------------	
def pltdm_ind_imf_1d(df, lsm, sfr, parlims, outfilename, fig_size, hide=False, bin_col='impf', data_col='losdm', given_ax=None, nobj=None, lsfr_lims=None, fortalk=False):
	#	Plot LoSDM vs impact factor for a given inclination range 
			
    df['bin'] = pd.cut(df[bin_col], bins=impbinegs, include_lowest=True)
    percentiles = [0.16, 0.25, 0.50, 0.75, 0.84] # Quantiles are 0 to 1 (so 16th percentile is 0.16)
    stats = df.groupby('bin')[data_col].quantile(percentiles).unstack()
    
    stats.columns = [f'p{int(p*100)}' for p in percentiles]
    for col in stats.columns: stats[col] *= 0.5
    stats = stats[np.isfinite(stats['p50'])]

    dmavg	= stats['p50']
    dmlower	= stats['p50'] - stats['p16']
    dmhier	= stats['p84'] - stats['p50']

    impx	= (impbinegs[:-1] + impbinegs[1:]) / 2.0

    #popt,pcov	= curve_fit(logbeselk, impx, np.log10(dmavg), p0=(2.0,100.0), maxfev=1000)
    popt,pcov	= curve_fit(logradialexp3, impx, np.log10(dmavg), p0=(10.0,100.0))
    perr 		= np.sqrt(np.diag(pcov))
    print(popt,perr)

    # ------------save the profile----------
    data_arr = np.array([impx, dmavg, dmlower, dmhier])
    np.save(outfilename + ".npy", data_arr)

    # ------------plot the profile----------
    if given_ax is None:
        fig 	= plt.figure(figsize=(1.2 * fig_size, fig_size))
        ax	 	= fig.add_axes([0.17,0.15,0.82,0.84])
    else:
        ax = given_ax
    
    ax.errorbar(impx, dmavg, yerr=[dmlower,dmhier],fmt='bo',lw=1,markersize=4,capsize=4)
    #ax.plot(impx, 10**logbeselk(impx,*popt),'k--',lw=1)
    ax.plot(impx, 10**logradialexp3(impx,*popt),'k--',lw=1)
    ax.plot(impx, 10**logradialexp3(impx,10.0**(0.61-0.53*(lsm-10)),10.0**(2.15+0.24*(lsm-10))),'r:',lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim([0.5, 3 * maxdmcol])
    ax.set_xlim([0.25 * impbinegs[1], 1.5 * impbinegs[-1]])
    ax.set_yticks(dm_ticks, dm_ticks)
    ax.set_ylabel("DM (pc cm$^{-3}$)")	
    ax.set_xlabel("Impact factor (kpc)")

    if lsfr_lims is None:
        ax.set_xticks(impbinegs[1:], impbinegs[1:])
        nobj_text = '' if nobj is None else f' ({nobj})'
        ax.text(x=0.4*impbinegs[1], y=300, s="%.2f < log ($M_* / M_{\odot}$) < %.2f%s"%(parlims[0],parlims[1], nobj_text))
        ax.text(x=0.4*impbinegs[1], y=1.6, s="log ($M_* / M_{\odot}$) = %.2f"%lsm)
        ax.text(x=0.4*impbinegs[1], y=0.8, s="SFR = %.2f $M_{\odot} yr^{-1}$"%sfr)
        ax.text(x=1.0*impbinegs[-4], y=150, s="$D_0$ = %d $\pm$ %d"%(popt[1],perr[1]))
        ax.text(x=1.0*impbinegs[-4], y=75, s="$r_0$ = %.1f $\pm$ %.1f"%(popt[0],perr[0]))
    else:
        ax.set_xticks(impbinegs[1::2], impbinegs[1::2])
        nobj_text = ''# if nobj is None else f' ({nobj})'
        ax.text(x=0.4*impbinegs[1], y=250, s="log $M_*$ $\in$ [%.1f, %.1f]%s"%(parlims[0],parlims[1], nobj_text))
        ax.text(x=0.4*impbinegs[1], y=1, s="log SFR $\in$ [%.1f, %.1f]"%(lsfr_lims[0],lsfr_lims[1]))

    if given_ax is None:
        figname = Path(outfilename + ".pdf")
        save_fig(fig, figname.parent, figname.name, fortalk=fortalk)
        if hide: plt.close()
        else: plt.show(block=False)

    return (popt, perr, ax)


#	----------------------------------------------------------------------------------------------------------
#	----------------------------------------------------------------------------------------------------------
def plt_dmpars_annotate(ax, ylabel, ylim, yticks, popt, perr, madex, yscale_log=True, coeff_label=['D', r'$\gamma$']):
	#   Annotate the DM0 or r0 vs stellar mass subplot

    if yscale_log: ax.set_yscale('log')
    ax.set_ylim(ylim)
    if yticks is not None:
        ax.set_yticks(yticks, yticks)
    ax.set_ylabel(ylabel)

    ax.set_xlim([8.51,11.49])
    ax.set_xlabel(r"log ($M_*/M_{\odot}$)")

    if len(popt) <= 2:
        ax.text(x=0.3, y=0.1, s="%s$_{10} = %.2f \pm %.2f$"%(coeff_label[0], popt[1],perr[1]), transform=ax.transAxes)
        ax.text(x=0.3, y=0.2, s="%s $= %.2f \pm %.2f$"%(coeff_label[1], popt[0],perr[0]), transform=ax.transAxes)
    
    ax.text(x=0.3, y=0.9, s="MAD = %.2f dex"%(madex), transform=ax.transAxes)

    return ax

#	----------------------------------------------------------------------------------------------------------
def plt_sfms(xdata_arr, ydata_arr, outfilename, fortalk=False):
    '''
    Mass vs SFR plot
    '''

    fillstyle_arr = ['none', 'full']

    fig, ax = plt.subplots(1, figsize=(6, 4), layout='constrained')

    for index in range(len(xdata_arr)):
        xdata = xdata_arr[index]
        ydata = ydata_arr[index]

        ax.plot(xdata, ydata, 'bo', fillstyle=fillstyle_arr[index])
    
    ax = annotate_axes(ax, r"log ($M_*/M_{\odot}$)", r'$\log SFR (M_{\odot}/yr)$', fontsize=8, fontfactor=1)
    
    # ------------save figure-------------
    figname = Path(outfilename + "_sfms.pdf")
    save_fig(fig, figname.parent, figname.name, fortalk=fortalk)

    return fig

#	----------------------------------------------------------------------------------------------------------
def plt_dmpars_fit_par(df, xcol, ycol, ax, ylabel, ylim, yticks, ycol2, ax2, ylabel2, ylim2, yticks2, fit_robust=True, outfilename=None, fortalk=False):
    # Fit DM0 or r0 vs log stellar mass

    # ------------now fitting D0-------------
    df_fit = df.copy()

    # ------------fit the parameter D0----------
    do_fit = True       
    while do_fit:
        popt,pcov	= np.polyfit(df_fit[xcol] - 10, np.log10(df_fit[ycol]), 1, cov=True)
        perr 		= np.sqrt(np.diag(pcov))
        dm0fit		= np.poly1d(popt)
        df_fit['devdex']		= np.log10(df_fit[ycol]) - dm0fit(df_fit[xcol] - 10)
        mad		= np.nanmedian(np.abs(df_fit['devdex']))

        df_fit['outlier_fl'] = np.abs(df_fit['devdex']) > scale_fit_thresh * mad # scale_fit_thresh is in globalpars.py
        
        if fit_robust and df_fit['outlier_fl'].any():
            unlucky = np.abs(df_fit['devdex']) == np.max(np.abs(df_fit['devdex']))
            df_fit = df_fit[~unlucky]
        else:
            do_fit = False
    
    df['devdex']		= np.log10(df[ycol]) - dm0fit(df[xcol] - 10)

    # ------------plot the parameter D0----------
    ax.plot(df[xcol], 10.0 ** dm0fit(df[xcol] - 10), 'k--')
    ax.errorbar(df[xcol], df[ycol], df['e' + ycol], fmt='bo', markersize=5, lw=0.5, capsize=2, fillstyle='none')
    ax.errorbar(df_fit[xcol], df_fit[ycol], df_fit['e' + ycol], fmt='bo', markersize=5, lw=0.5, capsize=2)

    ax = plt_dmpars_annotate(ax, ylabel, ylim, yticks, popt, perr, np.nanmedian(np.abs(df['devdex'])), coeff_label=['D', r'$\gamma$'])

    # ------------plot SFMS--------------------
    if outfilename is not None:
        mass_col, sfr_col = 'medlsm', 'medsfr'
        fig = plt_sfms([df[mass_col].values, df_fit[mass_col].values], [df[sfr_col].values, df_fit[sfr_col].values], outfilename, fortalk=fortalk)

    # ------------now fitting r0-------------
    popt2,pcov2	= np.polyfit(df_fit[xcol] - 10, np.log10(df_fit[ycol2]), 1, cov=True)
    perr2 		= np.sqrt(np.diag(pcov2))
    r0fit		= np.poly1d(popt2)
    
    df['devdex2']		= np.log10(df[ycol2]) - r0fit(df[xcol] - 10)

    # ------------plot the parameter r0----------
    ax2.plot(df[xcol], 10.0 ** r0fit(df[xcol] - 10), 'k--')
    ax2.errorbar(df[xcol], df[ycol2], df['e' + ycol2], fmt='bo', markersize=5, lw=0.5, capsize=2, fillstyle='none')
    ax2.errorbar(df_fit[xcol], df_fit[ycol2], df_fit['e' + ycol2], fmt='bo', markersize=5, lw=0.5, capsize=2)

    ax2 = plt_dmpars_annotate(ax2, ylabel2, ylim2, yticks2, popt2, perr2, np.nanmedian(np.abs(df['devdex2'])), coeff_label=['R', r'$\eta$'])

    return popt, perr, popt2, perr2

#	----------------------------------------------------------------------------------------------------------
def plt_dmpars_fit_multipar(df, xcol, x2col, ycol, ax, ylabel, ylim, yticks, ycol2, ax2, ylabel2, ylim2, yticks2, yscale_log=True, fit_robust=True, outfilename=None, fortalk=False):
    # Fit DM0 or r0 vs log stellar mass

    df_fit = df.copy()

    # ------------fit the param D0----------
    do_fit = True       
    while do_fit:
        lsmsfr_fit	= np.array([df_fit[xcol] - 10, np.log10(df_fit[x2col])]) 
        popt,pcov	= curve_fit(linearxy, lsmsfr_fit, np.log10(df_fit[ycol]), p0=(0.25,0.0,2.16))
        perr 		= np.sqrt(np.diag(pcov))
        df_fit['devdex']		= np.log10(df_fit[ycol]) - linearxy(lsmsfr_fit, *popt)
        mad         = np.nanmedian(np.abs(df_fit['devdex']))

        df_fit['outlier_fl'] = np.abs(df_fit['devdex']) > scale_fit_thresh * mad # scale_fit_thresh is in globalpars.py
        
        if fit_robust and df_fit['outlier_fl'].any():
            unlucky = np.abs(df_fit['devdex']) == np.max(np.abs(df_fit['devdex']))
            df_fit = df_fit[~unlucky]
        else:
            do_fit = False
    
    lsmsfr	= np.array([df[xcol] - 10, np.log10(df[x2col])]) 
    df['devdex']		= np.log10(df[ycol]) - linearxy(lsmsfr, *popt)
    
    # ------------plot the parameter D0----------
    ax.axhline(c='k',ls='--')
    ax.errorbar(df[xcol], df['devdex'], df['e' + ycol] / (np.log(10.0) * df[ycol]), fmt='bo', markersize=5, lw=0.5, capsize=2, fillstyle='none')
    ax.errorbar(df_fit[xcol], df_fit['devdex'], df_fit['e' + ycol] / (np.log(10.0) * df_fit[ycol]), fmt='bo', markersize=5, lw=0.5, capsize=2)
    
    ax = plt_dmpars_annotate(ax, ylabel, ylim, yticks, popt, perr, np.nanmedian(np.abs(df['devdex'])), yscale_log=yscale_log)

    # ------------plot SFMS--------------------
    if outfilename is not None:
        mass_col, sfr_col = 'medlsm', 'medsfr'
        fig = plt_sfms([df[mass_col].values, df_fit[mass_col].values], [df[sfr_col].values, df_fit[sfr_col].values], outfilename, fortalk=fortalk)

    # ------------now fitting r0-------------
    popt2,pcov2	= curve_fit(linearxy, lsmsfr_fit, np.log10(df_fit[ycol2]), p0=(-0.5,0.0,1.0))
    perr2 		= np.sqrt(np.diag(pcov2))
    df_fit['devdex2']		= np.log10(df_fit[ycol2]) - linearxy(lsmsfr_fit, *popt2)
    df['devdex2']		= np.log10(df[ycol2]) - linearxy(lsmsfr, *popt2)

    # ------------plot the parameter r0----------
    ax2.axhline(c='k',ls='--')
    ax2.errorbar(df[xcol], df['devdex2'], df['e' + ycol2] / (np.log(10.0) * df[ycol2]), fmt='bo', markersize=5, lw=0.5, capsize=2, fillstyle='none')
    ax2.errorbar(df_fit[xcol], df_fit['devdex2'], df_fit['e' + ycol2] / (np.log(10.0) * df_fit[ycol2]), fmt='bo', markersize=5, lw=0.5, capsize=2)

    ax2 = plt_dmpars_annotate(ax2, ylabel2, ylim2, yticks2, popt2, perr2, np.nanmedian(np.abs(df['devdex2'])), yscale_log=yscale_log)

    return popt, perr, popt2, perr2

#	----------------------------------------------------------------------------------------------------------
def plt_dmpars(df, outfilename, fig_size, xcol='medlsm', y1col='D0', y2col='r0', x2col='medsfr', fit_robust=True, fortalk=False):
	#	Plot LoSDM vs impact factor for a given inclination range
	
    fig 	= plt.figure(figsize=(2.4 * fig_size, fig_size))
    ax1	 	= fig.add_axes([0.08,0.15,0.42,0.83])
    ax2	    = fig.add_axes([0.57,0.15,0.42,0.83])

    popt, perr, popt2, perr2 = plt_dmpars_fit_par(df, xcol, y1col, ax1, r"$D_0\:(pc \: cm^{-3})$", [40,420], [50,100,200,400], 
                                  y2col, ax2, r"$r_0$ (kpc)", None, [1,2,4,8,16,32],
                                  fit_robust=fit_robust, outfilename=outfilename, fortalk=fortalk)
    
    # ------------save figure-------------
    figname = Path(outfilename + "_lsm.pdf")
    save_fig(fig, figname.parent, figname.name, fortalk=fortalk)

    # --------------Multiparameter fit--------------------------------
    fig 	= plt.figure(figsize=(2.4 * fig_size, fig_size))
    ax1	 	= fig.add_axes([0.09,0.15,0.40,0.83])
    ax2	 	= fig.add_axes([0.59,0.15,0.40,0.83])

    popt_d0, perr_d0, popt_r0, perr_r0 = plt_dmpars_fit_multipar(df, xcol, x2col, y1col, ax1, r"$\Delta \log D_0$", None, None, 
                                       y2col, ax2, r"$\Delta \log r_0$", None, None,
                                       yscale_log=False, fit_robust=fit_robust, outfilename=outfilename, fortalk=fortalk)

    # ------------write fit params-------------
    np.savetxt(f'{outfilename}_multifit_params.txt', np.array([popt_d0, perr_d0, popt_r0, perr_r0]), fmt='%.2f    %.2f    %.2f')
    print(f'Saved figures {outfilename}_multifit_params.txt')

    # ------------save figure-------------
    figname = Path(outfilename + "_lsmsfr.pdf")
    save_fig(fig, figname.parent, figname.name, fortalk=fortalk)

    return (0)

#	----------------------------------------------------------------------------------------------------------









