#!/usr/bin/env python3

"""

    Title :      make_lasso_phase_plot.py
    Notes :      Make a phase plot weighted by electron density using datashader for a given FOGGIE snapshot, 
                 and, optionally, interactively lasso regions over it to plot physical projection plots
                 Run this script in astroconda env
    Output :     Plot saved as png file
    Author :     Ayan Acharyya
    Started :    Nov 2026
    Examples :   run make_lasso_phase_plot.py --system ayan_local --halo 2878 --output RD0038 --upto_kpc 100 --use_cen_smoothed --weight edensity --selcol
                 run make_lasso_phase_plot.py --system ayan_local --halo 2878 --output RD0038 --upto_kpc 100 --use_cen_smoothed --weight edensity
"""
from foggie_header import *
from craft_utils import annotate_axes, save_fig
from matplotlib.widgets import LassoSelector
from matplotlib.widgets import SpanSelector

plt.rcParams['axes.linewidth'] = 1

start_time = datetime.now()

# ------------------------------------------------------------------
class SelectFromCollection:
    '''
    Object to interactively choose coordinate points in a given plot using the LassoSelector
    This is based on Raymond's foggie.angular_momentum.lasso_data_selection.ipynb
    '''
    def __init__(self, ax, xys, nbins, x_range, y_range, alpha_other=0.3):
        self.ax = ax
        self.alpha_other = alpha_other
        self.xys = xys

        self.lasso = LassoSelector(ax, lambda verts: self.onselect(verts, nbins), useblit=True)
        self.ind = []
        self.path = None  # store exact lasso polygon for fast cell selection

        mask = np.zeros((nbins, nbins))

        extent = x_range[0], x_range[1], y_range[0], y_range[1]
        text_xpos, text_ypos = x_range[0] + 0.5 * (x_range[1] - x_range[0]), y_range[0] + 0.95 * (y_range[1] - y_range[0])

        self.helpful_text = self.ax.text(text_xpos, text_ypos, '', ha='center', fontsize=args.fontsize)
        self.mask_plot = self.ax.imshow(mask, cmap=plt.cm.Greys_r, alpha=0, interpolation='bilinear', extent=extent, zorder=10, aspect='auto', vmin=0, vmax=1, origin='lower')
        self.blitman = BlitManager(self.ax.figure.canvas, [self.mask_plot, self.helpful_text])

    def onselect(self, verts, nbins):
        path = mpl_Path(verts)
        self.path = path
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        # -------to highlight the selected region-------------------
        if len(self.ind) > 0:
            mask = np.zeros((nbins, nbins))
            mask.ravel()[self.ind] = 1.
            self.mask_plot.set_data(mask)
            self.mask_plot.set_alpha(self.alpha_other)
            self.helpful_text.set_text('Press ENTER to proceed with this selection or select again')
        else:
            self.mask_plot.set_alpha(0.)
            self.helpful_text.set_text('No data selected; select points by dragging mouse in a loop')
        self.blitman.update()

    def disconnect(self):
        self.helpful_text.set_text('')
        self.mask_plot.set_alpha(0.)
        self.blitman.update()
        self.lasso.disconnect_events()

# ------------------------------------------------------------------
class SelectFromSlider:
    '''
    Object to interactively choose x axis range in a given histogram using the RangeSlider
   '''
    def __init__(self, ax, nbins, data_min, data_max, alpha_other=0.3):
        self.ax = ax
        self.alpha_other = alpha_other
        self.data_min = data_min
        self.data_max = data_max

        onselect = lambda min, max: self.update_span(min, max, nbins)
        self.slider = SpanSelector(ax, onselect, 'horizontal', onmove_callback=onselect, useblit=True, span_stays=True, rectprops=dict(alpha=0))
        self.val = [data_min, data_max]

        mask = np.zeros((nbins, nbins))
        extent = self.data_min, self.data_max, self.ax.get_ylim()[0], self.ax.get_ylim()[1]
        self.mask_plot = self.ax.imshow(mask, cmap=plt.cm.Greys_r, alpha=0, interpolation='none', extent=extent, zorder=10, aspect='auto', vmin=0, vmax=1)

        self.lower_limit_line = ax.axvline(self.data_min, color='k')
        self.upper_limit_line = ax.axvline(self.data_max, color='k')
        self.helpful_text = ax.text(ax.get_xlim()[0] + 0.5 * np.diff(ax.get_xlim())[0], ax.get_ylim()[0] + 0.95 * np.diff(ax.get_ylim())[0], 'Drag cursor to select a range in colors', ha='center', fontsize=args.fontsize)

        # ----------setup the Blit manager and key press event---------------------------
        self.blitman = BlitManager(self.ax.figure.canvas, [self.lower_limit_line, self.upper_limit_line, self.helpful_text])

    def update_span(self, min, max, nbins):
        '''
        Function to update the val passed to a callback by the RangeSlider
        '''
        self.val = [min, max]
        # -------to highlight the selected region-------------------
        if self.val[1] - self.val[0] > 0:
            mask = np.zeros((nbins, nbins))
            mask[:, int((self.val[0] - self.data_min) * nbins/(self.data_max - self.data_min)) + 1: int((self.val[1] - self.data_min) * nbins/(self.data_max - self.data_min)) + 1] = 1.
            self.mask_plot.set_data(mask)
            self.mask_plot.set_alpha(self.alpha_other)
            self.helpful_text.set_text('Press ENTER to proceed with this selection or select again')
        else:
            self.mask_plot.set_alpha(0)
            self.helpful_text.set_text('No range selected; select by dragging the cursor')
        # ------update the position of the vertical lines------
        self.lower_limit_line.set_xdata(self.val[0])
        self.upper_limit_line.set_xdata(self.val[1])
        self.blitman.update()

    def disconnect(self):
        self.helpful_text.set_text('')
        self.mask_plot.set_alpha(0)
        self.lower_limit_line.set_xdata(self.data_min)
        self.upper_limit_line.set_xdata(self.data_max)
        self.blitman.update()
        self.slider.disconnect_events()

# -------------------------------------------------------------------------------
def _get_proj_array(ds, field_tuple, axis, box_cut, center, width_kpc, unit, npix=512):
    '''
    Helper function to compute a 2D projected array directly from a yt data_source without saving a plot
    '''
    proj_obj = ds.proj(field_tuple, axis, data_source=box_cut, center=center)
    frb = proj_obj.to_frb((width_kpc, 'kpc'), (npix, npix))
    arr = np.array(frb[field_tuple].in_units(unit).ndarray_view(), dtype=np.float64)
    arr[arr <= 0] = np.nan
    return arr

# -------------------------------------------------------------------------------
def _add_ax_in_axes_frac(fig, host_ax, x0_frac, y0_frac, w_frac, h_frac):
    '''
    Create a plain Axes on `fig` positioned at [x0, y0, width, height]
    expressed as fractions of `host_ax`'s bounding box.
    '''
    p0_fig = fig.transFigure.inverted().transform(host_ax.transAxes.transform((x0_frac, y0_frac)))
    p1_fig = fig.transFigure.inverted().transform(host_ax.transAxes.transform((x0_frac + w_frac, y0_frac + h_frac)))
    return fig.add_axes([p0_fig[0], p0_fig[1], p1_fig[0] - p0_fig[0], p1_fig[1] - p0_fig[1]])

# -------------------------------------------------------------------------------
def _draw_selection_panel(dst_ax, df, args, selector=None, include_cbar=True):
    '''
    Draw the datashader phase plot and highlighted lasso mask directly onto dst_ax
    without saving or reading a PNG.
    '''
    if hasattr(args.colormap, 'colors'):
        color_key = [to_hex(item) for item in args.colormap.colors]
    else:
        color_key = args.colormap

    is_weighted = args.weight is not None and args.weight in df.columns
    if is_weighted:
        color_strategy = dsh.sum(args.col_to_color)
        clabel_prefix = 'Summed'
        vmin_plot, vmax_plot = None, None
    else:
        color_strategy = dsh.mean(args.col_to_color)
        clabel_prefix = 'Mean'
        vmin_plot, vmax_plot = args.cmin, args.cmax

    plot_df = df.copy()
    if islog_dict[args.col_to_color]:
        plot_df[args.col_to_color] = 10 ** plot_df[args.col_to_color]

    artist = dsshow(
        plot_df, dsh.Point(args.xcol, args.ycol), color_strategy,
        norm='log' if islog_dict[args.col_to_color] else 'linear',
        cmap=color_key, aspect='auto', ax=dst_ax,
        x_range=(args.xmin, args.xmax), y_range=(args.ymin, args.ymax),
        vmin=vmin_plot, vmax=vmax_plot,
    )

    xlabel = f'{label_dict[args.xcol]} [{unit_dict[args.xcol]}]'
    ylabel = f'{label_dict[args.ycol]} [{unit_dict[args.ycol]}]'
    if islog_dict[args.xcol]:
        xlabel = 'Log ' + xlabel
    if islog_dict[args.ycol]:
        ylabel = 'Log ' + ylabel

    dst_ax = annotate_axes(dst_ax, xlabel, ylabel, args=args)
    dst_ax.set_yticks([2, 4, 6, 8])

    if include_cbar:
        word_arr = f'{clabel_prefix} {label_dict[args.col_to_color]} [{unit_dict[args.col_to_color]}]'.split(' ')
        newline_at = int(np.ceil(len(word_arr) / 2))
        clabel = f'{" ".join(word_arr[:newline_at])}\n{" ".join(word_arr[newline_at:])}'

        # Set your tuned [x0, y0, width, height] fractions relative to dst_ax
        cax_xpos, cax_ypos, cax_width, cax_height = 0.4, 0.8, 0.5, 0.04
        cax = _add_ax_in_axes_frac(dst_ax.figure, dst_ax, cax_xpos, cax_ypos, cax_width, cax_height)
        cbar = plt.colorbar(artist, cax=cax, orientation='horizontal')
        if islog_dict[args.col_to_color]:
            cbar.ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10.0))
            cbar.ax.tick_params(labelsize=args.fontsize / args.fontfactor)
        else:
            cax.set_xticklabels(['%.0F' % index for index in cax.get_xticks()], fontsize=args.fontsize / args.fontfactor)
        cax.set_title(clabel, fontsize=args.fontsize / args.fontfactor, color='k', pad=2, loc='center')

    # Overlay the lasso highlight mask if available
    if selector is not None and hasattr(selector, 'mask_plot'):
        mask_data = selector.mask_plot.get_array()
        if mask_data is not None and np.any(mask_data > 0):
            dst_ax.imshow(mask_data, cmap=plt.cm.Greys_r, alpha=selector.alpha_other, interpolation='bilinear', extent=(args.xmin, args.xmax, args.ymin, args.ymax),
                zorder=10, aspect='auto', vmin=0, vmax=1, origin='lower')

# ------------------------------------------------------------------
def make_projections(box_cut, args, identifier, df=None, selector=None, range_selector=None, npix_proj=512):
    '''
    Function to directly build the single-row combined figure with selection panel(s),
    three compact line-of-sight density projections sharing y-axis and colorbar, and save once
    '''
    print('Preparing direct single-row combined projection figure..')
    bounds_dict.update(rad=(0, args.galrad))
    field_to_plot = 'density'
    field_tuple = field_dict[field_to_plot]
    proj_unit = projected_unit_dict[field_to_plot]
    zmin, zmax = projected_bounds_dict[field_to_plot][0], projected_bounds_dict[field_to_plot][1]

    proj_arr = ['x', 'y', 'z']
    proj_grids = {}
    for thisproj in proj_arr:
        proj_grids[thisproj] = _get_proj_array(
            box_cut.ds, field_tuple, thisproj, box_cut,
            box.center.in_units(kpc), 2 * args.galrad, npix=npix_proj, unit=proj_unit
        )

    # ------set up two subgrids in a single row: selection(s) on left, projections on right----------
    has_slider = bool(args.selcol) and (range_selector is not None)
    n_sel_panels = 2 if has_slider else 1

    fig = plt.figure(figsize=(3.2 * (n_sel_panels + len(proj_arr)) + 1.6, 3.8), num='combined_lasso_fig')
    outer_gs = gridspec.GridSpec(
        1, 2, figure=fig,
        width_ratios=[n_sel_panels * 1.0, len(proj_arr) * 1.0],
        wspace=0.2,
        left=0.05, right=0.89, top=0.92, bottom=0.18
    )

    sel_gs = gridspec.GridSpecFromSubplotSpec(1, n_sel_panels, subplot_spec=outer_gs[0], wspace=0.22)
    proj_gs = gridspec.GridSpecFromSubplotSpec(1, len(proj_arr), subplot_spec=outer_gs[1], wspace=0.06)

    # ------plot selection panel(s) directly in left subgrid----------
    sel_axes = [fig.add_subplot(sel_gs[0, i]) for i in range(n_sel_panels)]
    if df is not None:
        _draw_selection_panel(sel_axes[0], df, args, selector=selector, include_cbar=not has_slider)

    if has_slider and df is not None:
        color_quant = box_cut[field_dict[args.col_to_color]].in_units(unit_dict[args.col_to_color]).ndarray_view()
        if islog_dict[args.col_to_color]:
            color_quant = np.log10(color_quant)
            cmin = np.log10(bounds_dict[args.col_to_color][0])
            cmax = np.log10(bounds_dict[args.col_to_color][1])
        else:
            cmin, cmax = args.cmin, args.cmax
        Y, X = np.histogram(color_quant, bins=50, density=False)
        cmap_obj = plt.get_cmap(args.colormap) if isinstance(args.colormap, str) else args.colormap
        colors = [cmap_obj((x - cmin) / (cmax - cmin)) for x in X]
        sel_axes[1].bar(X[:-1], Y, color=colors, width=X[1] - X[0])
        sel_axes[1].axvspan(range_selector.val[0], range_selector.val[1], color='grey', alpha=range_selector.alpha_other)
        sel_axes[1].tick_params(axis='both', which='major', labelsize=args.fontsize)

    # ------plot x, y, z projections in right subgrid sharing y-axis and colorbar----------
    cmap_obj = cmap_dict[field_to_plot]
    if hasattr(cmap_obj, 'copy'):
        cmap_obj = cmap_obj.copy()
    try:
        cmap_obj.set_bad('k')
    except Exception:
        pass

    norm_obj = LogNorm(vmin=zmin, vmax=zmax) if islog_dict[field_to_plot] else Normalize(vmin=zmin, vmax=zmax)
    extent = [-args.galrad, args.galrad, -args.galrad, args.galrad]
    proj_ticks = np.array([-0.75, 0, 0.75]) * args.galrad

    first_proj_ax = None
    last_proj_ax = None
    last_im = None
    for idx, thisproj in enumerate(proj_arr):
        if first_proj_ax is None:
            ax_p = fig.add_subplot(proj_gs[0, idx])
            first_proj_ax = ax_p
        else:
            ax_p = fig.add_subplot(proj_gs[0, idx], sharey=first_proj_ax)

        last_im = ax_p.imshow(
            proj_grids[thisproj].T, origin='lower', extent=extent,
            cmap=cmap_obj, norm=norm_obj, aspect='equal'
        )
        ax_p.set_xticks(proj_ticks)
        ax_p.set_yticks(proj_ticks)
        ax_p.text(0.05, 0.93, f'Los = {thisproj}', color='white', transform=ax_p.transAxes,
                  fontsize=args.fontsize, ha='left', va='top')

        ax_p.set_xlabel('Offset (kpc)', fontsize=args.fontsize)
        ax_p.tick_params(axis='both', which='major', labelsize=args.fontsize)

        if idx == 0:
            ax_p.set_ylabel('Offset (kpc)', fontsize=args.fontsize)
        else:
            ax_p.tick_params(axis='y', which='both', labelleft=False)

        last_proj_ax = ax_p

    # ------single shared colorbar matched to projection axis height----------
    proj_pos = last_proj_ax.get_position()
    sel_pos = sel_axes[0].get_position()
    sel_axes[0].set_position([sel_pos.x0, proj_pos.y0, sel_pos.width, proj_pos.height])

    cbar_ax = fig.add_axes([0.905, proj_pos.y0, 0.015, proj_pos.height])
    cbar = fig.colorbar(last_im, cax=cbar_ax)
    cbar_label = f'{label_dict[field_to_plot]}\n[{proj_unit}]'
    if islog_dict[field_to_plot]:
        cbar_label = 'Log ' + cbar_label
    cbar.set_label(cbar_label, fontsize=args.fontsize)
    cbar.ax.tick_params(labelsize=args.fontsize)

    # ------save combined figure once----------
    selcol_text = f'_selcol_{args.col_to_color}' if args.selcol else ''
    outfile_rootname = f'{args.output}_{args.halo}_dsh_combined_projection_xyz{args.upto_text}_from_{args.ycolname}_vs_{args.xcolname}{selcol_text}_lasso{identifier}'
    saveplot(fig, args, outfile_rootname, outputdir=str(args.fig_dir) + '/')
    plt.show(block=False)

    return fig

# ---------------------------------------------------------------------------------
def on_second_key_press(event, box_cut, ax, args, range_selector, identifier, df=None, selector=None):
    if event.key == 'enter':
        start_time_on_press = time.time()
        x_lower_bound = range_selector.val[0]
        x_upper_bound = range_selector.val[1]
        print(f'\n\tSelected bounds {x_lower_bound}-{x_upper_bound}..')
        plt.close(ax.figure) ##

        if x_upper_bound > x_lower_bound:
            print('\tPreparing to cut regions further based on highlighted selection: %s from %.2F to %.2F %s' % (
                args.colorcol, x_lower_bound, x_upper_bound, unit_dict[args.colorcol]
            ))
            if islog_dict[args.colorcol]:
                x_lower_bound, x_upper_bound = 10 ** x_lower_bound, 10 ** x_upper_bound

            cut_crit = "((obj[{}] > obj.ds.quan({}, '{}')) & (obj[{}] < obj.ds.quan({}, '{}')))".format(
                field_dict[args.colorcol], x_lower_bound, unit_dict[args.colorcol],
                field_dict[args.colorcol], x_upper_bound, unit_dict[args.colorcol]
            )
            box_cut = box_cut.cut_region([cut_crit])
            prj = make_projections(box_cut, args, identifier, df=df, selector=selector, range_selector=range_selector)
        else:
            print('\tNo data point selected; therefore, not proceeding any further.')
        print('Second key press events completed in %s' % (timedelta(seconds=time.time() - start_time_on_press)))
        
# ---------------------------------------------------------------------------------
def convert_to_datashader_frame(data, data_min, data_max, npix_datashader, log_scale=False):
    '''
    Function to convert physical quantities to corresponding pixel values in the datashader image frame
    '''
    if log_scale:
        data, data_min, data_max = np.log10(data), np.log10(data_min), np.log10(data_max)
    data_in_pix = (data - data_min) * npix_datashader / (data_max - data_min)
    return data_in_pix

# ---------------------------------------------------------------------------------
def make_coordinate_axis(colname, data_min, data_max, ax, fontsize, npix_datashader=1000, dsh=True, log_scale=False, label=None):
    '''
    Function to make the coordinate axis
    Uses globally defined islog_dict and unit_dict
    '''
    if label is None:
        log_text = 'Log ' if islog_dict[colname] and not log_scale else ''
        ax.set_label_text(log_text + label_dict[colname] + ' (' + unit_dict[colname] + ')', fontsize=fontsize)
    else:
        ax.set_label_text(label, fontsize=fontsize)

    nticks = 50 if log_scale else 5
    ticks = np.linspace(data_min, data_max, nticks)
    if dsh: ax.set_ticks(convert_to_datashader_frame(ticks, data_min, data_max, npix_datashader, log_scale=log_scale))
    else: ax.set_ticks(ticks)

    nticklabels = 6 if log_scale else 5
    ticklabels = np.array(['     '] * nticks)
    if log_scale:
        label_arr = np.linspace(np.log10(data_min), np.log10(data_max), nticklabels + 1)
        dig = np.digitize(np.log10(ticks), label_arr)
        bins = [np.where(dig == i)[0][0] for i in np.unique(dig)]
        for bin in bins: ticklabels[bin] = '%.1F' % ticks[bin]
    else:
        ticklabel_every = int(nticks/nticklabels)
        ticklabels[:: ticklabel_every] = ['%.1F' % item for item in ticks[:: ticklabel_every]]

    ax.set_ticklabels(ticklabels, fontsize=fontsize)

    return ax

# -------------------------------------------------------------------
def make_histogram(box_cut, args, identifier, ncolbins, df=None, selector=None):
    '''
    Function to prepare histogram with interactive range_selector, once the 2D selection is finalised
    '''
    # ----------to plot the histogram of colors---------------------------
    fig, ax = plt.subplots(1, num='histogram_fig')
    color_quant = box_cut[field_dict[args.col_to_color]].in_units(unit_dict[args.col_to_color]).ndarray_view()
    if islog_dict[args.col_to_color]:
        color_quant = np.log10(color_quant)
        cmin = np.log10(bounds_dict[args.col_to_color][0])
        cmax = np.log10(bounds_dict[args.col_to_color][1])
    else:
        cmin, cmax = args.cmin, args.cmax

    Y, X = np.histogram(color_quant, bins=ncolbins, density=False)
    cmap_obj = plt.get_cmap(args.colormap) if isinstance(args.colormap, str) else args.colormap
    colors = [cmap_obj((x - cmin) / (cmax - cmin)) for x in X]

    h = ax.bar(X[:-1], Y, color=colors, width=X[1] - X[0])
    ax.xaxis = make_coordinate_axis(args.col_to_color, cmin, cmax, ax.xaxis, args.fontsize, dsh=False)
    ax.yaxis = make_coordinate_axis('PDF', 0, np.max(Y), ax.yaxis, args.fontsize, dsh=False, label='N cells')

    # ---------setup the range_selector on the histogram------------------
    range_selector = SelectFromSlider(ax, ncolbins, cmin, cmax)
    fig._range_selector = range_selector  # prevent garbage collection
    fig.canvas.mpl_connect("key_press_event", lambda event: on_second_key_press(event, box_cut, ax, args, range_selector, identifier, df=df, selector=selector))
    plt.show(block=False)

    return ax

# ------------------------------------------------------------------
def on_first_key_press(event, box, ax, args, selector, ncolbins=100, df=None):
    '''
    Function to accept a key press event and display the selected pixels and modify the axis title
    This is based on Raymond's foggie.angular_momentum.lasso_data_selection.ipynb
    '''
    if event.key == 'enter':
        start_time_on_press = time.time()
        nselection = len(selector.ind)

        if nselection > 0 and selector.path is not None:
            print(f'\n\tSelected {nselection} binned pixels; saving the selection plot..')
            identifier = random.randint(10, 99)
            plt.close(ax.figure) ##

            # ------to cut out the selected region using a single derived mask field----------
            print('\tPreparing to cut regions based on highlighted selection..')
            lasso_path = selector.path
            verts = lasso_path.vertices
            xmin_sel, ymin_sel = verts.min(axis=0)
            xmax_sel, ymax_sel = verts.max(axis=0)
            mask_field_name = f'lasso_mask_{identifier}'

            mask_cache = {}
            def _lasso_mask(field, data):
                key = (getattr(data, 'id', None), data.shape)
                if key[0] is not None and key in mask_cache:
                    return mask_cache[key]

                xvals = data[field_dict[args.xcol]].in_units(unit_dict[args.xcol]).ndarray_view()
                yvals = data[field_dict[args.ycol]].in_units(unit_dict[args.ycol]).ndarray_view()
                if islog_dict[args.xcol]: xvals = np.log10(xvals)
                if islog_dict[args.ycol]: yvals = np.log10(yvals)

                # fast bounding-box pre-filter before polygon check
                in_bbox = (xvals >= xmin_sel) & (xvals <= xmax_sel) & (yvals >= ymin_sel) & (yvals <= ymax_sel)
                mask = np.zeros(xvals.shape, dtype=np.float64)
                if np.any(in_bbox):
                    pts = np.column_stack((xvals[in_bbox], yvals[in_bbox]))
                    mask[in_bbox] = lasso_path.contains_points(pts).astype(np.float64)

                if key[0] is not None:
                    mask_cache[key] = mask
                return mask

            box.ds.add_field(('gas', mask_field_name), function=_lasso_mask, sampling_type='cell', units='', force_override=True)
            box_cut = box.cut_region([f"obj[('gas', '{mask_field_name}')] > 0.5"])

            if args.selcol:
                print('\tPreparing to make histogram for color selection..')
                ax_hist = make_histogram(box_cut, args, identifier, ncolbins=ncolbins, df=df, selector=selector)  # make histogram for selection in color-space
            else:
                print('\tNo color selection requested, so reparing to make projection plots directly..')
                prj = make_projections(box_cut, args, identifier, df=df, selector=selector)  # OR directly plot projections with all underlying color values
            selector.ind = []
        else:
            print('\tNo data point selected; therefore, not proceeding any further.')

        print('First key press events completed in %s' % (timedelta(seconds=time.time() - start_time_on_press)))

    elif event.key == 'escape':
        print('Exiting interactive mode..')
        selector.disconnect()

    else:
        print('You pressed ' + event.key + '; but you need to press ENTER to see your selection or ESC to exit interactive mode')

# -------------------------------------------------------------------------------------------------------------
def setup_interactive(box, ax, args, nselectionbins=100, ncolbins=50, df=None):
    '''
    Function to interactively select points from a live plot using LassoSelector
    This is based on Raymond's foggie.angular_momentum.lasso_data_selection.ipynb
    '''
    x_range, y_range = [args.xmin, args.xmax], [args.ymin, args.ymax]

    x = np.linspace(x_range[0], x_range[1], nselectionbins)
    y = np.linspace(y_range[0], y_range[1], nselectionbins)
    xgrid, ygrid = np.meshgrid(x, y)
    xys = np.array([xgrid.ravel(),ygrid.ravel()]).T
    selector = SelectFromCollection(ax, xys, nselectionbins, x_range, y_range)

    fig = plt.gcf()
    fig._selector = selector
    fig.canvas.mpl_connect("key_press_event", lambda event: on_first_key_press(event, box, ax, args, selector, ncolbins=ncolbins, df=df))

    selector.helpful_text.set_text('Select points by dragging mouse in a loop')
    fig.canvas.draw()
    selector.blitman.update()

    return fig

# -------------------------------------------------------------------------------
def get_df_from_ds(ds, args):
    '''
    Function to make a pandas dataframe from the yt dataset based on the given field list and color category
    The dataframe can be written to file for faster access in future
    This function is somewhat based on foggie.utils.prep_dataframe.prep_dataframe()
    :return: dataframe
    '''
    print('\nCreating datframe..')
    all_fields = [args.xcol, args.ycol, args.col_to_color]
    if args.weight is not None:
        all_fields += [args.weight]
    df = pd.DataFrame()
    for index,field in enumerate(all_fields):
        print('Doing property: ' + field + ', which is ' + str(index + 1) + ' of the ' + str(len(all_fields)) + ' fields..')
        arr = ds[field_dict[field]].in_units(unit_dict[field]).ndarray_view()
        df[field] = arr

    return df

# --------------------------------------------------------------------------------
def make_datashader_plot_mpl(df_input, args):
    '''
    Function to make data shader plot of y_field vs x_field, colored in bins of color_field
    This function is based on foggie.render.shade_maps.render_image()
    :return figure
    '''
    print(f'Making datashader plot..')
    df = df_input.copy()
    
    # --------to make the main datashader plot--------------------------
    fig, ax = plt.subplots(1, figsize=(7.5, 7), constrained_layout=True, num='datashader_fig')
    
    if hasattr(args.colormap, 'colors'): color_key = [to_hex(item) for item in args.colormap.colors]
    else: color_key = args.colormap
    
    is_weighted = args.weight is not None and args.weight in df.columns

    if is_weighted:
        color_strategy = dsh.sum(args.col_to_color)
        clabel_prefix = 'Summed'
        vmin_plot, vmax_plot = None, None
    else:
        color_strategy = dsh.mean(args.col_to_color)
        clabel_prefix = 'Mean'
        vmin_plot, vmax_plot = args.cmin, args.cmax

    if islog_dict[args.col_to_color]: df[args.col_to_color] = 10 ** df[args.col_to_color] # converting back to linear because they were already in log
    artist = dsshow(df, dsh.Point(args.xcol, args.ycol), color_strategy, 
                    norm='log' if islog_dict[args.col_to_color] else 'linear', cmap=color_key, aspect = 'auto', ax=ax,
                    x_range=(args.xmin, args.xmax), y_range=(args.ymin, args.ymax),
                    vmin=vmin_plot, vmax=vmax_plot,
                )

    # ------to make the axes-------------
    xlabel = f'{label_dict[args.xcol]} [{unit_dict[args.xcol]}]'
    ylabel = f'{label_dict[args.ycol]} [{unit_dict[args.ycol]}]'
    word_arr = f'{clabel_prefix} {label_dict[args.col_to_color]} [{unit_dict[args.col_to_color]}]'.split(' ')
    newline_at = int(np.ceil(len(word_arr) / 2))
    clabel = f'{" ".join(word_arr[:newline_at])}\n{" ".join(word_arr[newline_at:])}'
    if islog_dict[args.xcol]: xlabel = 'Log ' + xlabel
    if islog_dict[args.ycol]: ylabel = 'Log ' + ylabel

    ax = annotate_axes(ax, xlabel, ylabel, args=args)
    ax.set_yticks([2, 4, 6, 8])

    # ------to make the colorbar axis-------------
    cax_xpos, cax_ypos, cax_width, cax_height = 0.6, 0.8, 0.3, 0.035
    cax = fig.add_axes([cax_xpos, cax_ypos, cax_width, cax_height])
    cbar = plt.colorbar(artist, cax=cax, orientation='horizontal')
    if islog_dict[args.col_to_color]:
        cbar.ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10.0))
        cbar.ax.tick_params(labelsize=args.fontsize)
    else:
        cax.set_xticklabels(['%.0F' % index for index in cax.get_xticks()], fontsize=args.fontsize)
    fig.text(cax_xpos + cax_width / 2, cax_ypos + cax_height + 0.005, clabel, ha='center', va='bottom', fontsize=args.fontsize)

    # ---------to annotate and save the figure----------------------
    plt.show(block=False)
    return fig

# ----------to initialise dicts-----------------------
field_dict = {'rad': ('gas', 'radius_corrected'), 'metal': ('gas', 'metallicity'), 'vrad': ('gas', 'radial_velocity_corrected'), 'density':('gas', 'density'), 'temp':('gas', 'temperature'), 'edensity':('gas', 'El_number_density')}
unit_dict = {'rad': 'kpc', 'metal': r'Zsun', 'vrad': 'km/s', 'density':'g/cm**3', 'temp':'K', 'edensity':'cm**-3'}
label_dict = {'rad': 'Radius (kpc)', 'metal': 'Log Metallicity (Zsun)', 'vrad': 'Radial velocity (km/s)', 'density':'Gas density', 'temp':'Temperature', 'edensity':'Electron number density'}
bounds_dict = {'metal': (-2, 1), 'vrad': (-400, 400), 'density': (1e-30, 1e-22), 'temp': (1e1, 2e8), 'edensity': (1e-6, 2)}
colormap_dict = {'vrad': outflow_inflow_discrete_cmap, 'temp': temperature_discrete_cmap, 'density': density_discrete_cmap, 'metal': metal_discrete_cmap, 'edensity': 'terrain'}#edensity_discrete_cmap}
islog_dict = {'rad': False, 'metal': True, 'vrad': False, 'density': True, 'temp': True, 'edensity': True}

projected_unit_dict = defaultdict(lambda x: unit_dict[x], density='Msun/pc**2')
cmap_dict = {'density':density_color_map, 'metal':metal_color_map, 'temp':temperature_color_map, 'vrad':velocity_discrete_cmap} # for projection plots, if any
projected_bounds_dict = {'density': (5e-2, 1e3)}

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.system == "ayan_pleiades": args.code_dir = '/nobackupp19/aachary2/ayan_codes/foggie/foggie/'
    else: args.code_dir = '/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/'
    args.scale_factor = 2.4 # tune to desired readability when read into combined plot

    print(f'Lassoing snapshot {args.halo}:{args.output}..')

    # --------determining corresponding text suffixes and figname-------------
    args.fig_dir = Path(args.output_dir + 'plots/phase_plots_lasso')
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    args.fits_dir = Path(args.output_dir + 'data/')
    args.fits_dir.mkdir(parents=True, exist_ok=True)

    args.upto_text = '_upto%.1Fckpchinv' % args.upto_kpc if args.docomoving else '_upto%.1Fkpc' % args.upto_kpc

    # -------loading in snapshot-------------------
    halos_df_name = args.code_path + 'halo_infos/00' + args.halo + '/' + args.run + '/'
    halos_df_name += 'halo_cen_smoothed' if args.use_cen_smoothed else 'halo_c_v'
    ds, refine_box = load_sim(args, region='refine_box', do_filter_particles=False, disk_relative=False, halo_c_v_name=halos_df_name)

    # --------assigning additional keyword args-------------
    args.current_redshift = ds.current_redshift
    args.current_time = ds.current_time.in_units('Gyr').tolist()

    # ------tailoring the simulation box for individual snapshot analysis--------
    if args.docomoving: args.galrad = args.upto_kpc / (1 + args.current_redshift) / 0.695  # fit within a fixed comoving kpc h^-1, 0.695 is Hubble constant
    else: args.galrad = args.upto_kpc  # fit within a fixed physical kpc

    # ----------to initialise axes limits and colormaps-----------------------
    args.xcolname = 'log_' + args.xcol if islog_dict[args.xcol] else args.xcol
    args.ycolname = 'log_' + args.ycol if islog_dict[args.ycol] else args.ycol
    if args.weight is None: args.col_to_color = args.colorcol
    else: args.col_to_color = args.weight
    args.outfilename = Path(args.output_dir) / 'data' / f'{args.output}_{args.halo}_list_{args.xcol}_{args.ycol}_{args.col_to_color}{args.upto_text}.fits'

    bounds_dict.update(rad=(0, np.ceil(args.galrad * (1 + args.current_redshift)) if args.docomoving else args.galrad)) # if comoving is used, put plot limits to maximum radius, i.e. physical radius at z = 0
    args.xmin = np.log10(bounds_dict[args.xcol][0]) if islog_dict[args.xcol] else bounds_dict[args.xcol][0]
    args.xmax = np.log10(bounds_dict[args.xcol][1]) if islog_dict[args.xcol] else bounds_dict[args.xcol][1]
    args.ymin = np.log10(bounds_dict[args.ycol][0]) if islog_dict[args.ycol] else bounds_dict[args.ycol][0]
    args.ymax = np.log10(bounds_dict[args.ycol][1]) if islog_dict[args.ycol] else bounds_dict[args.ycol][1]
    args.cmin = bounds_dict[args.col_to_color][0]
    args.cmax = bounds_dict[args.col_to_color][1]
    args.colormap = colormap_dict[args.col_to_color]

    # -------extract the required box------------
    box_width = args.galrad * 2  # in kpc
    box_width_kpc = ds.arr(box_width, 'kpc')
    box_center = ds.halo_center_kpc
    box = ds.r[box_center[0] - box_width_kpc / 2.: box_center[0] + box_width_kpc / 2., 
            box_center[1] - box_width_kpc / 2.: box_center[1] + box_width_kpc / 2., 
            box_center[2] - box_width_kpc / 2.: box_center[2] + box_width_kpc / 2., ]
    
    # ----------to make dataframe from dataset-----------------------
    if not os.path.exists(args.outfilename) or args.clobber:
        df = get_df_from_ds(box, args)
        Table.from_pandas(df).write(args.outfilename, format='fits', overwrite=True)
        print(f'Saved dataframe in {args.outfilename}')
    else:
        print(f'Reading from existing {args.outfilename}')
        df = Table.read(args.outfilename, format='fits', memmap=True).to_pandas()

    # -----------prepping dataframe--------------
    print(f'Prepping dataframe..')
    for col in df.columns:
        if islog_dict[col]:
            df[col] = np.log10(df[col])

    args.col_text = 'weightby' if args.weight is not None and args.weight in df.columns else 'colorby'

    # ----------to make datashader plot-----------------------
    fig = make_datashader_plot_mpl(df, args)

    # --------to make interactive lasso------------
    print(f'\nThis plot is now in interactive mode..')
    fig = setup_interactive(box, fig.axes[0], args, nselectionbins=100, ncolbins=50, df=df)

    # -----------------------------------------------------------------------------------
    print('\nCompleted in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
    plt.show(block=False)
