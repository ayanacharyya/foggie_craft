##!/usr/bin/env python3

"""

    Title :      craft_utils
    Notes :      Contains various generic utility functions and classes used by the other scripts, including a 'master' function to parse args
    Author :     Ayan Acharyya
    Started :    31-03-26

"""
from craft_header import *
from globalpars import *

###################### Following routines are from auxfns.py ######################################
#	------------------------------------------------------------------------------------------------------
def neprofile (necube, dkpc, theta, phi, cenpx, radius):
	#	Calculate radial electron density profiles along theta, phi	
	
	nepro	=	np.zeros(len(radius), dtype=float)
	for i in range(0, len(radius)):
		nepro[i]= necube[ int(cenpx[0]+radius[i] * np.cos(theta) * np.cos(phi)), int(cenpx[1]+radius[i] * np.cos(theta) * np.sin(phi)), int(cenpx[2] + radius[i] * np.sin(theta))]
	
	return (nepro)

#	------------------------------------------------------------------------------------------------------
def inclinvec (vecomps, theta0, phi0):
#   Return inclination (in degrees) of a vector w.r.t. z axis
#	Arguments:	A 3-component vector, orientation of the galaxy

	dotpabs	= np.abs(vecomps[0] * np.sin(theta0) * np.cos(phi0) + vecomps[1] * np.sin(theta0) * np.sin(phi0) + vecomps[2] * np.cos(theta0))
	
	inclin  = np.rad2deg(np.arccos(dotpabs / np.sqrt(np.sum(vecomps**2))))

	return (inclin)

#	------------------------------------------------------------------------------------------------------
def impactfac (pt1, pt2):
#   Return the impact factor for a given LoS in units of "pixels"
#	Arguments:	Two fixed points on the LoS
    
    crosq	= (pt1[0]*pt2[1]-pt1[1]*pt2[0])**2 +(pt1[0]*pt2[2]-pt1[2]*pt2[0])**2 + (pt1[2]*pt2[1]-pt1[1]*pt2[2])**2
    modsq	= (pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2 + (pt1[2]-pt2[2])**2
    
    impactf = np.sqrt(crosq / modsq) 
    
    return (impactf)

#	------------------------------------------------------------------------------------------------------
def logradialexp (x, x0, a0):
#   Return a radial exponential
#	Arguments:	Radius, scale radius, normalization

	radexp	= np.log10(a0 * np.exp(-x/x0))

	return (radexp)

#	------------------------------------------------------------------------------------------------------
def logradialgaus (x, x0, a0):
#   Return a radial exponential
#	Arguments:	Radius, scale radius, normalization

	radexp	= np.log10(a0 * np.exp(-(x/x0)**2))

	return (radexp)

#	------------------------------------------------------------------------------------------------------
def logradialexp3 (x, x0, a0):
#   Return a radial exponential
#	Arguments:	Radius, scale radius, normalization

	radexp	= np.log10(a0 * np.exp(-(x/x0)**0.33))

	return (radexp)

#	------------------------------------------------------------------------------------------------------
def logbeselk (x, x0, a0):
#   Return a radial Bessel k
#	Arguments:	Radius, scale radius, normalization

	radexp	= np.log10(a0 * special.kn(0.5,x/x0))

	return (radexp)

#	------------------------------------------------------------------------------------------------------
def radialexp (x, x0, a0):
#   Return a radial exponential
#	Arguments:	Radius, scale radius, normalization

	radexp	= a0 * np.exp(-x/x0)

	return (radexp)

#	------------------------------------------------------------------------------------------------------
def radialexpower (x, x0, a0, power):
#   Return a radial powered exponential
#	Arguments:	Radius, normalization, exponent

	radexp	= np.log(a0) - ((x/x0)**power)

	return (radexp)

#	------------------------------------------------------------------------------------------------------
def distfrmajorax (theta0, phi0, pt1, pt2):
#   Return the projected distance of the LoS from the apparent major axis
#	Arguments:	LoS vector, Normal vector
    
	normvec		= np.array([ np.sin(theta0) * np.cos(phi0), np.sin(theta0) * np.sin(phi0), np.cos(theta0)])

	losvec		= pt2 - pt1

	majaxpt		= np.array([0.0, 0.0, 0.0])
	majaxvec	= np.cross(losvec, normvec)

	crpdct		= np.cross(majaxvec, losvec)
	mindist		= np.abs( np.dot(crpdct, (pt1 - majaxpt)) ) / np.sqrt(np.dot(crpdct, crpdct))
	
	return (mindist)

#	------------------------------------------------------------------------------------------------------
def linearxy (xy, a, b, c):
#	Retruns a linear function of two independent variables
# 	z	=	ax + by + c	

	z	= a*xy[0] + b*xy[1] + c 

	return (z)


###################### Following routines are from foggie_utils.py ######################################
# ----------------------------------------------------------------
def setup_plot_style():
    '''
    Function to set default style for all plots made in this project
    '''
    plt.rcParams['pdf.fonttype']	= 42
    plt.rcParams['ps.fonttype'] 	= 42
    plt.rcParams['savefig.dpi'] 	= 600
    plt.rcParams['font.family'] 	= 'sans-serif'
    plt.rcParams['font.size']		= 8
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['ytick.right'] = True
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True

# ----------------------------------------------------------------
def findReplace(directory, find, replace, filePattern):
    '''
    Function to replace string A with string B in ALL files with filenames of a certain pattern, within a given directory
    Modifies all matchinf files in place, therefore returns nothing
    Example use: findReplace("/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/gas_metallicity/", "from header impot *\nfrom util import *", "from header import *
from util import *", "*.py")
    Borrowed from: https://stackoverflow.com/questions/4205854/recursively-find-and-replace-string-in-text-files
    '''
    import fnmatch
    for path, dirs, files in os.walk(os.path.abspath(directory)):
        for filename in fnmatch.filter(files, filePattern):
            filepath = os.path.join(path, filename)
            with open(filepath) as f:
                s = f.read()
            s = s.replace(find, replace)
            with open(filepath, "w") as f:
                f.write(s)

# -----------------------------------------------------------------
def make_its_own_figure(ax, label, args):
    '''
    Function to take an already filled axis handle and turn it into its stand-alone figure
    Output: saved png
    '''
    fig, thisax = plt.subplots(figsize=(8, 8))
    fig.subplots_adjust(right=0.97, top=0.95, bottom=0.05, left=0.07)
    thisax = ax
    outfile_rootname = '%s_%s_%s%s%s%s%s.png' % (label, args.output, args.halo, args.Zgrad_den, args.upto_text, args.weightby_text, args.res_text)
    if args.do_all_sims: outfile_rootname = 'z=*_' + outfile_rootname[len(args.output) + 1:]
    figname = args.fig_dir + outfile_rootname.replace('*', '%.5F' % (args.current_redshift))
    fig.savefig(figname)
    myprint('Saved plot as ' + figname, args)

    return

# ---------------------------------------------------------------------------
def get_kpc_from_arc_at_redshift(arcseconds, redshift):
    '''
    Function to convert arcseconds on sky to physical kpc, at a given redshift
    '''
    cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
    d_A = cosmo.angular_diameter_distance(z=redshift)
    kpc = (d_A * arcseconds * u.arcsec).to(u.kpc, u.dimensionless_angles()).value # in kpc
    print('Converted resolution of %.2f arcseconds to %.2F kpc at target redshift of %.2f' %(arcseconds, kpc, redshift))
    return kpc

# -------------------------------------------------------------------------------------------
def myprint(text, args):
    '''
    Function to re-direct to print_mpi(), for backwards compatibility
    '''
    print_mpi(text, args)

# -------------------------------------------------------------------------------------------
def print_mpi(string, args):
    '''
    Function to print corresponding to each mpi thread
    '''
    comm = MPI.COMM_WORLD
    myprint_orig('[' + str(comm.rank) + '] {' + subprocess.check_output(['uname -n'],shell=True)[:-1].decode("utf-8") + '} ' + string + '\n', args)

# -------------------------------------------------------------------------------------------
def print_master(string, args):
    '''
    Function to print only if on the head node/thread
    '''
    comm = MPI.COMM_WORLD
    if comm.rank == 0: myprint_orig('[' + str(comm.rank) + '] ' + string + '\n', args)

# --------------------------------------------------------------------------------------------
def fix_time_format(text, keyword):
    '''
     Function to modify the way time is formatted in print statements
    '''
    arr = text.split(' ' + keyword)
    pre_time = ' '.join(arr[0].split(' ')[:-1])
    this_time = float(arr[0].split(' ')[-1])
    post_time = ' '.join(arr[1].split(' '))
    text = pre_time + ' %s' % (datetime.timedelta(minutes=this_time)) + post_time

    return text

# ------------------------------------------------------------------------
def insert_line_in_file(line, pos, filename, output=None):
    '''
    Function for nserting a line in a file
    '''
    f = open(filename, 'r')
    contents = f.readlines()
    f.close()

    if pos == -1: pos = len(contents)  # to append to end of file
    contents.insert(pos, line)

    if output is None: output = filename
    f = open(output, 'w')
    contents = ''.join(contents)
    f.write(contents)
    f.close()
    return

# -------------------------------------------------------------------------
def isfloat(str):
    '''
    Function to check if input is float
    '''

    try:
        float(str)
    except ValueError:
        return False
    return True

# -------------------------------------------------------------------------------
def num(s):
    '''
    Function to check if input is a number
    '''

    if s[-1].isdigit():
        return str(format(float(s), '0.2e'))
    else:
        return str(format(float(s[:-1]), '0.2e'))

# -----------------------------------------------------------------
def rebin(array, dimensions=None, scale=None):
    """ Return the array ``array`` to the new ``dimensions`` conserving flux the flux in the bins
    The sum of the array will remain the same

    >>> ar = numpy.array([
        [0,1,2],
        [1,2,3],
        [2,3,4]
        ])
    >>> rebin(ar, (2,2))
    array([
        [1.5, 4.5]
        [4.5, 7.5]
        ])
    Raises
    ------

    AssertionError
        If the totals of the input and result array don't agree, raise an error because computation may have gone wrong

    Reference
    =========
    +-+-+-+
    |1|2|3|
    +-+-+-+
    |4|5|6|
    +-+-+-+
    |7|8|9|
    +-+-+-+
    """
    if dimensions is not None:
        if isinstance(dimensions, float):
            dimensions = [int(dimensions)] * len(array.shape)
        elif isinstance(dimensions, int):
            dimensions = [dimensions] * len(array.shape)
        elif len(dimensions) != len(array.shape):
            raise RuntimeError('')
    elif scale is not None:
        if isinstance(scale, float) or isinstance(scale, int):
            dimensions = map(int, map(round, map(lambda x: x * scale, array.shape)))
        elif len(scale) != len(array.shape):
            raise RuntimeError('')
    else:
        raise RuntimeError('Incorrect parameters to rebin.\n\trebin(array, dimensions=(x,y))\n\trebin(array, scale=a')
    if np.shape(array) == dimensions: return array  # no rebinning actually needed
    import itertools
    # dY, dX = map(divmod, map(float, array.shape), dimensions)

    result = np.zeros(dimensions)
    for j, i in itertools.product(*map(range, array.shape)):
        (J, dj), (I, di) = divmod(j * dimensions[0], array.shape[0]), divmod(i * dimensions[1], array.shape[1])
        (J1, dj1), (I1, di1) = divmod(j + 1, array.shape[0] / float(dimensions[0])), divmod(i + 1,
                                                                                            array.shape[1] / float(
                                                                                                dimensions[1]))

        # Moving to new bin
        # Is this a discrete bin?
        dx, dy = 0, 0
        if (I1 - I == 0) | ((I1 - I == 1) & (di1 == 0)):
            dx = 1
        else:
            dx = 1 - di1
        if (J1 - J == 0) | ((J1 - J == 1) & (dj1 == 0)):
            dy = 1
        else:
            dy = 1 - dj1
        # Prevent it from allocating outide the array
        I_ = np.min([dimensions[1] - 1, I + 1])
        J_ = np.min([dimensions[0] - 1, J + 1])
        result[J, I] += array[j, i] * dx * dy
        result[J_, I] += array[j, i] * (1 - dy) * dx
        result[J, I_] += array[j, i] * dy * (1 - dx)
        result[J_, I_] += array[j, i] * (1 - dx) * (1 - dy)
    allowError = 0.1
    if array.sum() > 0: assert (array.sum() < result.sum() * (1 + allowError)) & (array.sum() > result.sum() * (1 - allowError))
    return result

# ---------------------------------------------------------------------------------------------
def pull_halo_redshift(args):
    '''
    Function to pull the current redshift of the halo (WITHOUT having to load the simulation), by matching with the corresponding halo_c_v file
    '''
    halo_cat_file = args.code_path + 'halo_infos/00' + args.halo + '/nref11c_nref9f/halo_c_v'
    #df = pd.read_csv(halo_cat_file, comment='#', sep='\s+|')
    #try: z = df.loc[df['name']==args.output, 'redshift'].values[0]
    # shifted from pandas to astropy because pandas runs in to weird error on pleiades
    df = Table.read(halo_cat_file, format='ascii')
    try: z = float(df['col2'][np.where(df['col3']==args.output)[0][0]])
    except IndexError: # if this snapshot is not yet there in halo_c_v file
        if args.halo == '4123' and args.output == 'RD0038': z = 0.15
        else: z = -99
    return z

# ------------------------------------------------------------------------
def setup_plots_for_talks():
    '''
    Function to setup plto themes etc for talks
    '''
    plt.style.use('cyberpunk')
    background_for_talks = 'cyberpunk'  # 'dark_background' #'Solarize_Light2' #
    plt.style.use(background_for_talks)
    new_foreground_color = '#FFF1D0'
    #plt.rcParams['grid.color'] = new_foreground_color
    plt.rcParams['text.color'] = new_foreground_color
    plt.rcParams['xtick.color'] = new_foreground_color
    plt.rcParams['ytick.color'] = new_foreground_color
    plt.rcParams['xtick.color'] = new_foreground_color
    plt.rcParams['axes.titlecolor'] = new_foreground_color
    plt.rcParams['axes.labelcolor'] = new_foreground_color
    plt.rcParams['axes.edgecolor'] = new_foreground_color
    plt.rcParams['figure.edgecolor'] = new_foreground_color
    plt.rcParams['savefig.edgecolor'] = new_foreground_color
    plt.rcParams['axes.linewidth'] = 2

    new_background_color = '#120000'
    plt.rcParams['axes.facecolor'] = new_background_color
    plt.rcParams['figure.facecolor'] = new_background_color
    plt.rcParams['savefig.facecolor'] = new_background_color
    plt.rcParams['grid.alpha'] = 0.5
    plt.rcParams['grid.linewidth'] = 0.3

# --------------------------------------------------------------------------------------------------------------------
def annotate_axes(ax, xlabel, ylabel, args=None, fontsize=10, fontfactor=1, label='', clabel='', hide_xaxis=False, hide_yaxis=False, hide_cbar=True, p=None, hide_cbar_ticks=False, cticks_integer=True, label_color='k', bbox=True, set_ticks=True):
    '''
    Annotates the axis of a given 2D image
    Returns the axis handle
    '''
    if args is not None: fontsize, fontfactor = args.fontsize, args.fontfactor
    ax.text(0.05, 0.9, label, c=label_color, fontsize=fontsize/fontfactor, ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9) if bbox else None, transform=ax.transAxes)

    if set_ticks:
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3, prune='both'))
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
         
    if hide_xaxis:
        ax.tick_params(axis='x', which='major', labelsize=fontsize, labelbottom=False)
    else:
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.tick_params(axis='x', which='major', labelsize=fontsize, labelbottom=True)

    if set_ticks:
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, prune='both'))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
         
    if hide_yaxis:
        ax.tick_params(axis='y', which='major', labelsize=fontsize, labelleft=False)
    else:
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.tick_params(axis='y', which='major', labelsize=fontsize, labelleft=True)

    if not hide_cbar and p is not None:
        cax = inset_axes(ax, width="5%", height="100%", loc='right', bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        cbar.set_label(clabel, fontsize=fontsize)

        cbar.locator = ticker.MaxNLocator(integer=cticks_integer, nbins=4)#, prune='both')
        cbar.update_ticks()
        if hide_cbar_ticks:
            cbar.ax.set_yticklabels([])
        else:
            cbar.ax.tick_params(labelsize=fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def save_fig(fig, fig_dir, figname, args, silent=False):
    '''
    Saves a given figure handle as a given output filename
    '''

    if args.fortalk:
        #mplcyberpunk.add_glow_effects()
        #try: mplcyberpunk.make_lines_glow()
        #except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig_dir.mkdir(exist_ok=True, parents=True)
    figname = fig_dir / figname
    fig.savefig(figname, transparent=args.fortalk)
    if not silent: print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    return

# --------------------------------------------------------------------------------------------------------------------
def make_colorbar_top(fig, axes, clabel, cmap, cmin, cmax, ncbins, fontsize, aspect=60):
    '''
    Creates a shared colorbar for the whole figure, at the top of the figure
    Returns figure handle
    '''
    norm = mplcolors.Normalize(vmin=cmin, vmax=cmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    cbar = fig.colorbar(sm, ax=axes, location='top', shrink=0.95, pad=0.01, aspect=aspect)
    cbar.set_label(clabel, fontsize=fontsize, labelpad=5)    
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.locator = ticker.MaxNLocator(integer=False, nbins=ncbins)#, prune='both')
    cbar.update_ticks()

    return fig

# --------------------------------------------------------------------------------------------------------------
def get_grid_size(n_total):
    '''
    Function to determine an appropriate number of nrows and ncols, given a total number of subplots
    Returns nrows, ncols
    '''
    # Start with the square root
    nrows = math.ceil(math.sqrt(n_total))
    ncols = math.ceil(n_total / nrows)
    
    return nrows, ncols

# --------------------------------------------------------------------------------------------------------------
def parse_args():
    '''
    Function to parse keyword arguments
    '''

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='''FRB DM analysis with FOGGIE simulations''')
    # ---- common args used widely over the full codebase ------------
    parser.add_argument('--system', metavar='system', type=str, action='store', default='ayan_local', help='Which system are you on? Default is Jase')
    parser.add_argument('--run', metavar='run', type=str, action='store', default='nref11c_nref9f', help='which run? default is natural')
    parser.add_argument('--pwd', dest='pwd', action='store_true', default=False, help='Just use the current working directory?, default is no')
    parser.add_argument('--halo', metavar='halo', type=str, action='store', default='8508', help='which halo? default is 8508 (Tempest)')
    parser.add_argument('--projection', metavar='projection', type=str, action='store', default='x', help='Which projection do you want to plot, i.e., which axis is your line of sight? Default is x')
    parser.add_argument('--do_all_halos', dest='do_all_halos', action='store_true', default=False, help='loop over all available halos (and all snapshots each halo has)?, default is no')
    parser.add_argument('--silent', dest='silent', action='store_true', default=False, help='Suppress all print statements?, default is no')
    parser.add_argument('--clobber', dest='clobber', action='store_true', default=False, help='overwrite existing outputs with same name?, default is no')
    parser.add_argument('--fontsize', metavar='fontsize', type=float, action='store', default=8, help='fontsize of plot labels, etc.; default is 20')
    parser.add_argument('--fontfactor', metavar='fontfactor', type=float, action='store', default=1.0, help='fontsize of plot labels, etc.; default is 1.2')
    parser.add_argument('--hide', dest='hide', action='store_true', default=False, help='Hide plots from screen, i.e. close them after creation?, default is no')
    parser.add_argument('--keep', dest='keep', action='store_true', default=False, help='keep previously displayed plots on screen?, default is no')
    parser.add_argument('--debug', dest='debug', action='store_true', default=False, help='run in debug mode (lots of print checks)?, default is no')
    parser.add_argument('--docomoving', dest='docomoving', action='store_true', default=False, help='consider the input upto_kpc as a comoving quantity?, default is no')

    # ------- args added for control on plot limits ------------------------------
    parser.add_argument('--xmin', metavar='xmin', type=float, action='store', default=None, help='minimum xaxis limit; default is None')
    parser.add_argument('--xmax', metavar='xmax', type=float, action='store', default=None, help='maximum xaxis limit; default is None')
    parser.add_argument('--ymin', metavar='ymin', type=float, action='store', default=None, help='minimum yaxis limit; default is None')
    parser.add_argument('--ymax', metavar='ymax', type=float, action='store', default=None, help='maximum yaxis limit; default is None')
    parser.add_argument('--ncolbins', metavar='ncolbins', type=int, action='store', default=None, help='number of bins in color space the data shader categories would be split across; default is None')
    parser.add_argument('--upto_kpc', metavar='upto_kpc', type=float, action='store', default=None, help='fit metallicity gradient out to what absolute kpc? default is None')
    parser.add_argument('--fortalk', dest='fortalk', action='store_true', default=False, help='Set plot labels, transparency etc for being used in a talk?, default is no')
    parser.add_argument('--forpaper', dest='forpaper', action='store_true', default=False, help='Set plot labels, transparency etc for being used in the paper?, default is no')
    parser.add_argument('--multi_panel', dest='multi_panel', action='store_true', default=False, help='Make multi-panel figure such that all subplots are in one figure?, default is no')

    # ------- args added for dmplot.py ------------------------------
    parser.add_argument('--mode', metavar='mode', type=str, action='store', default='lsmzsfr', help='which mode to run dmplot.py for? default is lsmzsfr')
    parser.add_argument('--rangekpc', metavar='rangekpc', type=float, action='store', default=200, help='Range (extent) in kpc; default is 200')
    parser.add_argument('--reskpc', metavar='reskpc', type=float, action='store', default=0.5, help='Resolution (cell size) in kpc; default is 0.5')
    parser.add_argument('--resfile_prefix', metavar='resfile_prefix', type=str, action='store', default='all_lsm', help='where to save the resulting data? default is defined later')

    # ------- args added for radialplot.py ------------------------------
    parser.add_argument('--quant', metavar='quant', type=str, action='store', default='electron', help='which quantity to make radial profile of (choose from electron or gas)? default is electron')

    # ------- args added for plots_for_frb_paper.py ------------------------------
    parser.add_argument('--inc_bin_edges', metavar='inc_bin_edges', type=str, action='store', default='0,90', help='Bin edges of inclination angles; Default is 0-90')
    parser.add_argument('--lsm_bin_edges', metavar='lsm_bin_edges', type=str, action='store', default='all', help='Bin edges of log stellar mass; Default is all mass bins')
    parser.add_argument('--lsfr_bin_edges', metavar='lsfr_bin_edges', type=str, action='store', default='-10,10', help='Bin edges of log SFR; Default is -10-10 i.e., all SFRs')
    parser.add_argument('--z_range', metavar='z_range', type=str, action='store', default='0,6', help='Range of redshifts; Default is 0-6')
    parser.add_argument('--set_ylin', dest='set_ylin', action='store_true', default=False, help='Set y-axis scale to linear?, default is no')
    parser.add_argument('--cmap', metavar='cmap', type=str, action='store', default='tab20', help='colormap to use; default is None')
    parser.add_argument('--fit_robust', dest='fit_robust', action='store_true', default=False, help='Fit in the robust-fit method?, default is no')

    parser.add_argument('--plot_dm_lsm', dest='plot_dm_lsm', action='store_true', default=False, help='Plot DM vs Impact factor for a given mass range? Default is no.')
    parser.add_argument('--plot_dm_all_lsm', dest='plot_dm_all_lsm', action='store_true', default=False, help='Plot DM vs Impact factor for ALL mass ranges? Default is no.')
    parser.add_argument('--plot_dm_fit', dest='plot_dm_fit', action='store_true', default=False, help='Plot DM0 and r0 vs log stellar mass? Default is no.')
    parser.add_argument('--make_latex_table', dest='make_latex_table', action='store_true', default=False, help='Convert the df_dmpars to a a latex table? Default is no.')

    # ------- wrap up and processing args ------------------------------
    args = parser.parse_args()
    
    # ----------------halo properties-------------
    args.halo_arr = [item for item in args.halo.split(',')]
    if len(args.halo_arr) == 1: args.halo = args.halo_arr[0]
    args.foggie_dir, args.output_dir, args.run_loc, args.code_path, args.trackname, args.haloname, args.spectra_dir, args.infofile = get_run_loc_etc(args)
    
    args.rangekpc = float(args.rangekpc)
    args.reskpc = float(args.reskpc)

    # --------------paths from globalpars.py-----------------
    docomoving_text = '_comoving' if args.docomoving else ''
    args.root_dir = Path(root_dir)
    args.los_dir = Path(f'{losdir}{docomoving_text}')
    args.data_dir = Path(f'{datadir}{docomoving_text}')
    args.plot_dir = Path(f'{plotdir}{docomoving_text}')
    
    args.fig_dir = args.plot_dir / 'plots_for_paper'
    args.fig_dir.mkdir(exist_ok=True, parents=True)

    args.resfile_prefix = str(args.data_dir / args.resfile_prefix)

    # ---------------parameter ranges---------------------
    args.inc_bin_edges = [float(item) for item in args.inc_bin_edges.split(',')]
    args.inc_bins = list(zip(args.inc_bin_edges, args.inc_bin_edges[1:]))
    
    args.z_range = [float(item) for item in args.z_range.split(',')]

    if args.lsm_bin_edges == 'all':
        args.lsm_bin_edges = all_lsm_bin_edges
    else:
        args.lsm_bin_edges = [float(item) for item in args.lsm_bin_edges.split(',')]
    args.lsm_bins = list(zip(args.lsm_bin_edges, args.lsm_bin_edges[1:]))

    if args.lsfr_bin_edges == 'all':
        args.lsfr_bin_edges = all_lsfr_bin_edges
    else:
        args.lsfr_bin_edges = [float(item) for item in args.lsfr_bin_edges.split(',')]
    args.lsfr_bins = list(zip(args.lsfr_bin_edges, args.lsfr_bin_edges[1:]))

    args.inc_range = args.inc_bins[0]
    args.lsfr_range = args.lsfr_bins[0]

    # ---------------aesthetics---------------------
    if args.fortalk:
        print(f'Setting up plots for talks..')
        setup_plots_for_talks()
        
    return args
