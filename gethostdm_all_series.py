#
#	Script for estimating FRB host galaxy DM from simulated electron density cubes
#
#								Originally by AB, August 2025
#                               Modified by AA, October 2026

#	--------------------------	Import modules	---------------------------
from craft_header import *
from craft_utils import *
setup_plot_style()
from globalpars import *
from nefns import *
from plotdm import *

start_time = datetime.now()

# -------------------------------------------------------------------------------------------
def print_instructions():
    '''
	Print instructions to terminal
	'''
    print("\n            You probably need some assistance here!\n")
    print("\n Arguments are       --- <mode> <nfixpts> <extent/ckpc> <scale/kpc> <filename/optional>\n")
    print(" Supported Modes are --- profile        (calculate electron density profiles)")
    print("                     --- plot_profile          (plot electron density profiles)")
    print("                     --- losdm          (calculate LoS DMs)")
    print("                     --- pltdm          (Plot LoS DMs)")
    print("                     --- dmscat         (Plot LoS DMs)")
    print("\n            Now let's try again!\n")
	
    return(0)

# -----main code-----------------
if __name__ == '__main__':
    #	--------------------------	Read inputs	-------------------------------
    if(len(sys.argv)<3):
        print_instructions()
        sys.exit()

    exmode		=	sys.argv[1]					#	What to do	
    nfixpts     =   int(sys.argv[2])            #   Number of fixed points on each face to simulate LoSs
    extent     =   float(sys.argv[3])             #   Extent of datacube on either side of the center, in units of kpc (to pick the right cubes from datadir)
    scalekpc	=	float(sys.argv[4])			#	Scale radius in kpc
    try: 
        filename	=	sys.argv[5]			#	Individual fits filename
    except:
        filename = ''
        pass

    #incranges	=	np.array([[0,20],[40,50],[80,90]])
    incranges	=	np.array([[item - dinc/2, item + dinc/2] for item in incvals])

    # --------determining file list-------------
    if filename != '':
        list_of_fits = glob.glob(datadir + f'{filename.replace(".fits", "")}.fits')
    else:
        list_of_fits = glob.glob(datadir + f'*El_number_density*{extent:.1f}kpc*{scalekpc:.1f}kpc*.fits') # all snapshots of this particular halo
    total_snaps = len(list_of_fits)

    # -------------loop over snapshots-----------------
    print(f'Operating on {total_snaps} snapshots')

    for index in range(total_snaps):
        start_time_this_snapshot = datetime.now()
        thisfile = Path(list_of_fits[index])
        fitsname = thisfile.stem

        if 'El' in fitsname: profsubdir = 'electron_density/'
        else: profsubdir = 'gas_density/'
        profdir = radialdir + profsubdir
        Path(profdir).mkdir(exist_ok=True, parents=True)
        profile_pkl_filename = profdir + fitsname + '_radprof.pkl'
        
        if fitsname[-3:] == str(nfixpts): fitsname = fitsname[:-4]
        this_sim = fitsname.split('_')[:2]
        print('Doing snapshot ' + this_sim[0] + ' of halo ' + this_sim[1] + ' which is ' + str(index + 1) + ' out of the total ' + str(total_snaps) + ' snapshots...')

        #	-------------------------	Load the fits file	---------------------------
        if exmode in ['losdm', 'profile']:
            if not (exmode == 'profile' and os.path.exists(profile_pkl_filename)):
                print("Reading "+fitsname)
                necub,dkpc,theta0,phi0	=	fitld(fitsname,3.2)
                print(f"Ne cube dimensions {necub.shape}")
                print(f"Spatial resolutions (kpc) {dkpc}")
                print(f"Orientation (deg) {np.rad2deg(theta0)},{np.rad2deg(phi0)}")

        #	-------------------------	Execute tasks	-------------------------------

        if (exmode=='profile'):
            if not os.path.exists(profile_pkl_filename):
                print("\nGenerating radial electron density profiles...\n")
                cubene	= neprofinc(necub,dkpc,1.0,theta0,phi0,1.0,1.0,1.0)
                
                with open(profile_pkl_filename, 'wb') as file_obj:
                    pkl.dump(cubene, file_obj) # dump the pickle file
            else:
                print(f"\nUsing existing radial electron density profile {profile_pkl_filename}\n")
            
            print("\nPlotting radial electron density profiles...\n")            
            with open(profile_pkl_filename, 'rb') as file_obj:
                cubene = pkl.load(file_obj) # load the pickle file
            
            prof_figdir = plotradial + profsubdir
            Path(prof_figdir).mkdir(exist_ok=True, parents=True)
            plot_nerad(cubene, incranges, prof_figdir + fitsname)
            if total_snaps > 10: plt.close('all')

        elif (exmode=='losdm'):
            print("\nEstimating LoS DMs...\n")
            losdms(fitsname,necub,dkpc,theta0,phi0,nfixpts,1.0,1.0,1.0, los_extent_kpc) # last argument is extent of shooting LoS (in kpc), the value is in globalpars.py

        elif (exmode=='pltdm'):
            print("\nPloting LoS DMs...\n")
            plotdms(fitsname,nfixpts,1.0,1.0,1.0,scalekpc)

        elif (exmode=='dmscat'):
            print("\nPloting LoS DMs...\n")
            plotdm2d(fitsname,nfixpts,1.0,1.0,1.0,scalekpc)

        else:
            print("\nHmm...What mode is that again...?\n")

        plt.show(block=False)
        print('This snapshots completed in %s' % timedelta(seconds=(datetime.now() - start_time_this_snapshot).seconds))

    # -----------------------------------------------------------------------------------
    print('Serially: time taken for ' + str(total_snaps) + ' snapshot was %s' % timedelta(seconds=(datetime.now() - start_time).seconds))