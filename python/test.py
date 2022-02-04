import rrtmg_lw_wrapper as lw
import numpy as np
import matplotlib.pyplot as plt

# Constants
cp = 1e3
g = 9.80
Lv = 2.5e6
Rv = 460e0
Rd = 287e0
e0 = 610e0
T0 = 273e0
Mair = 29e0
Mh2o = 18e0
Mco2 = 44e0
h = 6.626e-34
kB = 1.381e-23
c = 3.00e8
NUM_BANDS = 16


# Atmospheric profiles
def calc_qsat(T, p):
    esat = e0*np.exp(-Lv/Rv*(1/T - 1/T0))
    return esat*Rd/(p*Rv)

def calc_profiles(Ts, Gamma, zt, RH, 
        n = 1000, ps = 101325e0, ztop = 25e3):
    
    # Set height grid
    zf = np.linspace(0, ztop, n+1)
    zc = 0.5*(zf[1:] + zf[:-1])
    
    # Calculate temperatures
    Tt = Ts - Gamma*zt
    Tf = Ts - Gamma*zf
    Tf[zf >= zt] = Tt
    Tc = Ts - Gamma*zc
    Tc[zc >= zt] = Tt

    # Calculate pressures
    pt = ps*((Ts - Gamma*zt)/Ts)**(g/(Rd*Gamma))
    pf = ps*(Tf/Ts)**(g/(Rd*Gamma))
    pf[zf >= zt] = pt*np.exp(-g*(zf[zf >= zt] - zt)/(Rd*Tt))
    pc = ps*(Tc/Ts)**(g/(Rd*Gamma))
    pc[zc >= zt] = pt*np.exp(-g*(zc[zc >= zt] - zt)/(Rd*Tt))

    # Calculate humidity
    qc = RH*calc_qsat(Tc, pc)

    # Return profiles
    return (zc, zf, pc, pf, Tc, Tf, qc)

# Planck function
def calc_B(nu, T):
    return 2*h*c**2*nu**3/(np.exp(h*c*nu/(kB*T)) - 1)

# RRTMG interface
def rrtmg_zeros(shape):
    return np.zeros(shape, dtype = float, order = 'F')

def run_rrtmg_lw(Ts, Gamma, zt, pco2, pch4, RH, 
        n = 1000, ps = 101325e0, ztop = 25e3):

    # Generate profiles
    zc, zf, pc, pf, Tc, Tf, q = calc_profiles(
            Ts, Gamma, zt, RH, n = n, ps = ps, ztop = ztop)

    # Set inputs
    ncol = 1
    nlay = n
    icld = 0
    idrv = 1
    play = rrtmg_zeros((ncol,nlay))
    play[0,:] = pc/1e2
    plev = rrtmg_zeros((ncol,nlay+1))
    plev[0,:] = pf/1e2
    tlay = rrtmg_zeros((ncol,nlay))
    tlay[0,:] = Tc
    tlev = rrtmg_zeros((ncol,nlay+1))
    tlev[0,:] = Tf
    tsfc = rrtmg_zeros((ncol,))
    tsfc[0] = Ts
    h2ovmr = rrtmg_zeros((ncol,nlay))
    h2ovmr[0,:] = q*Rv/Rd
    o3vmr = rrtmg_zeros((ncol,nlay))
    co2vmr = rrtmg_zeros((ncol,nlay))
    co2vmr[0,:] = pco2
    ch4vmr = rrtmg_zeros((ncol,nlay))
    ch4vmr[0,:] = pch4
    n2ovmr = rrtmg_zeros((ncol,nlay))
    o2vmr = rrtmg_zeros((ncol,nlay))
    o2vmr[0,:] = 0.2
    cfc11vmr = rrtmg_zeros((ncol,nlay))
    cfc12vmr = rrtmg_zeros((ncol,nlay))
    cfc22vmr = rrtmg_zeros((ncol,nlay))
    ccl4vmr = rrtmg_zeros((ncol,nlay))
    emis = rrtmg_zeros((ncol,NUM_BANDS))
    emis[0,:] = 1
    inflglw = 2     # provide cloud water paths and droplet radii,
                    # treat water and ice clouds separately
    iceflglw = 1    # Ebert & Curry 1992
    liqflglw = 1    # radius-dependent absorption
    cldfr = rrtmg_zeros((ncol,nlay))
    cicewp = rrtmg_zeros((ncol,nlay))
    cliqwp = rrtmg_zeros((ncol,nlay))
    reice = rrtmg_zeros((ncol,nlay))
    reliq = rrtmg_zeros((ncol,nlay))
    taucld = rrtmg_zeros((NUM_BANDS,ncol,nlay))
    tauaer = rrtmg_zeros((ncol,nlay,NUM_BANDS))
    start_band = 1
    end_band = NUM_BANDS
    
    # Output arrays
    uflx = rrtmg_zeros((ncol,nlay+1))
    dflx = rrtmg_zeros((ncol,nlay+1))
    hr = rrtmg_zeros((ncol,nlay))
    uflxc = rrtmg_zeros((ncol,nlay+1))
    dflxc = rrtmg_zeros((ncol,nlay+1))
    hrc = rrtmg_zeros((ncol,nlay))
    wavenumber_range = rrtmg_zeros((2,))
    duflx_dt = rrtmg_zeros((ncol,nlay+1))
    duflxc_dt = rrtmg_zeros((ncol,nlay+1))

    # Run RRTM
    OLR = np.zeros((NUM_BANDS,))
    nu_low = np.zeros((NUM_BANDS,))
    nu_high = np.zeros((NUM_BANDS,))
    for iband in range(1,NUM_BANDS+1):
        lw.rrtmg_lw_wrapper(
            ncol, nlay, icld, idrv,                           
            play, plev, tlay, tlev, tsfc,                     
            h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,     
            cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis,      
            inflglw, iceflglw, liqflglw, cldfr,               
            taucld, cicewp, cliqwp, reice, reliq,             
            tauaer,         
            iband, iband,
            uflx, dflx, hr, uflxc, dflxc, hrc,  
            wavenumber_range,
            duflx_dt, duflxc_dt )
        OLR[iband-1] = uflx[0,-1]
        nu_low[iband-1] = wavenumber_range[0]
        nu_high[iband-1] = wavenumber_range[1]

    # Return results
    return OLR, nu_low, nu_high

# Initialize
lw.rrtmg_lw_ini_wrapper(cp)

# Compute spectrally-resolved OLR
# Use reference state from UChicago RRTM page
# http://climatemodels.uchicago.edu/rrtm/
OLR, nu_low, nu_high = run_rrtmg_lw(
    284.42, 6e-3, 15e3, 400e-6, 1.7e-6, 0.8)

# Create arrays for boxy plots
OLR_nu = OLR/(nu_high - nu_low)
OLR_nu_box = np.stack((OLR_nu, OLR_nu), axis = -1).ravel()
nu_box = np.stack((nu_low, nu_high), axis = -1).ravel()

# Plot
plt.rc('font', size = 8)
fig, axes = plt.subplots(figsize = (3.25, 3), nrows = 1, ncols = 1,
        constrained_layout = True, dpi = 200)
axes = np.array(axes)[np.newaxis,np.newaxis]
axes[0,0].plot(nu_box, OLR_nu_box, 'k-')
axes[0,0].set_title('Total OLR: %.1f W m$^{-2}$' % (np.sum(OLR)), fontsize = 8)
axes[0,0].set_xlabel(r'$\nu$ (cm $^{-1}$)')
axes[0,0].set_ylabel('OLR (W m$^{-2}$ cm)')

# Add reference Planck functions
axes[0,0].set_xlim([nu_box[0], nu_box[-1]])
nu = np.linspace(nu_box[0], nu_box[-1], 100)*1e2
axes[0,0].plot(nu/1e2, 1e2*np.pi*calc_B(nu, 280e0), 
        '-', color = 'gray', label = '280 K', zorder = -1)
axes[0,0].plot(nu/1e2, 1e2*np.pi*calc_B(nu, 250e0), 
        '--', color = 'gray', label = '240 K', zorder = -1)
axes[0,0].plot(nu/1e2, 1e2*np.pi*calc_B(nu, 200e0), 
        ':', color = 'gray', label = '200 K', zorder = -1)
axes[0,0].legend(loc = 'upper right', frameon = False)
axes[0,0].set_ylim([0, axes[0,0].get_ylim()[1]])

plt.savefig('test.pdf')
plt.show()
