import numpy as np 

from scipy.special import erf 

from setup import setup
from start import start

se = setup('./')
star = start(se)

# define stratification profile constants
N2s = 10**-4 # stratification in the surface layer
N2m = 1.5*10**-3 # stratification in the middle layer
zf = -4 # filament structure at x = 0
hf = 0.25 # vertical length scale for filament structure
N2f = 5*10**-3 # sets peak filament stratification
zm = -10 # base of the middle layer
N2p = 1.2*10**-2 # stratification peak
zp = -14 # stratification peak location
zt = (zm + zp)/2 # middle of the transition region
ht = (zm - zp)/8 # height scale of transition layer

def b_profile(z0):
    return (
        (N2s - N2m)*hf*np.log(1 + np.exp((z0-zf)/hf)) + N2f*hf*np.sqrt(np.pi/2)*erf((z0-zf)/np.sqrt(2)/hf)
        + (N2m - N2p)*ht*np.log(1 + np.exp((z0-zt)/ht)) + N2p*z0
    )

star.U[:,:,:] = 0
star.V[:,:,:] = 0
star.W[:,:,:] = 0
star.TH[0,:,:,:] = b_profile(star.GF[None,:,None])

star.write_start_file()