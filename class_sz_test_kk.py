import numpy as np
import matplotlib.pyplot as plt
# import healpy as hp
from scipy.interpolate import interp1d
import sys
from classy_sz import Class
from mcfit.transforms import *
from scipy import interpolate
import time 


#matplotlib.use('pdf')
font = {'size'   : 16, 'family':'STIXGeneral'}
plt.rcParams.update({
     "text.usetex": True,
     "font.family": "serif",
     "font.sans-serif": ['Computer Modern']})
plt.rc_context({'axes.autolimit_mode': 'round_numbers'})


def l_to_dl(lp):
    return lp*(lp+1.)/2./np.pi
cosmo_params = {
        'omega_b': 0.02242,
        'omega_cdm':  0.11933,
        'H0': 67.66, # use H0 because this is what is used by the emulators.
        'tau_reio': 0.0561,
        # 'sigma8': 0.81,
        # 'sigma8': 0.81,
        'ln10^{10}A_s': 3.0980,
        'n_s': 0.9665,

        'k_pivot': 0.05,
        'N_ncdm': 1,
        'N_ur': 2.0328,
        'm_ncdm': 0.06,
        # 'z_max_pk':20.

        'output': 'lens_lens_1h,lens_lens_2h',#,lens_lens_1h,lens_lens_2h',
        # 'ndim_redshifts':30, # this is the number of redshift points used to tabulate Pk from the emulators. If you dont need Pk, set this to 4. 
        # the parameter ndim_redshifts is critical ! the more z points you add, the more calls to emulators
        # 'skip_background_and_thermo':1,
        # 'skip_chi': 0,
        # 'skip_hubble':0,

        'ell_max': 20000.0,
        'ell_min': 2.0,
        'dlogell': 0.1,
        'dell': 0,
        'redshift_epsrel': 0.0001,
        'mass_epsabs': 1e-40,
        'mass_epsrel': 0.0001,
            
        'M_min':1e7,
        'M_max':1e+17,
        'z_min':1e-5,
        'z_max': 15.,
            
        # 'k_min_for_pk_class_sz' : 1e-4,
        # 'k_max_for_pk_class_sz' : 5e1,
        # 'k_per_decade_class_sz' : 20.,
        # 'P_k_max_h/Mpc' : 200.0,

            
        # 'ndim_masses' : 150, # important 128 is default ccl value
        'ndim_redshifts' : 50,
        'non_linear':'hmcode',
        # 'perturb_sampling_stepsize' : 0.005,
        # 'k_max_tau0_over_l_max':5.,
            
        'hm_consistency': 1,

        'use_pknl_in_2hterms' : 1,

}


# M = Class()
def compute_class_fast(M):
    
    M.set(cosmo_params)
    M.set({


    })
    M.compute_class_szfast()

# M = Class()
def compute_class_slow(M):
    
    M.set(cosmo_params)
    M.set({


    })
    M.compute()

M_fast = Class()
start = time.perf_counter()
compute_class_fast(M_fast)
end = time.perf_counter()
print('>>> class_szfast took %.3f s'%(end-start))


cl_kk_fast = M_fast.cl_kk
ell_fast = np.asarray(cl_kk_fast()['ell'])
fac_fast = ell_fast*(ell_fast+1.)/2./np.pi
cl_kk_1h_fast = np.asarray(cl_kk_fast()['1h'])
cl_kk_2h_fast = np.asarray(cl_kk_fast()['2h'])

# print(cl_kk_1h_fast,cl_kk_2h_fast)
# exit(0)

M_slow = Class()
start = time.perf_counter()
compute_class_slow(M_slow)
end = time.perf_counter()
print('>>> class_szslow took %.3f s'%(end-start))

cl_kk_slow = M_slow.cl_kk
ell_slow = np.asarray(cl_kk_slow()['ell'])
fac_slow = ell_slow*(ell_slow+1.)/2./np.pi
cl_kk_1h_slow = np.asarray(cl_kk_slow()['1h'])
cl_kk_2h_slow = np.asarray(cl_kk_slow()['2h'])







label_size = 17
title_size = 22
legend_size = 13
handle_length = 1.5
fig, (ax1) = plt.subplots(1,1,figsize=(10,5))

ax = ax1
ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
plt.setp(ax.get_xticklabels(), fontsize=label_size)
ax.grid( visible=True, which="both", alpha=0.1, linestyle='--')

ax.set_ylabel(r'$C_\ell$',size=title_size)
ax.set_xlabel(r'$\ell$',size=title_size)



ax.plot(ell_fast,cl_kk_1h_fast/fac_fast,ls='--',c='b',label=r'$\mathrm{1}$-$\mathrm{halo}$')
ax.plot(ell_fast,cl_kk_2h_fast/fac_fast,ls='-.',c='b',label=r'$\mathrm{2}$-$\mathrm{halo}$')
ax.plot(ell_fast,cl_kk_2h_fast/fac_fast+cl_kk_1h_fast/fac_fast,ls='-',c='k',label=r'$\mathrm{1+2}$-$\mathrm{halo}$')

ax.plot(ell_fast,cl_kk_1h_slow/fac_slow,ls='--',c='b',label=r'slow $\mathrm{1}$-$\mathrm{halo}$')
ax.plot(ell_fast,cl_kk_2h_slow/fac_slow,ls='-.',c='b',label=r'slow $\mathrm{2}$-$\mathrm{halo}$')
ax.plot(ell_fast,cl_kk_2h_slow/fac_slow+cl_kk_1h_slow/fac_slow,ls='--',c='r',label=r'slow $\mathrm{1+2}$-$\mathrm{halo}$')



ax.loglog()

# ax.set_ylim(1e1,1e5)
# ax.set_xlim(2e-3,1e1)

ax.legend(loc=1,frameon=True,framealpha=1,fontsize=11)

# ax.set_title(r'$z=%f$'%z)




fig.tight_layout()
plt.show()
# plt.savefig('figures/pkz.pdf')
