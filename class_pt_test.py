import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy_sz import Class
import os
import time




# the parameters needed for cosmology:
# use the last column of Planck 2018 (https://arxiv.org/pdf/1807.06209.pdf) Table 2
# TT,TE,EE+lowE+lensing+BAO
cosmo_params = {
'omega_b': 0.02242,
'omega_cdm':  0.11933,
'H0': 67.66, # use H0 because this is what is used by the emulators.
'tau_reio': 0.0561,
'ln10^{10}A_s': 3.047,
'n_s': 0.9665,

'k_pivot': 0.05,
'N_ncdm': 1,
'N_ur': 2.0328,
'm_ncdm': 0.06,

'non_linear': 'hmcode',
# 'l_max_scalars': 11000,
}

# a simple conversion from cl's to dl's
def l_to_dl(lp):
    return lp*(lp+1.)/2./np.pi

import classy_sz
classy_sz.__file__



z_pk = 0.5
cosmo = Class()
# cosmo.set({'A_s':2.089e-9,
#            'n_s':0.9649,
#            'tau_reio':0.052,
#            'omega_b':0.02237,
#            'omega_cdm':0.12,
#            'h':0.6736,
#            'YHe':0.2425,
#            'N_ur':2.0328,
#            'N_ncdm':1,
#            'm_ncdm':0.06,
#            'z_pk':z_pk
#           })
cosmo.set(cosmo_params)
# Set additional CLASS-PT settings
cosmo.set({'output':'mPk',
           'z_pk':z_pk,
           'non_linear':'PT',
           'IR resummation':'Yes',
           'Bias tracers':'Yes',
           'cb':'Yes', # use CDM+baryon spectra
           'RSD':'Yes',
           'AP':'Yes', # Alcock-Paczynski effect
           'Omfid':'0.31', # fiducial Omega_m
           'PNG':'No', # single-field inflation PNG
'nonlinear_pt_verbose':1,
'class_sz_verbose':1,
'skip_background_and_thermo': 0,
'skip_pkl': 0,
'skip_pknl': 1,
'skip_sigma8_and_der': 1,
'skip_sigma8_at_z': 1,
'skip_hubble':0,
'skip_chi':0,
'skip_cmb':1,
'ndim_redshifts':20
         })

start = time.time()
cosmo.compute_class_szfast()
end = time.time()
print('class_sz with pt took:',end-start)
# print('cosmo.pk(0.3,0.5)[0] = ',cosmo.pk(0.3,0.5)[0])
for i in range(10):
    print('cosmo.pk(0.3,0.5)[%d] = %.5e'%(i,cosmo.pk(0.3,0.5)[i]))

# cosmo.compute()
