import healpy as hp
import numpy as np

nmod = 100

# Make a cl import

map_I = hp.read_map('../data/COM_CMB_IQU-smica_1024_R2.02_full.fits', field=0)
map_Q = hp.read_map('../data/COM_CMB_IQU-smica_1024_R2.02_full.fits', field=1)
map_U = hp.read_map('../data/COM_CMB_IQU-smica_1024_R2.02_full.fits', field=2)

alm_I = hp.sphtfunc.map2alm(map_I, lmax=nmod)
alm_Q = hp.sphtfunc.map2alm(map_Q, lmax=nmod)
alm_U = hp.sphtfunc.map2alm(map_U, lmax=nmod)

hp.fitsfunc.write_alm('../data/alm_I_2048.fits', alm_I)
hp.fitsfunc.write_alm('../data/alm_Q_2048.fits', alm_Q)
hp.fitsfunc.write_alm('../data/alm_U_2048.fits', alm_U)

alm_norm_I = open('../data/alm_I_norm_2048.dat', 'w')
alm_norm_Q = open('../data/alm_Q_norm_2048.dat', 'w')
alm_norm_U = open('../data/alm_U_norm_2048.dat', 'w')

alm_file_I = hp.fitsfunc.read_alm('../data/alm_I_2048.fits')
alm_file_Q = hp.fitsfunc.read_alm('../data/alm_Q_2048.fits')
alm_file_U = hp.fitsfunc.read_alm('../data/alm_U_2048.fits')

lvals, mvals = hp.sphtfunc.Alm.getlm(nmod, np.arange(hp.sphtfunc.Alm.getsize(lmax)))

alm_norm__I = np.zeros((lmax + 1, lmax + 1), dtype=complex)
alm_norm__Q = np.zeros((lmax + 1, lmax + 1), dtype=complex)
alm_norm_U = np.zeros((lmax + 1, lmax + 1), dtype=complex)


for i in xrange(0, np.size(lvals)):
    m = mvals[i]
    l = lvals[i]

    # alm_norm_I[m][l] = alm_file_I[i]
    alm_norm_2_Q[m][l] = alm_file_Q[i]
    alm_norm_2_U[m][l] = alm_file_U[i]

for m in xrange(0, lmax+1):
    for l in xrange(0, lmax+1):
        # alm_norm_I.write(repr(real(alm_normal[m][l])) + '    ' + repr(m) + '    ' + repr(l) + '\n')
        alm_norm_Q.write(repr(m) + '  ' + repr(l) + '  ' + repr(alm_norm_2_Q[m][l].real) + '\n')
        alm_norm_U.write(repr(m) + '  ' + repr(l) + '  ' + repr(alm_norm_2_U[m][l].real) + '\n')