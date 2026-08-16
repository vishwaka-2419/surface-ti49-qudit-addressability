import numpy as np

KAPPA49_MHZ = -46.373841059602654
GN49 = 0.3146857142857143
MU_N_OVER_H_MHZ_PER_T = 7.622593285
I49 = 3.5

def comb_spacing_mhz(kappa_mhz=KAPPA49_MHZ, I=I49):
    return 3.0 * abs(kappa_mhz) / (2.0*I*(2.0*I-1.0))

def gamma_n_mhz_per_t(g_n=GN49):
    return g_n * MU_N_OVER_H_MHZ_PER_T

def period_t(kappa_mhz=KAPPA49_MHZ, g_n=GN49):
    return comb_spacing_mhz(kappa_mhz)/(2.0*gamma_n_mhz_per_t(g_n))

def ceiling_mhz(kappa_mhz=KAPPA49_MHZ):
    return comb_spacing_mhz(kappa_mhz)/2.0

def dmin_mhz(B_t, kappa_mhz=KAPPA49_MHZ, g_n=GN49, kmax=6):
    B=np.atleast_1d(np.asarray(B_t,dtype=float))
    S=comb_spacing_mhz(kappa_mhz); gam=gamma_n_mhz_per_t(g_n)
    ks=np.arange(-kmax,kmax+1,dtype=float)
    out=np.min(np.abs(S*ks[None,:]-2.0*gam*B[:,None]),axis=1)
    return out if np.ndim(B_t) else float(out[0])

def crossing_fields_t(n=4):
    p=period_t(); return np.arange(1,n+1,dtype=float)*p

def local_max_fields_t(n=4):
    p=period_t(); return (np.arange(n,dtype=float)+0.5)*p
