/** @file szpowerspectrum.h Documented includes for sz module. */
#ifndef __SZ__
#define __SZ__

#include "common.h"
#include "lensing.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_sf_expint.h"
#include "gsl/gsl_sf_lambert.h"
# include <fftw3.h>
// # include <gsl_integration.h>
// #include "fft.h"


#define _pk_at_z_1h_ ((ptsz->has_pk_at_z_1h == _TRUE_) && (index_md == ptsz->index_md_pk_at_z_1h))
#define _pk_at_z_2h_ ((ptsz->has_pk_at_z_2h == _TRUE_) && (index_md == ptsz->index_md_pk_at_z_2h))
#define _pk_gg_at_z_1h_ ((ptsz->has_pk_gg_at_z_1h == _TRUE_) && (index_md == ptsz->index_md_pk_gg_at_z_1h))
#define _pk_gg_at_z_2h_ ((ptsz->has_pk_gg_at_z_2h == _TRUE_) && (index_md == ptsz->index_md_pk_gg_at_z_2h))
#define _pk_bb_at_z_1h_ ((ptsz->has_pk_bb_at_z_1h == _TRUE_) && (index_md == ptsz->index_md_pk_bb_at_z_1h))
#define _pk_bb_at_z_2h_ ((ptsz->has_pk_bb_at_z_2h == _TRUE_) && (index_md == ptsz->index_md_pk_bb_at_z_2h))
#define _pk_b_at_z_2h_ ((ptsz->has_pk_b_at_z_2h == _TRUE_) && (index_md == ptsz->index_md_pk_b_at_z_2h))
#define _pk_em_at_z_1h_ ((ptsz->has_pk_em_at_z_1h == _TRUE_) && (index_md == ptsz->index_md_pk_em_at_z_1h))
#define _pk_em_at_z_2h_ ((ptsz->has_pk_em_at_z_2h == _TRUE_) && (index_md == ptsz->index_md_pk_em_at_z_2h))
#define _pk_HI_at_z_1h_ ((ptsz->has_pk_HI_at_z_1h == _TRUE_) && (index_md == ptsz->index_md_pk_HI_at_z_1h))
#define _pk_HI_at_z_2h_ ((ptsz->has_pk_HI_at_z_2h == _TRUE_) && (index_md == ptsz->index_md_pk_HI_at_z_2h))
#define _bk_at_z_1h_ ((ptsz->has_bk_at_z_1h == _TRUE_) && (index_md == ptsz->index_md_bk_at_z_1h))
#define _bk_at_z_2h_ ((ptsz->has_bk_at_z_2h == _TRUE_) && (index_md == ptsz->index_md_bk_at_z_2h))
#define _bk_at_z_3h_ ((ptsz->has_bk_at_z_3h == _TRUE_) && (index_md == ptsz->index_md_bk_at_z_3h))
#define _bk_ttg_at_z_1h_ ((ptsz->has_bk_ttg_at_z_1h == _TRUE_) && (index_md == ptsz->index_md_bk_ttg_at_z_1h))
#define _bk_ttg_at_z_2h_ ((ptsz->has_bk_ttg_at_z_2h == _TRUE_) && (index_md == ptsz->index_md_bk_ttg_at_z_2h))
#define _bk_ttg_at_z_3h_ ((ptsz->has_bk_ttg_at_z_3h == _TRUE_) && (index_md == ptsz->index_md_bk_ttg_at_z_3h))

//#define _bk_at_z_hf_ ((ptsz->has_bk_at_z_hf == _TRUE_) && (index_md == ptsz->index_md_bk_at_z_hf))
#define _mean_y_ ((ptsz->has_mean_y == _TRUE_) && (index_md == ptsz->index_md_mean_y))
#define _cib_monopole_ ((ptsz->has_cib_monopole == _TRUE_) && (index_md == ptsz->index_md_cib_monopole))
#define _cib_shotnoise_ ((ptsz->has_cib_shotnoise == _TRUE_) && (index_md == ptsz->index_md_cib_shotnoise))
#define _dcib0dz_ ((ptsz->has_dcib0dz == _TRUE_) && (index_md == ptsz->index_md_dcib0dz))
#define _dydz_ ((ptsz->has_dydz == _TRUE_) && (index_md == ptsz->index_md_dydz))
#define _hmf_ ((ptsz->has_hmf == _TRUE_) && (index_md == ptsz->index_md_hmf))
#define _tSZ_power_spectrum_ ((ptsz->has_sz_ps == _TRUE_) && (index_md == ptsz->index_md_sz_ps))
#define _trispectrum_ ((ptsz->has_sz_trispec == _TRUE_) && (index_md == ptsz->index_md_trispectrum))
#define _2halo_ ((ptsz->has_sz_2halo == _TRUE_) && (index_md == ptsz->index_md_2halo))
#define _te_y_y_ ((ptsz->has_sz_te_y_y == _TRUE_) && (index_md == ptsz->index_md_te_y_y))
#define _m_y_y_1h_ ((ptsz->has_sz_m_y_y_1h == _TRUE_) && (index_md == ptsz->index_md_m_y_y_1h))
#define _m_y_y_2h_ ((ptsz->has_sz_m_y_y_2h == _TRUE_) && (index_md == ptsz->index_md_m_y_y_2h))
#define _cov_Y_N_ ((ptsz->has_sz_cov_Y_N == _TRUE_) && (index_md == ptsz->index_md_cov_Y_N))
#define _cov_N_N_ ((ptsz->has_sz_cov_N_N == _TRUE_) && (index_md == ptsz->index_md_cov_N_N))
#define _cov_N_N_hsv_ ((ptsz->has_sz_cov_N_N_hsv == _TRUE_) && (index_md == ptsz->index_md_cov_N_N_hsv))
#define _cov_Y_Y_ssc_ ((ptsz->has_sz_cov_Y_Y_ssc == _TRUE_) && (index_md == ptsz->index_md_cov_Y_Y_ssc))
#define _cov_Y_N_next_order_ ((ptsz->has_sz_cov_Y_N_next_order == _TRUE_) && (index_md == ptsz->index_md_cov_Y_N_next_order))
#define _kSZ_kSZ_gal_1h_ ((ptsz->has_kSZ_kSZ_gal_1h == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gal_1h))
#define _kSZ_kSZ_gal_1h_fft_ ((ptsz->has_kSZ_kSZ_gal_1h_fft == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gal_1h_fft))
#define _kSZ_kSZ_gal_2h_fft_ ((ptsz->has_kSZ_kSZ_gal_2h_fft == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gal_2h_fft))
#define _kSZ_kSZ_gal_3h_fft_ ((ptsz->has_kSZ_kSZ_gal_3h_fft == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gal_3h_fft))
#define _kSZ_kSZ_gallens_1h_fft_ ((ptsz->has_kSZ_kSZ_gallens_1h_fft == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gallens_1h_fft))
#define _kSZ_kSZ_gallens_2h_fft_ ((ptsz->has_kSZ_kSZ_gallens_2h_fft == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gallens_2h_fft))
#define _kSZ_kSZ_gallens_3h_fft_ ((ptsz->has_kSZ_kSZ_gallens_3h_fft == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gallens_3h_fft))
#define _kSZ_kSZ_lens_1h_fft_ ((ptsz->has_kSZ_kSZ_lens_1h_fft == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_lens_1h_fft))
#define _kSZ_kSZ_lens_2h_fft_ ((ptsz->has_kSZ_kSZ_lens_2h_fft == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_lens_2h_fft))
#define _kSZ_kSZ_lens_3h_fft_ ((ptsz->has_kSZ_kSZ_lens_3h_fft == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_lens_3h_fft))
#define _kSZ_kSZ_gal_2h_ ((ptsz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gal_2h))
#define _kSZ_kSZ_gal_3h_ ((ptsz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gal_3h))
#define _kSZ_kSZ_gal_hf_ ((ptsz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gal_hf))
#define _kSZ_kSZ_gallens_hf_ ((ptsz->has_kSZ_kSZ_gallens_hf == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gallens_hf))
#define _kSZ_kSZ_lens_hf_ ((ptsz->has_kSZ_kSZ_lens_hf == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_lens_hf))
#define _kSZ_kSZ_lensmag_1halo_ ((ptsz->has_kSZ_kSZ_lensmag_1halo == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_lensmag_1halo))
#define _gal_gal_1h_ ((ptsz->has_gal_gal_1h == _TRUE_) && (index_md == ptsz->index_md_gal_gal_1h))
#define _gal_gal_2h_ ((ptsz->has_gal_gal_2h == _TRUE_) && (index_md == ptsz->index_md_gal_gal_2h))
#define _n5k_ ((ptsz->has_n5k == _TRUE_) && (index_md == ptsz->index_md_n5k))
#define _gal_gal_hf_ ((ptsz->has_gal_gal_hf == _TRUE_) && (index_md == ptsz->index_md_gal_gal_hf))
#define _gal_lens_2h_ ((ptsz->has_gal_lens_2h == _TRUE_) && (index_md == ptsz->index_md_gal_lens_2h))
#define _gal_lens_hf_ ((ptsz->has_gal_lens_hf == _TRUE_) && (index_md == ptsz->index_md_gal_lens_hf))
#define _gal_lens_1h_ ((ptsz->has_gal_lens_1h == _TRUE_) && (index_md == ptsz->index_md_gal_lens_1h))
#define _gal_lensmag_2h_ ((ptsz->has_gal_lensmag_2h == _TRUE_) && (index_md == ptsz->index_md_gal_lensmag_2h))
#define _gal_lensmag_1h_ ((ptsz->has_gal_lensmag_1h == _TRUE_) && (index_md == ptsz->index_md_gal_lensmag_1h))
#define _gal_gallens_2h_ ((ptsz->has_gal_gallens_2h == _TRUE_) && (index_md == ptsz->index_md_gal_gallens_2h))
#define _gal_gallens_1h_ ((ptsz->has_gal_gallens_1h == _TRUE_) && (index_md == ptsz->index_md_gal_gallens_1h))
#define _gallens_gallens_2h_ ((ptsz->has_gallens_gallens_2h == _TRUE_) && (index_md == ptsz->index_md_gallens_gallens_2h))
#define _gallens_gallens_1h_ ((ptsz->has_gallens_gallens_1h == _TRUE_) && (index_md == ptsz->index_md_gallens_gallens_1h))
#define _gallens_lens_2h_ ((ptsz->has_gallens_lens_2h == _TRUE_) && (index_md == ptsz->index_md_gallens_lens_2h))
#define _gallens_lens_1h_ ((ptsz->has_gallens_lens_1h == _TRUE_) && (index_md == ptsz->index_md_gallens_lens_1h))
#define _gal_lensmag_hf_ ((ptsz->has_gal_lensmag_hf == _TRUE_) && (index_md == ptsz->index_md_gal_lensmag_hf))
#define _tSZ_lensmag_2h_ ((ptsz->has_tSZ_lensmag_2h == _TRUE_) && (index_md == ptsz->index_md_tSZ_lensmag_2h))
#define _tSZ_lensmag_1h_ ((ptsz->has_tSZ_lensmag_1h == _TRUE_) && (index_md == ptsz->index_md_tSZ_lensmag_1h))
#define _lensmag_lensmag_hf_ ((ptsz->has_lensmag_lensmag_hf == _TRUE_) && (index_md == ptsz->index_md_lensmag_lensmag_hf))
#define _lensmag_lensmag_2h_ ((ptsz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == ptsz->index_md_lensmag_lensmag_2h))
#define _lensmag_lensmag_1h_ ((ptsz->has_lensmag_lensmag_1h == _TRUE_) && (index_md == ptsz->index_md_lensmag_lensmag_1h))
#define _lens_lensmag_2h_ ((ptsz->has_lens_lensmag_2h == _TRUE_) && (index_md == ptsz->index_md_lens_lensmag_2h))
#define _lens_lensmag_1h_ ((ptsz->has_lens_lensmag_1h == _TRUE_) && (index_md == ptsz->index_md_lens_lensmag_1h))
#define _lens_lensmag_hf_ ((ptsz->has_lens_lensmag_hf == _TRUE_) && (index_md == ptsz->index_md_lens_lensmag_hf))
#define _lens_lens_1h_ ((ptsz->has_lens_lens_1h == _TRUE_) && (index_md == ptsz->index_md_lens_lens_1h))
#define _lens_lens_2h_ ((ptsz->has_lens_lens_2h == _TRUE_) && (index_md == ptsz->index_md_lens_lens_2h))
#define _lens_lens_hf_ ((ptsz->has_lens_lens_hf == _TRUE_) && (index_md == ptsz->index_md_lens_lens_hf))
#define _tSZ_gal_1h_ ((ptsz->has_tSZ_gal_1h == _TRUE_) && (index_md == ptsz->index_md_tSZ_gal_1h))
#define _tSZ_gal_2h_ ((ptsz->has_tSZ_gal_2h == _TRUE_) && (index_md == ptsz->index_md_tSZ_gal_2h))
#define _tSZ_gallens_1h_ ((ptsz->has_tSZ_gallens_1h == _TRUE_) && (index_md == ptsz->index_md_tSZ_gallens_1h))
#define _tSZ_gallens_2h_ ((ptsz->has_tSZ_gallens_2h == _TRUE_) && (index_md == ptsz->index_md_tSZ_gallens_2h))
#define _tSZ_cib_1h_ ((ptsz->has_tSZ_cib_1h == _TRUE_) && (index_md == ptsz->index_md_tSZ_cib_1h))
#define _tSZ_cib_2h_ ((ptsz->has_tSZ_cib_2h == _TRUE_) && (index_md == ptsz->index_md_tSZ_cib_2h))
#define _gal_cib_1h_ ((ptsz->has_gal_cib_1h == _TRUE_) && (index_md == ptsz->index_md_gal_cib_1h))
#define _gal_cib_2h_ ((ptsz->has_gal_cib_2h == _TRUE_) && (index_md == ptsz->index_md_gal_cib_2h))
#define _gallens_cib_1h_ ((ptsz->has_gallens_cib_1h == _TRUE_) && (index_md == ptsz->index_md_gallens_cib_1h))
#define _gallens_cib_2h_ ((ptsz->has_gallens_cib_2h == _TRUE_) && (index_md == ptsz->index_md_gallens_cib_2h))
#define _ngal_ngal_1h_ ((ptsz->has_ngal_ngal_1h == _TRUE_) && (index_md == ptsz->index_md_ngal_ngal_1h))
#define _ngal_ngal_2h_ ((ptsz->has_ngal_ngal_2h == _TRUE_) && (index_md == ptsz->index_md_ngal_ngal_2h))
#define _ngal_ngal_hf_ ((ptsz->has_ngal_ngal_hf == _TRUE_) && (index_md == ptsz->index_md_ngal_ngal_hf))
#define _ngal_lens_1h_ ((ptsz->has_ngal_lens_1h == _TRUE_) && (index_md == ptsz->index_md_ngal_lens_1h))
#define _ngal_lens_2h_ ((ptsz->has_ngal_lens_2h == _TRUE_) && (index_md == ptsz->index_md_ngal_lens_2h))
#define _ngal_lens_hf_ ((ptsz->has_ngal_lens_hf == _TRUE_) && (index_md == ptsz->index_md_ngal_lens_hf))
#define _cib_cib_1h_ ((ptsz->has_cib_cib_1h == _TRUE_) && (index_md == ptsz->index_md_cib_cib_1h))
#define _cib_cib_2h_ ((ptsz->has_cib_cib_2h == _TRUE_) && (index_md == ptsz->index_md_cib_cib_2h))
#define _lens_cib_1h_ ((ptsz->has_lens_cib_1h == _TRUE_) && (index_md == ptsz->index_md_lens_cib_1h))
#define _lens_cib_2h_ ((ptsz->has_lens_cib_2h == _TRUE_) && (index_md == ptsz->index_md_lens_cib_2h))
#define _tSZ_lens_1h_ ((ptsz->has_tSZ_lens_1h == _TRUE_) && (index_md == ptsz->index_md_tSZ_lens_1h))
#define _tSZ_lens_2h_ ((ptsz->has_tSZ_lens_2h == _TRUE_) && (index_md == ptsz->index_md_tSZ_lens_2h))
#define _isw_lens_ ((ptsz->has_isw_lens == _TRUE_) && (index_md == ptsz->index_md_isw_lens))
#define _isw_tsz_ ((ptsz->has_isw_tsz == _TRUE_) && (index_md == ptsz->index_md_isw_tsz))
#define _isw_auto_ ((ptsz->has_isw_auto == _TRUE_) && (index_md == ptsz->index_md_isw_auto))
#define _dndlnM_ ((ptsz->has_dndlnM == _TRUE_) && (index_md == ptsz->index_md_dndlnM))
#define _szrates_ ((ptsz->has_sz_rates == _TRUE_) && (index_md == ptsz->index_md_szrates))
#define _kSZ_kSZ_1h_ ((ptsz->has_kSZ_kSZ_1h == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_1h))
#define _kSZ_kSZ_2h_ ((ptsz->has_kSZ_kSZ_2h == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_2h))
#define _kSZ_kSZ_tSZ_1h_ ((ptsz->has_kSZ_kSZ_tSZ_1h == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_tSZ_1h))
#define _kSZ_kSZ_tSZ_2h_ ((ptsz->has_kSZ_kSZ_tSZ_2h == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_tSZ_2h))
#define _kSZ_kSZ_tSZ_3h_ ((ptsz->has_kSZ_kSZ_tSZ_3h == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_tSZ_3h))
#define _tSZ_tSZ_tSZ_1halo_ ((ptsz->has_tSZ_tSZ_tSZ_1halo == _TRUE_) && (index_md == ptsz->index_md_tSZ_tSZ_tSZ_1halo))
#define _tSZ_tSZ_tSZ_2h_ ((ptsz->has_tSZ_tSZ_tSZ_2h == _TRUE_) && (index_md == ptsz->index_md_tSZ_tSZ_tSZ_2h))
#define _tSZ_tSZ_tSZ_3h_ ((ptsz->has_tSZ_tSZ_tSZ_3h == _TRUE_) && (index_md == ptsz->index_md_tSZ_tSZ_tSZ_3h))
//#define _tSZ_trispectrum_ ((ptsz->has_sz_trispec == _TRUE_))
//#define _tSZ_2halo_ ((ptsz->has_sz_2halo == _TRUE_))
//#define _tSZ_te_y_y_ ((ptsz->has_sz_te_y_y == _TRUE_))
// #define _cov_N_Cl_ ((ptsz->has_sz_cov_N_Cl == _TRUE_))






/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif


int class_sz_cosmo_init(struct background * pba,
                         struct thermo * pth,
                         struct perturbs * ppt,
                         struct nonlinear * pnl,
                         struct primordial * ppm,
                         struct spectra * psp,
                         struct lensing * ple,
                         struct tszspectrum * ptsz,
                         struct precision * ppr);

int szpowerspectrum_init(struct background * pba,
                         struct thermo * pth,
                         struct perturbs * ppt,
                         struct nonlinear * pnl,
                         struct primordial * ppm,
                         struct spectra * psp,
                         struct lensing * ple,
                         struct tszspectrum * ptsz,
                         struct precision * ppr);


  int szpowerspectrum_free(struct tszspectrum *ptsz);


  int compute_sz(struct background * pba,
                 struct nonlinear * pnl,
                 struct primordial * ppm,
                 struct perturbs * ppt,
                 struct tszspectrum * ptsz,
                 double * pvecback,
                 double * Pvectsz);



  //This evaluates the integrand which will be integrated
  //over M and then over z
  double integrand_at_m_and_z(double logM ,
                              double * pvecback,
                              double * pvectsz,
                              struct background * pba,
                              struct primordial * ppm,
                              struct nonlinear * pnl,
                              struct perturbs * ppt,
                              struct tszspectrum * ptsz);

  double delta_ell_lens_at_ell_and_z( double * pvecback,
                                  double * pvectsz,
                                  struct background * pba,
                                  struct primordial * ppm,
                                  struct nonlinear * pnl,
                                  struct tszspectrum * ptsz);

  double delta_ell_isw_at_ell_and_z( double * pvecback,
                                  double * pvectsz,
                                  struct background * pba,
                                  struct primordial * ppm,
                                  struct nonlinear * pnl,
                                  struct tszspectrum * ptsz);
  int do_mass_conversions(
                         double logM,
                         double z,
                         double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct tszspectrum * ptsz);
  int evaluate_HMF_at_logM_and_z(
                   double logM ,
                   double z,
                   double * pvecback,
                   double * pvectsz,
                   struct background * pba,
                   struct nonlinear * pnl,
                   struct tszspectrum * ptsz);

  int evaluate_completeness(double * pvecback,
                            double * pvectsz,
                            struct background * pba,
                            struct tszspectrum * ptsz);

  double evaluate_pressure_profile(double kl,
                                double * pvecback,
                                double * pvectsz,
                                struct background * pba,
                                struct tszspectrum * ptsz);

  int evaluate_tau_profile(double k,
                           double * pvecback,
                           double * pvectsz,
                           struct background * pba,
                           struct tszspectrum * ptsz);

  int evaluate_lensing_profile(double kl,
                               double m_delta,
                               double r_delta,
                               double c_delta,
                               double * pvecback,
                               double * pvectsz,
                               struct background * pba,
                               struct tszspectrum * ptsz);

double get_fstar_of_m_at_z(double m,
                           double z,
                           struct tszspectrum * ptsz);

double get_tau_profile_at_z_m_l(double z,
                                double m,
                                double k,
                                struct tszspectrum * ptsz,
                                struct background * pba);

double get_ksz_filter_at_l(double l,
                           struct tszspectrum * ptsz);

double get_M_min_of_z(double l,
                      struct tszspectrum * ptsz);

  int write_output_to_files_ell_indep_ints(struct nonlinear * pnl,
                                           struct background * pba,
                                           struct tszspectrum * ptsz);

  int write_output_to_files_cl(struct nonlinear * pnl,
                               struct background * pba,
                               struct primordial * ppm,
                               struct tszspectrum * ptsz);


  int show_preamble_messages(struct background * pba,
                             struct thermo * pth,
                             struct nonlinear * pnl,
                             struct primordial * ppm,
                             struct tszspectrum * ptsz);

  double gnu_tsz_of_nu_in_ghz(double nu_in_ghz,double Tcmb);

  int show_results(struct background * pba,
                   struct nonlinear * pnl,
                   struct primordial * ppm,
                   struct tszspectrum * ptsz);

  int select_multipole_array(struct tszspectrum * ptsz);

  int evaluate_halo_bias(double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct primordial * ppm,
                         struct nonlinear * pnl,
                         struct perturbs * ppt,
                         struct tszspectrum * ptsz);

  int evaluate_halo_bias_b2(double * pvecback,
                            double * pvectsz,
                            struct background * pba,
                            struct primordial * ppm,
                            struct nonlinear * pnl,
                            struct tszspectrum * ptsz);

 int evaluate_effective_galaxy_bias_ngal(int index_g,
                                         double * pvecback,
                                         double * pvectsz,
                                         struct background * pba,
                                         struct primordial * ppm,
                                         struct nonlinear * pnl,
                                         struct tszspectrum * ptsz);

 int evaluate_effective_galaxy_bias(double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct primordial * ppm,
                                    struct nonlinear * pnl,
                                    struct tszspectrum * ptsz);
double get_pk_lin_at_k_and_z(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct tszspectrum * ptsz);


double get2_pk_lin_at_k_and_z(//double * pvecback,//double * pvectsz,
  double * r,double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct tszspectrum * ptsz);
double get_pk_nonlin_at_k_and_z(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct tszspectrum * ptsz);

double get_pk_nonlin_at_k_and_z_fast(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct tszspectrum * ptsz);

 int evaluate_pk_at_ell_plus_one_half_over_chi(double * pvecback,
                                              double * pvectsz,
                                              struct background * pba,
                                              struct primordial * ppm,
                                              struct nonlinear * pnl,
                                              struct tszspectrum * ptsz);

 int evaluate_pk_at_ell_plus_one_half_over_chi_today(double * pvecback,
                                                      double * pvectsz,
                                                      struct background * pba,
                                                      struct primordial * ppm,
                                                      struct nonlinear * pnl,
                                                      struct tszspectrum * ptsz);


double evaluate_pk_halofit_over_pk_linear_at_ell_plus_one_half_over_chi(double * pvecback,
                                                                     double * pvectsz,
                                                                     struct background * pba,
                                                                     struct primordial * ppm,
                                                                     struct nonlinear * pnl,
                                                                     struct tszspectrum * ptsz);
int load_cl_ksz_template(struct tszspectrum * ptsz);

int load_nl_lensing_noise(struct tszspectrum * ptsz);


  int initialise_and_allocate_memory(struct tszspectrum * ptsz);


  int evaluate_temperature_mass_relation(double * pvecback,
                                         double * pvectsz,
                                         struct background * pba,
                                         struct tszspectrum * ptsz);

double evaluate_dlnMdeltadlnM(double logM,
                             double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct nonlinear * pnl,
                             struct tszspectrum * ptsz);

int evaluate_vrms2(double * pvecback,
                   double * pvectsz,
                   struct background * pba,
                   struct nonlinear * pnl,
                   struct tszspectrum * ptsz);


int evaluate_sigma2_hsv(double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct nonlinear * pnl,
                         struct tszspectrum * ptsz);

double integrand_mass(double xi, void *p);


int write_redshift_dependent_quantities(struct background * pba,
                                        struct tszspectrum * ptsz);


// int evaluate_tau_profile(double * pvecback,
//                         double * pvectsz,
//                         struct background * pba,
//                         struct tszspectrum * ptsz);



int tabulate_normalization_gas_density_profile(struct tszspectrum * ptsz,struct background * pba);

int tabulate_gas_pressure_profile_gNFW(struct background * pba,
                                   struct tszspectrum * ptsz);

int tabulate_gas_pressure_profile_B12(struct background * pba,
                                  struct tszspectrum * ptsz);

int tabulate_gas_pressure_profile_B12_fft(struct background * pba,
                                          struct tszspectrum * ptsz);

double evaluate_mean_galaxy_number_density_at_z(double z,
                                                struct tszspectrum * ptsz);

double evaluate_mean_galaxy_number_density_at_z_ngal(double z,
                                                     int index_g,
                                                     struct tszspectrum * ptsz);

double get_mean_galaxy_bias_at_z(double z,
                                 struct tszspectrum * ptsz);



double get_dyldzdlnm_at_l_z_and_m(double l,
                                  double z,
                                  double m,
                                  struct background * pba,
                                  struct nonlinear * pnl,
                                  struct tszspectrum * ptsz);

double HOD_mean_number_of_central_galaxies(double z,
                                           double M_halo,
                                           double M_min,
                                           double sigma_lnM,
                                           double f_cen,
                                           struct tszspectrum * ptsz,
                                           struct background * pba);

double HOD_mean_number_of_satellite_galaxies(double z,
                                             double M_halo,
                                             double Nc_mean,
                                             double M_min,
                                             double alpha_s,
                                             double M1_prime,
                                             struct tszspectrum * ptsz,
                                             struct background * pba);

double get_galaxy_profile_at_z_m_l_1h(double z,
                                      double m,
                                      double r_delta,
                                      double c_delta,
                                      double l,
                                      struct tszspectrum * ptsz,
                                      struct background * pba);


int evaluate_galaxy_profile_1h(double kl,
                               double m_delta,
                               double r_delta,
                               double c_delta,
                               double * pvecback,
                               double * pvectsz,
                               struct background * pba,
                               struct tszspectrum * ptsz);

int evaluate_galaxy_profile_ngal(double kl,
                               double m_delta,
                               double r_delta,
                               double c_delta,
                               double * pvecback,
                               double * pvectsz,
                               struct background * pba,
                               struct tszspectrum * ptsz);


int evaluate_galaxy_profile_2h(double kl,
                               double m_delta,
                               double r_delta,
                               double c_delta,
                               double * pvecback,
                               double * pvectsz,
                               struct background * pba,
                               struct tszspectrum * ptsz);

double get_truncated_nfw_profile_at_z_m_k_xout(//double * pvecback,
                                      double z,
                                      double m,
                                      double r_delta,
                                      double c_delta,
                                      double k,
                                      double xout,
                                      // double delta,
                                      struct background * pba,
                                      struct tszspectrum * ptsz);


double evaluate_truncated_nfw_profile(
                                   double z,
                                   double k,
                                   double r_delta,
                                   double c_delta,
                                   double xout);



int evaluate_c200m_D08(double * pvecback,
                        double * pvectsz,
                        struct background * pba,
                        struct tszspectrum * ptsz);

double get_c200c_at_m_and_z_D08(double M,
                                double z);

double get_c200c_at_m_and_z(//double * pvecback,
                        double m,
                        double z,
                        struct background * pba,
                        struct tszspectrum * ptsz);

double get_c500c_at_m_and_z(//double * pvecback,
                        double m,
                        double z,
                        struct background * pba,
                        struct tszspectrum * ptsz);

double get_galaxy_number_counts(double z,
                                struct tszspectrum * ptsz);


double get_f_of_sigma_at_m_and_z(double m,
                                 double z,
                                 struct background * pba,
                                 struct nonlinear * pnl,
                                 struct tszspectrum * ptsz);


double get_source_galaxy_number_counts(double z,
                                struct tszspectrum * ptsz);
double radial_kernel_W_galaxy_at_z( double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct tszspectrum * ptsz);

double radial_kernel_W_galaxy_ngal_at_z(  int index_g,
                                          double * pvecback,
                                          double * pvectsz,
                                          struct background * pba,
                                          struct tszspectrum * ptsz);

double radial_kernel_W_lensing_at_z(double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct primordial * ppm,
                                    struct nonlinear * pnl,
                                    struct tszspectrum * ptsz);

double radial_kernel_W_lensing_magnification_at_z(double * pvecback,
                                                  double * pvectsz,
                                                  struct background * pba,
                                                  struct primordial * ppm,
                                                  struct nonlinear * pnl,
                                                  struct tszspectrum * ptsz);

double radial_kernel_W_cmb_lensing_at_z(double z,
                                        double * pvectsz,
                                        struct background * pba,
                                        struct tszspectrum * ptsz);

double radial_kernel_W_galaxy_lensing_at_z(double z,
                                           // double * pvectsz,
                                           // struct background * pba,
                                           struct tszspectrum * ptsz);

double radial_kernel_W_galaxy_lensing_magnification_at_z(double z,
                                                         double * pvectsz,
                                                         struct background * pba,
                                                         struct tszspectrum * ptsz);


double evaluate_galaxy_number_counts( double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct tszspectrum * ptsz);

double evaluate_galaxy_number_counts_fdndz( double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct tszspectrum * ptsz);
double evaluate_unwise_m_min_cut(double z,
                                 int sample_id,
                                 struct tszspectrum * ptsz);


int evaluate_cib_profile(double m_delta,
                         double r_delta,
                         double c_delta,
                         double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct tszspectrum * ptsz);

double Luminosity_of_central_galaxies(double z,
                                      double  M_halo,
                                      double nu,
                                      double * pvectsz,
                                      struct tszspectrum * ptsz,
                                      struct background * pba);

double Luminosity_of_satellite_galaxies(double z,
                                        double  M_halo,
                                        double nu,
                                        struct tszspectrum * ptsz,
                                        struct background * pba);

double maniyar_cib_Mdot(double M, double z, struct tszspectrum * ptsz);
double evaluate_Sigma_cib(double M, struct tszspectrum * ptsz);
double evaluate_phi_cib(double z, struct tszspectrum * ptsz);
double evaluate_sed_cib(double z, double nu, struct tszspectrum * ptsz);
double evaluate_dust_temperature(double z, struct tszspectrum * ptsz);
double evaluate_galaxy_luminosity(double z, double M, double nu, struct tszspectrum * ptsz);
double subhalo_hmf_dndlnMs(double M_host,double M_sub,struct tszspectrum * ptsz);

double integrand_kSZ2_X_at_theta(double ell_prime, void *p);
double integrand_kSZ2_X(double theta, void *p);

double integrand_kSZ2_X_lensing_term_at_theta(double ell_prime, void *p);
double integrand_kSZ2_X_lensing_term(double theta, void *p);

int evaluate_matter_density_profile(
                             double k,
                             double r_delta,
                             double c_delta,
                             double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct tszspectrum * ptsz);
double get_matter_bispectrum_at_z_tree_level_PT(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);

double get_te_of_m500c_at_z_arnaud(double m, double z, struct background * pba,struct tszspectrum * ptsz);
double get_te_of_m500c_at_z_lee(double m, double z, struct background * pba,struct tszspectrum * ptsz);


int  evaluate_ttg_bispectrum_at_z_tree_level_PT(double * r,
                                                      double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);


int  evaluate_ttg_bispectrum_at_z_effective_approach(double * r,
                                                      double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);

double get_ttg_bispectrum_at_z_effective_approach(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);


double get_ttg_bispectrum_at_z_tree_level_PT(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);

double get_matter_bispectrum_at_z_effective_approach_smoothed(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);

double get_matter_bispectrum_at_z_effective_approach(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);
double get_matter_bispectrum_at_z_effective_approach_SC(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm);

double bispectrum_f2_kernel(double k,
                            double k_prime,
                            double k_prime_prime);
double bispectrum_f2_kernel_eff_SC(double k1,
                            double k2,
                            double k3,
                            double n1,
                            double n2,
                            double sig8_at_z,
                            double knl);

double bispectrum_f2_kernel_eff_a_SC(double k1,double n1,double sig8_at_z,double knl);
double bispectrum_f2_kernel_eff_b_SC(double k1,double n1,double knl);
double bispectrum_f2_kernel_eff_c_SC(double k1,double n1,double knl);
double bispectrum_f2_kernel_eff(double k1,
                            double k2,
                            double k3,
                            double n1,
                            double n2,
                            double sig8_at_z,
                            double knl);

double bispectrum_f2_kernel_eff_a(double k1,double n1,double sig8_at_z,double knl);
double bispectrum_f2_kernel_eff_b(double k1,double n1,double knl);
double bispectrum_f2_kernel_eff_c(double k1,double n1,double knl);
double bispectrum_f2_kernel_eff_Q3(double n1);

double get_rho_crit_at_z(double z_asked,
                         struct background * pba,
                         struct tszspectrum * ptsz);

double get_c200m_at_m_and_z(double M,
                            double z,
                            struct background * pba,
                            struct tszspectrum * ptsz);

double get_c200m_at_m_and_z_D08(double M,
                                double z);

double get_c200m_at_m_and_z_B13(double M,
                                double z,
                                struct background * pba,
                                struct tszspectrum * ptsz);

double get_c200c_at_m_and_z_B13(double M,
                                double z,
                                struct background * pba,
                                struct tszspectrum * ptsz);
double get_gas_profile_at_x_M_z_nfw_200c(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct tszspectrum * ptsz);

double get_rvir_of_m200c_at_z(
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct tszspectrum * ptsz);

double get_gas_profile_at_x_M_z_nfw_200m(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct tszspectrum * ptsz);

double get_gas_profile_at_x_M_z_bcm_200c(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct tszspectrum * ptsz);


double get_gas_profile_at_x_M_z_b16_200c(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         double c,
                                         double A_rho,
                                         double A_alpha,
                                         double A_beta,
                                         double alpha_m_rho0,
                                         double alpha_m_alpha,
                                         double alpha_m_beta,
                                         double alpha_z_rho0,
                                         double alpha_z_alpha,
                                         double alpha_z_beta,
                                         // break model param
		                                     double mcut,
		                                     double alphap_m_rho0,
                                         double alphap_m_alpha,
                                         double alphap_m_beta,
		                                     double alpha_c_rho0,
                                         double alpha_c_alpha,
                                         double alpha_c_beta,
                                         // end break model param
                                         double gamma,
                                         double xc,
                                         struct background * pba,
                                         struct tszspectrum * ptsz);


double get_second_order_bias_at_z_and_nu(double z,
                                         double nu,
                                         struct tszspectrum * ptsz,
                                         struct background * pba);

double get_first_order_bias_at_z_and_nu(double z,
                                         double nu,
                                         struct tszspectrum * ptsz);

double get_ng_bias_contribution_at_z_and_k(double z,
                                           double k,
                                           double bh,
                                           struct background * pba,
                                           struct perturbs * ppt,
                                           struct tszspectrum * ptsz);

double get_scale_dependent_bias_at_z_and_k(double z,
                                           double k,
                                           double bh,
                                           struct tszspectrum *ptsz);


double get_vrms2_at_z(double z,
                      struct tszspectrum * ptsz);

double get_dlnsigma_dlnR_at_z_and_m(double z,
                                    double m,
                                    struct tszspectrum * ptsz,
                                    struct background * pba);
double get_sigma_at_z_and_m(double z,
                            double m,
                            struct tszspectrum * ptsz,
                            struct background * pba);
double get_sigma8_at_z(double z,
                      struct tszspectrum * ptsz,
                      struct background * pba);
double get_nu_at_z_and_m(double z,
                         double m,
                         struct tszspectrum * ptsz,
                         struct background * pba);

// this is r_200c*P_200c
double get_1e6xdy_from_battaglia_pressure_at_x_z_and_m200c(double z,
                                                           double m,
                                                           double x,
                                                           struct background * pba,
                                                           struct tszspectrum * ptsz);


double get_1e6xdy_from_gnfw_pressure_at_x_z_and_m500c(double z,
                                                      double m,
                                                      double x,
                                                      double d,
                                                      struct background * pba,
                                                      struct tszspectrum * ptsz);

double get_pressure_P_over_P_delta_at_x_M_z_b12_200c(double x_asked,
                                                     double m_asked,
                                                     double z_asked,
                                                     double c_asked,
                                                     double A_P0,
                                                     double A_xc,
                                                     double A_beta,
                                                     double alpha_m_P0,
                                                     double alpha_m_xc,
                                                     double alpha_m_beta,
                                                     double alpha_z_P0,
                                                     double alpha_z_xc,
                                                     double alpha_z_beta,
                                  							     double mcut,
                                  							     double alphap_m_P0,
                                  							     double alphap_m_xc,
                                  							     double alphap_m_beta,
                                  							     double alpha_c_P0,
                                  							     double alpha_c_xc,
                                  							     double alpha_c_beta,
                                                     double alpha,
                                                     double gamma,
                                                     struct background * pba,
                                                     struct tszspectrum * tsz);

double get_pressure_P_over_P_delta_at_x_gnfw_500c(double x_asked,
                                                      double P0GNFW,
                                                      double alphaGNFW,
                                                      double betaGNFW,
                                                      double gammaGNFW,
                                                      double c500,
                                                      struct background * pba,
                                                      struct tszspectrum * tsz);


struct Parameters_for_integrand_gas_density_profile_2h{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  struct perturbs * ppt;
  double * pvecback;
  double * pvectsz;
  double z;
  double k;
};


struct Parameters_for_integrand_gas_pressure_profile_2h{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  struct perturbs * ppt;
  double * pvecback;
  double * pvectsz;
  double z;
  double k;
};


struct Parameters_for_integrand_kSZ2_X_at_theta{
struct nonlinear * pnl;
struct primordial * ppm;
struct tszspectrum * ptsz;
struct background * pba;
double * Pvecback;
double * Pvectsz;
double theta;
int index_ell_3;
double * b_l1_l2_l_1d;
double * ln_ell;
};



struct Parameters_for_integrand_kSZ2_X{
struct nonlinear * pnl;
struct primordial * ppm;
struct tszspectrum * ptsz;
struct background * pba;
double * Pvecback;
double * Pvectsz;
int index_ell_3;
double * b_l1_l2_l_1d;
double * ln_ell;
};


struct Parameters_for_integrand_kSZ2_X_lensing_term_at_theta{
struct nonlinear * pnl;
struct primordial * ppm;
struct tszspectrum * ptsz;
struct background * pba;
// double * Pvecback;
// double * Pvectsz;
double theta;
int index_ell;
double * integrand_l_lprime_phi;
double * ln_ellprime;
};



struct Parameters_for_integrand_kSZ2_X_lensing_term{
struct nonlinear * pnl;
struct primordial * ppm;
struct tszspectrum * ptsz;
struct background * pba;
// double * Pvecback;
// double * Pvectsz;
int index_ell;
double * integrand_l_lprime_phi;
double * ln_ellprime;
};





#ifdef __cplusplus
}
#endif

#endif
