/** @file nonlinear.h Documented includes for trg module */
#include <unistd.h>
#include <errno.h>
#include "primordial.h"

#include "nonlinear.h"
// #include "class_sz.h"
// #include "class_sz_tools.h"
// #include "Patterson.h"
// #include "r8lib.h"
// #include "fft.h"
# include <fftw3.h>
#include <complex.h>
#undef I

#ifndef __NONLINEAR_PT__
#define __NONLINEAR_PT__

/**
 * Maximum number of values of redshift at which the spectra will be
 * written in output files
 */

// #define _Z_PK_NUM_MAX_ 100


enum non_linear_method_pt {nlpt_none,nlpt_spt};
enum non_linear_irresumm_pt {irres_yes,irres_no};
enum non_linear_bias_pt {bias_yes,bias_no};
enum non_linear_rsd_pt {rsd_yes,rsd_no};
// enum non_linear_cb_pt {cb_yes,cb_no};
// enum non_linear_fftlogn {fftlog_fast,fftlog_norm};
enum non_linear_AP_effect_pt {AP_effect_yes,AP_effect_no};
enum non_linear_fNL_equil_ortho {fNL_equil_ortho_yes,fNL_equil_ortho_no};

//enum halofit_integral_type {halofit_integral_one, halofit_integral_two, halofit_integral_three};
//enum halofit_statement {ok, too_small};

/**
 * Structure containing all information on non-linear spectra.
 *
 * Once initialized by nonlinear_init(), contains a table for all two points correlation functions
 * and for all the ai,bj functions (containing the three points correlation functions), for each
 * time and wave-number.
 */


struct tszspectrum {



  #define __ALLOCATE_TSZ_PARAMETER__
  #include "class_sz_precisions.h"
  #undef __ALLOCATE_TSZ_PARAMETER__


  fftw_plan forward_plan, reverse_plan;
  fftw_plan forward_plan_counts_fft, reverse_plan_counts_fft;
  // int N_samp_fftw;


  int use_analytical_truncated_nfw;
  int use_hod; // Eq. 15 or 16 of KA20
  int unwise_galaxy_sample_id;
  int galaxy_sample;


  int no_b2;
  //double unwise_m_min_cut;

  int cib_nu0_norm;
  double sn_cutoff;

  double f_sky;
  double fsky_from_skyfracs;
  double szunbinned_loglike;
  double Omega_survey;

  double chi_star; //comoving distance to the surface of last scattering [Mpc/h]

  double hmf_int;
  double y_monopole;
  double * cib_monopole;
  double * cib_shotnoise;
  double * pk_at_z_1h;
  double * pk_at_z_2h;
  double * pk_gg_at_z_1h;
  double * pk_gg_at_z_2h;
  double * pk_bb_at_z_1h;
  double * pk_bb_at_z_2h;
  double * pk_b_at_z_2h;
  double * pk_em_at_z_1h;
  double * pk_em_at_z_2h;
  double * pk_HI_at_z_1h;
  double * pk_HI_at_z_2h;
  double * bk_at_z_1h;
  double * bk_at_z_2h;
  double * bk_at_z_3h;
  double * bk_ttg_at_z_1h;
  double * bk_ttg_at_z_2h;
  double * bk_ttg_at_z_3h;
  double * cl_sz_1h;
  double * cl_gal_gal_1h;
  double * cl_gal_gal_2h;
  double * cl_gal_gal_hf;
  double * cl_gal_lens_hf;
  double * cl_lens_lens_hf;
  double * cl_gal_lens_2h;
  double * cl_gal_lens_1h;
  double * cl_gal_lensmag_hf;
  double * cl_gal_lensmag_2h;
  double * cl_gal_gallens_1h;
  double * cl_gal_gallens_2h;
  double * cl_gallens_gallens_1h;
  double * cl_gallens_gallens_2h;
  double * cl_gallens_lens_1h;
  double * cl_gallens_lens_2h;
  double * thetas_arcmin;
  double * gamma_gal_gallens_1h;
  double * gamma_gal_gallens_2h;
  double * cl_gal_lensmag_1h;
  double * cl_tSZ_lensmag_2h;
  double * cl_tSZ_lensmag_1h;
  double * cl_lensmag_lensmag_hf;
  double * cl_lensmag_lensmag_2h;
  double * cl_lensmag_lensmag_1h;
  double * cl_lens_lensmag_hf;
  double * cl_lens_lensmag_2h;
  double * cl_lens_lensmag_1h;
  double * cl_lens_lens_1h;
  double * cl_lens_lens_2h;
  double * cl_tSZ_gal_1h;
  double * cl_tSZ_gal_2h;
  double * cl_tSZ_gallens_1h;
  double * cl_tSZ_gallens_2h;
  double *** cl_ngal_ngal_1h;
  double *** cl_ngal_ngal_2h;
  double *** cl_ngal_ngal_hf;
  double ** cl_ngal_lens_1h;
  double ** cl_ngal_lens_2h;
  double ** cl_ngal_lens_hf;
  double *** cl_cib_cib_1h;
  double *** cl_cib_cib_2h;
  double ** cl_tSZ_cib_1h;
  double ** cl_tSZ_cib_2h;
  double ** cl_gal_cib_1h;
  double ** cl_gal_cib_2h;
  double ** cl_gallens_cib_1h;
  double ** cl_gallens_cib_2h;
  double ** cl_lens_cib_1h;
  double ** cl_lens_cib_2h;
  double * cl_tSZ_lens_1h;
  double * cl_tSZ_lens_2h;
  double * szrate;
  double * cl_isw_lens;
  double * cl_isw_tsz;
  double * cl_isw_auto;
  double * cov_ll_kSZ_kSZ_gal;
  double * cl_t2t2f;
  double * cl_kSZ_kSZ_gal_lensing_term;
  double * cl_kSZ_kSZ_gal_1h;
  double * cl_kSZ_kSZ_gal_1h_fft;
  double * cl_kSZ_kSZ_gal_2h_fft;
  double * cl_kSZ_kSZ_gal_3h_fft;
  double * cl_kSZ_kSZ_gallens_1h_fft;
  double * cl_kSZ_kSZ_gallens_2h_fft;
  double * cl_kSZ_kSZ_gallens_3h_fft;
  double * cl_kSZ_kSZ_gallens_hf;
  double * cov_ll_kSZ_kSZ_gallens;
  double * cl_kSZ_kSZ_gallens_lensing_term;
  double * cl_kSZ_kSZ_lens_1h_fft;
  double * cl_kSZ_kSZ_lens_2h_fft;
  double * cl_kSZ_kSZ_lens_3h_fft;
  double * cl_kSZ_kSZ_lens_hf;
  double * cov_ll_kSZ_kSZ_lens;
  double * cl_kSZ_kSZ_lens_lensing_term;
  double * cl_kSZ_kSZ_gal_2h;
  double * cl_kSZ_kSZ_gal_3h;
  double * cl_kSZ_kSZ_gal_hf;
  double * cl_kSZ_kSZ_lensmag_1h;
  double * b_tSZ_tSZ_tSZ_1halo;
  double * b_tSZ_tSZ_tSZ_2h;
  double * b_tSZ_tSZ_tSZ_3h;
  double * cl_kSZ_kSZ_1h;
  double * cl_kSZ_kSZ_2h;
  double * b_kSZ_kSZ_tSZ_1h;
  double * b_kSZ_kSZ_tSZ_2h;
  double * b_kSZ_kSZ_tSZ_3h;
  double * cl_te_y_y;
  double * m_y_y_1h;
  double * m_y_y_2h;
  double ** tllprime_sz;
  double ** dndlnM_at_z_and_M;
  double * cl_sz_2h;
  double ** cov_N_cl;
  double ** r_N_cl;
  double ** cov_Y_N;
  double ** cov_Y_N_next_order;
  double ** cov_Y_Y_ssc;
  double * cov_N_N;
  double ** cov_N_N_hsv;
  double ** r_Y_N;
  double ** r_cl_clp;
  double ** trispectrum_ref;
  double * cov_cl_cl;
  double * sig_cl_squared_binned;

  int delta_def_galaxies;
  int delta_def_cib;
  int delta_def_matter_density;
  int delta_def_electron_pressure;
  int delta_def_electron_density;


  int delta_def_HI_pressure;
  int delta_def_HI_density;

  int bispec_conf_id;

  double M_min_ng_bar;
  double M_max_ng_bar;


  int index_d_tot;
  int index_phi;
  int index_psi;
  int number_of_titles;

  int need_m200m_to_m200c;
  int need_m200c_to_m200m;
  int need_m200m_to_m500c;
  int need_hmf;
  int need_sigma;
  int need_m200c_to_m500c;
  int need_m500c_to_m200c;

  int need_ng_bias;
  int nz_ng_bias;
  int nk_ng_bias;
  double * array_ln_1pz_ng_bias;
  double * array_ln_k_ng_bias;
  double * array_ln_ng_bias_at_z_and_k;

  double * array_ln_density_norm_at_z_and_m;

  int need_ksz_template;
  int need_tt_noise;
  int need_lensing_noise;


  int integrate_wrt_mvir;
  int integrate_wrt_m500c;
  int integrate_wrt_m200m;
  int integrate_wrt_m200c;

  int has_electron_pressure;
  int has_electron_density;
  int has_HI_density;
  int has_galaxy;
  int has_matter_density;
  int has_lensing;
  int has_cib;
  int has_dcib0dz;
  int has_dydz;
  int has_isw;

  int has_vir;
  int has_500c;
  int has_200m;
  int has_200c;


  int index_integrate_wrt_mvir;
  int index_integrate_wrt_m500c;
  int index_integrate_wrt_m200m;

  int index_has_electron_pressure;
  int index_has_electron_density;
  int index_has_HI_density;
  int index_has_galaxy;
  int index_has_matter_density;
  int index_has_lensing;
  int index_has_cib;
  int index_has_isw;

  int index_has_vir;
  int index_has_500c;
  int index_has_200m;
  int index_has_200c;

  int index_md;

  int has_sz_counts;
  int has_sz_counts_fft;

  int create_ref_trispectrum_for_cobaya;


  int use_m500c_in_ym_relation;
  //int has_sz_te_y_y;
  int has_sz_cov_N_Cl;

  int has_sz_cov_Y_N;
  int index_md_cov_Y_N;
  int index_integrand_id_cov_Y_N_first;
  int index_integrand_id_cov_Y_N_last;

  int has_sz_cov_Y_N_next_order;
  int index_md_cov_Y_N_next_order;
  int index_integrand_id_cov_Y_N_next_order_first;
  int index_integrand_id_cov_Y_N_next_order_last;

  int has_sz_cov_Y_Y_ssc;
  int index_md_cov_Y_Y_ssc;
  int index_integrand_id_cov_Y_Y_ssc_first;
  int index_integrand_id_cov_Y_Y_ssc_last;

  int has_sz_cov_N_N;
  int index_md_cov_N_N;
  int index_integrand_id_cov_N_N_first;
  int index_integrand_id_cov_N_N_last;

  int has_sz_cov_N_N_hsv;
  int index_md_cov_N_N_hsv;
  int index_integrand_id_cov_N_N_hsv_first;
  int index_integrand_id_cov_N_N_hsv_last;


  int has_sz_rates;
  int index_md_szrates;
  int index_integrand_id_szrates_first;
  int index_integrand_id_szrates_last;

  int has_hmf;
  int index_md_hmf;
  int index_integrand_id_hmf;

  int has_pk_bb_at_z_1h;
  int index_md_pk_bb_at_z_1h;
  int index_integrand_id_pk_bb_at_z_1h_first;
  int index_integrand_id_pk_bb_at_z_1h_last;

  int has_pk_bb_at_z_2h;
  int index_md_pk_bb_at_z_2h;
  int index_integrand_id_pk_bb_at_z_2h_first;
  int index_integrand_id_pk_bb_at_z_2h_last;


  int has_gas_pressure_profile_2h;
  int has_gas_density_profile_2h;

  int has_pk_b_at_z_2h;
  int index_md_pk_b_at_z_2h;
  int index_integrand_id_pk_b_at_z_2h_first;
  int index_integrand_id_pk_b_at_z_2h_last;

  int has_pk_em_at_z_1h;
  int index_md_pk_em_at_z_1h;
  int index_integrand_id_pk_em_at_z_1h_first;
  int index_integrand_id_pk_em_at_z_1h_last;

  int has_pk_em_at_z_2h;
  int index_md_pk_em_at_z_2h;
  int index_integrand_id_pk_em_at_z_2h_first;
  int index_integrand_id_pk_em_at_z_2h_last;


  int has_pk_HI_at_z_1h;
  int index_md_pk_HI_at_z_1h;
  int index_integrand_id_pk_HI_at_z_1h_first;
  int index_integrand_id_pk_HI_at_z_1h_last;

  int has_pk_HI_at_z_2h;
  int index_md_pk_HI_at_z_2h;
  int index_integrand_id_pk_HI_at_z_2h_first;
  int index_integrand_id_pk_HI_at_z_2h_last;



  int has_pk_gg_at_z_1h;
  int index_md_pk_gg_at_z_1h;
  int index_integrand_id_pk_gg_at_z_1h_first;
  int index_integrand_id_pk_gg_at_z_1h_last;

  int has_pk_gg_at_z_2h;
  int index_md_pk_gg_at_z_2h;
  int index_integrand_id_pk_gg_at_z_2h_first;
  int index_integrand_id_pk_gg_at_z_2h_last;

  int has_pk_at_z_1h;
  int index_md_pk_at_z_1h;
  int index_integrand_id_pk_at_z_1h_first;
  int index_integrand_id_pk_at_z_1h_last;

  int has_pk_at_z_2h;
  int index_md_pk_at_z_2h;
  int index_integrand_id_pk_at_z_2h_first;
  int index_integrand_id_pk_at_z_2h_last;

  int has_bk_at_z_1h;
  int index_md_bk_at_z_1h;
  int index_integrand_id_bk_at_z_1h_first;
  int index_integrand_id_bk_at_z_1h_last;

  int has_bk_at_z_2h;
  int index_md_bk_at_z_2h;
  int index_integrand_id_bk_at_z_2h_first;
  int index_integrand_id_bk_at_z_2h_last;

  int has_bk_at_z_3h;
  int index_md_bk_at_z_3h;
  int index_integrand_id_bk_at_z_3h_first;
  int index_integrand_id_bk_at_z_3h_last;

  int has_bk_ttg_at_z_1h;
  int index_md_bk_ttg_at_z_1h;
  int index_integrand_id_bk_ttg_at_z_1h_first;
  int index_integrand_id_bk_ttg_at_z_1h_last;

  int has_bk_ttg_at_z_2h;
  int index_md_bk_ttg_at_z_2h;
  int index_integrand_id_bk_ttg_at_z_2h_first;
  int index_integrand_id_bk_ttg_at_z_2h_last;

  int has_bk_ttg_at_z_3h;
  int index_md_bk_ttg_at_z_3h;
  int index_integrand_id_bk_ttg_at_z_3h_first;
  int index_integrand_id_bk_ttg_at_z_3h_last;



  int has_bk_at_z_hf;
  int has_bk_ttg_at_z_hf;

  int has_mean_galaxy_bias;
  int has_ng_in_bh;
  double fNL;

  int has_mean_y;
  int index_md_mean_y;
  int index_integrand_id_mean_y;

  int has_sz_ps;
  int index_md_sz_ps;
  int index_integrand_id_sz_ps_first;
  int index_integrand_id_sz_ps_last;

  int has_sz_2halo;
  int index_md_2halo;
  int index_integrand_id_sz_ps_2halo_first;
  int index_integrand_id_sz_ps_2halo_last;

  int has_sz_te_y_y;
  int index_md_te_y_y;
  int index_integrand_id_sz_ps_te_y_y_first;
  int index_integrand_id_sz_ps_te_y_y_last;

  int has_sz_m_y_y_1h;
  int index_md_m_y_y_1h;
  int index_integrand_id_sz_ps_m_y_y_1h_first;
  int index_integrand_id_sz_ps_m_y_y_1h_last;

  int has_sz_m_y_y_2h;
  int index_md_m_y_y_2h;
  int index_integrand_id_sz_ps_m_y_y_2h_first;
  int index_integrand_id_sz_ps_m_y_y_2h_last;

  int has_kSZ_kSZ_gal_1h;
  int index_md_kSZ_kSZ_gal_1h;
  int index_integrand_id_kSZ_kSZ_gal_1h_first;
  int index_integrand_id_kSZ_kSZ_gal_1h_last;

  int has_kSZ_kSZ_gal_1h_fft;
  int index_md_kSZ_kSZ_gal_1h_fft;
  int index_integrand_id_kSZ_kSZ_gal_1h_fft_first;
  int index_integrand_id_kSZ_kSZ_gal_1h_fft_last;

  int has_kSZ_kSZ_gal_2h_fft;
  int index_md_kSZ_kSZ_gal_2h_fft;
  int index_integrand_id_kSZ_kSZ_gal_2h_fft_first;
  int index_integrand_id_kSZ_kSZ_gal_2h_fft_last;

  int has_kSZ_kSZ_gal_3h_fft;
  int index_md_kSZ_kSZ_gal_3h_fft;
  int index_integrand_id_kSZ_kSZ_gal_3h_fft_first;
  int index_integrand_id_kSZ_kSZ_gal_3h_fft_last;


  int has_kSZ_kSZ_gal_2h;
  int index_md_kSZ_kSZ_gal_2h;
  int index_integrand_id_kSZ_kSZ_gal_2h_first;
  int index_integrand_id_kSZ_kSZ_gal_2h_last;

  int has_kSZ_kSZ_gal_3h;
  int index_md_kSZ_kSZ_gal_3h;
  int index_integrand_id_kSZ_kSZ_gal_3h_first;
  int index_integrand_id_kSZ_kSZ_gal_3h_last;

  int has_kSZ_kSZ_gal_covmat;
  int has_kSZ_kSZ_gallens_covmat;
  int has_kSZ_kSZ_lens_covmat;
  int has_kSZ_kSZ_gal_lensing_term;
  int has_kSZ_kSZ_gallens_lensing_term;
  int has_kSZ_kSZ_lens_lensing_term;

  int has_kSZ_kSZ_gal_hf;
  int index_md_kSZ_kSZ_gal_hf;
  int index_integrand_id_kSZ_kSZ_gal_hf_first;
  int index_integrand_id_kSZ_kSZ_gal_hf_last;

  int has_kSZ_kSZ_gallens_1h_fft;
  int index_md_kSZ_kSZ_gallens_1h_fft;
  int index_integrand_id_kSZ_kSZ_gallens_1h_fft_first;
  int index_integrand_id_kSZ_kSZ_gallens_1h_fft_last;

  int has_kSZ_kSZ_gallens_2h_fft;
  int index_md_kSZ_kSZ_gallens_2h_fft;
  int index_integrand_id_kSZ_kSZ_gallens_2h_fft_first;
  int index_integrand_id_kSZ_kSZ_gallens_2h_fft_last;

  int has_kSZ_kSZ_gallens_3h_fft;
  int index_md_kSZ_kSZ_gallens_3h_fft;
  int index_integrand_id_kSZ_kSZ_gallens_3h_fft_first;
  int index_integrand_id_kSZ_kSZ_gallens_3h_fft_last;

  int has_kSZ_kSZ_gallens_hf;
  int index_md_kSZ_kSZ_gallens_hf;
  int index_integrand_id_kSZ_kSZ_gallens_hf_first;
  int index_integrand_id_kSZ_kSZ_gallens_hf_last;

  int has_kSZ_kSZ_lens_1h_fft;
  int index_md_kSZ_kSZ_lens_1h_fft;
  int index_integrand_id_kSZ_kSZ_lens_1h_fft_first;
  int index_integrand_id_kSZ_kSZ_lens_1h_fft_last;

  int has_kSZ_kSZ_lens_2h_fft;
  int index_md_kSZ_kSZ_lens_2h_fft;
  int index_integrand_id_kSZ_kSZ_lens_2h_fft_first;
  int index_integrand_id_kSZ_kSZ_lens_2h_fft_last;

  int has_kSZ_kSZ_lens_3h_fft;
  int index_md_kSZ_kSZ_lens_3h_fft;
  int index_integrand_id_kSZ_kSZ_lens_3h_fft_first;
  int index_integrand_id_kSZ_kSZ_lens_3h_fft_last;

  int has_kSZ_kSZ_lens_hf;
  int index_md_kSZ_kSZ_lens_hf;
  int index_integrand_id_kSZ_kSZ_lens_hf_first;
  int index_integrand_id_kSZ_kSZ_lens_hf_last;


  int has_kSZ_kSZ_lensmag_1halo;
  int index_md_kSZ_kSZ_lensmag_1halo;
  int index_integrand_id_kSZ_kSZ_lensmag_1halo_first;
  int index_integrand_id_kSZ_kSZ_lensmag_1halo_last;




  int has_kSZ_kSZ_1h;
  int index_md_kSZ_kSZ_1h;
  int index_integrand_id_kSZ_kSZ_1h_first;
  int index_integrand_id_kSZ_kSZ_1h_last;

  int has_kSZ_kSZ_2h;
  int index_md_kSZ_kSZ_2h;
  int index_integrand_id_kSZ_kSZ_2h_first;
  int index_integrand_id_kSZ_kSZ_2h_last;


  int has_kSZ_kSZ_tSZ_1h;
  int index_md_kSZ_kSZ_tSZ_1h;
  int index_integrand_id_kSZ_kSZ_tSZ_1h_first;
  int index_integrand_id_kSZ_kSZ_tSZ_1h_last;

  int has_kSZ_kSZ_tSZ_2h;
  int index_md_kSZ_kSZ_tSZ_2h;
  int index_integrand_id_kSZ_kSZ_tSZ_2h_first;
  int index_integrand_id_kSZ_kSZ_tSZ_2h_last;

  int has_kSZ_kSZ_tSZ_3h;
  int index_md_kSZ_kSZ_tSZ_3h;
  int index_integrand_id_kSZ_kSZ_tSZ_3h_first;
  int index_integrand_id_kSZ_kSZ_tSZ_3h_last;

  int has_tSZ_tSZ_tSZ_1halo;
  int index_md_tSZ_tSZ_tSZ_1halo;
  int index_integrand_id_tSZ_tSZ_tSZ_1halo_first;
  int index_integrand_id_tSZ_tSZ_tSZ_1halo_last;

  int has_tSZ_tSZ_tSZ_2h;
  int index_md_tSZ_tSZ_tSZ_2h;
  int index_integrand_id_tSZ_tSZ_tSZ_2h_first;
  int index_integrand_id_tSZ_tSZ_tSZ_2h_last;

  int has_tSZ_tSZ_tSZ_3h;
  int index_md_tSZ_tSZ_tSZ_3h;
  int index_integrand_id_tSZ_tSZ_tSZ_3h_first;
  int index_integrand_id_tSZ_tSZ_tSZ_3h_last;



  int has_tSZ_lens_1h;
  int index_md_tSZ_lens_1h;
  int index_integrand_id_tSZ_lens_1h_first;
  int index_integrand_id_tSZ_lens_1h_last;

  int has_gal_gal_1h;
  int index_md_gal_gal_1h;
  int index_integrand_id_gal_gal_1h_first;
  int index_integrand_id_gal_gal_1h_last;

  int has_gal_gal_2h;
  int index_md_gal_gal_2h;
  int index_integrand_id_gal_gal_2h_first;
  int index_integrand_id_gal_gal_2h_last;

  int has_gal_gal_hf;
  int index_md_gal_gal_hf;
  int index_integrand_id_gal_gal_hf_first;
  int index_integrand_id_gal_gal_hf_last;

  int has_n5k;
  int index_md_n5k;

  int has_gal_lens_hf;
  int index_md_gal_lens_hf;
  int index_integrand_id_gal_lens_hf_first;
  int index_integrand_id_gal_lens_hf_last;

  int has_lens_lens_hf;
  int index_md_lens_lens_hf;
  int index_integrand_id_lens_lens_hf_first;
  int index_integrand_id_lens_lens_hf_last;


  int has_gal_lens_2h;
  int index_md_gal_lens_2h;
  int index_integrand_id_gal_lens_2h_first;
  int index_integrand_id_gal_lens_2h_last;

  int has_gal_lens_1h;
  int index_md_gal_lens_1h;
  int index_integrand_id_gal_lens_1h_first;
  int index_integrand_id_gal_lens_1h_last;

  int has_gal_lensmag_hf;
  int index_md_gal_lensmag_hf;
  int index_integrand_id_gal_lensmag_hf_first;
  int index_integrand_id_gal_lensmag_hf_last;

  int has_gal_lensmag_2h;
  int index_md_gal_lensmag_2h;
  int index_integrand_id_gal_lensmag_2h_first;
  int index_integrand_id_gal_lensmag_2h_last;

  int has_gal_lensmag_1h;
  int index_md_gal_lensmag_1h;
  int index_integrand_id_gal_lensmag_1h_first;
  int index_integrand_id_gal_lensmag_1h_last;

  int convert_cls_to_gamma;
  int has_gal_gallens_2h;
  int index_md_gal_gallens_2h;
  int index_integrand_id_gal_gallens_2h_first;
  int index_integrand_id_gal_gallens_2h_last;

  int has_gal_gallens_1h;
  int index_md_gal_gallens_1h;
  int index_integrand_id_gal_gallens_1h_first;
  int index_integrand_id_gal_gallens_1h_last;


  int has_gallens_gallens_2h;
  int index_md_gallens_gallens_2h;
  int index_integrand_id_gallens_gallens_2h_first;
  int index_integrand_id_gallens_gallens_2h_last;

  int has_gallens_gallens_1h;
  int index_md_gallens_gallens_1h;
  int index_integrand_id_gallens_gallens_1h_first;
  int index_integrand_id_gallens_gallens_1h_last;


  int has_gallens_lens_2h;
  int index_md_gallens_lens_2h;
  int index_integrand_id_gallens_lens_2h_first;
  int index_integrand_id_gallens_lens_2h_last;

  int has_gallens_lens_1h;
  int index_md_gallens_lens_1h;
  int index_integrand_id_gallens_lens_1h_first;
  int index_integrand_id_gallens_lens_1h_last;




  int has_tSZ_lensmag_2h;
  int index_md_tSZ_lensmag_2h;
  int index_integrand_id_tSZ_lensmag_2h_first;
  int index_integrand_id_tSZ_lensmag_2h_last;

  int has_tSZ_lensmag_1h;
  int index_md_tSZ_lensmag_1h;
  int index_integrand_id_tSZ_lensmag_1h_first;
  int index_integrand_id_tSZ_lensmag_1h_last;


  int has_lensmag_lensmag_hf;
  int index_md_lensmag_lensmag_hf;
  int index_integrand_id_lensmag_lensmag_hf_first;
  int index_integrand_id_lensmag_lensmag_hf_last;

  int has_lensmag_lensmag_2h;
  int index_md_lensmag_lensmag_2h;
  int index_integrand_id_lensmag_lensmag_2h_first;
  int index_integrand_id_lensmag_lensmag_2h_last;

  int has_lensmag_lensmag_1h;
  int index_md_lensmag_lensmag_1h;
  int index_integrand_id_lensmag_lensmag_1h_first;
  int index_integrand_id_lensmag_lensmag_1h_last;


  int has_lens_lensmag_hf;
  int index_md_lens_lensmag_hf;
  int index_integrand_id_lens_lensmag_hf_first;
  int index_integrand_id_lens_lensmag_hf_last;


  int has_lens_lensmag_2h;
  int index_md_lens_lensmag_2h;
  int index_integrand_id_lens_lensmag_2h_first;
  int index_integrand_id_lens_lensmag_2h_last;

  int has_lens_lensmag_1h;
  int index_md_lens_lensmag_1h;
  int index_integrand_id_lens_lensmag_1h_first;
  int index_integrand_id_lens_lensmag_1h_last;

  int has_lens_lens_1h;
  int index_md_lens_lens_1h;
  int index_integrand_id_lens_lens_1h_first;
  int index_integrand_id_lens_lens_1h_last;

  int has_lens_lens_2h;
  int index_md_lens_lens_2h;
  int index_integrand_id_lens_lens_2h_first;
  int index_integrand_id_lens_lens_2h_last;

  int has_tSZ_gal_1h;
  int index_md_tSZ_gal_1h;
  int index_integrand_id_tSZ_gal_1h_first;
  int index_integrand_id_tSZ_gal_1h_last;

  int has_tSZ_gal_2h;
  int index_md_tSZ_gal_2h;
  int index_integrand_id_tSZ_gal_2h_first;
  int index_integrand_id_tSZ_gal_2h_last;

  int has_tSZ_gallens_1h;
  int index_md_tSZ_gallens_1h;
  int index_integrand_id_tSZ_gallens_1h_first;
  int index_integrand_id_tSZ_gallens_1h_last;

  int has_tSZ_gallens_2h;
  int index_md_tSZ_gallens_2h;
  int index_integrand_id_tSZ_gallens_2h_first;
  int index_integrand_id_tSZ_gallens_2h_last;

  int has_gallens_cib_1h;
  int index_md_gallens_cib_1h;
  int index_integrand_id_gallens_cib_1h_first;
  int index_integrand_id_gallens_cib_1h_last;

  int has_gallens_cib_2h;
  int index_md_gallens_cib_2h;
  int index_integrand_id_gallens_cib_2h_first;
  int index_integrand_id_gallens_cib_2h_last;

  int has_gal_cib_1h;
  int index_md_gal_cib_1h;
  int index_integrand_id_gal_cib_1h_first;
  int index_integrand_id_gal_cib_1h_last;

  int has_gal_cib_2h;
  int index_md_gal_cib_2h;
  int index_integrand_id_gal_cib_2h_first;
  int index_integrand_id_gal_cib_2h_last;

  int has_tSZ_cib_1h;
  int index_md_tSZ_cib_1h;
  int index_integrand_id_tSZ_cib_1h_first;
  int index_integrand_id_tSZ_cib_1h_last;

  int has_tSZ_cib_2h;
  int index_md_tSZ_cib_2h;
  int index_integrand_id_tSZ_cib_2h_first;
  int index_integrand_id_tSZ_cib_2h_last;

  int has_lens_cib_1h;
  int index_md_lens_cib_1h;
  int index_integrand_id_lens_cib_1h_first;
  int index_integrand_id_lens_cib_1h_last;

  int has_lens_cib_2h;
  int index_md_lens_cib_2h;
  int index_integrand_id_lens_cib_2h_first;
  int index_integrand_id_lens_cib_2h_last;

  int has_cib_monopole;
  int index_md_cib_monopole;
  int index_integrand_id_cib_monopole_first;
  int index_integrand_id_cib_monopole_last;

  int has_cib_shotnoise;
  int index_md_cib_shotnoise;
  int index_integrand_id_cib_shotnoise_first;
  int index_integrand_id_cib_shotnoise_last;

  int has_ngal_ngal_1h;
  int index_md_ngal_ngal_1h;
  int index_integrand_id_ngal_ngal_1h_first;
  int index_integrand_id_ngal_ngal_1h_last;

  int has_ngal_ngal_2h;
  int index_md_ngal_ngal_2h;
  int index_integrand_id_ngal_ngal_2h_first;
  int index_integrand_id_ngal_ngal_2h_last;

  int has_ngal_ngal_hf;
  int index_md_ngal_ngal_hf;
  int index_integrand_id_ngal_ngal_hf_first;
  int index_integrand_id_ngal_ngal_hf_last;

  int has_ngal_lens_1h;
  int index_md_ngal_lens_1h;
  int index_integrand_id_ngal_lens_1h_first;
  int index_integrand_id_ngal_lens_1h_last;

  int has_ngal_lens_2h;
  int index_md_ngal_lens_2h;
  int index_integrand_id_ngal_lens_2h_first;
  int index_integrand_id_ngal_lens_2h_last;

  int has_ngal_lens_hf;
  int index_md_ngal_lens_hf;
  int index_integrand_id_ngal_lens_hf_first;
  int index_integrand_id_ngal_lens_hf_last;


  int has_cib_cib_1h;
  int index_md_cib_cib_1h;
  int index_integrand_id_cib_cib_1h_first;
  int index_integrand_id_cib_cib_1h_last;

  int has_cib_cib_2h;
  int index_md_cib_cib_2h;
  int index_integrand_id_cib_cib_2h_first;
  int index_integrand_id_cib_cib_2h_last;

  int has_tSZ_lens_2h;
  int index_md_tSZ_lens_2h;
  int index_integrand_id_tSZ_lens_2h_first;
  int index_integrand_id_tSZ_lens_2h_last;

  int has_isw_lens;
  int index_md_isw_lens;
  int index_integrand_id_isw_lens_first;
  int index_integrand_id_isw_lens_last;

  int has_isw_tsz;
  int index_md_isw_tsz;
  int index_integrand_id_isw_tsz_first;
  int index_integrand_id_isw_tsz_last;

  int has_isw_auto;
  int index_md_isw_auto;
  int index_integrand_id_isw_auto_first;
  int index_integrand_id_isw_auto_last;


  int has_dndlnM;
  int index_md_dndlnM;
  int index_integrand_id_dndlnM_first;
  int index_integrand_id_dndlnM_last;

  int has_sz_trispec;
  //int index_md_sz_trispec;
  int index_integrand_id_trispectrum_first;
  int index_integrand_id_trispectrum_last;
  int index_md_trispectrum;


  int number_of_integrands;
  int index_integrand;
  int index_integrand_te_y_y;
  int index_integrand_2halo_term;

  int index_integrand_trispectrum_first; //for trispectrum
  int index_integrand_trispectrum_last;  //for trispectrum

  int index_integrand_cov_N_cl_first;
  int index_integrand_cov_N_cl_last;


  int index_integrand_N_for_cov_N_cl_first;
  int index_integrand_N_for_cov_N_cl_last;


  int index_integrand_id;

  int number_of_integrals_per_thread;

  int index_integrands_first;
  int index_integrands_last;

  int index_md_dcib0dz;
  int index_md_dydz;








  //double  pk;

  // FileName root; /**< root for all file names */
  // FileName path_to_class; /**< root for all file names */
  // FileName append_name_cobaya_ref;
  // FileName path_to_ref_trispectrum_for_cobaya;
  // FileName full_path_to_noise_curve_for_y_y;
  //FileName full_path_to_dndz_gal;

 /* vector of all SZ quantities function of redshift*/

  int  tsz_size;

  int  index_flag_cov_N_cl;
  int  index_Rho_crit;
  int  index_Delta_c;
  int  index_rVIR;
  int  index_cVIR;
  int  index_c200m;
  int  index_r200m;
  int  index_mVIR;
  int  index_m500;
  int  index_r500;
  int  index_l500;
  int  index_ls;
  int  index_rs;
  int  index_m200;
  int  index_m180m;
  int  index_m200m;
  int  index_m1600m;
  int  index_m500c;
  int  index_mass_for_hmf;
  int  index_mass_for_galaxies;
  int  index_mass_for_cib;
  int  index_mass_for_matter_density;
  int  index_mass_for_electron_pressure;
  int  index_mass_for_electron_density;
  int  index_mass_for_HI_pressure;
  int  index_mass_for_HI_density;
  int  index_concentration_for_galaxies;
  int  index_concentration_for_cib;
  int  index_concentration_for_matter_density;
  int  index_concentration_for_electron_pressure;
  int  index_concentration_for_electron_density;
  int  index_concentration_for_HI_pressure;
  int  index_concentration_for_HI_density;
  int  index_radius_for_galaxies;
  int  index_radius_for_cib;
  int  index_radius_for_matter_density;
  int  index_radius_for_electron_pressure;
  int  index_radius_for_electron_density;
  int  index_radius_for_HI_pressure;
  int  index_radius_for_HI_density;
  int  index_r500c;
  int  index_Rh;
  int  index_mf;
  int  index_dlognudlogRh;
  int  index_lognu;
  int  index_dlogSigma2dlogRh;
  int  index_dndlogRh;
  int  index_logSigma2;
  int  index_z;
  int  index_c200c;
  int  index_m200c;
  int  index_l200c;
  int  index_characteristic_multipole_for_nfw_profile;
  int  index_r200c;
  int  index_multipole;
  int  index_szrate;
  int  index_multipole_prime;
  int  index_mass_bin_1;
  int  index_mass_bin_2;
  int  index_multipole_1;
  int  index_multipole_2;
  int  index_multipole_3;
  int  index_redshift_for_dndlnM;
  int  index_mass_for_dndlnM;
  int  index_multipole_for_pressure_profile;
  int  index_pressure_profile;
  int  index_multipole_for_tau_profile;
  int  index_multipole_for_nfw_profile;
  int  index_tau_profile;
  int  index_lensing_profile;
  int  index_multipole_for_lensing_profile;
  int  index_completeness;
  int  index_te_of_m;
  int  index_volume;
  int  index_chi2; // in [Mpc/h]^2
  int  index_dgdz; // d(D/a)/dz = D(1-f)
  int  index_lensing_Sigma_crit;
  int  index_vrms2;
  int  index_pk_for_halo_bias;
  int  index_dlnMdeltadlnM;
  int  index_part_id_cov_hsv;

  int  index_mean_y;
  int  index_hmf;

  int index_sigma2_hsv;

  int  index_halo_bias;
  int  index_halo_bias_b2;
  int  index_k_value_for_halo_bias;

  int index_phi_galaxy_counts;
  int index_mean_galaxy_number_density;
  int index_mean_galaxy_bias;
  int index_c500c;
  int index_multipole_for_galaxy_profile;
  int index_multipole_for_truncated_nfw_profile;
  int index_galaxy_profile;

  int index_ngal_for_galaxy_profile;
  int index_ngal_prime_for_galaxy_profile;


  int index_multipole_for_cib_profile;
  int index_frequency_for_cib_profile;
  int index_frequency_prime_for_cib_profile;
  int index_cib_profile;

  int index_W_lensmag;

  int index_W_gallens_sources;

  int index_k_for_pk_hm;
  int index_density_profile;

  int index_multipole_for_pk;

  int index_ell_1;
  int index_ell_2;
  int index_ell_3;

  //////////////

  int index_integral;
  int index_integral_te_y_y;
  int index_integral_2halo_term;

  int index_integral_trispectrum_first;
  int index_integral_trispectrum_last;

  int index_integral_cov_N_cl_first;
  int index_integral_cov_N_cl_last;

  int index_integral_N_for_cov_N_cl_first;
  int index_integral_N_for_cov_N_cl_last;


  int  index_integrals_over_m_first;
  int  index_integrals_over_m_last;

  int  index_integrals_over_z_first;
  int  index_integrals_over_z_last;



  int  index_integral_over_m;
  int  index_integral_te_y_y_over_m;
  int  index_integral_2halo_term_over_m;
  int  index_integral_trispectrum_first_over_m;
  int  index_integral_trispectrum_last_over_m;
  int  index_integral_cov_N_cl_first_over_m;
  int  index_integral_cov_N_cl_last_over_m;
  int  index_integral_N_for_cov_N_cl_first_over_m;
  int  index_integral_N_for_cov_N_cl_last_over_m;




  //mass bins for covariance between cluster counts and power spectrum
  int nbins_M;
  double * M_bins;
  double dlogM;
  double * cov_Y_N_mass_bin_edges;

double * szcounts_fft_qobs;
double * szcounts_fft_z;
double * szcounts_fft_sigmayobs;
int ** szcounts_fft_index_zsig;
// double ** szcounts_fft_rates_at_z_sigy_qobs;
// double * szcounts_fft_nexpected_dndzdqgt;
double * szcounts_fft_dndzdq;
double * szcounts_fft_nexpected_qobs;
int ** szcounts_fft_index_zq;
int ** szcounts_fft_index_zq_final;
double ** szcounts_fft_qmconv_all_patches;

double szcounts_ntot;

  //HOD
  double M_min_HOD;
  double M_min_HOD_cib;
  double M0_HOD;
  double sigma_log10M_HOD;
  double alpha_s_HOD;
  double M1_prime_HOD;

  double * M_min_HOD_ngal;
  double * M0_HOD_ngal;
  double * sigma_log10M_HOD_ngal;
  double * alpha_s_HOD_ngal;
  double * M1_prime_HOD_ngal;


  double rho_y_gal;

  int M0_Mmin_flag;

  double M_min_HOD_mass_factor_unwise;
  double x_out_truncated_nfw_profile;
  double x_out_truncated_nfw_profile_electrons;
  double x_out_truncated_density_profile;
  double M_min_HOD_satellite_mass_factor_unwise;
  double M1_prime_HOD_factor;
  double cvir_tau_profile_factor;


  double effective_galaxy_bias;
  double * effective_galaxy_bias_ngal;

  int use_bg_eff_in_ksz2g_eff;

  int hm_consistency;

  int use_class_sz_fast_mode;
  double * array_lnk;
  double * array_pknl_at_z_and_k;
  double * array_pkl_at_z_and_k;

  // int cszfast_pk_grid_nk;
  // int cszfast_pk_grid_nz;

  int check_consistency_conditions;

  // noise curve for cov(y,y)

  int include_noise_cov_y_y;

  //units for tSZ spectrum
  double exponent_unit;

  //completeness
  double theta_bin_min;
  double theta_bin_max;
  int nthetas;
  double * thetas;

  double *skyfracs;
  int nskyfracs;

  int Ny;
  int Nth;
  double * erfs_2d_to_1d_th_array;
  double * erfs_2d_to_1d_y_array;

  double ** ylims;
  double * sky_averaged_ylims;

  //SZ catalog
  double * szcat_z;
  double * szcat_snr;
  int  szcat_size;

  double shape_noise_siggamma2;
  double ns_gal_per_arcmin2;
  double cl_gal_gal_A_sn;

  int experiment;
  //SO completeness
  double * SO_thetas;
  double * SO_Qfit;
  int  SO_Q_size;

  double * SO_RMS;
  double * SO_skyfrac;
  int  SO_RMS_size;

  double csat_over_cdm;
  //INPUT PARAMETERS
  int nlSZ;
  int n_ell_independent_integrals;
  int n_frequencies_for_cib;

  double * l_unwise_filter;
  double * f_unwise_filter;
  int unwise_filter_size;

  double * M_min_of_z_z;
  double * M_min_of_z_M_min;
  int M_min_of_z_size;

  double * nl_lensing_noise;
  double * l_lensing_noise;
  int lensing_noise_size;

  double * l_ksz_template;
  double * cl_ksz_template;
  int ksz_template_size;

  int damping_1h_term;
  double kstar_damping_1h_term_Mpc; // inverse Mpc

  /*Redshift limits for the integration*/
  double z1SZ;
  double z2SZ;

  double z1SZ_dndlnM;
  double z2SZ_dndlnM;

  double M1SZ_dndlnM;
  double M2SZ_dndlnM;


  double y_min;
  double y_max;


  double theta_ej_bcm;
  double delta_bcm;
  double gamma_bcm;
  double eta_star_bcm;
  double log10Mc_bcm;
  double mu_bcm;
  double nu_log10Mc_bcm;
  // double xxx_bcdm;

  // int n_y_y_to_m;
  // int n_z_y_to_m;
  // int n_mass_y_to_m;
  double * array_y_to_m_y;
  double * array_y_to_m_at_z_y;
  double * array_y_to_m_redshift;

  // double z1SZ_L_sat;
  // double z2SZ_L_sat;
  //
  // double M1SZ_L_sat;
  // double M2SZ_L_sat;
  //
  // double epsabs_L_sat;
  // double epsrel_L_sat;

  double mass_epsrel_cluster_counts;
  double mass_epsabs_cluster_counts;

  double redshift_epsrel_cluster_counts;
  double redshift_epsabs_cluster_counts;

  double dlnM_cluster_count_completeness_grid;
  double dz_cluster_count_completeness_grid_low_z;
  double dz_cluster_count_completeness_grid_mid_z;
  double dz_cluster_count_completeness_grid_high_z;

  double lnymin;
  double lnymax;
  double dlny;

  double cluster_count_completeness_grid_z_cutoff_low;
  double cluster_count_completeness_grid_z_cutoff_mid;

  int n_z_W_lensmag;
  int n_z_W_gallens_sources;

  /*Array size*/
  int n_arraySZ;//number of z in the interpolation
  int n_arraySZ_for_integral;//number of z in the integration

  int n_k;
  int n_z_dndlnM;
  int n_m_dndlnM;

  int compute_ksz2ksz2;

  // int n_z_L_sat;
  // int n_m_L_sat;
  // int n_nu_L_sat;

  int N_redshift_dndlnM;
  int N_mass_dndlnM;

  //mass limits: h^-1 Msun
  double M1SZ;
  double M2SZ;

  double delta_alpha;
  double alpha_p;

  double alpha_b;
  double Ap;
  int mass_dependent_bias;

  //Planck pressure profile
  double P0GNFW;
  double c500;
  double gammaGNFW;
  double alphaGNFW;
  double betaGNFW;

  double ln_x_size_for_pp;
  double * ln_x_for_pp;

  double x_size_for_pp;
  double * x_for_pp;

  int use_websky_m200m_to_m200c_conversion;
  //Battaglia pressure profile
  double alpha_B12;
  double gamma_B12;
  double P0_B12;
  double xc_B12;
  double beta_B12;

  double alpha_m_P0_B12;
  double alpha_m_xc_B12;
  double alpha_m_beta_B12;

  double alpha_z_P0_B12;
  double alpha_z_xc_B12;
  double alpha_z_beta_B12;


  // B.H.
  double mcut_B12;
  double c_B12;
  double c_B16;
  double alphap_m_P0_B12;
  double alphap_m_xc_B12;
  double alphap_m_beta_B12;

  double alpha_c_P0_B12;
  double alpha_c_xc_B12;
  double alpha_c_beta_B12;


    // B.H.
  double mcut;
  double alphap_m_rho0;
  double alphap_m_alpha;
  double alphap_m_beta;

  double alpha_c_rho0;
  double alpha_c_alpha;
  double alpha_c_beta;

  // Battaglia density profile:
  double A_rho0;
  double A_alpha;
  double A_beta;

  double alpha_m_rho0;
  double alpha_m_alpha;
  double alpha_m_beta;

  double alpha_z_rho0;
  double alpha_z_alpha;
  double alpha_z_beta;

  double gamma_B16;
  double xc_B16;


// JCH
//double precision :: f_free=0.85d0 !for kSZ calculations, fraction of free electrons w.r.t. total
//double precision :: mu_e=1.14d0 !mean molecular weight per electron, for primordial composition

  double f_free;
  double mu_e;
  double f_b_gas;

  /*Pressure profile is considered between x_in and x_out*/
  double x_inSZ;
  double x_outSZ;

  double HSEbias;

  /*For the computation of sigma2*/
  int  ndimSZ;
  double logR1SZ; // 0.0034Mpc/h, 1.8e4  solar mass
  double logR2SZ; // 54.9Mpc/h, 7.5e16 solar mass
  double delta_cSZ;



  /*Multplicity function Tinker 2010*/

  double alphaSZ;
  double beta0SZ;
  double gamma0SZ;

  double phi0SZ;
  double eta0SZ;
  int T10_alpha_fixed;


  /*Multplicity function Bocquet 2015*/

  double Ap0;
  double a0;
  double b0;
  double c0;

  int pk_nonlinear_for_vrms2;

  int MF;
  //1:Tinker 2010 (T10)
  //2:Bocquet 2015 (B15)
  //3:Jenkins 2001 (J01)
  //4:Tinker 2008 (T08)
  //5:Tinker 2008 interpolated @ M500 (T08@M500)
  int SHMF;

  //Precision Parameters For qromb_sz_integrand
  int K;
  double EPS;
  double JMAX;


  //Precision Parameters For qromb_sz_sigma
  int K_sigma;
  double EPS_sigma;
  double JMAX_sigma;

  ////////////////////////
  //integration method and parameters (mass)
  int integration_method_mass;

  double redshift_epsrel;
  double redshift_epsabs;

  double mass_epsrel;
  double mass_epsabs;

  double pressure_profile_epsabs;
  double pressure_profile_epsrel;
  double nu_y_dist_GHz;


  int * galaxy_samples_list;
  int galaxy_samples_list_num;
  int ngal_dim;

  int cib_frequency_list_num;
  int cib_dim;
  double * cib_frequency_list;
  double * cib_Snu_cutoff_list_in_mJy;

  int id_nu_cib_to_save;
  int id_nu_prime_cib_to_save;

  double ystar_ym;
  double alpha_ym;
  double sigmaM_ym;
  double beta_ym;
  double A_ym;
  double B_ym;
  double C_ym;
  double m_pivot_ym;

  double alpha_theta;
  int y_m_relation;
  double thetastar;

  int use_maniyar_cib_model;
  double maniyar_cib_tau;
  double maniyar_cib_zc;
  double maniyar_cib_etamax;
  double maniyar_cib_fsub;

  //BB: added for class_sz
  int ln_k_size_for_tSZ;
  double k_per_decade_for_tSZ;
  double k_min_for_pk_in_tSZ;
  double k_max_for_pk_in_tSZ;
  double * ln_k_for_tSZ;


int nsteps_m;
int nsteps_z;

double * steps_z;
double * steps_m;

  // Table 1  of MM20
  double alpha_cib; //redshift evolution of dust temperature
  double T0_cib; // dust temperature today
  double beta_cib; // emissivity index of sed
  double gamma_cib; // Power law index of SED at high frequency
  double delta_cib; // Redshift evolution of L − M normalisation
  double m_eff_cib; // Most efficient halo mass in Msun/h
  double L0_cib; // Normalisation of L − M relation
  double sigma2_LM_cib; // Size of of halo masses sourcing CIB emission
  int has_cib_flux_cut;
  double z_obs_cib;
  double z_plateau_cib;
  double M_min_subhalo_in_Msun;
  int use_nc_1_for_all_halos_cib_HOD;
  int use_redshift_dependent_M_min;

  double nfw_profile_epsabs;
  double nfw_profile_epsrel;

  int patterson_show_neval;

  int number_of_mass_bins; //for trapezoidal rule
  ////////////////////////

  ////////////////////////
  //integration method and parameters (pressure profile)
  int integration_method_pressure_profile;

  //Foreground parameters
  double A_cib, A_rs, A_ir, A_cn;

  double * k_for_pk_hm;
  double dlnk_for_pk_hm;
  double k_min_for_pk_hm;
  double k_max_for_pk_hm;
  int n_k_for_pk_hm;
  double z_for_pk_hm;

  //Cl spectrum
  double * ell;
  double * ell_plc;
  double * ell_plc_no_low_ell;
  double * ell_plc_low;
  double * ell_mock;
  double * ell_trispectrum;
  double * x_gauss;
  double * w_gauss;

  double * frequencies_for_cib;

  double * ell_kSZ2_gal_multipole_grid;
  int N_kSZ2_gal_multipole_grid;


  double * theta_kSZ2_gal_theta_grid;
  int N_kSZ2_gal_theta_grid;


  double dlogell;
  double dell;
  double ell_min_mock;
  double ell_max_mock;


  double freq_max;
  double freq_min;
  double dlogfreq;
  double dfreq;

  double Tcmb_gNU_at_150GHz;
  double Tcmb_gNU;

  double Rho_crit_0;
  double D_0;
  double D_z1SZ;
  double Omega_m_0;
  double Omega_r_0;
  double Omega_ncdm_0;
  double Omega0_b;
  double Omega0_cdm;
  double bispectrum_lambda_k2;
  double bispectrum_lambda_k3;

  double Sigma8OmegaM_SZ;
  double sigma8_Pcb;

  short has_knl;
  short has_nl_index;
  short has_vrms2;
  short has_sigma2_hsv;

  short has_tszspectrum;  //do we need tSZ spectrum? */
  short sz_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */
  short write_sz;  //do we need SZ quatitiies vs redshift? */

  int use_planck_binned_proba;
  double bin_z_min_cluster_counts;
  double bin_z_max_cluster_counts;
  double bin_dz_cluster_counts;
  int apply_relativistic_correction_to_y_m;
  double bin_dlog10_snr;
  double bin_dlog10_snr_last_bin;
  double log10_snr_min;
  double log10_snr_max;

  double x_out_truncated_nfw_profile_satellite_galaxies;
  double * x_out_truncated_nfw_profile_satellite_galaxies_ngal;

  double f_cen_HOD;
  double * f_cen_HOD_ngal;
  double Delta_z_lens;
  double Delta_z_source;

  short has_completeness_for_ps_SZ;
  short has_completeness;
  short which_ps_sz;
  double H0_in_class_units;
  double sky_area_deg2;
  int ell_sz;
  // Figure 7 of KS02 -> KS02
  // Planck 2015 effective multipoles -> P15
  // SZFASTDKS -> DKS

  // halo occupation distribution
  int hod_model;


  int concentration_parameter;
  //Duffy 2008: D08
  //Seljak 2000: S00

  int pressure_profile;
  //Planck 2013 (P13)
  //Arnaud et al 2010 (A10)
  //Custom. GNFW

  int tau_profile;
  int tau_profile_mode;

  int HMF_prescription_NCDM;
  int effective_temperature;
  int temperature_mass_relation;
  int mean_y;

  double * PP_lnx;
  double * PP_lnI;
  double * PP_d2lnI;

  int PP_lnx_size;

  double * RNFW_lnx;
  double * RNFW_lnI;

  int RNFW_lnx_size;


  double * T10_ln1pz;
  double * T10_lnalpha;
  int T10_lnalpha_size;

  double * normalized_source_dndz_z;
  double * normalized_source_dndz_phig;

  int normalized_source_dndz_size;

  double * normalized_dndz_z;
  double * normalized_dndz_phig;

  double ** normalized_dndz_ngal_z;
  double ** normalized_dndz_ngal_phig;

  int normalized_dndz_size;
  int * normalized_dndz_ngal_size;

  double * normalized_fdndz_z;
  double * normalized_fdndz_phig;

  int normalized_fdndz_size;

  double * normalized_cosmos_dndz_z;
  double * normalized_cosmos_dndz_phig;

  int normalized_cosmos_dndz_size;

  double * unbinned_nl_yy_ell;
  double * unbinned_nl_yy_n_ell;
  int nl_yy_is_binned;
  int unbinned_nl_yy_size;

  double * unbinned_nl_tt_ell;
  double * unbinned_nl_tt_n_ell;
  int unbinned_nl_tt_size;

  int truncate_wrt_rvir;

  int no_tt_noise_in_kSZ2X_cov;

  double * CM_redshift;
  double * CM_logM;

  int CM_redshift_size;
  int CM_logM_size;
  double * CM_logC;


  // double * array_profile_2h_ln_1pz;
  double * array_profile_ln_rho_2h_at_k_and_z;
  double * array_profile_rho_2h_at_r_and_z;


  double * array_m_m200m_to_m200c;
  double * array_ln_1pz_m200m_to_m200c;
  double * array_m200m_to_m200c_at_z_and_M;


  double * array_m_m200c_to_m200m;
  double * array_ln_1pz_m200c_to_m200m;
  double * array_m200c_to_m200m_at_z_and_M;


  double * array_m_m200m_to_m500c;
  double * array_ln_1pz_m200m_to_m500c;
  double * array_m200m_to_m500c_at_z_and_M;

  double * array_m_m200c_to_m500c;
  double * array_ln_1pz_m200c_to_m500c;
  double * array_m200c_to_m500c_at_z_and_M;


  double * array_m_m500c_to_m200c;
  double * array_ln_1pz_m500c_to_m200c;
  double * array_m500c_to_m200c_at_z_and_M;

  double ** array_pressure_profile_ln_p_at_lnk_lnm_z;
  double * array_pressure_profile_ln_k;
  double * array_pressure_profile_2h_ln_k;
  double * array_pressure_profile_ln_r;
  double * array_pressure_profile_ln_m;
  double * array_pressure_profile_ln_1pz;
  double * array_pressure_profile_ln_pressure_2h_at_k_and_z;
  double * array_pressure_profile_pressure_2h_at_r_and_z;

  double ** array_profile_ln_rho_at_lnk_lnM_z;
  double * array_profile_ln_r;
  double * array_profile_ln_k;
  double * array_profile_ln_m;
  double * array_profile_ln_1pz;

  // int array_profile_ln_PgNFW_at_lnl_over_ls_size; defined in class_sz_precisions.h
  double * array_profile_ln_l_over_ls;
  double * array_profile_ln_PgNFW_at_lnl_over_ls;

  double * dndlnM_array_z;
  double * dndlnM_array_m;

  double * array_m_dndlnM;
  double * array_z_dndlnM;
  double * array_dndlnM_at_z_and_M;

  double * array_m_L_sat;
  double * array_z_L_sat;
  double * array_nu_L_sat;

  double ** array_L_sat_at_M_z_nu;
  double ** array_L_sat_at_z_and_M_at_nu;
  //double * array_L_sat_at_z_and_M_at_nu_prime;



  double * array_z_W_lensmag;
  double * array_W_lensmag;

  double * array_z_W_gallens_sources;
  double * array_W_gallens_sources;

  double * array_redshift;
  double * array_radius;
  // double * array_k;
  double * array_nl_index_at_z_and_k;
  double * array_nl_index_at_z_and_k_no_wiggles;
  double * array_sigma_at_z_and_R;
  double * array_dsigma2dR_at_z_and_R;

  double * array_knl_at_z;
  double * array_vrms2_at_z;
  double * array_sigma2_hsv_at_z;

  double * array_mean_galaxy_number_density;
  double ** array_mean_galaxy_number_density_ngal;

  double * array_mean_galaxy_bias;

  // int n_z_hmf_counter_terms;
  int hm_consistency_counter_terms_done;
  double * array_redshift_hmf_counter_terms;
  double * array_hmf_counter_terms_nmin;
  double * array_hmf_counter_terms_b1min;
  double * array_hmf_counter_terms_b2min;
  ErrorMsg error_message; /**< zone for writing error messages */


  double * array_n5k_F1_F;
  double * array_n5k_F1_k;
  int * array_n5k_F1_l;

  double * n5k_pk_z;
  double * n5k_pk_k;
  double * n5k_pk_pk;
  int n5k_pk_z_size;
  int n5k_pk_k_size;


  double * cib_Snu_z;
  double * cib_Snu_nu;
  double * cib_Snu_snu;
  int cib_Snu_z_size;
  int cib_Snu_nu_size;

  double * n5k_cl_K1_K1;
  double * n5k_cl_K1_chi;
  int n5k_cl_K1_size;

  double * n5k_z_of_chi_z;
  double * n5k_z_of_chi_chi;
  int n5k_z_of_chi_size;


  double * array_psi_b2t_redshift;
  double * array_psi_b2t_multipole;
  double * array_psi_b2t_psi;

  double * array_psi_b2g_redshift;
  double * array_psi_b2g_multipole;
  double * array_psi_b2g_psi;

  double * array_psi_b1g_redshift;
  double * array_psi_b1g_multipole;
  double * array_psi_b1g_psi;

  double * array_psi_b2kg_redshift;
  double * array_psi_b2kg_multipole;
  double * array_psi_b2kg_psi;

  double * array_psi_b1kg_redshift;
  double * array_psi_b1kg_multipole;
  double * array_psi_b1kg_psi;



  double * array_psi_b1t_redshift;
  double * array_psi_b1t_multipole;
  double * array_psi_b1t_psi;

  double * array_dcib0dz_nu;
  double * array_dcib0dz_redshift;
  double * array_dcib0dz_at_z_nu;


  double * array_m_to_xout_mass;
  double * array_m_to_xout_redshift;
  double * array_m_to_xout_at_z_m;


  double * array_dydz_redshift;
  double * array_dydz_at_z;


  double * array_psi_b1gt_redshift;
  double * array_psi_b1gt_multipole;
  double ** array_psi_b1gt_psi;

  double * array_psi_b1kgt_redshift;
  double * array_psi_b1kgt_multipole;
  double ** array_psi_b1kgt_psi;

  // int n_z_psi_b1g;
  // int n_l_psi_b1g;


};

double get_pk_lin_at_k_and_z_fast(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct tszspectrum * ptsz);

struct nonlinear_pt {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these
      parameters and the content of the 'precision', 'background',
      'thermo', 'primordial' and 'spectra' structures) */

  //@{

int use_class_sz_fast_mode;

enum non_linear_method_pt method; /**< method for computing non-linear corrections (none, Halogit, etc.) */
enum non_linear_irresumm_pt irres;
enum non_linear_bias_pt bias;
enum non_linear_rsd_pt rsd;
// enum non_linear_cb_pt cb;
// enum non_linear_fftlogn norm;
enum non_linear_AP_effect_pt AP_effect;
enum non_linear_fNL_equil_ortho fNL_equil_ortho_switch;

  //@}

  /** @name - table non-linear corrections for matter density, sqrt(P_NL(k,z)/P_NL(k,z)) */

  //@{

    // M.I. Global collable variables have to be defined here (they are called with prefix pnlpt->)

    int z_pk_num;
    double z_pk[_Z_PK_NUM_MAX_];

  int k_size;      /**< k_size = total number of k values */
  int ln_k_size;   /**< k_size = total number of k values */
  double * k;      /**< k[index_k] = list of k values */
  double * ln_k;      /**< k[index_k] = list of k values */
  int tau_size;    /**< tau_size = number of values */
  double * tau;    /**< tau[index_tau] = list of time values */
  double * ln_tau;    /**< ln_tau is absolutely different from psp->ln_tau!! It is a simple log of tau with size 609*/

  int index_md_scalars,ic_size,tp_size;
  double * dd_sources_tp_delta_m;
  double * dd_sources_tp_delta_cb;
  double * sources_tp_delta_m;
  double * sources_tp_delta_cb;

short fast_output;
short cb;

//  double * pk_nl_out;

  //@{

  int k_size_cmb;  /**< k_size_cmb[index_md] number of k values used
                        for CMB calculations, requiring a fine
                        sampling in k-space */

  int k_size_cl;  /**< k_size_cl[index_md] number of k values used
                       for non-CMB \f$ C_l \f$ calculations, requiring a coarse
                       sampling in k-space. */

  double k_min;     /**< minimum value (over all modes) */
  double k_max;     /**< maximum value (over all modes) */

  //@}

double OmfidAP;

double * M13_oneline;
double * M22_oneline;
double * M22basic_oneline;
double * IFG2_oneline;
//double * M12_oneline;    //GC!


double * M12_oneline;

    //GC: ORTHOGONAL -- start

double * M12_oneline_ortho;

    //GC: ORTHOGONAL -- finish

/**/
double * M12_oneline_matter_multipoles_vv0_f2;
double * M12_oneline_matter_multipoles_vv0_f3;
double * M12_oneline_matter_multipoles_vd0_f1;
double * M12_oneline_matter_multipoles_vd0_f2;
double * M12_oneline_matter_multipoles_dd0_f0;
double * M12_oneline_matter_multipoles_dd0_f1;
double * M12_oneline_matter_multipoles_vv2_f3;
double * M12_oneline_matter_multipoles_vd2_f2;
double * M12_oneline_matter_multipoles_vv4_f3;
double * M12_oneline_matter_multipoles_vd4_f2;
/**/

    //GC: ORTHOGONAL -- start

    double * M12_oneline_matter_multipoles_vv0_f2_ortho;
    double * M12_oneline_matter_multipoles_vv0_f3_ortho;
    double * M12_oneline_matter_multipoles_vd0_f1_ortho;
    double * M12_oneline_matter_multipoles_vd0_f2_ortho;
    double * M12_oneline_matter_multipoles_dd0_f0_ortho;
    double * M12_oneline_matter_multipoles_dd0_f1_ortho;
    double * M12_oneline_matter_multipoles_vv2_f3_ortho;
    double * M12_oneline_matter_multipoles_vd2_f2_ortho;
    double * M12_oneline_matter_multipoles_vv4_f3_ortho;
    double * M12_oneline_matter_multipoles_vd4_f2_ortho;

    //GC: ORTHOGONAL -- finish


double * M12_oneline_matter_mu_powers_vd2_f1;
double * M12_oneline_matter_mu_powers_vd2_f2;
double * M12_oneline_matter_mu_powers_dd2_f1;
double * M12_oneline_matter_mu_powers_vv4_f2;
double * M12_oneline_matter_mu_powers_vd4_f2;
double * M12_oneline_matter_mu_powers_vv6_f3;


    //GC: ORTHOGONAL -- start

double * M12_oneline_matter_mu_powers_vd2_f1_ortho;
double * M12_oneline_matter_mu_powers_vd2_f2_ortho;
double * M12_oneline_matter_mu_powers_dd2_f1_ortho;
double * M12_oneline_matter_mu_powers_vv4_f2_ortho;
double * M12_oneline_matter_mu_powers_vd4_f2_ortho;
double * M12_oneline_matter_mu_powers_vv6_f3_ortho;


    //GC: ORTHOGONAL -- finish


/**/
double * M12_oneline_bias_real_space_b2;
double * M12_oneline_bias_real_space_bG2;
/**/


    //GC: ORTHOGONAL -- start


    double * M12_oneline_bias_real_space_b2_ortho;
    double * M12_oneline_bias_real_space_bG2_ortho;


    //GC: ORTHOGONAL -- finish


double * M12_oneline_bias_multipoles_b2_vv0_f1;
double * M12_oneline_bias_multipoles_bG2_vv0_f1;

    //GC: ORTHOGONAL -- finish

    double * M12_oneline_bias_multipoles_b2_vv0_f1_ortho;
    double * M12_oneline_bias_multipoles_bG2_vv0_f1_ortho;

    //GC: ORTHOGONAL -- finish


    //GC --> go now to the complexified ones...

    //GC --> THESE ALL REMAIN...

double complex * M13_oneline_complex;
double complex * M22_oneline_complex;
double complex * M22basic_oneline_complex;
double complex * IFG2_oneline_complex;

double complex * M22_oneline_0_vv_complex;
double complex * M13_0_vv_oneline_complex;

double complex * M22_oneline_0_vd_complex;
double complex * M13_0_vd_oneline_complex;

double complex * M22_oneline_0_dd_complex;
double complex * M13_0_dd_oneline_complex;

    double complex * M22_oneline_2_vv_complex;
    double complex * M13_2_vv_oneline_complex;

    double complex * M22_oneline_4_vv_complex;
    double complex * M13_4_vv_oneline_complex;

    double complex * M22_oneline_2_vd_complex;
    double complex * M13_2_vd_oneline_complex;

    double complex * M22_oneline_4_vd_complex;
    double complex * M13_4_vd_oneline_complex;
    double complex * M22_oneline_4_dd_complex;

    double complex * M22_oneline_2_dd_complex;
    double complex * M13_2_dd_oneline_complex;

double complex * M22_0_b1b2_oneline_complex;
double complex * M22_0_b2_oneline_complex;
double complex * M22_0_b1bG2_oneline_complex;
    double complex * M22_0_bG2_oneline_complex;


    double complex * M22_2_b1b2_oneline_complex;
    double complex * M22_2_b2_oneline_complex;
    double complex * M22_2_b1bG2_oneline_complex;
    double complex * M22_2_bG2_oneline_complex;

    double complex * M22_4_b2_oneline_complex;
    double complex * M22_4_bG2_oneline_complex;

double complex * M_Id2;
double complex * M_IG2;
double complex * M_Id2G2;
double complex * M_IG2G2;

double complex * M22_oneline_mu2_vd_complex;
double complex * M22_oneline_mu2_dd_complex;
    double complex * M22_oneline_mu4_vv_complex;
    double complex * M22_oneline_mu4_vd_complex;

    double complex * M22_oneline_mu4_dd_complex;
    double complex * M22_oneline_mu6_vv_complex;
    double complex * M22_oneline_mu6_vd_complex;
    double complex * M22_oneline_mu8_complex;

double complex * M13_mu2_dd_oneline_complex;
    double complex * M13_mu2_vd_oneline_complex;
    double complex * M13_mu4_vv_oneline_complex;
    double complex * M13_mu4_vd_oneline_complex;
    double complex * M13_mu6_oneline_complex;


    //GC!


    double complex * M12_oneline_complex;

    //GC: ORTHOGONAL -- start

    double complex * M12_oneline_complex_ortho;

    //GC: ORTHOGONAL -- finish



    double complex * M12_oneline_complex_matter_multipoles_vv0_f2;
    double complex * M12_oneline_complex_matter_multipoles_vv0_f3;
    double complex * M12_oneline_complex_matter_multipoles_vd0_f1;
    double complex * M12_oneline_complex_matter_multipoles_vd0_f2;
    double complex * M12_oneline_complex_matter_multipoles_dd0_f0;
    double complex * M12_oneline_complex_matter_multipoles_dd0_f1;
    double complex * M12_oneline_complex_matter_multipoles_vv2_f3;
    double complex * M12_oneline_complex_matter_multipoles_vd2_f2;
    double complex * M12_oneline_complex_matter_multipoles_vv4_f3;
    double complex * M12_oneline_complex_matter_multipoles_vd4_f2;


    //GC: ORTHOGONAL -- start


    double complex * M12_oneline_complex_matter_multipoles_vv0_f2_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vv0_f3_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vd0_f1_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vd0_f2_ortho;
    double complex * M12_oneline_complex_matter_multipoles_dd0_f0_ortho;
    double complex * M12_oneline_complex_matter_multipoles_dd0_f1_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vv2_f3_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vd2_f2_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vv4_f3_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vd4_f2_ortho;


    //GC: ORTHOGONAL -- finish



    double complex * M12_oneline_complex_matter_mu_powers_vd2_f1;
    double complex * M12_oneline_complex_matter_mu_powers_vd2_f2;
    double complex * M12_oneline_complex_matter_mu_powers_dd2_f1;
    double complex * M12_oneline_complex_matter_mu_powers_vv4_f2;
    double complex * M12_oneline_complex_matter_mu_powers_vd4_f2;
    double complex * M12_oneline_complex_matter_mu_powers_vv6_f3;


    //GC: ORTHOGONAL -- start


    double complex * M12_oneline_complex_matter_mu_powers_vd2_f1_ortho;
    double complex * M12_oneline_complex_matter_mu_powers_vd2_f2_ortho;
    double complex * M12_oneline_complex_matter_mu_powers_dd2_f1_ortho;
    double complex * M12_oneline_complex_matter_mu_powers_vv4_f2_ortho;
    double complex * M12_oneline_complex_matter_mu_powers_vd4_f2_ortho;
    double complex * M12_oneline_complex_matter_mu_powers_vv6_f3_ortho;


    //GC: ORTHOGONAL -- finish



    double complex * M12_oneline_complex_bias_real_space_b2;
    double complex * M12_oneline_complex_bias_real_space_bG2;


    //GC: ORTHOGONAL -- start

    double complex * M12_oneline_complex_bias_real_space_b2_ortho;
    double complex * M12_oneline_complex_bias_real_space_bG2_ortho;

    //GC: ORTHOGONAL -- finish



    double complex * M12_oneline_complex_bias_multipoles_b2_vv0_f1;
    double complex * M12_oneline_complex_bias_multipoles_bG2_vv0_f1;


    //GC: ORTHOGONAL -- start

    double complex * M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho;
    double complex * M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho;


    //GC: ORTHOGONAL -- finish



    double complex * M12_oneline_0_vv_complex;
    double complex * M12_oneline_0_vd_complex;
    double complex * M12_oneline_0_dd_complex;
    double complex * M12_oneline_2_vv_complex;
    double complex * M12_oneline_2_vd_complex;
    double complex * M12_oneline_2_dd_complex;
    double complex * M12_oneline_4_vv_complex;
    double complex * M12_oneline_4_vd_complex;
    double complex * M12_oneline_mu2_vd_complex;
    double complex * M12_oneline_mu2_dd_complex;
    double complex * M12_oneline_mu4_vv_complex;
    double complex * M12_oneline_mu4_vd_complex;
    double complex * M12_oneline_mu6_vv_complex;


    //GC: ORTHOGONAL -- start

    double complex * M12_oneline_0_vv_complex_ortho;
    double complex * M12_oneline_0_vd_complex_ortho;
    double complex * M12_oneline_0_dd_complex_ortho;
    double complex * M12_oneline_2_vv_complex_ortho;
    double complex * M12_oneline_2_vd_complex_ortho;
    double complex * M12_oneline_2_dd_complex_ortho;
    double complex * M12_oneline_4_vv_complex_ortho;
    double complex * M12_oneline_4_vd_complex_ortho;
    double complex * M12_oneline_mu2_vd_complex_ortho;
    double complex * M12_oneline_mu2_dd_complex_ortho;
    double complex * M12_oneline_mu4_vv_complex_ortho;
    double complex * M12_oneline_mu4_vd_complex_ortho;
    double complex * M12_oneline_mu6_vv_complex_ortho;

    //GC: ORTHOGONAL -- finish



    double complex * M12_2_bG2_oneline_complex;
    double complex * M12_2_b2_oneline_complex;
    double complex * M12_0_bG2_oneline_complex;
    double complex * M12_0_b1bG2_oneline_complex;
    double complex * M12_0_b2_oneline_complex;
    double complex * M12_0_b1b2_oneline_complex;


    //GC: ORTHOGONAL -- start

    double complex * M12_2_bG2_oneline_complex_ortho;
    double complex * M12_2_b2_oneline_complex_ortho;
    double complex * M12_0_bG2_oneline_complex_ortho;
    double complex * M12_0_b1bG2_oneline_complex_ortho;
    double complex * M12_0_b2_oneline_complex_ortho;
    double complex * M12_0_b1b2_oneline_complex_ortho;


    //GC: ORTHOGONAL -- finish


    double complex * M_fNLd2;
    double complex * M_fNLG2;


    //GC: ORTHOGONAL -- start

    double complex * M_fNLd2_ortho;
    double complex * M_fNLG2_ortho;

    //GC: ORTHOGONAL -- finish


    //GC!


double * ln_pk_nl;     /*For classy[i_z*pnlpt->k_size+index_k]*/
double * ln_pk_Id2d2;
    double * ln_pk_Id2d2_2;
    double * ln_pk_Id2d2_4;
double * ln_pk_Id2;
double * ln_pk_IG2;
double * ln_pk_Id2G2;
    double * ln_pk_Id2G2_2;
    double * ln_pk_Id2G2_4;
double * ln_pk_IG2G2;
    double * ln_pk_IG2G2_2;
    double * ln_pk_IG2G2_4;
double * ln_pk_IFG2;
    double * ln_pk_IFG2_0b1;
    double * ln_pk_IFG2_0;
    double * ln_pk_IFG2_2;
double * ln_pk_CTR;
    double * ln_pk_CTR_0;
    double * ln_pk_CTR_2;
    double * ln_pk_CTR_4;
double * ln_pk_Tree;
    double * ln_pk_Tree_0_vv;
    double * ln_pk_Tree_0_vd;
    double * ln_pk_Tree_0_dd;
    double * ln_pk_Tree_2_vv;
    double * ln_pk_Tree_2_vd;
    double * ln_pk_Tree_4_vv;

double * ln_pk_0_vv;
double * ln_pk_0_vd;
double * ln_pk_0_dd;

    double * ln_pk_2_vv;
    double * ln_pk_2_vd;
    double * ln_pk_2_dd;

    double * ln_pk_4_vv;
    double * ln_pk_4_vd;

    double * ln_pk_4_dd;

double * ln_pk_0_b1b2;
double * ln_pk_0_b1bG2;
double * ln_pk_0_b2;
    double * ln_pk_0_bG2;

    double * ln_pk_2_b1b2;
    double * ln_pk_2_b1bG2;
    double * ln_pk_2_b2;
    double * ln_pk_2_bG2;

    double * ln_pk_4_b2;
    double * ln_pk_4_bG2;
    double * ln_pk_4_b1b2;
    double * ln_pk_4_b1bG2;

    //GC!

    double * ln_pk_fNL_0_vv;
    double * ln_pk_fNL_0_vd;
    double * ln_pk_fNL_0_dd;

    double * ln_pk_fNL_2_vv;
    double * ln_pk_fNL_2_vd;
    double * ln_pk_fNL_2_dd;

    double * ln_pk_fNL_4_vv;
    double * ln_pk_fNL_4_vd;
    double * ln_pk_fNL_4_dd;


    //GC: ORTHOGONAL -- start


    double * ln_pk_fNL_0_vv_ortho;
    double * ln_pk_fNL_0_vd_ortho;
    double * ln_pk_fNL_0_dd_ortho;

    double * ln_pk_fNL_2_vv_ortho;
    double * ln_pk_fNL_2_vd_ortho;
    double * ln_pk_fNL_2_dd_ortho;

    double * ln_pk_fNL_4_vv_ortho;
    double * ln_pk_fNL_4_vd_ortho;
    double * ln_pk_fNL_4_dd_ortho;


    //GC: ORTHOGONAL -- finish



    double * ln_pk12_0_b1b2;
    double * ln_pk12_0_b2;
    double * ln_pk12_0_b1bG2;
    double * ln_pk12_0_bG2;
    double * ln_pk12_2_b1b2;
    double * ln_pk12_2_b2;
    double * ln_pk12_2_b1bG2;
    double * ln_pk12_2_bG2;
    double * ln_pk12_4_b1b2;
    double * ln_pk12_4_b2;
    double * ln_pk12_4_b1bG2;
    double * ln_pk12_4_bG2;



    //GC: ORTHOGONAL -- start


    double * ln_pk12_0_b1b2_ortho;
    double * ln_pk12_0_b2_ortho;
    double * ln_pk12_0_b1bG2_ortho;
    double * ln_pk12_0_bG2_ortho;
    double * ln_pk12_2_b1b2_ortho;
    double * ln_pk12_2_b2_ortho;
    double * ln_pk12_2_b1bG2_ortho;
    double * ln_pk12_2_bG2_ortho;
    double * ln_pk12_4_b1b2_ortho;
    double * ln_pk12_4_b2_ortho;
    double * ln_pk12_4_b1bG2_ortho;
    double * ln_pk12_4_bG2_ortho;


    //GC: ORTHOGONAL -- finish



    double * ln_pk_nl_fNL;
    double * ln_pk_fNLd2;
    double * ln_pk_fNLG2;


    //GC: ORTHOGONAL -- start


    double * ln_pk_nl_fNL_ortho;
    double * ln_pk_fNLd2_ortho;
    double * ln_pk_fNLG2_ortho;


    //GC: ORTHOGONAL -- finish


    //GC!

    double * growthf;   /*For RSD effect[i_z]*/
    double * hratio_array;   /*For RSD effect[i_z]*/
    double * Dratio_array;   /*For RSD effect[i_z]*/

    double * gauss;
    double * gauss_x;
    double * gauss_w;

    // double *koff;
    // double *Poff;
    // double *offtab;

/*< nl_corr_density[index_tau * ppt->k_size + index_k] */
    double * nl_corr_density;
    /*
    double * nl_corr_density;
    double * nl_corr_Id2d2;
    double * nl_corr_IG2;
    double * nl_corr_Id2;
    double * nl_corr_Id2G2;
    double * nl_corr_IG2G2;
    double * nl_corr_IFG2;
    double * nl_corr_CTR;
    double * nl_corr_Tree;
    */


    char input_pk[500];
    int replace_pk;
    int replace_background;
    int no_wiggle;
    int wiggle_only;
    double alpha_rs;
    double replace_Hz_value;
    double replace_DAz_value;
    double replace_Dz_value;
    double replace_fz_value;

//  double * k_nl;  /**< wavenumber at which non-linear corrections become important, defined differently by different non_linear_method's */
  int index_tau_min_nl; /**< index of smallest value of tau at which nonlinear corrections have been computed (so, for tau<tau_min_nl, the array nl_corr_density only contains some factors 1 */

  //@}

  /** @name - technical parameters */

  //@{

  short nonlinear_pt_verbose;  	/**< amount of information written in standard output */

  ErrorMsg error_message; 	/**< zone for writing error messages */

  //@}
};

/********************************************************************************/


extern void zspmv_(char*, int*, double complex*, double complex*, double complex*, int*, double complex*, double complex*, int*);
extern double complex zdotu_(int*, double complex*, int*, double complex*, int*);

/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif


  int nonlinear_pt_init(
                     struct precision *ppr,
                     struct background *pba,
                     struct thermo *pth,
                     struct perturbs *ppt,
                     struct primordial *ppm,
                     struct nonlinear_pt *pnlpt,
                     struct nonlinear * pnl,
                     struct tszspectrum * ptsz
                     );

  int perturb_get_k_list_nl(
                         struct precision * ppr,
                         struct background * pba,
                         struct thermo * pth,
                         struct perturbs * ppt,
                         struct nonlinear_pt *pnlpt
                         );

  int nonlinear_pt_free(
                     struct nonlinear_pt *pnlpt
                     );

  int nonlinear_pt_pk_l(
                     struct background *pba,
                     struct perturbs *ppt,
                     struct primordial *ppm,
                     struct nonlinear_pt *pnlpt,
                     struct nonlinear * pnl,
                     struct tszspectrum * ptsz,
                     int index_tau,
                     double *pk_l,
                     double *lnk,
                     double *lnpk,
                     double *ddlnpk); //,
                     //double *tk_l,
                     //double *lntk,
                     //double *ddlntk); //GC

//GC!

//GC: the function that extracts the primordial PS to get the transfer function {\cal M}(k)...

    int nonlinear_pt_pPRIMk_l(
                              //struct background *pba,
                              struct perturbs *ppt,
                              struct primordial *ppm,
                              struct nonlinear_pt *pnlpt,
                              double *lnk,
                              double *pPRIMk_l,
                              double *lnpPRIMk,
                              double *ddlnpPRIMk);

//GC!


    int nonlinear_pt_loop(
                          struct precision *ppr,
                          struct background *pba,
                          struct primordial *ppm,
                          struct thermo *pth,
                          struct nonlinear_pt *pnlpt,
                          double tau,
                          double f,
                          double hratio,
                          double Dratio,
                          double *pk_l,
                          double *pPRIMk_l, //GC! Never used...
                          double *pk_l_0_vv,
                          double *pk_l_0_vd,
                          double *pk_l_0_dd,
                          double *pk_l_2_vv,
                          double *pk_l_2_vd,
                          double *pk_l_2_dd,
                          double *pk_l_4_vv,
                          double *pk_l_4_vd,
                          double *pk_l_4_dd,
                          double *pk_nl,
                          double *pk_nl_fNL, //GC!
                          double *pk_fNLd2, //GC!
                          double *pk_fNLG2, //GC!
                          //GC: ORTHOGONAL -- start
                          double *pk_nl_fNL_ortho, //GC!
                          double *pk_fNLd2_ortho, //GC!
                          double *pk_fNLG2_ortho, //GC!
                          //GC: ORTHOGONAL -- finish
                          double *pk_Id2d2,
                          double *pk_Id2d2_2,
                          double *pk_Id2d2_4,
                          double *pk_l_0_b1b2,
                          double *pk_l_0_b2,
                          double *pk_l_0_b1bG2,
                          double *pk_l_0_bG2,
                          double *pk_l_2_b1b2,
                          double *pk_l_2_b2,
                          double *pk_l_2_b1bG2,
                          double *pk_l_2_bG2,
                          double *pk_l_4_b2,
                          double *pk_l_4_bG2,
                          double *pk_l_4_b1b2,
                          double *pk_l_4_b1bG2,
                          double *pk_Id2,
                          double *pk_IG2,
                          double *pk_Id2G2,
                          double *pk_Id2G2_2,
                          double *pk_Id2G2_4,
                          double *pk_IG2G2,
                          double *pk_IG2G2_2,
                          double *pk_IG2G2_4,
                          double *pk_IFG2,
                          double *pk_IFG2_0b1,
                          double *pk_IFG2_0,
                          double *pk_IFG2_2,
                          double *pk_CTR,
                          double *pk_CTR_0,
                          double *pk_CTR_2,
                          double *pk_CTR_4,
                          double *pk_Tree,
                          double *pk_Tree_0_vv,
                          double *pk_Tree_0_vd,
                          double *pk_Tree_0_dd,
                          double *pk_Tree_2_vv,
                          double *pk_Tree_2_vd,
                          double *pk_Tree_4_vv,
                          //GC!
                          double *pk_l_fNL_0_vv,
                          double *pk_l_fNL_0_vd,
                          double *pk_l_fNL_0_dd,
                          double *pk_l_fNL_2_vv,
                          double *pk_l_fNL_2_vd,
                          double *pk_l_fNL_2_dd,
                          double *pk_l_fNL_4_vv,
                          double *pk_l_fNL_4_vd,
                          double *pk_l_fNL_4_dd,
                          double *pk12_l_0_b1b2,
                          double *pk12_l_0_b2,
                          double *pk12_l_0_b1bG2,
                          double *pk12_l_0_bG2,
                          double *pk12_l_2_b1b2,
                          double *pk12_l_2_b2,
                          double *pk12_l_2_b1bG2,
                          double *pk12_l_2_bG2,
                          double *pk12_l_4_b1b2,
                          double *pk12_l_4_b2,
                          double *pk12_l_4_b1bG2,
                          double *pk12_l_4_bG2,
                          //GC: ORTHOGONAL -- start
                          double *pk_l_fNL_0_vv_ortho,
                          double *pk_l_fNL_0_vd_ortho,
                          double *pk_l_fNL_0_dd_ortho,
                          double *pk_l_fNL_2_vv_ortho,
                          double *pk_l_fNL_2_vd_ortho,
                          double *pk_l_fNL_2_dd_ortho,
                          double *pk_l_fNL_4_vv_ortho,
                          double *pk_l_fNL_4_vd_ortho,
                          double *pk_l_fNL_4_dd_ortho,
                          double *pk12_l_0_b1b2_ortho,
                          double *pk12_l_0_b2_ortho,
                          double *pk12_l_0_b1bG2_ortho,
                          double *pk12_l_0_bG2_ortho,
                          double *pk12_l_2_b1b2_ortho,
                          double *pk12_l_2_b2_ortho,
                          double *pk12_l_2_b1bG2_ortho,
                          double *pk12_l_2_bG2_ortho,
                          double *pk12_l_4_b1b2_ortho,
                          double *pk12_l_4_b2_ortho,
                          double *pk12_l_4_b1bG2_ortho,
                          double *pk12_l_4_bG2_ortho,
                          //GC: ORTHOGONAL -- finish
                          //GC!
                          double *lnk_l,
                          double *lnpk_l,
                          double *lnpPRIMk_l //GC!
                          );


    /**
     * Function definitions of the FFT tool, used by transfer.c and spectra.c
     * For more information, see fft.c
     */

#ifndef FFT_DEFINED
#define FFT_DEFINED
    /**
     * Compute FFT of two real inputs and store into two complex outputs of sizes N
     * Assumes arrays are allocated and of size N
     */
    void FFT_real(double* input_real_1, double* input_real_2, double* output_real_1, double* output_imag_1, double* output_real_2, double* output_imag_2, int N);
    /**
     * Compute FFT of two real inputs and store into two complex outputs of sizes N
     * Assumes arrays are allocated and of size N
     * Only gives up nonredundant part for real FFT which have c_(N-n)=c_n
     */
    void FFT_real_short(double* input_real_1, double* input_real_2, double* output_real_1, double* output_imag_1, double* output_real_2, double* output_imag_2, int N);
    /**
     * Compute FFT of single complex input and stores into single complex output of size N
     * Assumes arrays are allocated and of size N
     */
    void FFT(double* input_real, double* input_imag, double* output_real, double* output_imag, int N, int stepsize);
    void DCT_real(double* input_real,double* input_imag,double* output_real,double* output_imag,int N);
#endif


#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
/* @endcond */
