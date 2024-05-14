#pragma once

#include <matpack.h>

#include <functional>

#include "lin_alg.h"
#include "matpack_view.h"
#include "sorted_grid.h"

namespace disort {
struct BDRF {
  using func = std::function<void(ExhaustiveMatrixView,
                                  const ExhaustiveConstVectorView&,
                                  const ExhaustiveConstVectorView&)>;
  func f{[](auto, auto&, auto&) {
    throw std::runtime_error("BDRF function not set");
  }};

  void operator()(ExhaustiveMatrixView x,
                  const ExhaustiveConstVectorView& a,
                  const ExhaustiveConstVectorView& b) const {
    ARTS_ASSERT(x.nrows() == a.size());
    ARTS_ASSERT(x.ncols() == b.size());
    f(x, a, b);
  }

  Matrix operator()(const Vector& a, const Vector& b) const {
    Matrix x(a.size(), b.size());
    f(x, a, b);
    return x;
  }
};

struct mathscr_v_data {
  Matrix src;
  Matrix G;
  solve_workdata solve_work;
  Vector k1, k2;
  Vector cvec;

  mathscr_v_data(const Index Nk = 0, const Index Nc = 0)
      : src(Nk, Nc), G(Nk, Nk), solve_work(Nk), k1(Nk), k2(Nk), cvec(Nc) {}
  mathscr_v_data(const mathscr_v_data&) = default;
  mathscr_v_data(mathscr_v_data&&) = default;
  mathscr_v_data& operator=(const mathscr_v_data&) = default;
  mathscr_v_data& operator=(mathscr_v_data&&) = default;

  void resize(const Index Nk, const Index Nc) {
    src.resize(Nk, Nc);
    G.resize(Nk, Nk);
    solve_work.resize(Nk);
    k1.resize(Nk);
    k2.resize(Nk);
    cvec.resize(Nc);
  }
};

struct u_data {
  Matrix exponent;
  Matrix um;
  mathscr_v_data src;
  Vector intensities;
};

struct u0_data {
  mathscr_v_data src;
  Vector exponent;
  Vector u0;
};

struct tms_data {
  Vector nu;
  Vector TMS;
  Matrix mathscr_B;
  Matrix contribution_from_other_layers_pos;
  Matrix contribution_from_other_layers_neg;
};

struct flux_data {
  Vector exponent;
  mathscr_v_data src;
  Vector u0_pos;
  Vector u0_neg;
};

/** The main data structure for the DISORT algorithm
 * 
 * This pre-allocates all the memory needed for the computations
 *
 * Note to cite the original author of the algorithm when this is presented:
 * https://github.com/LDEO-CREW/Pythonic-DISORT
 *
 * Also note to present the original author of the new Legendre Gauss algorithm:
 * https://doi.org/10.1137/140954969
 */
class main_data {
  Index NLayers{0};
  Index NQuad{0};
  Index NLeg{0};
  Index NFourier{0};
  Index N{0};
  Index Nscoeffs{0};
  Index NLeg_all{0};
  Index NBDRF{0};
  bool has_source_poly{false};
  bool is_multilayer{false};
  bool has_beam_source{false};

  //! User inputs
  AscendingGrid tau_arr{};                 // [NLayers]
  Vector omega_arr{};                      // [NLayers]
  Vector f_arr{};                          // [NLayers]
  Matrix source_poly_coeffs{};             // [NLayers, Nscoeffs]
  Matrix Leg_coeffs_all{};                 // [NLayers, NLeg_all]
  Matrix b_pos{};                          // [NFourier, N]
  Matrix b_neg{};                          // [NFourier, N]
  std::vector<BDRF> brdf_fourier_modes{};  // [NBDRF]
  Numeric mu0{};
  Numeric I0{};
  Numeric phi0{};

  //! Derived values
  Vector scale_tau{};                   // [NLayers]
  Vector scaled_omega_arr{};            // [NLayers]
  Vector scaled_tau_arr_with_0{};       // [NLayers + 1]
  Vector mu_arr{};                      // [NQuad]
  Vector inv_mu_arr{};                  // [NQuad]
  Vector W{};                           // [N]
  Vector Leg_coeffs_residue_avg{};      // [NLeg_all]
  Matrix weighted_scaled_Leg_coeffs{};  // [NLayers, NLeg]
  Matrix weighted_Leg_coeffs_all{};     // [NLayers, NLeg_all]
  Tensor4 GC_collect{};                 // [NFourier, NLayers, NQuad, NQuad]
  Tensor4 G_collect{};                  // [NFourier, NLayers, NQuad, NQuad]
  Tensor3 K_collect{};                  // [NFourier, NLayers, NQuad]
  Tensor3 B_collect{};                  // [NFourier, NLayers, NQuad]
  Numeric I0_orig{};
  Numeric f_avg{};
  Numeric omega_avg{};
  Numeric scaled_mu0{};

  //! Internal compute data
  Index n{};                            // NQuad * NLayers;
  Vector RHS{};                         // [n]
  Vector jvec{};                        // [NQuad]
  Vector fac{};                         // [NLeg]
  Vector weighted_asso_Leg_coeffs_l{};  // [NLeg]
  Vector asso_leg_term_mu0{};           // [NLeg]
  Vector X_temp{};                      // [NLeg]
  Vector mathscr_X_pos{};               // [N]
  Vector E_Lm1L{};                      // [N]
  Vector E_lm1l{};                      // [N]
  Vector E_llp1{};                      // [N]
  Vector BDRF_RHS_contribution{};       // [N]
  Matrix Gml{};                         // [NQuad, NQuad]
  Matrix BDRF_LHS{};                    // [N, NQuad]
  Matrix R{};                           // [N, N]
  Matrix mathscr_D_neg{};               // [N, N]
  Matrix D_pos{};                       // [N, N]
  Matrix D_neg{};                       // [N, N]
  Matrix apb{};                         // [N, N]
  Matrix amb{};                         // [N, N]
  Matrix sqr{};                         // [N, N]
  Matrix asso_leg_term_pos{};           // [N, NLeg]
  Matrix asso_leg_term_neg{};           // [N, NLeg]
  Matrix D_temp{};                      // [N, NLeg]

  //! [NQuad / 2]
  solve_workdata solve_work{};

  //! [4 * N] + [N, N]
  diagonalize_workdata diag_work{};

  //! [n, 9 * N - 2] + [n / 2]
  matpack::band_matrix LHSB{};

  //! [NQuad, Nscoeffs] + [NQuad, NQuad] + 3 * [Nquad] + [Nscoeffs]
  mathscr_v_data comp_data{};

 public:
  main_data() = default;
  main_data(const main_data&) = default;
  main_data(main_data&&) = default;
  main_data& operator=(const main_data&) = default;
  main_data& operator=(main_data&&) = default;

  main_data(const Index NLayers,
            const Index NQuad,
            const Index NLeg,
            const Index NFourier,
            const Index Nscoeffs,
            const Index NLeg_all,
            const Index NBDRF);

  main_data(const Index NQuad,
            const Index NLeg,
            const Index NFourier,
            AscendingGrid tau_arr,
            Vector omega_arr,
            Matrix Leg_coeffs_all,
            Matrix b_pos,
            Matrix b_neg,
            Vector f_arr,
            Matrix source_poly_coeffs,
            std::vector<BDRF> brdf_fourier_modes,
            Numeric mu0,
            Numeric I0,
            Numeric phi0);

  /** Get the index of the tau value closest to the given tau
    *
    * Throws if tau is out-of-bounds
    *
    * Safe for parallel use
    *
    * @param tau The point-wise optical thickness
    * @return Index of the tau value closest to the given tau
    */
  [[nodiscard]] Index tau_index(const Numeric tau) const;

  /** Get the TMS correction factor
    *
    * Called by u_corr, exists for testing purposes
    *
    * Safe for parallel use if data is unique for each call
    *
    * @param data Compute data, the main output is in data.TMS
    * @param tau The point-wise optical thickness
    * @param phi The azimuthal angle of observation [0, 2 * pi)
    */
  void TMS(tms_data& data, const Numeric tau, const Numeric phi) const;

  /** Get the IMS correction factor
    *
    * Called by u_corr, exists for testing purposes
    *
    * Safe for parallel use if ims is unique for each call
    *
    * @param oms Compute data
    * @param tau The point-wise optical thickness
    * @param phi The azimuthal angle of observation [0, 2 * pi)
    */
  void IMS(Vector& ims, const Numeric tau, const Numeric phi) const;

  /** Spectral radiance at a given tau and phi
    *
    * The computations are done at the quadrature points
    *
    * Safe for parallel use if u_data is unique per thread
    * 
    * @param data Compute data, the spectral radiance is in data.intensities
    * @param tau The point-wise optical thickness 
    * @param phi The azimuthal angle of observation [0, 2 * pi)
    */
  void u(u_data& data, const Numeric tau, const Numeric phi) const;

  /** Spectral radiance at a given tau, only for the 0th Fourier mode
    *
    * The computations are done at the quadrature points
    *
    * Safe for parallel use if u0_data is unique per thread
    * 
    * @param data Compute data, the spectral radiance is in data.u0
    * @param tau The point-wise optical thickness 
    */
  void u0(u0_data& data, const Numeric tau) const;

  /** Spectral radiance at a given tau and phi
    *
    * The computations are done at the quadrature points
    *
    * Both the IMS and TMS corrections are applied
    *
    * Safe for parallel use if u_data, ims, and tms_data are unique per thread
    * 
    * @param u_data Compute data, the corrected spectral radiance is in data.intensities
    * @param ims The IMS correction compute data
    * @param tms_data The TMS correction compute data
    * @param tau The point-wise optical thickness 
    * @param phi The azimuthal angle of observation [0, 2 * pi)
    */
  void u_corr(u_data& u_data,
              Vector& ims,
              tms_data& tms_data,
              const Numeric tau,
              const Numeric phi) const;

  /** Compute the upward flux at a given tau
    *
    * Safe for parallel use if flux_data is unique per thread
    * 
    * @param tau The point-wise optical thickness 
    * @return Numeric value of the upward flux
    */
  [[nodiscard]] Numeric flux_up(flux_data&, const Numeric tau) const;

  void gridded_u(ExhaustiveTensor3View, const Vector& phi) const;
  void gridded_flux(ExhaustiveVectorView up,
                    ExhaustiveVectorView down,
                    ExhaustiveVectorView down_direct) const;

  void ungridded_u(ExhaustiveTensor3View out,
                   const AscendingGrid& tau,
                   const Vector& phi) const;
  void ungridded_flux(ExhaustiveVectorView flux_up,
                      ExhaustiveVectorView flux_do,
                      ExhaustiveVectorView flux_dd,
                      const AscendingGrid& tau) const;

  /** Compute the downward flux at a given tau
    *
    * Safe for parallel use if flux_data is unique per thread
    * 
    * @param tau The point-wise optical thickness 
    * @return std::pair<Numeric, Numeric> Diffuse and direct downward flux, respectively
    */
  [[nodiscard]] std::pair<Numeric, Numeric> flux_down(flux_data&,
                                                      const Numeric tau) const;

  /** Computes the IMS correction factors
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    * 
    * Depends on:
    * - omega_arr
    * - tau_arr
    * - f_arr
    * - mu0
    * - Leg_coeffs_all
    *
    * Modifies:
    * - scaled_mu0
    * - Leg_coeffs_residue_avg
    * - omega_avg
    * - f_avg
    */
  void set_ims_factors();

  /** Set the scaled objects
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    *
    * Depends on:
    * - omega_arr
    * - tau_arr
    * - f_arr
    * - Leg_coeffs_all
    *
    * Modifies:
    * - scale_tau
    * - scaled_tau_arr_with_0
    * - weighted_scaled_Leg_coeffs
    * - scaled_omega_arr
    */
  void set_scales();

  /** Diagonalize the system of equations
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    * 
    * Depends on:
    * - weighted_scaled_Leg_coeffs 
    * - scaled_omega_arr 
    * - mu0 
    * - I0
    *
    * Modifies:
    * - G_collect 
    * - K_collect 
    * - B_collect 
    */
  void diagonalize();

  /** Solves the system of equations
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    * 
    * Depends on:
    * - G_collect 
    * - K_collect 
    * - B_collect 
    * - source_poly_coeffs 
    * - b_pos 
    * - b_neg 
    * - tau_arr 
    * - scaled_tau_arr_with_0 
    * - brdf_fourier_modes 
    * - mu0 
    * - I0 
    *
    * Modifies:
    * - GC_collect 
    */
  void solve_for_coefs();

  /** Set the weighted Leg coeffs all object
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    * 
    * Depends on:
    * - Leg_coeffs_all
    *
    * Modifies:
    * - weighted_Leg_coeffs_all
    */
  void set_weighted_Leg_coeffs_all();

  /** Set the beam origin
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    * 
    * Depends on:
    * - b_pos
    *
    * Modifies:
    * - I0_orig
    * - I0
    * 
    * @param I0 The new beam intensity
    */
  void set_beam_source(const Numeric I0);

  //! Checks the input sizes
  void check_input_size() const;

  //! Checks the input values
  void check_input_value() const;

  /** Updates all the internal data
    *
    * This function should be called after changing any of the input values
    * to ensure that the internal data is consistent.
    *
    * Additionally, this calls check_input_value to ensure that the input
    * values are valid.
    *
    * Not safe for parallel use.
    * 
    * @param I0 The new beam intensity if it should be changed, otherwise -1
    */
  void update_all(const Numeric I0 = -1);

  //! The angles of quadrature - NQuad
  [[nodiscard]] auto&& mu() const { return mu_arr; }

  //! The weights of quadrature - N
  [[nodiscard]] auto&& weights() const { return W; }

  //! The optical thicknesses grid - NLayers
  [[nodiscard]] auto&& tau() const { return tau_arr; }

  //! The single scattering albedo - NLayers
  [[nodiscard]] auto&& omega() const { return omega_arr; }

  //! The fractional scattering into the peak - NLayers or 0
  [[nodiscard]] auto&& f() const { return f_arr; }

  //! Polynomial coefficients of isotropic internal sources - NLayers x Nscoeffs or 0 x 0
  [[nodiscard]] auto&& source_poly() const { return source_poly_coeffs; }

  //! Legendre coefficients of the scattering phase function (unweighted) - NLayers x NLeg_all
  [[nodiscard]] auto&& all_legendre_coeffs() const { return Leg_coeffs_all; }

  //! Positive Fourier coefficients of the beam source - 1 x 1 or N x NFourier
  [[nodiscard]] auto&& positive_boundary() const { return b_pos; }

  //! Negative Fourier coefficients of the beam source - 1 x 1 or N x NFourier
  [[nodiscard]] auto&& negative_boundary() const { return b_neg; }

  //! Fourier modes of the BRDF - NBDRF
  [[nodiscard]] auto&& brdf_modes() const { return brdf_fourier_modes; }

  //! The cosine of the beam source zenith angle
  [[nodiscard]] auto&& solar_zenith() const { return mu0; }

  //! The intensity of the beam source
  [[nodiscard]] auto&& beam_source() const { return I0; }

  //! The azimuthal angle of the beam source
  [[nodiscard]] auto&& beam_azimuth() const { return phi0; }

  //! The optical thicknesses grid - NLayers
  [[nodiscard]] auto tau() { return AscendingGridView{tau_arr}; }

  //! The single scattering albedo - NLayers
  [[nodiscard]] auto omega() { return ExhaustiveVectorView{omega_arr}; }

  //! The fractional scattering into the peak - NLayers or 0
  [[nodiscard]] auto f() { return ExhaustiveVectorView{f_arr}; }

  //! Polynomial coefficients of isotropic internal sources - NLayers x Nscoeffs or 0 x 0
  [[nodiscard]] auto source_poly() {
    return ExhaustiveMatrixView{source_poly_coeffs};
  }

  //! Legendre coefficients of the scattering phase function (unweighted) - NLayers x NLeg_all
  [[nodiscard]] auto all_legendre_coeffs() {
    return ExhaustiveMatrixView{Leg_coeffs_all};
  }

  //! Positive Fourier coefficients of the beam source - 1 x 1 or N x NFourier
  [[nodiscard]] auto positive_boundary() { return ExhaustiveMatrixView{b_pos}; }

  //! Negative Fourier coefficients of the beam source - 1 x 1 or N x NFourier
  [[nodiscard]] auto negative_boundary() { return ExhaustiveMatrixView{b_neg}; }

  //! Fourier modes of the BRDF - NBDRF
  [[nodiscard]] auto brdf_modes() { return std::span{brdf_fourier_modes}; }

  //! The cosine of the beam source zenith angle
  [[nodiscard]] Numeric& solar_zenith() { return mu0; }

  /** The intensity of the beam source
    * 
    * This function is disabled, use set_beam_source instead, or if
    * it is part of a larger update, use update_all.
    *
    * If you really want the beam source intensity, you can use the
    * const version of this function by first const-casting your object.
    */
  [[nodiscard]] Numeric& beam_source() = delete;

  //! The azimuthal angle of the beam source
  [[nodiscard]] Numeric& beam_azimuth() { return phi0; }

  friend std::ostream& operator<<(std::ostream& os, const main_data&) {
    return os << "Not implemented, output stream operator\n";
  }
};
}  // namespace disort
