#include <disort-test.h>

#include "artstime.h"
#include "disort.h"
#include "matpack_data.h"
#include "matpack_math.h"

void test_11a_1layer() try {
  const AscendingGrid tau_arr{8.};
  const Vector omega_arr{0.999999};
  const Index NQuad = 16;
  const Matrix Leg_coeffs_all{
      Vector{1.00000000e+00, 7.50000000e-01, 5.62500000e-01, 4.21875000e-01,
             3.16406250e-01, 2.37304688e-01, 1.77978516e-01, 1.33483887e-01,
             1.00112915e-01, 7.50846863e-02, 5.63135147e-02, 4.22351360e-02,
             3.16763520e-02, 2.37572640e-02, 1.78179480e-02, 1.33634610e-02,
             1.00225958e-02, 7.51694682e-03, 5.63771011e-03, 4.22828259e-03,
             3.17121194e-03, 2.37840895e-03, 1.78380672e-03, 1.33785504e-03,
             1.00339128e-03, 7.52543458e-04, 5.64407594e-04, 4.23305695e-04,
             3.17479271e-04, 2.38109454e-04, 1.78582090e-04, 1.33936568e-04}
          .reshape(tau_arr.nelem(), 32)};

  const Numeric mu0 = 0.6;
  const Numeric I0 = Constant::pi / mu0;
  const Numeric phi0 = 0.9 * Constant::pi;
  Matrix b_neg(NQuad, NQuad / 2, 0);
  b_neg[0] = 1;
  Matrix b_pos(NQuad, NQuad / 2, 0);
  b_pos[0] = 1;
  const std::vector<disort::BDRF> BDRF_Fourier_modes{
      disort::BDRF{[](auto c, auto&, auto&) { c = 1; }}};
  const Matrix s_poly_coeffs{
      Vector{172311.79936609, -102511.4417051}.reshape(tau_arr.nelem(), 2)};
  const Vector f_arr{Leg_coeffs_all(joker, NQuad)};

  // Optional (unused)
  const Index NLeg = NQuad;
  const Index NFourier = NQuad;

  const disort::main_data dis(NQuad,
                              NLeg,
                              NFourier,
                              tau_arr,
                              omega_arr,
                              Leg_coeffs_all,
                              b_pos,
                              b_neg,
                              f_arr,
                              s_poly_coeffs,
                              BDRF_Fourier_modes,
                              mu0,
                              I0,
                              phi0);

  const Vector taus{Vector{
      0.06354625877794251,
      0.6354625877794251,
      6.354625877794252,
  }
                        .reshape_as(3)};

  const Vector phis{Vector{
      0.0,
      1.5705463267948965,
      3.141092653589793,
      4.71163898038469,
      6.282185307179586,
  }
                        .reshape_as(5)};

  const Tensor3 u{Vector{
      -542007159849.03156, -664949102438.8613,  -839637987333.5054,
      -1035805674770.3422, -1231938613823.8433, -1407681199351.6865,
      -1543638602703.844,  -1624366141667.5176, -452072110490.451,
      -176463523131.0405,  -64325918641.56033,  -27040816484.165386,
      -13603684906.620405, -8201190412.133552,  -5846358435.193959,
      -4869771665.426168,  -790242606575.9447,  -881438825714.7039,
      -1026942140062.2853, -1204154478852.7078, -1389390648390.349,
      -1559098349248.1858, -1691599509314.7083, -1770464094682.6465,
      -744691625289.5955,  -646223522460.6571,  -464904941223.2196,
      -284313241449.4057,  -170356771941.81042, -110255227734.72589,
      -80607006071.84906,  -67597179969.75537,  -2265179672683.9697,
      -2348927474365.532,  -2490386663339.6855, -2670429930683.9067,
      -2841162286327.642,  -2964480785488.3647, -3038400674926.433,
      -3074416823835.482,  -2224855371700.144,  -2142340847702.23,
      -2006946607264.6624, -1837837321118.4453, -1658028213256.2122,
      -1491973286930.4834, -1362420547412.6404, -1285715077306.0273,
      -542007159848.76764, -664949102438.5804,  -839637987333.2377,
      -1035805674770.1227, -1231938613823.6829, -1407681199351.5806,
      -1543638602703.7812, -1624366141667.4907, -452072110490.2173,
      -176463523130.93448, -64325918641.509026, -27040816484.13253,
      -13603684906.594194, -8201190412.11028,   -5846358435.174582,
      -4869771665.414905,  -790242606575.6283,  -881438825714.4117,
      -1026942140062.0413, -1204154478852.5247, -1389390648390.2236,
      -1559098349248.1057, -1691599509314.662,  -1770464094682.6272,
      -744691625289.2698,  -646223522460.3213,  -464904941222.9079,
      -284313241449.1437,  -170356771941.588,   -110255227734.53209,
      -80607006071.6926,   -67597179969.66686,  -2265179672683.9546,
      -2348927474365.52,   -2490386663339.677,  -2670429930683.901,
      -2841162286327.6387, -2964480785488.363,  -3038400674926.4326,
      -3074416823835.4814, -2224855371700.1284, -2142340847702.2112,
      -2006946607264.637,  -1837837321118.4104, -1658028213256.1633,
      -1491973286930.4246, -1362420547412.5842, -1285715077305.9954,
      -542007159847.95764, -664949102437.8752,  -839637987332.7219,
      -1035805674769.7993, -1231938613823.5005, -1407681199351.485,
      -1543638602703.7366, -1624366141667.475,  -452072110489.41174,
      -176463523130.47278, -64325918641.15577,  -27040816483.72276,
      -13603684906.123302, -8201190411.859676,  -5846358435.104522,
      -4869771665.4016905, -790242606574.8704,  -881438825713.8123,
      -1026942140061.6403, -1204154478852.287,  -1389390648390.0933,
      -1559098349248.0383, -1691599509314.6304, -1770464094682.6157,
      -744691625288.4229,  -646223522459.271,   -464904941221.52246,
      -284313241447.2442,  -170356771939.27344, -110255227733.15106,
      -80607006071.24594,  -67597179969.57239,  -2265179672683.9453,
      -2348927474365.5127, -2490386663339.671,  -2670429930683.898,
      -2841162286327.637,  -2964480785488.362,  -3038400674926.432,
      -3074416823835.481,  -2224855371700.1177, -2142340847702.1982,
      -2006946607264.6184, -1837837321118.3809, -1658028213256.117,
      -1491973286930.3638, -1362420547412.535,  -1285715077305.974,
      -542007159848.94543, -664949102438.7667,  -839637987333.4106,
      -1035805674770.2592, -1231938613823.7783, -1407681199351.641,
      -1543638602703.8157, -1624366141667.5042, -452072110490.3762,
      -176463523131.00784, -64325918641.54545,  -27040816484.156708,
      -13603684906.61382,  -8201190412.127581,  -5846358435.188248,
      -4869771665.421796,  -790242606575.8364,  -881438825714.6008,
      -1026942140062.1954, -1204154478852.6357, -1389390648390.2966,
      -1559098349248.1511, -1691599509314.6865, -1770464094682.6367,
      -744691625289.4858,  -646223522460.5479,  -464904941223.12286,
      -284313241449.3294,  -170356771941.74823, -110255227734.67093,
      -80607006071.79932,  -67597179969.71972,  -2265179672683.9624,
      -2348927474365.527,  -2490386663339.6816, -2670429930683.904,
      -2841162286327.6406, -2964480785488.364,  -3038400674926.4326,
      -3074416823835.4814, -2224855371700.137,  -2142340847702.221,
      -2006946607264.6506, -1837837321118.4294, -1658028213256.1912,
      -1491973286930.4578, -1362420547412.6155, -1285715077306.012,
      -542007159849.0316,  -664949102438.8616,  -839637987333.5055,
      -1035805674770.3422, -1231938613823.843,  -1407681199351.6863,
      -1543638602703.844,  -1624366141667.5176, -452072110490.45105,
      -176463523131.04053, -64325918641.56033,  -27040816484.165386,
      -13603684906.62041,  -8201190412.133553,  -5846358435.193962,
      -4869771665.426171,  -790242606575.9447,  -881438825714.7039,
      -1026942140062.2854, -1204154478852.7078, -1389390648390.349,
      -1559098349248.1858, -1691599509314.7083, -1770464094682.6465,
      -744691625289.5955,  -646223522460.6572,  -464904941223.21954,
      -284313241449.40576, -170356771941.81036, -110255227734.72594,
      -80607006071.84904,  -67597179969.75538,  -2265179672683.9697,
      -2348927474365.532,  -2490386663339.6855, -2670429930683.9067,
      -2841162286327.642,  -2964480785488.3647, -3038400674926.433,
      -3074416823835.482,  -2224855371700.144,  -2142340847702.23,
      -2006946607264.6624, -1837837321118.4453, -1658028213256.2122,
      -1491973286930.4834, -1362420547412.6404, -1285715077306.0273,
  }
                      .reshape_as(5, 3, 16)};

  const Matrix u0{Vector{
      -542007159848.68823, -664949102438.5311,  -839637987333.224,
      -1035805674770.1328, -1231938613823.7017, -1407681199351.5989,
      -1543638602703.7942, -1624366141667.497,  -452072110490.1273,
      -176463523130.873,   -64325918641.45187,  -27040816484.0495,
      -13603684906.495117, -8201190412.062683,  -5846358435.166687,
      -4869771665.41626,   -790242606575.5796,  -881438825714.389,
      -1026942140062.0438, -1204154478852.54,   -1389390648390.241,
      -1559098349248.12,   -1691599509314.6719, -1770464094682.6316,
      -744691625289.2051,  -646223522460.2166,  -464904941222.7234,
      -284313241448.8089,  -170356771941.14282, -110255227734.29614,
      -80607006071.65399,  -67597179969.6792,   -2265179672683.958,
      -2348927474365.523,  -2490386663339.6787, -2670429930683.9023,
      -2841162286327.6396, -2964480785488.3633, -3038400674926.4326,
      -3074416823835.4814, -2224855371700.132,  -2142340847702.215,
      -2006946607264.642,  -1837837321118.4165, -1658028213256.1711,
      -1491973286930.4329, -1362420547412.594,  -1285715077306.0022,
  }
                      .reshape_as(3, 16)};

  const Vector flux_down_diffuse{Vector{
      -63531167081.666534,
      -560958120277.2043,
      -4985647758594.26,
  }
                                     .reshape_as(3)};

  const Vector flux_down_direct{Vector{
      14.796267664990843,
      5.704076430616008,
      0.00041353941227710225,
  }
                                    .reshape_as(3)};

  const Vector flux_up{Vector{
      -4092716310876.2905,
      -4590149266091.762,
      -9014935283001.768,
  }
                           .reshape_as(3)};

  //flat_print(u, compute_u(dis, taus, phis, true) );
  //const auto [flux_up_, flux_down_diffuse_, flux_down_direct_] =  compute_flux(dis, taus);
  //flat_print(flux_up, flux_up_);

  compare("test_11a-1layer",
          dis,
          taus,
          phis,
          u,
          u0,
          flux_down_diffuse,
          flux_down_direct,
          flux_up,
          true);
} catch (std::exception& e) {
  throw std::runtime_error(var_string("Error in test-11a-1layer:\n", e.what()));
}

void test_11a_multilayer() try {
  const AscendingGrid tau_arr{
      0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.};
  const Vector omega_arr(tau_arr.size(), 0.999999);
  const Index NQuad = 16;
  Matrix Leg_coeffs_all(tau_arr.size(), 32);
  for (auto&& v : Leg_coeffs_all)
    v = {1.00000000e+00, 7.50000000e-01, 5.62500000e-01, 4.21875000e-01,
         3.16406250e-01, 2.37304688e-01, 1.77978516e-01, 1.33483887e-01,
         1.00112915e-01, 7.50846863e-02, 5.63135147e-02, 4.22351360e-02,
         3.16763520e-02, 2.37572640e-02, 1.78179480e-02, 1.33634610e-02,
         1.00225958e-02, 7.51694682e-03, 5.63771011e-03, 4.22828259e-03,
         3.17121194e-03, 2.37840895e-03, 1.78380672e-03, 1.33785504e-03,
         1.00339128e-03, 7.52543458e-04, 5.64407594e-04, 4.23305695e-04,
         3.17479271e-04, 2.38109454e-04, 1.78582090e-04, 1.33936568e-04};

  const Numeric mu0 = 0.6;
  const Numeric I0 = Constant::pi / mu0;
  const Numeric phi0 = 0.9 * Constant::pi;
  Matrix b_neg(NQuad, NQuad / 2, 0);
  b_neg[0] = 1;
  Matrix b_pos(NQuad, NQuad / 2, 0);
  b_pos[0] = 1;
  const std::vector<disort::BDRF> BDRF_Fourier_modes{
      disort::BDRF{[](auto c, auto&, auto&) { c = 1; }}};
  Matrix s_poly_coeffs(tau_arr.size(), 2);
  for (auto&& v : s_poly_coeffs) v = {172311.79936609, -102511.4417051};
  const Vector f_arr{Leg_coeffs_all(joker, NQuad)};

  // Optional (unused)
  const Index NLeg = NQuad;
  const Index NFourier = NQuad;

  const disort::main_data dis(NQuad,
                              NLeg,
                              NFourier,
                              tau_arr,
                              omega_arr,
                              Leg_coeffs_all,
                              b_pos,
                              b_neg,
                              f_arr,
                              s_poly_coeffs,
                              BDRF_Fourier_modes,
                              mu0,
                              I0,
                              phi0);

  const Vector taus{Vector{
      0.06354625877794251,
      0.6354625877794251,
      6.354625877794252,
  }
                        .reshape_as(3)};

  const Vector phis{Vector{
      0.0,
      1.5705463267948965,
      3.141092653589793,
      4.71163898038469,
      6.282185307179586,
  }
                        .reshape_as(5)};

  const Tensor3 u{Vector{
      -542007159849.03156, -664949102438.8613,  -839637987333.5054,
      -1035805674770.3422, -1231938613823.8433, -1407681199351.6865,
      -1543638602703.844,  -1624366141667.5176, -452072110490.451,
      -176463523131.0405,  -64325918641.56033,  -27040816484.165386,
      -13603684906.620405, -8201190412.133552,  -5846358435.193959,
      -4869771665.426168,  -790242606575.9447,  -881438825714.7039,
      -1026942140062.2853, -1204154478852.7078, -1389390648390.349,
      -1559098349248.1858, -1691599509314.7083, -1770464094682.6465,
      -744691625289.5955,  -646223522460.6571,  -464904941223.2196,
      -284313241449.4057,  -170356771941.81042, -110255227734.72589,
      -80607006071.84906,  -67597179969.75537,  -2265179672683.9697,
      -2348927474365.532,  -2490386663339.6855, -2670429930683.9067,
      -2841162286327.642,  -2964480785488.3647, -3038400674926.433,
      -3074416823835.482,  -2224855371700.144,  -2142340847702.23,
      -2006946607264.6624, -1837837321118.4453, -1658028213256.2122,
      -1491973286930.4834, -1362420547412.6404, -1285715077306.0273,
      -542007159848.76764, -664949102438.5804,  -839637987333.2377,
      -1035805674770.1227, -1231938613823.6829, -1407681199351.5806,
      -1543638602703.7812, -1624366141667.4907, -452072110490.2173,
      -176463523130.93448, -64325918641.509026, -27040816484.13253,
      -13603684906.594194, -8201190412.11028,   -5846358435.174582,
      -4869771665.414905,  -790242606575.6283,  -881438825714.4117,
      -1026942140062.0413, -1204154478852.5247, -1389390648390.2236,
      -1559098349248.1057, -1691599509314.662,  -1770464094682.6272,
      -744691625289.2698,  -646223522460.3213,  -464904941222.9079,
      -284313241449.1437,  -170356771941.588,   -110255227734.53209,
      -80607006071.6926,   -67597179969.66686,  -2265179672683.9546,
      -2348927474365.52,   -2490386663339.677,  -2670429930683.901,
      -2841162286327.6387, -2964480785488.363,  -3038400674926.4326,
      -3074416823835.4814, -2224855371700.1284, -2142340847702.2112,
      -2006946607264.637,  -1837837321118.4104, -1658028213256.1633,
      -1491973286930.4246, -1362420547412.5842, -1285715077305.9954,
      -542007159847.95764, -664949102437.8752,  -839637987332.7219,
      -1035805674769.7993, -1231938613823.5005, -1407681199351.485,
      -1543638602703.7366, -1624366141667.475,  -452072110489.41174,
      -176463523130.47278, -64325918641.15577,  -27040816483.72276,
      -13603684906.123302, -8201190411.859676,  -5846358435.104522,
      -4869771665.4016905, -790242606574.8704,  -881438825713.8123,
      -1026942140061.6403, -1204154478852.287,  -1389390648390.0933,
      -1559098349248.0383, -1691599509314.6304, -1770464094682.6157,
      -744691625288.4229,  -646223522459.271,   -464904941221.52246,
      -284313241447.2442,  -170356771939.27344, -110255227733.15106,
      -80607006071.24594,  -67597179969.57239,  -2265179672683.9453,
      -2348927474365.5127, -2490386663339.671,  -2670429930683.898,
      -2841162286327.637,  -2964480785488.362,  -3038400674926.432,
      -3074416823835.481,  -2224855371700.1177, -2142340847702.1982,
      -2006946607264.6184, -1837837321118.3809, -1658028213256.117,
      -1491973286930.3638, -1362420547412.535,  -1285715077305.974,
      -542007159848.94543, -664949102438.7667,  -839637987333.4106,
      -1035805674770.2592, -1231938613823.7783, -1407681199351.641,
      -1543638602703.8157, -1624366141667.5042, -452072110490.3762,
      -176463523131.00784, -64325918641.54545,  -27040816484.156708,
      -13603684906.61382,  -8201190412.127581,  -5846358435.188248,
      -4869771665.421796,  -790242606575.8364,  -881438825714.6008,
      -1026942140062.1954, -1204154478852.6357, -1389390648390.2966,
      -1559098349248.1511, -1691599509314.6865, -1770464094682.6367,
      -744691625289.4858,  -646223522460.5479,  -464904941223.12286,
      -284313241449.3294,  -170356771941.74823, -110255227734.67093,
      -80607006071.79932,  -67597179969.71972,  -2265179672683.9624,
      -2348927474365.527,  -2490386663339.6816, -2670429930683.904,
      -2841162286327.6406, -2964480785488.364,  -3038400674926.4326,
      -3074416823835.4814, -2224855371700.137,  -2142340847702.221,
      -2006946607264.6506, -1837837321118.4294, -1658028213256.1912,
      -1491973286930.4578, -1362420547412.6155, -1285715077306.012,
      -542007159849.0316,  -664949102438.8616,  -839637987333.5055,
      -1035805674770.3422, -1231938613823.843,  -1407681199351.6863,
      -1543638602703.844,  -1624366141667.5176, -452072110490.45105,
      -176463523131.04053, -64325918641.56033,  -27040816484.165386,
      -13603684906.62041,  -8201190412.133553,  -5846358435.193962,
      -4869771665.426171,  -790242606575.9447,  -881438825714.7039,
      -1026942140062.2854, -1204154478852.7078, -1389390648390.349,
      -1559098349248.1858, -1691599509314.7083, -1770464094682.6465,
      -744691625289.5955,  -646223522460.6572,  -464904941223.21954,
      -284313241449.40576, -170356771941.81036, -110255227734.72594,
      -80607006071.84904,  -67597179969.75538,  -2265179672683.9697,
      -2348927474365.532,  -2490386663339.6855, -2670429930683.9067,
      -2841162286327.642,  -2964480785488.3647, -3038400674926.433,
      -3074416823835.482,  -2224855371700.144,  -2142340847702.23,
      -2006946607264.6624, -1837837321118.4453, -1658028213256.2122,
      -1491973286930.4834, -1362420547412.6404, -1285715077306.0273,
  }
                      .reshape_as(5, 3, 16)};

  const Matrix u0{Vector{
      -542007159848.704,   -664949102438.5885,  -839637987333.277,
      -1035805674770.2255, -1231938613823.7686, -1407681199351.6738,
      -1543638602703.9316, -1624366141667.4453, -452072110490.16077,
      -176463523130.90643, -64325918641.467865, -27040816484.100464,
      -13603684906.434326, -8201190412.101501,  -5846358435.129089,
      -4869771665.520386,  -790242606575.5867,  -881438825714.4045,
      -1026942140062.0756, -1204154478852.571,  -1389390648390.2861,
      -1559098349248.1511, -1691599509314.7944, -1770464094682.6614,
      -744691625289.238,   -646223522460.2239,  -464904941222.7064,
      -284313241448.8547,  -170356771941.12115, -110255227734.29156,
      -80607006071.65015,  -67597179969.75745,  -2265179672683.9707,
      -2348927474365.5474, -2490386663339.7075, -2670429930683.9473,
      -2841162286327.6465, -2964480785488.459,  -3038400674926.45,
      -3074416823835.5444, -2224855371700.18,   -2142340847702.2676,
      -2006946607264.6577, -1837837321118.4775, -1658028213256.171,
      -1491973286930.4614, -1362420547412.6282, -1285715077306.039,
  }
                      .reshape_as(3, 16)};

  const Vector flux_down_diffuse{Vector{
      -63531167081.666534,
      -560958120277.2043,
      -4985647758594.26,
  }
                                     .reshape_as(3)};

  const Vector flux_down_direct{Vector{
      14.796267664990843,
      5.704076430616008,
      0.00041353941227710225,
  }
                                    .reshape_as(3)};

  const Vector flux_up{Vector{
      -4092716310876.2905,
      -4590149266091.762,
      -9014935283001.768,
  }
                           .reshape_as(3)};

  //flat_print(u, compute_u(dis, taus, phis, true) );
  //  const auto [flux_up_, flux_down_diffuse_, flux_down_direct_] =  compute_flux(dis, taus);
  //flat_print(flux_down_diffuse, flux_down_diffuse_);

  compare("test_11a-multilayer",
          dis,
          taus,
          phis,
          u,
          u0,
          flux_down_diffuse,
          flux_down_direct,
          flux_up,
          true);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("Error in test-11a-multilayer:\n", e.what()));
}

int main() try {
  std::cout << std::setprecision(16);
  test_11a_1layer();
  test_11a_multilayer();
} catch (std::exception& e) {
  std::cerr << "Error in main:\n" << e.what() << '\n';
  return EXIT_FAILURE;
}