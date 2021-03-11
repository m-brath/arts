/* Copyright (C) 2019

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA.
*/

/*!
  \file   newdoit.cc
  \author Manfred Brath manfred.brath@uni-hamburg.de
  \date   06/09/2019
  
  \brief  This file will contain the internal functions
          to calculate the radiative transfer
          inside the cloudbox using the newDoit method.
  
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include "newdoit.h"
#include "agenda_class.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpackVII.h"
#include "messages.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "propagationmatrix.h"
#include "rte.h"
#include "sorting.h"
#include "special_interp.h"
#include "xml_io.h"

extern const Numeric PI;
extern const Numeric RAD2DEG;
extern const Numeric DEG2RAD;
extern const Numeric SPEED_OF_LIGHT;

void RTDomainScatteringProperties::CalcScatteredField(
    Tensor6& cloudbox_scat_field,
    //WS Input:
    const Tensor6& cloudbox_field_mono,
    const Index& atmosphere_dim,
    const Verbosity& verbosity) {
  if (atmosphere_dim == 1) {
    CalcScatteredField1D(cloudbox_scat_field, cloudbox_field_mono, verbosity);

  } else {
    CalcScatteredField3D(cloudbox_scat_field, cloudbox_field_mono, verbosity);
  }
}

void RTDomainScatteringProperties::CalcScatteredField1D(
    // Output
    Tensor6& cloudbox_scat_field,
    // Input:
    const ConstTensor6View& cloudbox_field_mono,
    const Verbosity& verbosity) {
  CREATE_OUT3;

  // Number of zenith angles.
  const Index Npza = mPdirIdx0.nelem();
  const Index Niza = mScatZaGrid.nelem();

  // Get stokes dimension from *cloudbox_scat_field*:
  const Index stokes_dim = cloudbox_field_mono.ncols();

  // and pressure dimension
  const Index Np = cloudbox_field_mono.nvitrines();

  //Prepare intensity field interpolated on equidistant grid.
  Matrix cloudbox_field_int(Niza, stokes_dim, 0);

  //Prepare product field
  Matrix product_field(Niza, stokes_dim, 0);

  for (Index i_p = 0; i_p < Np; i_p++) {
    // Interpolate intensity field:
    // Only one interpolation is required. We have to interpolate the
    // intensity field on the equidistant grid. There is no further interpolation
    // needed, because we evaluate the scattering matrix directly for the propagation/
    // outgoing direction, which is the direction of the actual radiation field.
    for (Index i = 0; i < stokes_dim; i++) {
      interp(cloudbox_field_int(joker, i),
             mItw(joker, 0, joker),
             cloudbox_field_mono(i_p, 0, 0, joker, 0, i),
             mGpZaI);
    }

    //There is only loop over zenith angle grid; no azimuth angle grid.
    for (Index i_pza = 0; i_pza < Npza; i_pza++) {
      product_field = 0;

      //as we did the azimuth integration already during scattering data preparation
      //we just need to do the zenith integration
      for (Index i_iza = 0; i_iza < Niza; i_iza++) {
        // Multiplication of scattering matrix with incoming
        // intensity field.

        for (Index i = 0; i < stokes_dim; i++) {
          for (Index j = 0; j < stokes_dim; j++) {
            product_field(i_iza, i) +=
                mScatteringMatrix(i_p, 0, 0, i_pza, i_iza, i, j) *
                cloudbox_field_int(i_iza, j);
          }
        }
      }

      out3 << "Compute integral. \n";
      for (Index i = 0; i < stokes_dim; i++) {
        cloudbox_scat_field(i_p, 0, 0, i_pza, 0, i) =
            AngIntegrate_trapezoid(product_field(joker, i), mScatZaGrid) / 2 /
            PI;
      }
    }
  }

  out3 << "DONE calculating scattering field 1d. \n";
}

void RTDomainScatteringProperties::CalcScatteredField3D(
    // Output
    Tensor6& cloudbox_scat_field,
    // Input:
    const ConstTensor6View& cloudbox_field_mono,
    const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  // Number of incoming zenith angles.
  const Index Niza = mScatZaGrid.nelem();

  // Number of incoming azimuth angles.
  const Index Niaa = mScatAaGrid.nelem();

  // Get stokes dimension from *cloudbox_scat_field*:
  const Index stokes_dim = cloudbox_field_mono.ncols();

  const Index Np = cloudbox_field_mono.nvitrines();
  const Index Nlat = cloudbox_field_mono.nshelves();
  const Index Nlon = cloudbox_field_mono.nbooks();

  //Prepare intensity field interpolated on equidistant grid.
  Tensor3 cloudbox_field_int(Niza, Niaa, stokes_dim, 0);

  //Grid stepsize of zenith and azimuth angle grid, these are needed for the
  //integration function.
  Vector grid_stepsize(2);
  grid_stepsize[0] = 180. / (Numeric)(Niza - 1);
  grid_stepsize[1] = 360. / (Numeric)(Niaa - 1);

  //
  Tensor3 product_field(Niza, Niaa, stokes_dim, 0);

  for (Index i_p = 0; i_p < Np; i_p++) {
    for (Index i_lat = 0; i_lat < Nlat; i_lat++) {
      for (Index i_lon = 0; i_lon < Nlon; i_lon++) {
        // Interpolate intensity field:
        // Only one interpolation is required. We have to interpolate the
        // intensity field on the equidistant grid. There is no further interpolation
        // needed, because we evaluate the scattering matrix directly for the propagation/
        // outgoing direction, which is the direction of the actual radiation field.
        for (Index i = 0; i < stokes_dim; i++) {
          interp(cloudbox_field_int(joker, joker, i),
                 mItw,
                 cloudbox_field_mono(i_p, i_lat, i_lon, joker, joker, i),
                 mGpZaI,
                 mGpAaI);
        }

        //loop over outgoing directions
        for (Index i_prop = 0; i_prop < mPdirIdx0.nelem(); i_prop++) {
          //loop over incoming directions
          for (Index i_inc = 0; i_inc < mIdirIdx0.nelem(); i_inc++) {
            //Calculate the product of scattering matrix and incoming radiation field
            for (Index i = 0; i < stokes_dim; i++) {
              for (Index j = 0; j < stokes_dim; j++) {
                product_field(mIdirIdx0[i_inc], mIdirIdx1[i_inc], i) +=
                    mScatteringMatrix(i_p,
                                      i_lon,
                                      i_lat,
                                      mIdirIdx0[i_inc],
                                      mIdirIdx1[i_inc],
                                      i,
                                      j) *
                    cloudbox_field_int(mIdirIdx0[i_inc], mIdirIdx1[i_inc], j);
              }
            }
          }

          //Integrate over incoming zenith and azimuth direction
          out3 << "Compute integral. \n";
          for (Index i = 0; i < stokes_dim; i++) {
            cloudbox_scat_field(
                i_p, i_lat, i_lon, mPdirIdx0[i_prop], mPdirIdx1[i_prop], i) =
                AngIntegrate_trapezoid_opti(product_field(joker, joker, i),
                                            mScatZaGrid,
                                            mScatAaGrid,
                                            grid_stepsize);
          }
        }
      }
    }
  }
}

void Initialize_cloudbox_field(  //Output
    Tensor7& cloudbox_field,
    // WS Input
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Vector& f_grid,
    const Vector& rt_za_grid,
    const Vector& rt_aa_grid,
    const ArrayOfIndex& cloudbox_limits) {
  //------------- end of checks ---------------------------------------

  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index Nrza = rt_za_grid.nelem();
  const Index Ns = stokes_dim;

  // Resize and initialize radiation field in the cloudbox
  if (atmosphere_dim == 1) {
    cloudbox_field.resize(Nf, Np_cloud, 1, 1, Nrza, 1, Ns);
  } else if (atmosphere_dim == 3) {
    const Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
    const Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
    const Index Nraa = rt_aa_grid.nelem();

    cloudbox_field.resize(Nf, Np_cloud, Nlat_cloud, Nlon_cloud, Nrza, Nraa, Ns);
  }

  cloudbox_field = NAN;
}

void SetAngularGrids(Vector& za_grid,
                     Vector& aa_grid,
                     Vector& scat_za_grid,
                     Vector& scat_aa_grid,
                     // Keywords:
                     const Index& N_za_grid,
                     const Index& N_aa_grid,
                     const Index& N_scat_za_grid,
                     const Index& N_scat_aa_grid,
                     const String& za_grid_type) {
  // Azimuth angle grid
  if (N_aa_grid > 1)
    nlinspace(aa_grid, 0, 360, N_aa_grid);
  else if (N_aa_grid < 1) {
    ostringstream os;
    os << "N_aa_grid must be > 0 (even for 1D).";
    throw std::runtime_error(os.str());
  } else {
    aa_grid.resize(1);
    aa_grid[0] = 0.;
  }

  // Azimuth angle grid (scattering)
  if (N_scat_aa_grid > 1)
    nlinspace(scat_aa_grid, 0, 360, N_scat_aa_grid);
  else if (N_scat_aa_grid < 1) {
    ostringstream os;
    os << "N_scat_aa_grid must be > 0 (even for 1D).";
    throw std::runtime_error(os.str());
  } else {
    scat_aa_grid.resize(1);
    scat_aa_grid[0] = 0.;
  }

  // Make scattering zenith angle
  nlinspace(scat_za_grid, 0, 180, N_scat_za_grid);

  //calculate zenith angle grid
  za_grid.resize(N_za_grid);
  za_grid = 0.;
  if (za_grid_type == "linear") {
    nlinspace(za_grid, 0, 180, N_za_grid);

  } else if (za_grid_type == "linear_mu") {
    Vector x;
    nlinspace(x, -1, 1, N_za_grid);

    //allocate
    Vector za_grid_temp;
    za_grid_temp.resize(x.nelem());

    for (Index i = 0; i < N_za_grid; i++) {
      za_grid_temp[i] = acos(x[i]) / DEG2RAD;
    }

    //#sort weights and theta in increasing direction of scat_za_grid
    za_grid = za_grid_temp[Range(x.nelem() - 1, x.nelem(), -1)];

    //dirty hack to be assure that 90 deg is exactly 90 deg.
    if (N_za_grid % 2){
      za_grid[int(N_za_grid/2)]=90.;
    }


  } else {
    ostringstream os;
    os << "The selected grid type is not implemented";
    throw std::runtime_error(os.str());
  }

  //be sure that the first and the last angle are within the closed interval
  //between 0 and 180 deg, because ARTS is picky if the angles are due to numerics
  // slightly below and above,respectively.
  if (za_grid[0] < 0) {
    za_grid[0] = 0.;
  }

  if (scat_za_grid[0] < 0) {
    scat_za_grid[0] = 0.;
  }

  if (za_grid[za_grid.nelem() - 1] > 180) {
    za_grid[za_grid.nelem() - 1] = 180.;
  }

  if (scat_za_grid[scat_za_grid.nelem() - 1] > 180) {
    scat_za_grid[scat_za_grid.nelem() - 1] = 180.;
  }
}

void SetClearsky_cloudbox(Tensor7& cloudbox_field,
                          const Vector& f_grid,
                          const Vector& p_grid,
                          const Vector& lat_grid,
                          const Vector& lon_grid,
                          const ArrayOfIndex& cloudbox_limits,
                          const Index& atmosphere_dim,
                          const Verbosity& verbosity) {
  CREATE_OUT2;

  out2
      << "  Interpolate boundary clearsky field to obtain the initial field.\n";

  // Initial field only needs to be calculated from clearsky field for the
  // first frequency. For the next frequencies the solution field from the
  // previous frequencies is used.
  if (atmosphere_dim == 1) {
    const Index nf = f_grid.nelem() ? cloudbox_field.nlibraries() : 1;

    for (Index f_index = 0; f_index < nf; f_index++) {
      Index N_za = cloudbox_field.npages();
      Index N_aa = cloudbox_field.nrows();
      Index N_i = cloudbox_field.ncols();

      //1. interpolation - pressure grid

      /*the old grid is having only two elements, corresponding to the
             cloudbox_limits and the new grid have elements corresponding to
             all grid points inside the cloudbox plus the cloud_box_limits*/

      ArrayOfGridPos p_gp((cloudbox_limits[1] - cloudbox_limits[0]) + 1);

      p2gridpos(p_gp,
                p_grid[Range(cloudbox_limits[0],
                             2,
                             (cloudbox_limits[1] - cloudbox_limits[0]))],
                p_grid[Range(cloudbox_limits[0],
                             (cloudbox_limits[1] - cloudbox_limits[0]) + 1)]);



      Matrix itw((cloudbox_limits[1] - cloudbox_limits[0]) + 1, 2);
      interpweights(itw, p_gp);

      Tensor6 scat_i_p(2, 1, 1, N_za, 1, N_i);
      scat_i_p(0, joker, joker, joker, joker, joker) =
          cloudbox_field(f_index, 0, joker, joker, joker, joker, joker);
      scat_i_p(1, joker, joker, joker, joker, joker) =
          cloudbox_field(f_index,
                       cloudbox_field.nvitrines() - 1,
                       joker,
                       joker,
                       joker,
                       joker,
                       joker);

      for (Index za_index = 0; za_index < N_za; ++za_index) {
        for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
          for (Index i = 0; i < N_i; ++i) {
            VectorView target_field = cloudbox_field(
                f_index, Range(joker), 0, 0, za_index, aa_index, i);

            ConstVectorView source_field =
                scat_i_p(Range(joker), 0, 0, za_index, aa_index, i);

            interp(target_field, itw, source_field, p_gp);
          }
        }
      }
    }
  } else if (atmosphere_dim == 3) {
    if (f_grid.nelem() == 0)
      throw runtime_error(
          "Error in cloudbox_fieldSetClearsky: For 3D "
          "all_frequencies option is not implemented \n");

    for (Index f_index = 0; f_index < cloudbox_field.nvitrines(); f_index++) {
      Index N_p = cloudbox_field.nvitrines();
      Index N_lat = cloudbox_field.nshelves();
      Index N_lon = cloudbox_field.nbooks();
      Index N_za = cloudbox_field.npages();
      Index N_aa = cloudbox_field.nrows();
      Index N_i = cloudbox_field.ncols();

      Tensor6 scat_i_p(2, N_lat, N_lon, N_za, N_aa, N_i);
      scat_i_p(0, joker, joker, joker, joker, joker) =
          cloudbox_field(f_index, 0, joker, joker, joker, joker, joker);
      scat_i_p(1, joker, joker, joker, joker, joker) =
          cloudbox_field(f_index, N_p - 1, joker, joker, joker, joker, joker);

      Tensor6 scat_i_lat(N_p, 2, N_lon, N_za, N_aa, N_i);
      scat_i_lat(joker, 0, joker, joker, joker, joker) =
          cloudbox_field(f_index, joker, 0, joker, joker, joker, joker);
      scat_i_lat(joker, 1, joker, joker, joker, joker) =
          cloudbox_field(f_index, joker, N_lat - 1, joker, joker, joker, joker);

      Tensor6 scat_i_lon(N_p, N_lat, 2, N_za, N_aa, N_i);
      scat_i_lon(joker, joker, 0, joker, joker, joker) =
          cloudbox_field(f_index, joker, joker, 0, joker, joker, joker);
      scat_i_lon(joker, joker, 1, joker, joker, joker) =
          cloudbox_field(f_index, joker, joker, N_lon - 1, joker, joker, joker);

      //1. interpolation - pressure grid, latitude grid and longitude grid

      ArrayOfGridPos p_gp((cloudbox_limits[1] - cloudbox_limits[0]) + 1);
      ArrayOfGridPos lat_gp((cloudbox_limits[3] - cloudbox_limits[2]) + 1);
      ArrayOfGridPos lon_gp((cloudbox_limits[5] - cloudbox_limits[4]) + 1);

      /*the old grid is having only two elements, corresponding to the
             cloudbox_limits and the new grid have elements corresponding to
             all grid points inside the cloudbox plus the cloud_box_limits*/

      p2gridpos(p_gp,
                p_grid[Range(cloudbox_limits[0],
                             2,
                             (cloudbox_limits[1] - cloudbox_limits[0]))],
                p_grid[Range(cloudbox_limits[0],
                             (cloudbox_limits[1] - cloudbox_limits[0]) + 1)]);
      gridpos(lat_gp,
              lat_grid[Range(cloudbox_limits[2],
                             2,
                             (cloudbox_limits[3] - cloudbox_limits[2]))],
              lat_grid[Range(cloudbox_limits[2],
                             (cloudbox_limits[3] - cloudbox_limits[2]) + 1)]);
      gridpos(lon_gp,
              lon_grid[Range(cloudbox_limits[4],
                             2,
                             (cloudbox_limits[5] - cloudbox_limits[4]))],
              lon_grid[Range(cloudbox_limits[4],
                             (cloudbox_limits[5] - cloudbox_limits[4]) + 1)]);

      //interpolation weights corresponding to pressure, latitude and
      //longitude grids.

      Matrix itw_p((cloudbox_limits[1] - cloudbox_limits[0]) + 1, 2);
      Matrix itw_lat((cloudbox_limits[3] - cloudbox_limits[2]) + 1, 2);
      Matrix itw_lon((cloudbox_limits[5] - cloudbox_limits[4]) + 1, 2);

      interpweights(itw_p, p_gp);
      interpweights(itw_lat, lat_gp);
      interpweights(itw_lon, lon_gp);

      // interpolation - pressure grid
      for (Index lat_index = 0;
           lat_index <= (cloudbox_limits[3] - cloudbox_limits[2]);
           ++lat_index) {
        for (Index lon_index = 0;
             lon_index <= (cloudbox_limits[5] - cloudbox_limits[4]);
             ++lon_index) {
          for (Index za_index = 0; za_index < N_za; ++za_index) {
            for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
              for (Index i = 0; i < N_i; ++i) {
                VectorView target_field = cloudbox_field(f_index,
                                                       Range(joker),
                                                       lat_index,
                                                       lon_index,
                                                       za_index,
                                                       aa_index,
                                                       i);

                ConstVectorView source_field = scat_i_p(
                    Range(joker), lat_index, lon_index, za_index, aa_index, i);

                interp(target_field, itw_p, source_field, p_gp);
              }
            }
          }
        }
      }
      //interpolation latitude
      for (Index p_index = 0;
           p_index <= (cloudbox_limits[1] - cloudbox_limits[0]);
           ++p_index) {
        for (Index lon_index = 0;
             lon_index <= (cloudbox_limits[5] - cloudbox_limits[4]);
             ++lon_index) {
          for (Index za_index = 0; za_index < N_za; ++za_index) {
            for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
              for (Index i = 0; i < N_i; ++i) {
                VectorView target_field = cloudbox_field(f_index,
                                                       p_index,
                                                       Range(joker),
                                                       lon_index,
                                                       za_index,
                                                       aa_index,
                                                       i);

                ConstVectorView source_field = scat_i_lat(
                    p_index, Range(joker), lon_index, za_index, aa_index, i);

                interp(target_field, itw_lat, source_field, lat_gp);
              }
            }
          }
        }
      }
      //interpolation -longitude
      for (Index p_index = 0;
           p_index <= (cloudbox_limits[1] - cloudbox_limits[0]);
           ++p_index) {
        for (Index lat_index = 0;
             lat_index <= (cloudbox_limits[3] - cloudbox_limits[2]);
             ++lat_index) {
          for (Index za_index = 0; za_index < N_za; ++za_index) {
            for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
              for (Index i = 0; i < N_i; ++i) {
                VectorView target_field = cloudbox_field(f_index,
                                                       p_index,
                                                       lat_index,
                                                       Range(joker),
                                                       za_index,
                                                       aa_index,
                                                       i);

                ConstVectorView source_field = scat_i_lon(
                    p_index, lat_index, Range(joker), za_index, aa_index, i);

                interp(target_field, itw_lon, source_field, lon_gp);
              }
            }
          }
        }
      }  //end of interpolation
    }    // end of frequency loop
  }      //ends atmosphere_dim = 3
}

void GetIncomingRadiation(Workspace& ws,
                          Tensor7& cloudbox_field,
                          const Agenda& iy_main_agenda,
                          const Index& atmosphere_dim,
                          const Vector& lat_grid,
                          const Vector& lon_grid,
                          const Tensor3& z_field,
                          const EnergyLevelMap& nlte_field,
                          const ArrayOfIndex& cloudbox_limits,
                          const Vector& f_grid,
                          const Vector& za_grid,
                          const Vector& aa_grid,
                          const Verbosity& verbosity) {
  CREATE_OUT0;

  // iy_unit hard.coded to "1" here
  const String iy_unit = "1";

  //Warning threshold if ratio of two radiances
  const Numeric maxratio = 100.;

  Index Nf = f_grid.nelem();
  Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  Index Nza = za_grid.nelem();
  Matrix iy;
  Ppath ppath;

  if (atmosphere_dim == 1) {
    //Define the variables for position and direction.
    Vector los(1), pos(1);

    //--- Get cloudbox_field at lower and upper boundary
    //    (boundary=0: lower, boundary=1: upper)
    for (Index boundary = 0; boundary <= 1; boundary++) {
      const Index boundary_index = boundary ? cloudbox_field.nvitrines() - 1 : 0;
      pos[0] = z_field(cloudbox_limits[boundary], 0, 0);

      // doing the first angle separately for allowing dy between 2 angles
      // in the loop
      los[0] = za_grid[0];
      get_iy(ws,
             iy,
             0,
             f_grid,
             nlte_field,
             pos,
             los,
             Vector(0),
             iy_unit,
             iy_main_agenda);
      cloudbox_field(joker, boundary_index, 0, 0, 0, 0, joker) = iy;

      for (Index za_index = 1; za_index < Nza; za_index++) {
        los[0] = za_grid[za_index];

        get_iy(ws,
               iy,
               0,
               f_grid,
               nlte_field,
               pos,
               los,
               Vector(0),
               iy_unit,
               iy_main_agenda);

        cloudbox_field(joker, boundary_index, 0, 0, za_index, 0, joker) = iy;

        for (Index fi = 0; fi < Nf; fi++) {
          if (cloudbox_field(fi, boundary_index, 0, 0, za_index - 1, 0, 0) /
                      cloudbox_field(fi, boundary_index, 0, 0, za_index, 0, 0) >
                  maxratio ||
              cloudbox_field(fi, boundary_index, 0, 0, za_index - 1, 0, 0) /
                      cloudbox_field(fi, boundary_index, 0, 0, za_index, 0, 0) <
                  1 / maxratio) {
            out0 << "ERROR: Radiance difference between "
                 << "interpolation points is too large (factor " << maxratio
                 << ")\n"
                 << "to safely interpolate. This might be due to "
                 << "za_grid being too coarse or the radiance field "
                 << "being a step-like function.\n"
                 << "Happens at boundary " << boundary_index
                 << " between zenith angles " << za_grid[za_index - 1]
                 << " and " << za_grid[za_index] << "deg\n"
                 << "for frequency #" << fi << ", where radiances are "
                 << cloudbox_field(fi, boundary_index, 0, 0, za_index - 1, 0, 0)
                 << " and "
                 << cloudbox_field(fi, boundary_index, 0, 0, za_index, 0, 0)
                 << " W/(sr m2 Hz).";
          }
        }
      }
    }
  }

  //--- atmosphere_dim = 3: --------------------------------------------------
  else {
    Index Naa = aa_grid.nelem();

    if (aa_grid[0] != 0. || aa_grid[Naa - 1] != 360.)
      throw runtime_error(
          "*aa_grid* must include 0 and 360 degrees as endpoints.");

    Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
    Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;

    // Convert aa_grid to "sensor coordinates"
    // (-180° < azimuth angle < 180°)
    //
    Vector sensor_aa_grid(Naa);
    for (Index i = 0; i < Naa; i++) sensor_aa_grid[i] = aa_grid[i] - 180;

    // Define the variables for position and direction.
    Vector los(2), pos(3);

    //--- Get cloudbox_field at lower and upper boundary
    //    (boundary=0: lower, boundary=1: upper)
    for (Index boundary = 0; boundary <= 1; boundary++) {
      const Index boundary_index = boundary ? cloudbox_field.nvitrines() - 1 : 0;
      for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++) {
        for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++) {
          pos[2] = lon_grid[lon_index + cloudbox_limits[4]];
          pos[1] = lat_grid[lat_index + cloudbox_limits[2]];
          pos[0] = z_field(cloudbox_limits[boundary],
                           lat_index + cloudbox_limits[2],
                           lon_index + cloudbox_limits[4]);

          for (Index za_index = 0; za_index < Nza; za_index++) {
            for (Index scat_aa_index = 0; scat_aa_index < Naa;
                 scat_aa_index++) {
              los[0] = za_grid[za_index];
              los[1] = sensor_aa_grid[scat_aa_index];

              // For end points of za_index (0 & 180deg), we
              // only need to perform calculations for one scat_aa
              // and set the others to same value
              if ((za_index != 0 && za_index != (Nza - 1)) ||
                  scat_aa_index == 0) {
                get_iy(ws,
                       iy,
                       0,
                       f_grid,
                       nlte_field,
                       pos,
                       los,
                       Vector(0),
                       iy_unit,
                       iy_main_agenda);
              }

              cloudbox_field(joker,
                           boundary_index,
                           lat_index,
                           lon_index,
                           za_index,
                           scat_aa_index,
                           joker) = iy;
            }
          }
        }
      }
    }

    //--- Get scat_i_lat (2nd and 3rd boundary)
    for (Index boundary = 0; boundary <= 1; boundary++) {
      const Index boundary_index = boundary ? cloudbox_field.nshelves() - 1 : 0;
      for (Index p_index = 0; p_index < Np_cloud; p_index++) {
        for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++) {
          pos[2] = lon_grid[lon_index + cloudbox_limits[4]];
          pos[1] = lat_grid[cloudbox_limits[boundary + 2]];
          pos[0] = z_field(p_index + cloudbox_limits[0],
                           cloudbox_limits[boundary + 2],
                           lon_index + cloudbox_limits[4]);

          for (Index za_index = 0; za_index < Nza; za_index++) {
            for (Index scat_aa_index = 0; scat_aa_index < Naa;
                 scat_aa_index++) {
              los[0] = za_grid[za_index];
              los[1] = sensor_aa_grid[scat_aa_index];

              // For end points of za_index, we need only to
              // perform calculations for first scat_aa
              if ((za_index != 0 && za_index != (Nza - 1)) ||
                  scat_aa_index == 0) {
                get_iy(ws,
                       iy,
                       0,
                       f_grid,
                       nlte_field,
                       pos,
                       los,
                       Vector(0),
                       iy_unit,
                       iy_main_agenda);
              }

              cloudbox_field(joker,
                           p_index,
                           boundary_index,
                           lon_index,
                           za_index,
                           scat_aa_index,
                           joker) = iy;
            }
          }
        }
      }
    }

    //--- Get scat_i_lon (1st and 2nd boundary):
    for (Index boundary = 0; boundary <= 1; boundary++) {
      const Index boundary_index = boundary ? cloudbox_field.nbooks() - 1 : 0;
      for (Index p_index = 0; p_index < Np_cloud; p_index++) {
        for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++) {
          pos[2] = lon_grid[cloudbox_limits[boundary + 4]];
          pos[1] = lat_grid[lat_index + cloudbox_limits[2]];
          pos[0] = z_field(p_index + cloudbox_limits[0],
                           lat_index + cloudbox_limits[2],
                           cloudbox_limits[boundary + 4]);

          for (Index za_index = 0; za_index < Nza; za_index++) {
            for (Index scat_aa_index = 0; scat_aa_index < Naa;
                 scat_aa_index++) {
              los[0] = za_grid[za_index];
              los[1] = sensor_aa_grid[scat_aa_index];

              // For end points of za_index, we need only to
              // perform calculations for first scat_aa
              if ((za_index != 0 && za_index != (Nza - 1)) ||
                  scat_aa_index == 0) {
                get_iy(ws,
                       iy,
                       0,
                       f_grid,
                       nlte_field,
                       pos,
                       los,
                       Vector(0),
                       iy_unit,
                       iy_main_agenda);
              }

              cloudbox_field(joker,
                           p_index,
                           lat_index,
                           boundary_index,
                           za_index,
                           scat_aa_index,
                           joker) = iy;
            }
          }
        }
      }
    }
  }  // End atmosphere_dim = 3.
}

void LimitInputGridsAndFieldsToCloudbox(Vector& p_grid_cldbx,
                                        Vector& lat_grid_cldbx,
                                        Vector& lon_grid_cldbx,
                                        Tensor3& t_field_cldbx,
                                        Tensor3& z_field_cldbx,
                                        Tensor4& vmr_field_cldbx,
                                        const ConstVectorView p_grid,
                                        const ConstVectorView lat_grid,
                                        const ConstVectorView lon_grid,
                                        const ConstTensor3View t_field,
                                        const ConstTensor3View z_field,
                                        const ConstTensor4View vmr_field,
                                        const ArrayOfIndex& cloudbox_limits,
                                        const Verbosity& verbosity) {
  CREATE_OUT3;

  //Limit input fields and spatial grids to cloudbox
  const Index Np_cldbx = cloudbox_limits[1] - cloudbox_limits[0] + 1;

  const Index Idx0_lat_cldbx =
      cloudbox_limits.nelem() > 2 ? cloudbox_limits[2] : 0;
  const Index Nlat_cldbx = cloudbox_limits.nelem() > 2
                               ? (cloudbox_limits[3] - cloudbox_limits[2] + 1)
                               : 1;
  const Index Idx0_lon_cldbx =
      cloudbox_limits.nelem() > 2 ? cloudbox_limits[4] : 0;
  const Index Nlon_cldbx = cloudbox_limits.nelem() > 2
                               ? (cloudbox_limits[5] - cloudbox_limits[4] + 1)
                               : 1;

  ostringstream os;
  os << "limit fields to cloudbox: index are set \n"
     << "Idx0_lat_cldbx = " << Idx0_lat_cldbx << " \n"
     << "Nlat_cldbx = " << Nlat_cldbx << "\n"
     << "Idx0_lon_cldbx = " << Idx0_lon_cldbx << " \n"
     << "Nlon_cldbx = " << Nlon_cldbx << "\n";
  out3 << os.str();

  p_grid_cldbx = p_grid[Range(cloudbox_limits[0], Np_cldbx)];
  if (cloudbox_limits.nelem() > 2) {
    lat_grid_cldbx = lat_grid[Range(Idx0_lat_cldbx, Nlat_cldbx)];
    lon_grid_cldbx = lon_grid[Range(Idx0_lon_cldbx, Nlon_cldbx)];
  } else {
    lat_grid_cldbx = lat_grid;
    lon_grid_cldbx = lat_grid;
  }


  t_field_cldbx = t_field(Range(cloudbox_limits[0], Np_cldbx),
                          Range(Idx0_lat_cldbx, Nlat_cldbx),
                          Range(Idx0_lon_cldbx, Nlon_cldbx));


  z_field_cldbx = z_field(Range(cloudbox_limits[0], Np_cldbx),
                          Range(Idx0_lat_cldbx, Nlat_cldbx),
                          Range(Idx0_lon_cldbx, Nlon_cldbx));


  if (vmr_field.nbooks()) {
    vmr_field_cldbx = vmr_field(joker,
                                Range(cloudbox_limits[0], Np_cldbx),
                                Range(Idx0_lat_cldbx, Nlat_cldbx),
                                Range(Idx0_lon_cldbx, Nlon_cldbx));

  }
}

void NewDoitMonoCalc(Workspace& ws,
    //Input and Output:
                     Tensor6& cloudbox_field_mono,
                     Tensor3& gas_extinction,
                     Tensor6& extinction_matrix,
                     Tensor5& absorption_vector,
                     Tensor7& scattering_matrix,
                     Index& convergence_flag,
                     Index& iteration_counter,
                     const ArrayOfIndex& cloudbox_limits,
                     const Agenda& propmat_clearsky_agenda,
                     const Agenda& surface_rtprop_agenda,
                     const Agenda& ppath_step_agenda,
                     const Index& atmosphere_dim,
                     const Index& stokes_dim,
                     const Tensor4& pnd_field,
                     const Tensor3& t_field,
                     const Tensor3& z_field,
                     const Tensor4& vmr_field,
                     const Matrix& z_surface,
                     const Vector& p_grid,
                     const Vector& lat_grid,
                     const Vector& lon_grid,
                     const Vector& za_grid,
                     const Vector& aa_grid,
                     const Vector& scat_za_grid,
                     const Vector& scat_aa_grid,
                     const Numeric& f_mono,
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     const Index& t_interp_order,
                     const String& iy_unit,
                     const Vector& refellipsoid,
                     const Vector& epsilon,
                     const Index& max_num_iterations,
                     const Numeric& tau_max,
                     const Index& accelerated,
                     const Index& ForwardCorrectionFlag,
                     const Numeric& ppath_lmax,
                     const Numeric& ppath_lraytrace,
                     const Verbosity& verbosity)

{

  CREATE_OUT0;
  ostringstream os, os1, os2;


  //-------Scattering related quantities----------------------------------------

  // set up main rt domain
  RTDomain MainRTDomain(p_grid,
                        lat_grid,
                        lon_grid,
                        za_grid,
                        aa_grid,
                        atmosphere_dim);

  //Set cloudbox_field of MainRTDomain
  MainRTDomain.set_CloudboxField(cloudbox_field_mono);

  RTDomainScatteringProperties MainDomainScatteringProperties;

  // calculate the bulk particle extinction, absorption and scattering
  // properties
  CalcDomainExtAbsScatProperties(MainDomainScatteringProperties,
                                 MainRTDomain,
                                 t_field,
                                 scat_data,
                                 pnd_field,
                                 stokes_dim,
                                 atmosphere_dim,
                                 za_grid,
                                 aa_grid,
                                 scat_za_grid,
                                 scat_aa_grid,
                                 t_interp_order,
                                 ForwardCorrectionFlag,
                                 verbosity);

  //-------Surface related quantities-------------------------------------------

  Matrix surface_skin_t;
  Tensor4 surface_los;
  Tensor6 surface_reflection_matrix;
  Tensor5 surface_emission;

  //calculate surf_optpropCalc_doit
  if (atmosphere_dim == 1) {
    if (z_surface(0, 0) >= z_field(0, 0, 0)) {
      CalcSurfaceProperties(ws,
                            surface_skin_t,
                            surface_los,
                            surface_reflection_matrix,
                            surface_emission,
                            surface_rtprop_agenda,
                            f_mono,
                            za_grid,
                            aa_grid,
                            lat_grid,
                            lon_grid,
                            atmosphere_dim,
                            stokes_dim,
                            z_surface);
    }
  }


  CalcGasExtinctionField(ws,
                         gas_extinction,
                         p_grid,
                         lat_grid,
                         lon_grid,
                         t_field,
                         vmr_field,
                         propmat_clearsky_agenda,
                         f_mono);

  //TODO: Include refinement here. This includes the check where to refine and
  // the prepration of the subdomains.

  Tensor4 DlogK;
  Tensor4 tau;
  ArrayOfIndex refine;

  CheckForRefinement(
      refine,
      DlogK,
      tau,
      MainRTDomain.get_ExtinctionMatrix(),  //(Np,Nlat,Nlon,Nza,Naa,nst,nst)
      gas_extinction,
      p_grid,
      lat_grid,
      lon_grid,
      scat_za_grid,
      z_field,
      refellipsoid,
      atmosphere_dim,
      1.,
      tau_max);

  //TODO: Remove this! For the adaptive this is not needed anymore.
//  Tensor3 p_path_maxlength;
//  if (tau_max > 0) {
//    //To estimate the maximum ppath length, we need gas and particle extinction,
//    //so we calculate the gas extinction on the cloud box pressure grid first.
//    //The actual gas extinction which will be used for the RT calculation will
//    //calculated additionally.
//
//    //calculate local ppath_lmax
//    CalcPropagationPathMaxLength(
//        p_path_maxlength,
//        MainRTDomain.get_ExtinctionMatrix(),  //(Np,Nlat,Nlon,ndir,nst,nst)
//        gas_extinction,
//        p_grid,
//        lat_grid,
//        lon_grid,
//        za_grid,
//        ppath_lmax,
//        ppath_lraytrace,
//        tau_max);
//  }




  ArrayOfMatrix VmrArray;

  //Estimate the ppath for main rt domain
  if (atmosphere_dim ==1) {
      EstimatePPathElements1D(ws,
                              MainRTDomain.get_PressureArray(),
                              MainRTDomain.get_TemperatureArray(),
                              VmrArray,
                              MainRTDomain.get_InterpWeightsArray(),
                              MainRTDomain.get_InterpWeightsAngleArray(),
                              MainRTDomain.get_GposPArray(),
                              MainRTDomain.get_GposZenithArray(),
                              MainRTDomain.get_LStepArray(),
                              MainRTDomain.get_MaxLimbIndex(),
                              cloudbox_limits,
                              za_grid,
                              ppath_step_agenda,
                              ppath_lmax,
                              ppath_lraytrace,
//                              p_path_maxlength,
                              p_grid,
                              z_field,
                              t_field,
                              vmr_field,
                              refellipsoid,
                              Vector(1, f_mono),
                              verbosity);

      //TODO: Add here the ppath estimation for the subdomains
      //Estimate the ppath for the subdomains and save them ArrayOf...
      //Each array element is one subdomain.
      //Ask Olli for the best strategy


  } else {
    throw runtime_error("3d is not implemented yet! Sorry");
//    EstimatePPathElements3D(ws,
//                            PressureArray,
//                            TemperatureArray,
//                            VmrArray,
//                            InterpWeightsArray,
//                            InterpWeightsZenithArray,
//                            GposArray,
//                            GposZenithArray,
//                            LstepArray,
//                            MaxLimbIndex,
//                            cloudbox_limits,
//                            za_grid,
//                            ppath_step_agenda,
//                            ppath_lmax,
//                            ppath_lraytrace,
//                            p_path_maxlength,
//                            p_grid,
//                            z_field,
//                            t_field,
//                            vmr_field,
//                            refellipsoid,
//                            Vector(1, f_mono),
//                            verbosity);
  }



  //calculate gas extinction on p-path points
  MainRTDomain.get_GasExtinctionArray().resize(MainRTDomain.get_PressureArray().nelem());

  for (Index i = 0; i < MainRTDomain.get_PressureArray().nelem(); i++) {
    Vector gas_extinction_temp;
    CalcGasExtinction(ws,
                      MainRTDomain.get_GasExtinctionArray()[i],
                      MainRTDomain.get_PressureArray()[i],
                      MainRTDomain.get_TemperatureArray()[i],
                      VmrArray[i],
                      propmat_clearsky_agenda,
                      f_mono);

  }

  //TODO: Add here the Gasextinction of the subdomains

  os << "gas absorption calculated \n";
  out0 << os.str();
  os.clear();

  //TODO:add subdomain variables. Ask Olli if there is a better way to include
  // the subdomain.



  //run new doit
  RunNewDoit(
             MainRTDomain,
             convergence_flag,
             iteration_counter,
             MainDomainScatteringProperties,
             surface_reflection_matrix,
             surface_emission,
             cloudbox_limits,
             atmosphere_dim,
             f_mono,
             iy_unit,
             epsilon,
             max_num_iterations,
             accelerated,
             verbosity);

  cloudbox_field_mono=MainRTDomain.get_CloudboxField();
  extinction_matrix=MainRTDomain.get_ExtinctionMatrix()(joker,joker,joker,joker,0,joker,joker);
  absorption_vector=MainRTDomain.get_AbsorptionVector()(joker,joker,joker,joker,0,joker);
  scattering_matrix=MainDomainScatteringProperties.get_ScatteringMatrix();

}

void CalcGasExtinction(Workspace& ws,
                       Vector& gas_extinct,
                       const ConstVectorView& p_grid,
                       const ConstVectorView& t_vector,
                       const ConstMatrixView& vmr_matrix,
                       const Agenda& propmat_clearsky_agenda,
                       const ConstVectorView& f_mono) {
  const Index Np = p_grid.nelem();

  // Initialization
  gas_extinct.resize(Np);
  gas_extinct = 0.;

  // Local variables to be used in agendas

  PropagationMatrix propmat_clearsky_local;

  const EnergyLevelMap rtp_temperature_nlte_local_dummy;

  // Calculate layer averaged gaseous extinction
  for (Index ip = 0; ip < Np; ip++) {
    const Vector rtp_mag_dummy(3, 0);
    const Vector ppath_los_dummy;

    StokesVector nlte_dummy;
    ArrayOfPropagationMatrix partial_dummy;
    ArrayOfStokesVector partial_source_dummy;
    propmat_clearsky_agendaExecute(ws,
                                   propmat_clearsky_local,
                                   nlte_dummy,
                                   partial_dummy,
                                   partial_source_dummy,
                                   ArrayOfRetrievalQuantity(0),
                                   (Vector)f_mono,  // monochromatic calculation
                                   rtp_mag_dummy,
                                   ppath_los_dummy,
                                   p_grid[ip],
                                   t_vector[ip],
                                   rtp_temperature_nlte_local_dummy,
                                   vmr_matrix(joker, ip),
                                   propmat_clearsky_agenda);

    //Assuming non-polarized light and only one frequency
    //TODO: Check if polarization is needed for absorption

    gas_extinct[ip] += propmat_clearsky_local.Kjj(0,0)[0];


  }
}

void CalcGasExtinctionField(Workspace& ws,
                            Tensor3& gas_extinct,
                            const ConstVectorView& p_grid,
                            const ConstVectorView& lat_grid,
                            const ConstVectorView& lon_grid,
                            const ConstTensor3View& t_field,
                            const ConstTensor4View& vmr_field,
                            const Agenda& propmat_clearsky_agenda,
                            const ConstVectorView& f_mono) {

  const Index Np = p_grid.nelem();
  const Index Nlat = lat_grid.nelem() > 0 ? lat_grid.nelem() : 1;
  const Index Nlon = lon_grid.nelem() > 0 ? lon_grid.nelem() : 1;

  // Initialization
  gas_extinct.resize(Np, Nlat, Nlon);
  gas_extinct = 0.;

  for (Index ilat = 0; ilat < Nlat; ilat++) {
    for (Index ilon = 0; ilon < Nlon; ilon++) {
      Vector gas_extinc_temp;
      CalcGasExtinction(ws,
                        gas_extinc_temp,
                        p_grid,
                        t_field(joker, ilat, ilon),
                        vmr_field(joker, joker, ilat, ilon),
                        propmat_clearsky_agenda,
                        f_mono);

      gas_extinct(joker, ilat, ilon) = gas_extinc_temp;
    }
  }
}

void CalcDomainExtAbsScatProperties(
    RTDomainScatteringProperties& DomainScatteringProperties,
    RTDomain& Domain,
    const Tensor3& t_field,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const ConstTensor4View& pnd_field,
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Vector& scat_za_grid,
    const Vector& scat_aa_grid,
    const Index& t_interp_order,
    const Index& ForwardCorrectionFlag,
    const Verbosity& verbosity){


  // some resizing
  Domain.get_ExtinctionMatrix().resize(t_field.npages(),
                                         t_field.nrows(),
                                         t_field.ncols(),
                                         za_grid.nelem(),
                                         aa_grid.nelem(),
                                         stokes_dim,
                                         stokes_dim);

  Domain.get_AbsorptionVector().resize(t_field.npages(),
                                         t_field.nrows(),
                                         t_field.ncols(),
                                         za_grid.nelem(),
                                       aa_grid.nelem(),
                                         stokes_dim);

  DomainScatteringProperties.set_ScatZaGrid(scat_za_grid);
  DomainScatteringProperties.set_ScatAaGrid(scat_aa_grid);


  //calculate particle optical properties
  //It assumens at least azimuthally randomly oriented particles
  CalcParticleOpticalProperties(Domain.get_ExtinctionMatrix(),
                                Domain.get_AbsorptionVector(),
                                scat_data,
                                za_grid,
                                pnd_field,
                                t_field,
                                stokes_dim);

  //calculate sca_optpropCalc_doit
  CalcScatteringProperties(DomainScatteringProperties.get_ScatteringMatrix(),
                           DomainScatteringProperties.get_IdirIdx0(),
                           DomainScatteringProperties.get_IdirIdx1(),
                           DomainScatteringProperties.get_PdirIdx0(),
                           DomainScatteringProperties.get_PdirIdx1(),
                           DomainScatteringProperties.get_GpZaI(),
                           DomainScatteringProperties.get_GpAaI(),
                           DomainScatteringProperties.get_Itw(),
                           t_field,
                           scat_data,
                           pnd_field,
                           stokes_dim,
                           atmosphere_dim,
                           za_grid,
                           aa_grid,
                           DomainScatteringProperties.get_ScatZaGrid(),
                           DomainScatteringProperties.get_ScatAaGrid(),
                           t_interp_order,
                           verbosity);

  //forward scattering correction
  if (ForwardCorrectionFlag) {
    ForwardScatteringCorrection(
        DomainScatteringProperties.get_ScatteringMatrix(),
        Domain.get_ExtinctionMatrix(),  //(Np,Nlat,Nlon,Nza,Naa,nst,nst)
        Domain.get_AbsorptionVector(),  //(Np,Nlat,Nlon,Nza,Naa,nst)
        DomainScatteringProperties.get_IdirIdx0(),
        DomainScatteringProperties.get_IdirIdx1(),
        DomainScatteringProperties.get_PdirIdx0(),
        DomainScatteringProperties.get_PdirIdx1(),
        atmosphere_dim,
        za_grid,
        aa_grid,
        DomainScatteringProperties.get_ScatZaGrid(),
        DomainScatteringProperties.get_ScatAaGrid());
  }





};



void CalcParticleOpticalProperties(Tensor7& extinction_matrix,//(Np,Nlat,Nlon,Nza,Naa,nst,nst)
                    Tensor6& absorption_vector,//(Np,Nlat,Nlon,Nza,Naa,nst)
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    const Vector& za_grid,
                    const ConstTensor4View& pnd_field,
                    const ConstTensor3View& t_field,
                    const Index& stokes_dim) {
  // Initialization
  extinction_matrix = 0.;
  absorption_vector = 0.;

  const Index Np = pnd_field.npages();
  const Index Nlat = pnd_field.nrows();
  const Index Nlon = pnd_field.ncols();


  assert(absorption_vector.nvitrines() == Np);
  assert(extinction_matrix.nlibraries() == Np);

  // preparing input data
  Matrix dir_array(za_grid.nelem(), 2, 0.);
  dir_array(joker, 0) = za_grid;

  // making output containers
  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor5 ext_mat_ssbulk;
  ArrayOfTensor4 abs_vec_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor5 ext_mat_bulk_ii;
  Tensor4 abs_vec_bulk_ii;
  Index ptype_bulk_ii;


  for (Index ilat = 0; ilat < Nlat; ilat++) {
    for (Index ilon = 0; ilon < Nlon; ilon++) {

      //Calculate the optical properties for each scattering element
      opt_prop_NScatElems(ext_mat_Nse,
                          abs_vec_Nse,
                          ptypes_Nse,
                          t_ok,
                          scat_data,
                          stokes_dim,
                          t_field(joker,ilat,ilon),
                          dir_array,
                          0);

      //Calculate the optical properties for each scattering species
      opt_prop_ScatSpecBulk(ext_mat_ssbulk,
                            abs_vec_ssbulk,
                            ptype_ssbulk,
                            ext_mat_Nse,
                            abs_vec_Nse,
                            ptypes_Nse,
                            pnd_field(joker, joker, ilat, ilon),
                            t_ok);

      //Calculate the bulk optical properties
      opt_prop_Bulk(ext_mat_bulk_ii, //(nf,nT,ndir,nst,nst)
                    abs_vec_bulk_ii, //(nf,nT,ndir,nst)
                    ptype_bulk_ii,
                    ext_mat_ssbulk,
                    abs_vec_ssbulk,
                    ptype_ssbulk);

      extinction_matrix(joker, ilat, ilon, joker, 0, joker, joker) =
          ext_mat_bulk_ii(0, joker, joker, joker, joker);

      absorption_vector(joker, ilat, ilon, joker, 0, joker) =
          abs_vec_bulk_ii(0, joker, joker, joker);
    }
  }
}

void CalcScatteringProperties(  //Output
    Tensor7& scattering_matrix,
    ArrayOfIndex& idir_idx0,
    ArrayOfIndex& idir_idx1,
    ArrayOfIndex& pdir_idx0,
    ArrayOfIndex& pdir_idx1,
    ArrayOfGridPos& gp_za_i,
    ArrayOfGridPos& gp_aa_i,
    Tensor3& itw,
    const Tensor3& t_field,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const ConstTensor4View& pnd_field,
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Vector& scat_za_grid,
    const Vector& scat_aa_grid,
    const Index& t_interp_order,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  const Index Np = pnd_field.npages();
  const Index Nlat = pnd_field.nrows();
  const Index Nlon = pnd_field.ncols();

  // preparing input data
  Matrix idir_array;

  //Flatten incoming directions
  FlattenMeshGrid(idir_array, scat_za_grid, scat_aa_grid);

  //Flatten outgoing directions
  Matrix pdir_array;
  if (atmosphere_dim == 1) {
    pdir_array.resize(za_grid.nelem(), 2);
    pdir_array(joker, 0) = za_grid;
    pdir_array(joker, 1) = 0.;

  } else {
    FlattenMeshGrid(pdir_array, za_grid, aa_grid);
  }

  // making output containers
  ArrayOfArrayOfTensor6 sca_mat_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor6 sca_mat_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor6 scat_mat_bulk_ii;
  Index ptype_bulk_ii;

  //We also prepare the gridpos-array, the interpolation weights and the index-
  //arrays for the flattened meshgrid of the angular grids. They all are
  //constant, so we calculate them only once. They are all needed
  //for the calculation of the scattering field within the actual solver.

  //grid pos array for zenith angle interpolation.
  gp_za_i.resize(scat_za_grid.nelem());
  gridpos(gp_za_i, za_grid, scat_za_grid);

  if (atmosphere_dim == 1) {
    scattering_matrix.resize(Np,
                             Nlat,
                             Nlon,
                             pdir_array.nrows(),
                             scat_za_grid.nelem(),
                             stokes_dim,
                             stokes_dim);

    //index arrays for flattened meshgrid of the angular grids, also needed
    //for the calculation of the scattering field later on.
    FlattenMeshGridIndex(idir_idx0, idir_idx1, scat_za_grid, Vector(1, 0.));
    FlattenMeshGridIndex(pdir_idx0, pdir_idx1, za_grid, Vector(1, 0.));

    //Interpolation weights, for 1D we need only interpolation weight for the
    //zenith angle interpolation.
    Matrix itw_za_i(scat_za_grid.nelem(), 2);
    itw.resize(scat_za_grid.nelem(), 1, 2);
    itw = 0;

    interpweights(itw_za_i, gp_za_i);
    itw(joker, 0, joker) = itw_za_i;

  } else {
    scattering_matrix.resize(Np,
                             Nlat,
                             Nlon,
                             pdir_array.nrows(),
                             idir_array.nrows(),
                             stokes_dim,
                             stokes_dim);

    FlattenMeshGridIndex(idir_idx0, idir_idx1, scat_za_grid, scat_aa_grid);
    FlattenMeshGridIndex(pdir_idx0, pdir_idx1, za_grid, aa_grid);

    //grid pos array for azimuth angle interpolation
    gp_aa_i.resize(scat_aa_grid.nelem());
    gridpos(gp_aa_i, aa_grid, scat_aa_grid);

    //Interpolation weights
    itw.resize(scat_za_grid.nelem(), scat_aa_grid.nelem(), 4);
    interpweights(itw, gp_za_i, gp_aa_i);
  }
  scattering_matrix = 0.;

  ostringstream os;
  os << "calculating ensemble scattering matrix \n";
  out2 << os.str();
  os.clear();

  for (Index ilat = 0; ilat < Nlat; ilat++) {
    for (Index ilon = 0; ilon < Nlon; ilon++) {
      pha_mat_NScatElems(  //Output
          sca_mat_Nse,     // [nss][nse](nf,nT,npdir,nidir,nst,nst)
          ptypes_Nse,
          t_ok,
          //Input
          scat_data,
          stokes_dim,
          t_field(joker, ilat, ilon),
          pdir_array,
          idir_array,
          0,
          t_interp_order);

      pha_mat_ScatSpecBulk(  //Output
          sca_mat_ssbulk,    // [nss](nf,nT,npdir,nidir,nst,nst)
          ptype_ssbulk,
          //Input
          sca_mat_Nse,  // [nss][nse](nf,nT,npdir,nidir,nst,nst)
          ptypes_Nse,
          pnd_field(joker, joker, ilat, ilon),
          t_ok);

      pha_mat_Bulk(          //Output
          scat_mat_bulk_ii,  // (nf,nT,npdir,nidir,nst,nst)
          ptype_bulk_ii,
          //Input
          sca_mat_ssbulk,  // [nss](nf,nT,npdir,nidir,nst,nst)
          ptype_ssbulk);

      //For 1D atmosphere we can do already the azimuth integration here
      //as for 1D atmosphere, there is no azimuth dependency of the incoming
      //radiation
      if (atmosphere_dim == 1) {
        os << "As we have a 1D atmosphere, we can already integrate over \n"
           << "incoming azimuth. \n";
        out2 << os.str();
        os.clear();

        Tensor5 scat_mat_bulk_ii_temp(1,
                                      Np,
                                      za_grid.nelem(),
                                      stokes_dim,
                                      stokes_dim);
        Tensor6 scat_mat_bulk_ii_int(1,
                                     Np,
                                     za_grid.nelem(),
                                     scat_za_grid.nelem(),
                                     stokes_dim,
                                     stokes_dim,0.);

        Index idx_i = 0;
        Index idx_i1 = 0;
        Numeric delta_phi = 360. * DEG2RAD / Numeric((scat_aa_grid.nelem() - 1));

        //integrate over incoming azimuth
        for (Index i_za = 0; i_za < scat_za_grid.nelem(); i_za++) {
          for (Index i_aa = 0; i_aa < (scat_aa_grid.nelem() - 1); i_aa++) {
            //flattened indices
            idx_i = i_aa * scat_za_grid.nelem() + i_za;
            idx_i1 = (i_aa + 1) * scat_za_grid.nelem() + i_za;

            scat_mat_bulk_ii_temp =
                scat_mat_bulk_ii(joker, joker, joker, idx_i, joker, joker);

            scat_mat_bulk_ii_temp +=
                scat_mat_bulk_ii(joker, joker, joker, idx_i1, joker, joker);

            scat_mat_bulk_ii_int(joker, joker, joker, i_za, joker, joker) +=
                scat_mat_bulk_ii_temp;
          }
        }
        scat_mat_bulk_ii_int *= delta_phi / 2.;

        scattering_matrix(joker, ilat, ilon, joker, joker, joker, joker) =
            scat_mat_bulk_ii_int(0, joker, joker, joker, joker, joker);

      } else {
        scattering_matrix(joker, ilat, ilon, joker, joker, joker, joker) =
            scat_mat_bulk_ii(0, joker, joker, joker, joker, joker);
      }
    }
  }
}

void ForwardScatteringCorrection(  //Output
    Tensor7& scattering_matrix,    //(Np,Nlat,Nlon,npdir,nidir,nst,nst)
    Tensor7& extinction_matrix,    //(Np,Nlat,Nlon,Nza,Naa,nst,nst)
    Tensor6& absorption_vector,    //(Np,Nlat,Nlon,Nza,Naa,nst)
    const ArrayOfIndex& idir_idx0,
    const ArrayOfIndex& idir_idx1,
    const ArrayOfIndex& pdir_idx0,
    const ArrayOfIndex& pdir_idx1,
    const Index& atmosphere_dim,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Vector& scat_za_grid,
    const Vector& scat_aa_grid) {

  Index Np = scattering_matrix.nlibraries();
  Index Nst = scattering_matrix.ncols();

  if (atmosphere_dim == 1) {


    Vector Zij(1);
    Vector EvalPoints_Za(1);
    Vector scat_za_temp(scat_za_grid.nelem()-1);
    Vector Z_temp(scat_za_grid.nelem()-1);
    Numeric Ki1;

    Index counter;

    for (Index i_p = 0; i_p < Np; i_p++) {
      for (Index i_pz = 0; i_pz <= za_grid.nelem() - 1; i_pz++) {


        if (extinction_matrix(i_p, 0, 0, i_pz,0, 0, 0) > 0) {
          // Find nearest incidence angle to forward propagation direction
          Index i_sz = 0;
          while (i_sz < idir_idx0.nelem() - 1 &&
                 abs(scat_za_grid[i_sz] - za_grid[i_pz]) >
                     abs(scat_za_grid[i_sz + 1] - za_grid[i_pz])) {
            i_sz++;
          }

          EvalPoints_Za[0]=scat_za_grid[i_sz]; //nearest to forward direction

          //Take the scattering grid without the nearest one to the forward direction
          counter=-1;
          for (Index i_tz = 0; i_tz <= scat_za_grid.nelem()-1; i_tz++){
            if (i_tz!=i_sz){
              counter++;
              scat_za_temp[counter]=scat_za_grid[i_tz];
            }
          }

          //Prepare interpolation, we do a fixed Lagrangian interpolation of order 3
          const FixedLagrangeInterpolation<3> lag(0,EvalPoints_Za[0],scat_za_temp);
          const auto iwlag = interpweights(lag);

          for (Index i_sti = 0; i_sti < Nst; i_sti++) {
            for (Index i_stj = 0; i_stj < Nst; i_stj++) {

              //Take the scattering data without the nearest one to the forward direction
              counter=-1;
              for (Index i_tz = 0; i_tz <= scat_za_grid.nelem()-1; i_tz++){
                if (i_tz!=i_sz){
                  counter++;
                  Z_temp[counter]=scattering_matrix(i_p, 0, 0, i_pz, i_tz, i_sti, i_stj);
                }
              }

              //Remove forward scattering peak and replace it with interpolated
              // value. Interpolate the direction within the scattering grid nearest to
              // the forward direction.

              Zij=interp(Z_temp, iwlag, lag);

              //Interpolated value on scattering matrix
              scattering_matrix(i_p, 0, 0, i_pz, i_sz, i_sti, i_stj)=Zij[0];

            }
          }

          // The four extinction matrix elements based on energy conservation
          // The other elements are not changed.

          for (Index i_sti = 0; i_sti < Nst; i_sti++) {
            Ki1 = AngIntegrate_trapezoid(
                      scattering_matrix(i_p, 0, 0, i_pz, joker, i_sti, 0),
                      scat_za_grid) /
                      2 / PI +
                  absorption_vector(i_p, 0, 0, i_pz,0 , i_sti);


            extinction_matrix(i_p, 0, 0, i_pz, 0, i_sti, 0)=Ki1;
            if (i_sti>0) {
              extinction_matrix(i_p, 0, 0, i_pz, 0, 0, i_sti) = Ki1;

              extinction_matrix(i_p, 0, 0, i_pz, 0, i_sti, i_sti) = extinction_matrix(i_p, 0, 0, i_pz, 0, 0, 0);
            }
          }
        }
      }
    }
  } else {
    throw runtime_error("3d is not implemented yet! Sorry");
  }
}

void CalcSurfaceProperties(Workspace& ws,
                      //Output
                      Matrix& surface_skin_t,
                      Tensor4& surface_los,
                      Tensor6& surface_reflection_matrix,
                      Tensor5& surface_emission,

                      //Input
                      const Agenda& surface_rtprop_agenda,
                      const ConstVectorView& f_grid,
                      const ConstVectorView& za_grid,
                      const ConstVectorView& aa_grid,
                      const Vector& lat_grid,
                      const Vector& lon_grid,
                      const Index& atmosphere_dim,
                      const Index& stokes_dim,
                      const Matrix& z_surface) {

  chk_not_empty("surface_rtprop_agenda", surface_rtprop_agenda);
  const Index Nza = za_grid.nelem();
  const Index Nlat = lat_grid.nelem() > 0 ? lat_grid.nelem() : 1;
  const Index Nlon = lon_grid.nelem() > 0 ? lon_grid.nelem() : 1;

  Index Naa = 1;
  Vector rte_pos(1, 0);
  Vector rte_los(1, 0);

  // rte_pos and rte_los are of different size than for 1d atm.
  if (atmosphere_dim == 3) {
    rte_pos.resize(3);
    rte_pos = 0.;
    rte_pos.resize(2);
    rte_pos = 0.;

    Naa = aa_grid.nelem();
  }

  //Allocate
  Numeric surface_skin_t_ii;
  Matrix surface_los_ii;
  Tensor4 surface_rmatrix_ii;
  Matrix surface_emission_ii;
  Index los_idx;

  //prepare output container
  surface_skin_t.resize(Nlat, Nlon);
  surface_los.resize(Nlat, Nlon, Naa * Nza, 2);
  surface_reflection_matrix.resize(Nlat, Nlon, Naa * Nza, 1, stokes_dim, stokes_dim);
  surface_emission.resize(Nlat, Nlon,Naa * Nza, 1, stokes_dim);
  surface_skin_t = 0.;
  surface_los = NAN;
  surface_reflection_matrix = NAN;
  surface_emission = 0;

  // loop over geographic coordinates
  for (Index ilat = 0; ilat < Nlat; ilat++) {
    for (Index ilon = 0; ilon < Nlon; ilon++) {

      rte_pos[0] = z_surface(ilat, ilon);
      // for a 3d atmosphere rte_pos has 3 components
      if (atmosphere_dim == 3) {
        rte_pos[1] = lat_grid[ilat];
        rte_pos[2] = lon_grid[ilon];
      }

      los_idx = 0;
      for (Index iaa = 0; iaa < Naa; iaa++) {
        for (Index iza = 0; iza < Nza; iza++) {
          if (za_grid[iza] > 90) {
            rte_los[0] = za_grid[iza];


            if (atmosphere_dim == 3) {
              rte_los[1] = aa_grid[iaa];
            }

            // calculate the local surface properties
            surface_rtprop_agendaExecute(ws,
                                         surface_skin_t_ii,
                                         surface_emission_ii,
                                         surface_los_ii,
                                         surface_rmatrix_ii,
                                         f_grid,
                                         rte_pos,
                                         rte_los,
                                         surface_rtprop_agenda);

            // surface_skin_t is independent of rte_los, so we need it only for
            // first rte_los
            if (iaa == 0 && iza == 0) {
              surface_skin_t(ilat, ilon) = surface_skin_t_ii;
            }

            //Store it
            surface_reflection_matrix(ilat, ilon, los_idx, 0, joker, joker) =
                surface_rmatrix_ii(0, 0, joker, joker);
            surface_emission(ilat, ilon, los_idx, 0, joker) =
                surface_emission_ii(0, joker);

            surface_los(ilat, ilon, los_idx, 0) = surface_los_ii(0, 0);
            //for 3d atmosphere sensor_los has two angles.
            if (atmosphere_dim == 3) {
              surface_los(ilat, ilon, los_idx, 1) = surface_los_ii(0, 1);
            }
          }
          los_idx++;
        }
      }
    }
  }
}

void CalcPropagationPathMaxLength(
    Tensor3& p_path_maxlength,
    const Tensor7View& extinction_matrix,  //(Np,Nlat,Nlon,Nza,Naa,nst,nst)
    const ConstTensor3View& gas_extinct,
    const ConstVectorView& p_grid,
    const ConstVectorView& lat_grid,
    const ConstVectorView& lon_grid,
    const Vector& scat_za_grid,
    const Numeric& ppath_lmax,
    const Numeric& ppath_lraytrace,
    const Numeric& tau_max) {


  const Index Np = p_grid.nelem();
  const Index Nlat = lat_grid.nelem() > 0 ? lat_grid.nelem() : 1;
  const Index Nlon = lon_grid.nelem() > 0 ? lon_grid.nelem() : 1;
  const Index Ndir = extinction_matrix.npages();
  Numeric ext_mat_elem;


  //prepare output container
  p_path_maxlength.resize(Np, Nlat, Nlon);
  p_path_maxlength = min(ppath_lmax,ppath_lraytrace);

  for (Index ip = 0; ip < Np; ip++) {
    for (Index ilat = 0; ilat < Nlat; ilat++) {
      for (Index ilon = 0; ilon < Nlon; ilon++) {
        ext_mat_elem = 0.;

        // calculate average over zenith direction,
        // for azimuthal random orientation, this is not fully correct,
        // but since this is more an rough estimation we assume, that the
        // optical thickness is only a function of the position.
        for (Index idir = 0; idir < Ndir; idir++) {
          ext_mat_elem += (extinction_matrix(ip, ilat, ilon, idir, 0, 0, 0) *
                           sin(scat_za_grid[idir] * DEG2RAD));
        }
        ext_mat_elem /= Numeric(Ndir);

        // maximum propagation path length for each grid point.
        const Numeric maxlength=tau_max / (ext_mat_elem + gas_extinct(ip, ilat, ilon));

        if (maxlength < p_path_maxlength(ip, ilat, ilon)){
          p_path_maxlength(ip, ilat, ilon) = maxlength;
        }
      }
    }
  }
}

void CheckForRefinement(
    ArrayOfIndex& refine,
    Tensor4& DlogK,
    Tensor4& tau,
    const Tensor7View& extinction_matrix,  //(Np,Nlat,Nlon,Nza,Naa,nst,nst)
    const ConstTensor3View& gas_extinct,
    const ConstVectorView& p_grid,
    const ConstVectorView& lat_grid,
    const ConstVectorView& lon_grid,
    const Vector& scat_za_grid,
    const Tensor3& z_field,
    const Vector& refellipsoid,
    const Index& atmosphere_dim,
    const Numeric& refine_crit,
    const Numeric& tau_crit) {

  const Index Np = p_grid.nelem();
  const Index Nlat = lat_grid.nelem() > 0 ? lat_grid.nelem() : 1;
  const Index Nlon = lon_grid.nelem() > 0 ? lon_grid.nelem() : 1;
  const Index Nza = extinction_matrix.nbooks();
  const Index Naa = extinction_matrix.npages();
  Numeric ext_mat_elem;

  //allocate
  Tensor3 ExtMat(Np, Nlat, Nlon,0);

  // get the average total extinction at each position
  for (Index ip = 0; ip < Np; ip++) {
    for (Index ilat = 0; ilat < Nlat; ilat++) {
      for (Index ilon = 0; ilon < Nlon; ilon++) {
        ext_mat_elem = 0.;

        // calculate average over all direction,
        // for azimuthal random and arbitrary orientation, this is not fully correct,
        // but this is more an rough estimation.
        for (Index iza = 0; iza < Nza; iza++) {
          for (Index iaa = 0; iaa < Naa; iaa++) {
            ext_mat_elem += (extinction_matrix(ip, ilat, ilon, iza, iaa, 0, 0) *
                             sin(scat_za_grid[iza] * DEG2RAD));
          }
        }

        ext_mat_elem /= Numeric(Nza * Naa);

        // Total extinction
        ext_mat_elem += gas_extinct(ip, ilat, ilon);

        ExtMat(ip, ilat, ilon) = ext_mat_elem;
      }
    }
  }

  //allocate
  DlogK.resize(Np - 1, Nlat, Nlon, atmosphere_dim);
  tau.resize(Np - 1, Nlat, Nlon, atmosphere_dim);

  // estimate the change of the the extinction in terms of log difference between two
  // adjacent cell(layers)
  if (atmosphere_dim == 1) {
    refine.resize(Np - 1);

    for (Index ip = 0; ip < Np - 1; ip++) {
      // calculate log10 extinction difference
      if (ExtMat(ip + 1, 0, 0) > 0 && ExtMat(ip, 0, 0) > 0) {
        DlogK(ip, 0, 0, 0) =
            abs(log10(ExtMat(ip + 1, 0, 0)) - log10(ExtMat(ip, 0, 0)));
      } else {
        DlogK(ip, 0, 0, 0) = -1;
      }

      // calculate optical thickness
      tau(ip, 0, 0, 0) = ( ExtMat(ip + 1, 0, 0) + ExtMat(ip, 0, 0)) / 2. *
          abs(z_field(ip + 1, 0, 0) - z_field(ip, 0, 0));


      if (DlogK(ip, 0, 0, 0) > refine_crit && tau(ip, 0, 0, 0) > tau_crit) {
        refine[ip] = 1;
      } else {
        refine[ip] = 0;
      }
    }

  } else {
    refine.resize((Np - 1) * (Nlat - 1) * (Nlon - 1));

    throw runtime_error("3d is not implemented yet! Sorry");
  }
}

void PrepareSubdomains(
    ArrayOfRTDomain& Subdomains,
    ArrayOfRTDomainScatteringProperties& SubdomainsScatteringProperties,
    const RTDomain MainDomain,
    const RTDomainScatteringProperties MainDomainScatteringProperties,
    const Tensor3& t_field,
    const Tensor3& z_field,
    const Tensor4& vmr_field,
    const Tensor4& pnd_field,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& t_interp_order,
    const Index& ForwardCorrectionFlag,
    const ArrayOfIndex& refine,
    const Tensor4& DlogK,
    const Index& N_decade,
    const Verbosity& verbosity){

  if (MainDomain.get_AtmosphereDim() == 1) {
    PrepareSubdomains1D(Subdomains,
                        SubdomainsScatteringProperties,
                        MainDomain,
                        MainDomainScatteringProperties,
                        t_field,
                        z_field,
                        vmr_field,
                        pnd_field,
                        scat_data,
                        t_interp_order,
                        ForwardCorrectionFlag,
                        refine,
                        DlogK,
                        N_decade,
                        verbosity);
  } else {
    throw runtime_error("3d is not implemented yet! Sorry");
  }
}

void PrepareSubdomains1D(
    ArrayOfRTDomain& Subdomains,
    ArrayOfRTDomainScatteringProperties& SubdomainsScatteringProperties,
    const RTDomain MainDomain,
    const RTDomainScatteringProperties MainDomainScatteringProperties,
    const Tensor3& t_field,
    const Tensor3& z_field,
    const Tensor4& vmr_field,
    const Tensor4& pnd_field,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& t_interp_order,
    const Index& ForwardCorrectionFlag,
    const ArrayOfIndex& refine,
    const Tensor4& DlogK,
    const Index& N_decade,
    const Verbosity& verbosity) {
  const Index Nlayers = refine.nelem();

  Index Nz;
  Tensor3 z_field_sub;
  Tensor3 t_field_sub;
  Tensor4 vmr_field_sub;
  Tensor4 pnd_field_sub;
  Vector p_grid_sub;

  Vector dump;
  Numeric dz_sub;
  Matrix itw;
  ArrayOfGridPos gridpositions_z;

  Subdomains.resize(refine.nelem());
  SubdomainsScatteringProperties.resize(refine.nelem());

  //number of zenith angles of propagation direction
  Index N_za = MainDomain.get_za_grid().nelem();

  //number of azimuth angles of propagation direction
  Index N_aa = MainDomain.get_aa_grid().nelem();

  //Stokes dimension
  Index N_st = MainDomain.get_CloudboxField().ncols();

  //Number of incidence direction of scattering matrix. For 1D atmosphere the
  //number of propagation directions is the same as for the cloudbox field. For
  //3D it would be N_za*N_aa
  Index N_idir = MainDomainScatteringProperties.get_ScatteringMatrix().npages();


  for (Index il = 0; il < Nlayers; il++) {
    if (refine[il]) {
      // Get the number of sublayers, we want to have for the subdomain
      Nz = Index(ceil(DlogK(il, 0, 0, 0) * Numeric(N_decade)) + 1);

      // allocate subdomain fields
      z_field_sub.resize(Nz + 1, 1, 1);
      t_field_sub.resize(Nz + 1, 1, 1);
      vmr_field_sub.resize(vmr_field.nbooks(), Nz + 1, 1, 1);
      pnd_field_sub.resize(pnd_field.nbooks(), Nz + 1, 1, 1);

      p_grid_sub.resize(Nz + 1);

      // define z_field of subdomain
      nlinspace(dump, z_field(il , 0, 0), z_field(il + 1, 0, 0), Nz + 1);
      z_field_sub(joker, 0, 0) = dump;

      // set up interpolation
      gridpositions_z.resize(z_field_sub.npages());
      gridpos(gridpositions_z, z_field(joker, 0, 0), z_field_sub(joker, 0, 0));

      // store interpolation weights:
      itw.resize(gridpositions_z.nelem(), 2);
      interpweights(itw, gridpositions_z);

      // Now interpolate the p_grid to the subdamain
      itw2p(p_grid_sub, MainDomain.get_p_grid(), gridpositions_z, itw);

      // init sub rt domain
      Subdomains[il] = RTDomain(p_grid_sub,
                                MainDomain.get_lat_grid(),
                                MainDomain.get_lon_grid(),
                                MainDomain.get_za_grid(),
                                MainDomain.get_aa_grid(),
                                MainDomain.get_AtmosphereDim());



      //allocate
      Subdomains[il].get_CloudboxField().resize(Nz, 1, 1, N_za, 1, N_st);
      Subdomains[il].get_CloudboxScatteringField().resize(
          Nz, 1, 1, N_za, 1, N_st);


      //interpolate temperature field
      interp(
          t_field_sub(joker, 0, 0), itw, t_field(joker, 0, 0), gridpositions_z);

      // interpolate vmr field
      for (Index i_v = 0; i_v < vmr_field.nbooks(); i_v++) {
        interp(vmr_field_sub(i_v, joker, 0, 0),
               itw,
               vmr_field(i_v, joker, 0, 0),
               gridpositions_z);
      }

      // interpolate pnd field
      for (Index i_pnd = 0; i_pnd < pnd_field.nbooks(); i_pnd++) {
        interp(pnd_field_sub(i_pnd, joker, 0, 0),
               itw,
               pnd_field(i_pnd, joker, 0, 0),
               gridpositions_z);
      }


      CalcDomainExtAbsScatProperties(SubdomainsScatteringProperties[il],
                                     Subdomains[il],
                                     t_field,
                                     scat_data,
                                     pnd_field,
                                     MainDomain.get_CloudboxField().ncols(),
                                     MainDomain.get_AtmosphereDim(),
                                     MainDomain.get_za_grid(),
                                     MainDomain.get_aa_grid(),
                                     MainDomainScatteringProperties.get_ScatZaGrid(),
                                     MainDomainScatteringProperties.get_ScatAaGrid(),
                                     t_interp_order,
                                     ForwardCorrectionFlag,
                                     verbosity);



    }
  }
}

void RunNewDoit(  //Input and Output:
    RTDomain& MainDomain,
    Index& convergence_flag,
    Index& iteration_counter,
    RTDomainScatteringProperties& MainDomainScatteringProperties,
    //Input
    const ConstTensor6View& surface_reflection_matrix,
    const ConstTensor5View& surface_emission,
    const ArrayOfIndex& cloudbox_limits,
    const Index& atmosphere_dim,
    const Numeric& f_mono,
    const String& iy_unit,
    const Vector& epsilon,
    const Index& max_num_iterations,
    const Index& accelerated,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  Tensor6 cloudbox_field_mono_old;

  //TODO: move this into NewDoitMonoCalc
  MainDomain.get_CloudboxScatteringField().resize(
      MainDomain.get_CloudboxField().nvitrines(),
      MainDomain.get_CloudboxField().nshelves(),
      MainDomain.get_CloudboxField().nbooks(),
      MainDomain.get_CloudboxField().npages(),
      MainDomain.get_CloudboxField().nrows(),
      MainDomain.get_CloudboxField().ncols());

  convergence_flag = 0;
  iteration_counter = 0;

  // Array to save the last iteration steps
  ArrayOfTensor6 acceleration_input;
  if (accelerated) {
    acceleration_input.resize(4);
  }
  while (convergence_flag == 0) {

    // store cloudbox field, because we need it for the convergence check
    cloudbox_field_mono_old = MainDomain.get_CloudboxField();

    // Calculate the scattered field.
    out2 << "  Calculate scattering field. \n";
    MainDomainScatteringProperties.CalcScatteredField(MainDomain.get_CloudboxScatteringField(),
                       MainDomain.get_CloudboxField(),
                       atmosphere_dim,
                       verbosity);

    // Update cloudbox_field.
    out2 << "  Execute doit_rte_agenda. \n";
    Vector f_grid(1, f_mono);
    UpdateSpectralRadianceField(
                                MainDomain,
                                surface_reflection_matrix,
                                surface_emission,
                                cloudbox_limits,
                                atmosphere_dim,
                                f_grid,
                                verbosity);

    //Convergence test.
    CheckConvergence(convergence_flag,
                     iteration_counter,
                     MainDomain.get_CloudboxField(),
                     cloudbox_field_mono_old,
                     f_mono,
                     epsilon,
                     max_num_iterations,
                     iy_unit,
                     verbosity);

    // Convergence Acceleration, if wished.
    //TODO: I think this acceleration work only for 1D
    if (accelerated > 0 && convergence_flag == 0) {
      acceleration_input[(iteration_counter - 1) % 4] = MainDomain.get_CloudboxField();
      // NG - Acceleration
      if (iteration_counter % 4 == 0) {
        cloudbox_field_ngAcceleration(
            MainDomain.get_CloudboxField(), acceleration_input, accelerated, verbosity);
      }
    }
  }  //end of while loop, convergence is reached.
}

void UpdateSpectralRadianceField(//Input and Output:
                                 RTDomain& Domain,
                                 //Input:
                                 const ConstTensor6View& surface_reflection_matrix,
                                 const ConstTensor5View& surface_emission,
                                 const ArrayOfIndex& cloudbox_limits,
                                 const Index& atmosphere_dim,
                                 const Vector& f_grid,
                                 const Verbosity& verbosity) {
  if (atmosphere_dim == 1) {
    UpdateSpectralRadianceField1D(Domain,
                                  surface_reflection_matrix,
                                  surface_emission,
                                  cloudbox_limits,
                                  f_grid,
                                  verbosity);
  } else if (atmosphere_dim == 3) {
//    UpdateSpectralRadianceField3D(ws,
//                                  cloudbox_field_mono,
//                                  cloudbox_scat_field,
//                                  gas_extinction,
//                                  extinction_matrix,
//                                  absorption_vector,
//                                  cloudbox_limits,
//                                  za_grid,
//                                  aa_grid,
//                                  ppath_step_agenda,
//                                  ppath_lmax,
//                                  ppath_lraytrace,
//                                  p_path_maxlength,
//                                  p_grid,
//                                  z_field,
//                                  refellipsoid,
//                                  t_field,
//                                  f_grid,
//                                  verbosity);
  }
}

void UpdateSpectralRadianceField1D(
    //Input and Output:
    RTDomain& Domain,
    const ConstTensor6View& surface_reflection_matrix,
    const ConstTensor5View& surface_emission,
    const ArrayOfIndex& cloudbox_limits,
    const Vector& f_grid,
    const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  out2
      << "  UpdateSpectralRadianceField1D: Radiative transfer calculation in cloudbox\n";
  out2 << "  ------------------------------------------------------------- \n";

  // Number of zenith angles.
  const Index N_za = Domain.get_CloudboxField().npages();
  const Index N_p = Domain.get_CloudboxField().nvitrines();
  const Index stokes_dim = Domain.get_CloudboxField().ncols();

  // Epsilon for additional limb iterations
  Vector epsilon(4);
  epsilon[0] = 0.1;
  epsilon[1] = 0.01;
  epsilon[2] = 0.01;
  epsilon[3] = 0.01;

  Matrix cloudbox_field_limb;


  Index MaxUpwardAngleIndex=int(Numeric(N_za)/2.-0.5);

  //Loop over all directions, defined by za_grid
  for (Index i_za = 0; i_za < N_za; i_za++) {

    //======================================================================
    // Radiative transfer inside the cloudbox
    //=====================================================================

    // Sequential update for uplooking angles
    if (i_za <= MaxUpwardAngleIndex ) {
      // Loop over all positions inside the cloud box defined by the
      // cloudbox_limits excluding the upper boundary. For uplooking
      // directions, we start from cloudbox_limits[1]-1 and go down
      // to cloudbox_limits[0] to do a sequential update of the
      // radiation field
      for (Index i_p = N_p - 2; i_p >= 0; i_p--) {

        const Index idx = subscript2index(i_za, i_p, N_za);

        UpdateCloudPropagationPath1D(Domain,
                                     i_p,
                                     i_za,
                                     idx,
                                     f_grid,
                                     surface_reflection_matrix,
                                     surface_emission,
                                     cloudbox_limits,
                                     verbosity);
      }
    } else if (i_za > Domain.get_MaxLimbIndex()) {
      //
      // Sequential updating for downlooking angles
      //
      for (Index i_p = 0 + 1; i_p <= N_p-1; i_p++) {

        const Index idx = subscript2index(i_za, i_p, N_za);

        UpdateCloudPropagationPath1D(Domain,
                                     i_p,
                                     i_za,
                                     idx,
                                     f_grid,
                                     surface_reflection_matrix,
                                     surface_emission,
                                     cloudbox_limits,
                                     verbosity);
      }  // Close loop over p_grid (inside cloudbox).
    }    // end if downlooking.

    //
    // Limb looking:
    // We have to include a special case here, as we may miss the endpoints
    // when the intersection point is at the same level as the aactual point.
    // To be save we loop over the full cloudbox. Inside the function
    // cloud_ppath_update1D it is checked whether the intersection point is
    // inside the cloudbox or not.

    else {
      bool conv_flag = false;
      Index limb_it = 0;
      while (!conv_flag && limb_it < 10) {
        limb_it++;
        cloudbox_field_limb = Domain.get_CloudboxField()(joker, 0, 0, i_za, 0, joker);
        for (Index i_p = 0; i_p <= N_p-1; i_p++) {
          // For this case the cloudbox goes down to the surface and we
          // look downwards. These cases are outside the cloudbox and
          // not needed. Switch is included here, as ppath_step_agenda
          // gives an error for such cases.

          if (i_p != 0) {
            const Index idx = subscript2index(i_za, i_p, N_za);

            UpdateCloudPropagationPath1D(Domain,
                                         i_p,
                                         i_za,
                                         idx,
                                         f_grid,
                                         surface_reflection_matrix,
                                         surface_emission,
                                         cloudbox_limits,
                                         verbosity);
          }
        }

        conv_flag = true;
        for (Index i_p = 0; conv_flag && i_p < Domain.get_CloudboxField().nvitrines();
             i_p++) {
          for (Index stokes_index = 0; conv_flag && stokes_index < stokes_dim;
               stokes_index++) {
            Numeric diff = Domain.get_CloudboxField()(i_p, 0, 0, i_za, 0, stokes_index) -
                           cloudbox_field_limb(i_p, stokes_index);

            // If the absolute difference of the components
            // is larger than the pre-defined values, continue with
            // another iteration
            Numeric diff_bt = invrayjean(diff, f_grid[0]);
            if (abs(diff_bt) > epsilon[stokes_index]) {
              out2 << "Limb BT difference: " << diff_bt << " in stokes dim "
                   << stokes_index << "\n";

              conv_flag = false;
            }
          }
        }
      }
      out2 << "Limb iterations: " << limb_it << "\n";
    }
  }  // Closes loop over za_grid.
}

void UpdateSpectralRadianceField3D(Workspace& ws,
                                   // WS Input and Output:
                                   Tensor6& cloudbox_field_mono,
                                   Tensor6& cloudbox_scat_field,
                                   const ConstTensor3View& gas_extinction,
                                   const ConstTensor6View& extinction_matrix,
                                   const ConstTensor5View& absorption_vector,
                                   // WS Input:
                                   const ArrayOfIndex& cloudbox_limits,
                                   const Vector& za_grid,
                                   const Vector& aa_grid,
                                   // Propagation path calculation:
                                   const Agenda& ppath_step_agenda,
                                   const Numeric& ppath_lmax,
                                   const Numeric& ppath_lraytrace,
                                   const Tensor3& p_path_maxlength,
                                   const Vector& p_grid,
                                   const Tensor3& z_field,
                                   const Vector& refellipsoid,
                                   // Calculate thermal emission:
                                   const Tensor3& t_field,
                                   const Vector& f_grid,
                                   const Verbosity& verbosity) {

}

void UpdateCloudPropagationPath1D(
    RTDomain& Domain,
    const Index& p_index,
    const Index& za_index,
    const Index& idx,
    const ConstVectorView& f_grid,
    const ConstTensor6View& surface_reflection_matrix,
    const ConstTensor5View& surface_emission,
    const ArrayOfIndex& cloudbox_limits,
    const Verbosity& verbosity) {

  Index Npath = Domain.get_PressureArray()[idx].nelem();
  Index Ncloud = Domain.get_CloudboxScatteringField().nvitrines();

  // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.
  if ((0 <= Domain.get_GposPArray()[idx][1].idx &&
       Ncloud > Domain.get_GposPArray()[idx][1].idx) ||
      (Ncloud == Domain.get_GposPArray()[idx][1].idx &&
       abs(Domain.get_GposPArray()[idx][1].fd[0]) < 1e-6)) {
    // Stokes dimension
    const Index stokes_dim = Domain.get_CloudboxField().ncols();


    // Initialize variables for interpolated values
    Tensor3 ext_mat_int(stokes_dim, stokes_dim, Npath, 0.);
    Matrix abs_vec_int(stokes_dim, Npath, 0.);
    Matrix sca_vec_int(stokes_dim, Npath, 0.);
    Matrix cloudbox_field_mono_int(stokes_dim, Npath, 0.);


    InterpolateOnPropagation1D(
                         ext_mat_int,
                         abs_vec_int,
                         sca_vec_int,
                         cloudbox_field_mono_int,
                         Domain.get_ExtinctionMatrix()(joker, joker, joker, za_index, 0, joker, joker),
                         Domain.get_AbsorptionVector()(joker, joker, joker, za_index, 0, joker),
                         Domain.get_CloudboxScatteringField(),
                         Domain.get_CloudboxField(),
                         Domain.get_GposPArray()[idx],
                         Domain.get_GposZenithArray()[idx],
                         Domain.get_InterpWeightsArray()[idx],
                         Domain.get_InterpWeightsAngleArray()[idx],
                         verbosity);

    // ppath_what_background(ppath_step) tells the
    // radiative background.  More information in the
    // function get_iy_of_background.
    // if there is no background we proceed the RT
//    Index bkgr = ppath_what_background(ppath_step);

    // Radiative transfer from one layer to the next, starting
    // at the intersection with the next layer and propagating
    // to the considered point.
    RTStepInCloudNoBackground(Domain.get_CloudboxField(),
                              Domain.get_LStepArray()[idx],
                              Domain.get_TemperatureArray()[idx],
                              Domain.get_PressureArray()[idx],
                              Domain.get_GasExtinctionArray()[idx],
                              ext_mat_int,
                              abs_vec_int,
                              sca_vec_int,
                              cloudbox_field_mono_int,
                              cloudbox_limits,
                              f_grid,
                              p_index,
                              0,
                              0,
                              za_index,
                              0,
                              verbosity);

    // bkgr=2 indicates that the background is the surface
//    if (bkgr == 2) {
//       cout << "hit surface "<< ppath_step.gp_p << endl;
//      cloud_RT_surface(ws,
//                       cloudbox_field_mono,
//                       surface_rtprop_agenda,
//                       f_grid,
//                       stokes_dim,
//                       ppath_step,
//                       cloudbox_limits,
//                       za_grid,
//                       za_index);
//    }

  }  //end if inside cloudbox
}

void InterpolateOnPropagation1D(  //Output
    Tensor3View& ext_mat_int,
    MatrixView& abs_vec_int,
    MatrixView& sca_vec_int,
    MatrixView& cloudbox_field_mono_int,
    const ConstTensor5View& ext_mat_field,
    const ConstTensor4View& abs_vec_field,
    const ConstTensor6View& cloudbox_scat_field,
    const ConstTensor6View& cloudbox_field_mono,
    const ArrayOfGridPos& cloud_gp_p,
    const ArrayOfGridPos& cloud_gp_za,
    const MatrixView& itw,
    const MatrixView& itw_za,
    const Verbosity& verbosity) {

  CREATE_OUT3;

  // Stokes dimension
  const Index stokes_dim = cloudbox_field_mono.ncols();

  // For the interpolation, the precalculated weights and grid positions are
  // used.

  //Interpolate the quantities on the ppath_step
  for (Index i = 0; i < stokes_dim; i++) {

    // Interpolation of extinction matrix
    // Extinction matrix requires a second loop
    // over stokes_dim
    out3 << "Interpolate ext_mat:\n";
    for (Index j = 0; j < stokes_dim; j++) {

      interp(ext_mat_int(i, j, joker),
             itw,
             ext_mat_field(joker, 0, 0, i, j),
             cloud_gp_p);
    }

    // Interpolation of absorption vector
    out3 << "Interpolate abs_vec:\n";
    interp(
        abs_vec_int(i, joker), itw, abs_vec_field(joker, 0, 0, i), cloud_gp_p);


    out3 << "Interpolate cloudbox_scat_field and cloudbox_field_mono:\n";
    // Interpolate scattered field
    interp(sca_vec_int(i, joker),
             itw_za,
             cloudbox_scat_field(joker, 0, 0, joker, 0, i),
             cloud_gp_p,
             cloud_gp_za);

    // Interpolate radiation field
    interp(cloudbox_field_mono_int(i, joker),
             itw_za,
             cloudbox_field_mono(joker, 0, 0, joker, 0, i),
             cloud_gp_p,
             cloud_gp_za);

  }
}

void RTStepInCloudNoBackground(Tensor6View cloudbox_field_mono,
                               const ConstVectorView& lstep_ppath,
                               const ConstVectorView& temperature_ppath,
                               const ConstVectorView& pressure_ppath,
                               const ConstVectorView& gas_extinction_ppath,
                               const ConstTensor3View& ext_mat_int,
                               const ConstMatrixView& abs_vec_int,
                               const ConstMatrixView& sca_vec_int,
                               const ConstMatrixView& cloudbox_field_mono_int,
                               const ArrayOfIndex& cloudbox_limits,
                               const ConstVectorView& f_grid,
                               const Index& p_index,
                               const Index& lat_index,
                               const Index& lon_index,
                               const Index& za_index,
                               const Index& aa_index,
                               const Verbosity& verbosity) {
  CREATE_OUT3;

  const Index stokes_dim = cloudbox_field_mono.ncols();
  const Index atmosphere_dim = cloudbox_limits.nelem() / 2;
  const Index Nppath = pressure_ppath.nelem();

  Vector sca_vec_av(stokes_dim, 0);
  Vector stokes_vec(stokes_dim, 0.);
  Vector rtp_temperature_nlte_dummy(0);

  // Two propmat_clearsky to average between
  PropagationMatrix cur_propmat_clearsky(1, stokes_dim);

  PropagationMatrix prev_propmat_clearsky(1, stokes_dim);

  PropagationMatrix ext_mat_local;
  StokesVector abs_vec_local;
  Matrix matrix_tmp(stokes_dim, stokes_dim);
  Vector vector_tmp(stokes_dim);

  // Incoming stokes vector
  stokes_vec = cloudbox_field_mono_int(joker, Nppath - 1);

  for (Index k = Nppath - 1; k >= 0; k--) {
    // Save propmat_clearsky from previous level
    std::swap(cur_propmat_clearsky, prev_propmat_clearsky);

    //Set current propmat clearsky
    cur_propmat_clearsky.Kjj() = gas_extinction_ppath[k];

    // Skip any further calculations for the first point.
    // We need values at two ppath points before we can average.
    if (k == Nppath - 1) {
      continue;
    }

    // Average prev_propmat_clearsky with cur_propmat_clearsky
    prev_propmat_clearsky += cur_propmat_clearsky;
    prev_propmat_clearsky *= 0.5;

    opt_prop_sum_propmat_clearsky(
        ext_mat_local, abs_vec_local, prev_propmat_clearsky);

    for (Index i = 0; i < stokes_dim; i++) {
      // Averaging of sca_vec:
      sca_vec_av[i] = 0.5 * (sca_vec_int(i, k) + sca_vec_int(i, k + 1));
    }

    // Add average particle absorption to abs_vec.
    abs_vec_local.AddAverageAtPosition(abs_vec_int(joker, k),
                                       abs_vec_int(joker, k + 1));

    // Add average particle extinction to ext_mat.
    ext_mat_local.AddAverageAtPosition(ext_mat_int(joker, joker, k),
                                       ext_mat_int(joker, joker, k + 1));

    // Frequency
    Numeric f = f_grid[0];

    // Calculate Planck function
    Numeric rte_planck_value =
        planck(f, 0.5 * (temperature_ppath[k] + temperature_ppath[k + 1]));

    // Length of the path between the two layers.
    Numeric lstep = lstep_ppath[k];

    // Some messages:
    if (out3.sufficient_priority()) {
      vector_tmp=abs_vec_local.VectorAtPosition();
      ext_mat_local.MatrixAtPosition(matrix_tmp);
      out3 << "-----------------------------------------\n";
      out3 << "Input for radiative transfer step \n"
           << "calculation inside"
           << " the cloudbox:"
           << "\n";
      out3 << "Stokes vector at intersection point: \n" << stokes_vec << "\n";
      out3 << "lstep: ..." << lstep << "\n";
      out3 << "------------------------------------------\n";
      out3 << "Averaged coefficients: \n";
      out3 << "Planck function: " << rte_planck_value << "\n";
      out3 << "Scattering vector: " << sca_vec_av << "\n";
      out3 << "Absorption vector: " << vector_tmp << "\n";
      out3 << "Extinction matrix: " << matrix_tmp << "\n";

      assert(!is_singular(matrix_tmp));
    }

    // Radiative transfer step calculation. The Stokes vector
    // is updated until the considered point is reached.
    RadiativeTransferStep(stokes_vec,
                          Matrix(stokes_dim, stokes_dim),
                          ext_mat_local,
                          abs_vec_local,
                          sca_vec_av,
                          lstep,
                          rte_planck_value);

  }  // End of loop over a ppath step.

  // Assign calculated Stokes Vector to cloudbox_field_mono.
  if (atmosphere_dim == 1)
    cloudbox_field_mono(p_index, 0, 0, za_index, 0, joker) = stokes_vec;
  else if (atmosphere_dim == 3)
    cloudbox_field_mono(p_index,
                        lat_index - cloudbox_limits[2],
                        lon_index - cloudbox_limits[4],
                        za_index,
                        aa_index,
                        joker) = stokes_vec;
}

void RadiativeTransferStep(  //Output and Input:
    //TODO: Remove trans_mat, it is not used within doit.
    VectorView stokes_vec,
    MatrixView trans_mat,
    //Input
    const PropagationMatrix& ext_mat_av,
    const StokesVector& abs_vec_av,
    const ConstVectorView& sca_vec_av,
    const Numeric& lstep,
    const Numeric& rtp_planck_value,
    const bool& trans_is_precalc) {
  //Stokes dimension:
  Index stokes_dim = stokes_vec.nelem();

  //Test sizes
  assert(ext_mat_av.NumberOfFrequencies() == 1 and
         abs_vec_av.NumberOfFrequencies() == 1);

  //Check inputs:
  assert(is_size(trans_mat, 1, stokes_dim, stokes_dim));
  assert(stokes_dim == ext_mat_av.StokesDimensions() and
         stokes_dim == abs_vec_av.StokesDimensions());
  assert(is_size(sca_vec_av, stokes_dim));
  assert(rtp_planck_value >= 0);
  assert(lstep >= 0);
  //assert (not ext_mat_av.AnySingular());  This is asserted at a later time in this version...

  // Check, if only the first component of abs_vec is non-zero:
  const bool unpol_abs_vec = abs_vec_av.IsUnpolarized(0);

  bool unpol_sca_vec = true;

  for (Index i = 1; i < stokes_dim; i++)
    if (sca_vec_av[i] != 0) unpol_sca_vec = false;

  // Calculate transmission by general function, if not precalculated
  Index extmat_case = 0;
  if (!trans_is_precalc) {
    compute_transmission_matrix_from_averaged_matrix_at_frequency(
        trans_mat, lstep, ext_mat_av, 0);
  }

  //--- Scalar case: ---------------------------------------------------------
  if (stokes_dim == 1) {
    stokes_vec[0] = stokes_vec[0] * trans_mat(0, 0) +
                    (abs_vec_av.Kjj()[0] * rtp_planck_value + sca_vec_av[0]) /
                    ext_mat_av.Kjj()[0] * (1 - trans_mat(0, 0));
  }

    //--- Vector case: ---------------------------------------------------------

    // We have here two cases, diagonal or non-diagonal ext_mat_gas
    // For diagonal ext_mat_gas, we expect abs_vec_gas to only have a
    // non-zero value in position 1.

    //- Unpolarised
  else if (extmat_case == 1 && unpol_abs_vec && unpol_sca_vec) {
    const Numeric invK = 1.0 / ext_mat_av.Kjj()[0];
    // Stokes dim 1
    stokes_vec[0] = stokes_vec[0] * trans_mat(0, 0) +
                    (abs_vec_av.Kjj()[0] * rtp_planck_value + sca_vec_av[0]) *
                    invK * (1 - trans_mat(0, 0));

    // Stokes dims > 1
    for (Index i = 1; i < stokes_dim; i++) {
      stokes_vec[i] = stokes_vec[i] * trans_mat(i, i) +
                      sca_vec_av[i] * invK * (1 - trans_mat(i, i));
    }
  }

    //- General case
  else {
    Matrix invK(stokes_dim, stokes_dim);
    ext_mat_av.MatrixInverseAtPosition(invK);

    Vector source = abs_vec_av.VectorAtPosition();
    source *= rtp_planck_value;

    for (Index i = 0; i < stokes_dim; i++)
      source[i] += sca_vec_av[i];  // b = abs_vec * B + sca_vec

    // solve K^(-1)*b = x
    Vector x(stokes_dim);
    mult(x, invK, source);

    Vector term1(stokes_dim);
    Vector term2(stokes_dim);

    Matrix ImT(stokes_dim, stokes_dim);
    id_mat(ImT);
    ImT -= trans_mat;
    mult(term2, ImT, x);  // term2: second term of the solution of the RTE with
    //fixed scattered field

    // term1: first term of solution of the RTE with fixed scattered field
    mult(term1, trans_mat, stokes_vec);

    for (Index i = 0; i < stokes_dim; i++)
      stokes_vec[i] = term1[i] + term2[i];  // Compute the new Stokes Vector
  }
}

void NewRTStepInCloudNoBackground(
    Tensor6View cloudbox_field_mono,
    const ConstVectorView& lstep_ppath,
    const ConstVectorView& temperature_ppath,
    const ConstVectorView& pressure_ppath,
    const ConstVectorView& gas_extinction_ppath,
    const ConstTensor3View& ext_mat_int,
    const ConstMatrixView& abs_vec_int,
    const ConstMatrixView& sca_vec_int,
    const ConstMatrixView& cloudbox_field_mono_int,
    const ArrayOfIndex& cloudbox_limits,
    const ConstVectorView& f_grid,
    const Index& p_index,
    const Index& lat_index,
    const Index& lon_index,
    const Index& za_index,
    const Index& aa_index,
    const Verbosity& verbosity) {
  const Index stokes_dim = cloudbox_field_mono.ncols();
  const Index atmosphere_dim = cloudbox_limits.nelem() / 2;
  const Index Nppath = pressure_ppath.nelem();

  Vector sca_vec_av(stokes_dim, 0);
  Vector stokes_vec(stokes_dim, 0.);
  Vector rtp_temperature_nlte_dummy(0);

  // Two propmat_clearsky to average between
  PropagationMatrix propmat_clearsky1(1, stokes_dim);

  PropagationMatrix ext_mat_01;
  PropagationMatrix ext_mat_0;
  PropagationMatrix ext_mat_1;
  StokesVector abs_vec_1;
  Vector sca_vec_1(stokes_dim, 0);
  Vector source_0(stokes_dim, 0);
  Vector source_1(stokes_dim, 0);
  Vector j_0(stokes_dim, 0);
  Vector j_1(stokes_dim, 0);
  Vector j_01(stokes_dim, 0);
  Vector Qmax(stokes_dim, 0);
  Matrix trans_mat_01(stokes_dim, stokes_dim, 0);
  Numeric tau_01 = 0.;
  Numeric rte_planck_value_1 = 0.;
  Numeric lstep = 0.;

  Matrix matrix_tmp(stokes_dim, stokes_dim);
  Vector vector_tmp(stokes_dim);

  // Incoming stokes vector
  stokes_vec = cloudbox_field_mono_int(joker, Nppath - 1);

  // Frequency
  Numeric f = f_grid[0];


  for (Index k = Nppath - 1; k >= 0; k--) {
    std::swap(source_1, source_0);
    std::swap(j_1, j_0);
    std::swap(ext_mat_1, ext_mat_0);


    //Set current propmat clearsky
    propmat_clearsky1.Kjj() = gas_extinction_ppath[k];

    opt_prop_sum_propmat_clearsky(ext_mat_1, abs_vec_1, propmat_clearsky1);

    //Scattered intensity
    for (Index i = 0; i < stokes_dim; i++) {
      sca_vec_1[i] = sca_vec_int(i, k);
    }

    // Add particle absorption to abs_vec.
    abs_vec_1 += (StokesVector)abs_vec_int(joker, k);

    // Add particle extinction to ext_mat.
    PropagationMatrix PM(ext_mat_int(joker, joker, k));
    ext_mat_1 += PM;

    rte_planck_value_1 = planck(f, temperature_ppath[k]);
    CalcSourceForRTStep(
        source_1, j_1, ext_mat_1, abs_vec_1, sca_vec_1, rte_planck_value_1);

    if (k == Nppath - 1) {
      continue;
    }



    //mean extinction between point 2 and  point 1
    ext_mat_01 = ext_mat_0;
    ext_mat_01 += ext_mat_1;
    ext_mat_01 *= 0.5;

    lstep = lstep_ppath[k];

    //calc optical thickness and transmission
    compute_transmission_matrix_from_averaged_matrix_at_frequency(
        trans_mat_01, lstep, ext_mat_01, 0);

    tau_01 = ext_mat_01.Kjj(0, 0)[0] * lstep;

    Matrix T(stokes_dim, stokes_dim);
    id_mat(T);
    T *= tau_01;

    j_01 = j_0;
    j_01 += j_1;
    j_01 *= 0.5;

    Qmax=j_01;
    Qmax*=lstep;

    RTSolver2ndOrder(stokes_vec,
                     source_0,
                     source_1,
                     trans_mat_01,
                     tau_01,
                     Qmax);


  }

  // Assign calculated Stokes Vector to cloudbox_field_mono.
  if (atmosphere_dim == 1)
    cloudbox_field_mono(p_index, 0, 0, za_index, 0, joker) = stokes_vec;
  else if (atmosphere_dim == 3)
    cloudbox_field_mono(p_index,
                        lat_index - cloudbox_limits[2],
                        lon_index - cloudbox_limits[4],
                        za_index,
                        aa_index,
                        joker) = stokes_vec;
}

void NewRTStepInCloudNoBackground2(
    Tensor6View cloudbox_field_mono,
    const ConstVectorView& lstep_ppath,
    const ConstVectorView& temperature_ppath,
    const ConstVectorView& pressure_ppath,
    const ConstVectorView& gas_extinction_ppath,
    const ConstTensor3View& ext_mat_int,
    const ConstMatrixView& abs_vec_int,
    const ConstMatrixView& sca_vec_int,
    const ConstMatrixView& cloudbox_field_mono_int,
    const ArrayOfIndex& cloudbox_limits,
    const ConstVectorView& f_grid,
    const Index& p_index,
    const Index& lat_index,
    const Index& lon_index,
    const Index& za_index,
    const Index& aa_index,
    const Verbosity& verbosity) {
  const Index stokes_dim = cloudbox_field_mono.ncols();
  const Index atmosphere_dim = cloudbox_limits.nelem() / 2;
  const Index Nppath = pressure_ppath.nelem();

  Vector sca_vec_av(stokes_dim, 0);
  Vector stokes_vec(stokes_dim, 0.);
  Vector rtp_temperature_nlte_dummy(0);

  // Two propmat_clearsky to average between
  PropagationMatrix propmat_clearsky2(1, stokes_dim);

  PropagationMatrix ext_mat_01;
  PropagationMatrix ext_mat_12;
  PropagationMatrix ext_mat_0;
  PropagationMatrix ext_mat_1;
  PropagationMatrix ext_mat_2;
  StokesVector abs_vec_2;
  Vector sca_vec_2(stokes_dim, 0);
  Vector source_0(stokes_dim, 0);
  Vector source_1(stokes_dim, 0);
  Vector source_2(stokes_dim, 0);
  Vector j_0(stokes_dim, 0);
  Vector j_1(stokes_dim, 0);
  Vector j_2(stokes_dim, 0);
  Vector j_01(stokes_dim, 0);
  Vector Qmax(stokes_dim, 0);
  Matrix trans_mat_01(stokes_dim, stokes_dim, 0);
  Matrix trans_mat_12(stokes_dim, stokes_dim, 0);
  Numeric tau_01 = 0.;
  Numeric tau_12 = 0.;
  Numeric rte_planck_value_2 = 0.;
  Numeric lstep = 0.;

  Matrix matrix_tmp(stokes_dim, stokes_dim);
  Vector vector_tmp(stokes_dim);

  // Incoming stokes vector
  stokes_vec = cloudbox_field_mono_int(joker, Nppath - 1);

  // Frequency
  Numeric f = f_grid[0];


  for (Index k = Nppath - 1; k >= 0; k--) {
    std::swap(source_1, source_0);
    std::swap(source_2, source_1);
    std::swap(j_1, j_0);
    std::swap(j_2, j_1);

    std::swap(ext_mat_1, ext_mat_0);
    std::swap(ext_mat_2, ext_mat_1);


    //Set current propmat clearsky
    propmat_clearsky2.Kjj() = gas_extinction_ppath[k];

    opt_prop_sum_propmat_clearsky(ext_mat_2, abs_vec_2, propmat_clearsky2);

    //Scattered intensity
    for (Index i = 0; i < stokes_dim; i++) {
      sca_vec_2[i] = sca_vec_int(i, k);
    }

    // Add particle absorption to abs_vec.
    abs_vec_2 += (StokesVector)abs_vec_int(joker, k);

    // Add particle extinction to ext_mat.
    PropagationMatrix PM(ext_mat_int(joker, joker, k));
    ext_mat_2 += PM;

    rte_planck_value_2 = planck(f, temperature_ppath[k]);
    CalcSourceForRTStep(
        source_2, j_2, ext_mat_2, abs_vec_2, sca_vec_2, rte_planck_value_2);

    if (k >= Nppath - 2) {
      continue;
    }



    //mean extinction between point 0 and  point 1
    ext_mat_01 = ext_mat_0;
    ext_mat_01 += ext_mat_1;
    ext_mat_01 *= 0.5;

    ext_mat_12 = ext_mat_1;
    ext_mat_12 += ext_mat_2;
    ext_mat_12 *= 0.5;

    lstep = lstep_ppath[k];

    //calc optical thickness and transmission
    compute_transmission_matrix_from_averaged_matrix_at_frequency(
        trans_mat_01, lstep, ext_mat_01, 0);



    tau_01 = ext_mat_01.Kjj(0, 0)[0] * lstep;
    tau_12 = ext_mat_12.Kjj(0, 0)[0] * lstep;

    j_01 = j_0;
    j_01 += j_1;
    j_01 *= 0.5;

    Qmax=j_01;
    Qmax*=lstep;

    RTSolver3rdOrder(stokes_vec,
                     source_0,
                     source_1,
                     source_2,
                     trans_mat_01,
                     tau_01,
                     tau_12,
                     Qmax);

    if (k==0){

      compute_transmission_matrix_from_averaged_matrix_at_frequency(
          trans_mat_12, lstep, ext_mat_12, 0);

      j_01 = j_1;
      j_01 += j_2;
      j_01 *= 0.5;

      Qmax=j_01;
      Qmax*=lstep;

      RTSolver2ndOrder(stokes_vec,
                       source_1,
                       source_2,
                       trans_mat_12,
                       tau_12,
                       Qmax);

    }


  }

  // Assign calculated Stokes Vector to cloudbox_field_mono.
  if (atmosphere_dim == 1)
    cloudbox_field_mono(p_index, 0, 0, za_index, 0, joker) = stokes_vec;
  else if (atmosphere_dim == 3)
    cloudbox_field_mono(p_index,
                        lat_index - cloudbox_limits[2],
                        lon_index - cloudbox_limits[4],
                        za_index,
                        aa_index,
                        joker) = stokes_vec;
}

void ShortCharacteristicsRT_step3rdOrder(
    Tensor6View cloudbox_field_mono,
    const ConstVectorView& lstep_ppath,
    const ConstVectorView& temperature_ppath,
    const ConstVectorView& gas_extinction_ppath,
    const ConstTensor3View& ext_mat_int,
    const ConstMatrixView& abs_vec_int,
    const ConstMatrixView& sca_vec_int,
    const ConstMatrixView& cloudbox_field_mono_int,
    const ConstVectorView& lstep_ppath_ip1,
    const ConstVectorView& temperature_ppath_ip1,
    const ConstVectorView& gas_extinction_ppath_ip1,
    const ConstTensor3View& ext_mat_int_ip1,
    const ConstMatrixView& abs_vec_int_ip1,
    const ConstMatrixView& sca_vec_int_ip1,
    const ArrayOfIndex& cloudbox_limits,
    const ConstVectorView& f_grid,
    const Index& p_index,
    const Index& lat_index,
    const Index& lon_index,
    const Index& za_index,
    const Index& aa_index,
    const Verbosity& verbosity) {
  const Index stokes_dim = cloudbox_field_mono.ncols();
  const Index atmosphere_dim = cloudbox_limits.nelem() / 2;
  const Index Nppath = temperature_ppath.nelem();

  Vector sca_vec_av(stokes_dim, 0);
  Vector stokes_vec(stokes_dim, 0.);
  Vector rtp_temperature_nlte_dummy(0);

  // Two propmat_clearsky to average between
  PropagationMatrix propmat_clearsky(1, stokes_dim);

  PropagationMatrix ext_mat_01;
  PropagationMatrix ext_mat_12;
  PropagationMatrix ext_mat_0;
  PropagationMatrix ext_mat_1;
  PropagationMatrix ext_mat_2;
  StokesVector abs_vec_0;
  StokesVector abs_vec_1;
  StokesVector abs_vec_2;
  Vector sca_vec_0(stokes_dim, 0);
  Vector sca_vec_1(stokes_dim, 0);
  Vector sca_vec_2(stokes_dim, 0);
  Vector source_0(stokes_dim, 0);
  Vector source_1(stokes_dim, 0);
  Vector source_2(stokes_dim, 0);
  Vector j_0(stokes_dim, 0);
  Vector j_1(stokes_dim, 0);
  Vector j_2(stokes_dim, 0);
  Vector j_01(stokes_dim, 0);
  Vector Qmax(stokes_dim, 0);
  Matrix trans_mat_01(stokes_dim, stokes_dim, 0);
  Matrix trans_mat_12(stokes_dim, stokes_dim, 0);
  Numeric tau_01 = 0.;
  Numeric tau_12 = 0.;
  Numeric rte_planck_value_0 = 0.;
  Numeric rte_planck_value_1 = 0.;
  Numeric rte_planck_value_2 = 0.;
  Numeric lstep01 = 0.;
  Numeric lstep12 = 0.;

  Matrix matrix_tmp(stokes_dim, stokes_dim);
  Vector vector_tmp(stokes_dim);

  // Incoming stokes vector
  stokes_vec = cloudbox_field_mono_int(joker, Nppath - 1);

  // Frequency
  Numeric f = f_grid[0];

  //Set current propmat clearsky
  propmat_clearsky.Kjj() = gas_extinction_ppath[Nppath - 1];
  opt_prop_sum_propmat_clearsky(ext_mat_0, abs_vec_0, propmat_clearsky);

  propmat_clearsky.Kjj() = gas_extinction_ppath[0];
  opt_prop_sum_propmat_clearsky(ext_mat_1, abs_vec_1, propmat_clearsky);

  propmat_clearsky.Kjj() = gas_extinction_ppath_ip1[0];
  opt_prop_sum_propmat_clearsky(ext_mat_2, abs_vec_2, propmat_clearsky);

  //Scattered intensity
  for (Index i = 0; i < stokes_dim; i++) {
    sca_vec_0[i] = sca_vec_int(i, Nppath-1);
    sca_vec_1[i] = sca_vec_int(i, 0);
    sca_vec_2[i] = sca_vec_int_ip1(i, 0);
  }

  // Add particle absorption to abs_vec.
  abs_vec_0 += (StokesVector)abs_vec_int(joker, Nppath-1);
  abs_vec_1 += (StokesVector)abs_vec_int(joker, 0);
  abs_vec_2 += (StokesVector)abs_vec_int_ip1(joker, 0);

  // Add particle extinction to ext_mat.
  PropagationMatrix PM0(ext_mat_int(joker, joker, Nppath-1));
  ext_mat_0 += PM0;

  PropagationMatrix PM1(ext_mat_int(joker, joker, 0));
  ext_mat_1 += PM1;

  PropagationMatrix PM2(ext_mat_int_ip1(joker, joker, 0));
  ext_mat_2 += PM2;

  rte_planck_value_0 = planck(f, temperature_ppath[Nppath-1]);
  rte_planck_value_1 = planck(f, temperature_ppath[0]);
  rte_planck_value_2 = planck(f, temperature_ppath_ip1[0]);

  CalcSourceForRTStep(
      source_0, j_0, ext_mat_0, abs_vec_0, sca_vec_0, rte_planck_value_0);

  CalcSourceForRTStep(
      source_1, j_1, ext_mat_1, abs_vec_1, sca_vec_1, rte_planck_value_1);

  CalcSourceForRTStep(
      source_2, j_2, ext_mat_2, abs_vec_2, sca_vec_2, rte_planck_value_2);

  //mean extinction between point 0 and  point 1
  ext_mat_01 = ext_mat_0;
  ext_mat_01 += ext_mat_1;
  ext_mat_01 *= 0.5;

  ext_mat_12 = ext_mat_1;
  ext_mat_12 += ext_mat_2;
  ext_mat_12 *= 0.5;

  lstep01 = lstep_ppath[0];
  lstep12 = lstep_ppath_ip1[0];

  //calc optical thickness and transmission
  compute_transmission_matrix_from_averaged_matrix_at_frequency(
      trans_mat_01, lstep01, ext_mat_01, 0);

  tau_01 = ext_mat_01.Kjj(0, 0)[0] * lstep01;
  tau_12 = ext_mat_12.Kjj(0, 0)[0] * lstep12;

  j_01 = j_0;
  j_01 += j_1;
  j_01 *= 0.5;

  Qmax=j_01;
  Qmax*=lstep01;

  RTSolver3rdOrder(stokes_vec,
                   source_0,
                   source_1,
                   source_2,
                   trans_mat_01,
                   tau_01,
                   tau_12,
                   Qmax);


  // Assign calculated Stokes Vector to cloudbox_field_mono.
  if (atmosphere_dim == 1)
    cloudbox_field_mono(p_index, 0, 0, za_index, 0, joker) = stokes_vec;
  else if (atmosphere_dim == 3)
    cloudbox_field_mono(p_index,
                        lat_index - cloudbox_limits[2],
                        lon_index - cloudbox_limits[4],
                        za_index,
                        aa_index,
                        joker) = stokes_vec;
}


void CalcSourceForRTStep(  //Output:
    VectorView source,
    VectorView j, //emission density
    //Input
    const PropagationMatrix& ext_mat,
    const StokesVector& abs_vec,
    const ConstVectorView& sca_vec,
    const Numeric& rtp_planck_value) {
  //Stokes dimension:
  Index stokes_dim = sca_vec.nelem();

  //--- Scalar case: ---------------------------------------------------------
  if (stokes_dim == 1) {

    j[0]= (abs_vec.Kjj()[0] * rtp_planck_value + sca_vec[0]);

    source[0] = j[0] / ext_mat.Kjj()[0];
  }
    //- General case
  else {
    Matrix invK(stokes_dim, stokes_dim);
    ext_mat.MatrixInverseAtPosition(invK);

    j = abs_vec.VectorAtPosition();
    j *= rtp_planck_value;

    for (Index i = 0; i < stokes_dim; i++)
      j[i] += sca_vec[i];  // b = abs_vec * B + sca_vec

    // solve K^(-1)*b = x
    mult(source, invK, j);
  }
}

void RTSolver2ndOrder(  //in and out:
    VectorView stokes_vec,
    ConstVectorView& source_0,
    ConstVectorView& source_1,
    ConstMatrixView& trans_mat_01,
    const Numeric& tau_01,
    ConstVectorView& Qmax) {
  Index stokes_dim = source_0.nelem();

  Vector Q2nd(stokes_dim, 0);
  Numeric Q;

  Numeric e0 = 1 - exp(-tau_01);
  Numeric e1 = tau_01 - e0;

  Numeric u = e0 -e1 / tau_01 ;
  Numeric v = e1 /tau_01;

  //    Q3rd=u*source_0+v*source_1+w*source_2;
  for (Index j = 0; j < stokes_dim; j++) {
    Q2nd[j] = u * source_0[j] + v * source_1[j] ;
  }

  //--- Scalar case: ---------------------------------------------------------
  if (stokes_dim == 1) {
    Q = std::min(Q2nd[0], Qmax[0]);
    Q = std::max(Q, 0.);

    stokes_vec *= trans_mat_01(0, 0);
    stokes_vec += Q;
    //      stokes_vec[0] = stokes_vec[0]*trans_mat_01 (0, 0)+Q;

  } else {
    Vector term1(stokes_dim);
    Vector term2(stokes_dim);

    // term1: Attenuated Input
    mult(term1, trans_mat_01, stokes_vec);

    for (Index j = 0; j < stokes_dim; j++) {
      Q = std::min(Q2nd[j], Qmax[j]);
      Q = std::max(Q, 0.);

      stokes_vec[j] = term1[j] + Q;  // Compute the new Stokes Vector
    }
  }
}

void RTSolver3rdOrder(  //in and out:
    VectorView stokes_vec,
    ConstVectorView& source_0,
    ConstVectorView& source_1,
    ConstVectorView& source_2,
    ConstMatrixView& trans_mat_01,
    const Numeric& tau_01,
    const Numeric& tau_12,
    ConstVectorView& Qmax) {
  Index stokes_dim = source_0.nelem();

  Vector Q3rd(stokes_dim, 0);
  Numeric Q;

  Numeric e0 = 1 - exp(-tau_01);
  Numeric e1 = tau_01 - e0;
  Numeric e2 = tau_01 * tau_01 - 2 * e1;

  Numeric u =
      e0 + (e2 - (2. * tau_01 + tau_12) * e1) / (tau_01 * (tau_01 + tau_12));
  Numeric v = ((tau_01 + tau_12) * e1 - e2) / (tau_01 * tau_12);
  Numeric w = (e2 - tau_01 * e1) / (tau_12 * (tau_01 + tau_12));

  //    Q3rd=u*source_0+v*source_1+w*source_2;
  for (Index j = 0; j < stokes_dim; j++) {
    Q3rd[j] = u * source_0[j] + v * source_1[j] + w * source_2[j];
  }

  //--- Scalar case: ---------------------------------------------------------
  if (stokes_dim == 1) {
    Q = std::min(Q3rd[0], Qmax[0]);
    Q = std::max(Q, 0.);

    stokes_vec *= trans_mat_01(0, 0);
    stokes_vec += Q;
    //      stokes_vec[0] = stokes_vec[0]*trans_mat_01 (0, 0)+Q;

  } else {
    Vector term1(stokes_dim);
    Vector term2(stokes_dim);

    // term1: Attenuated Input
    mult(term1, trans_mat_01, stokes_vec);

    for (Index j = 0; j < stokes_dim; j++) {
      Q = std::min(Q3rd[j], Qmax[j]);
      Q = std::max(Q, 0.);

      stokes_vec[j] = term1[j] + Q;  // Compute the new Stokes Vector
    }
  }
}


void CheckConvergence(
    Index& convergence_flag,
    Index& iteration_counter,
    Tensor6& cloudbox_field_mono,
    const Tensor6& cloudbox_field_mono_old,
    const Numeric& f_mono,
    const Vector& epsilon,
    const Index& max_iterations,
    const String& iy_unit,
    const Verbosity& verbosity) {


  CREATE_OUT1;
  CREATE_OUT2;

  const Index N_p = cloudbox_field_mono.nvitrines();
  const Index N_lat = cloudbox_field_mono.nshelves();
  const Index N_lon = cloudbox_field_mono.nbooks();
  const Index N_za = cloudbox_field_mono.npages();
  const Index N_aa = cloudbox_field_mono.nrows();
  const Index stokes_dim = cloudbox_field_mono.ncols();

  bool Tb_unit = (iy_unit != "1");


  iteration_counter += 1;
  out2 << "  Number of DOIT iteration: " << iteration_counter << "\n";


  if (iteration_counter > max_iterations) {
    ostringstream out;
    out << "At frequency " << f_mono/1e9 << " GHz \n"
        << "method does not converge (number of iterations \n"
        << "is > " << max_iterations << "). "
        << " number density \n"
        << "is too large or the numerical setup for the DOIT \n"
        << "calculation is not correct. In case of limb \n"
        << "simulations please make sure that you use an \n"
        << "optimized zenith angle grid. \n"
        << "*cloudbox_field* might be wrong.\n";
    out1 << "Warning in DOIT calculation (output set to NaN):\n" << out.str();
//    cloudbox_field_mono = NAN;
    convergence_flag = 1;

  } else {
    Numeric diff;

    for (Index i_p = 0; i_p < N_p; i_p++) {
      for (Index i_lat = 0; i_lat < N_lat; i_lat++) {
        for (Index i_lon = 0; i_lon < N_lon; i_lon++) {
          for (Index i_za = 0; i_za < N_za; i_za++) {
            for (Index i_aa = 0; i_aa < N_aa; i_aa++) {
              for (Index i_s = 0; i_s < stokes_dim; i_s++) {
                diff =
                    cloudbox_field_mono(i_p, i_lat, i_lon, i_za, i_aa, i_s) -
                    cloudbox_field_mono_old(i_p, i_lat, i_lon, i_za, i_aa, i_s);

                // If the absolute difference of the components
                // is larger than the pre-defined values, return
                // to *cloudbox_fieldIterate* and do next iteration
                if (Tb_unit) {
                  diff = invrayjean(diff, f_mono);
                  if (abs(diff) > epsilon[i_s]) {
                    out1 << "BT difference: " << diff << " in stokes dim "
                         << i_s << "\n";
                    return;
                  }
                } else {
                  if (abs(diff) > epsilon[i_s]) {
                    out1 << "difference: " << diff << "\n";
                    return;
                  }
                }
              }  // End loop stokes_dom.
            }    // End loop scat_aa_grid.
          }      // End loop scat_za_grid.
        }        // End loop lon_grid.
      }          // End loop lat_grid.
    }            // End p_grid.
    convergence_flag = 1;
  }
}


void EstimatePPathElements1D(
    Workspace& ws,
    ArrayOfVector& PressureArray,
    ArrayOfVector& TemperatureArray,
    ArrayOfMatrix& VmrArray,
    ArrayOfMatrix& InterpWeightsArray,
    ArrayOfMatrix& InterpWeightsZenithArray,
    ArrayOfArrayOfGridPos& GposArray,
    ArrayOfArrayOfGridPos& GposZenithArray,
    ArrayOfVector& LstepArray,
    Index& MaxLimbIndex,
    const ArrayOfIndex& cloudbox_limits,
    const Vector& za_grid,
    const Agenda& ppath_step_agenda,
    const Numeric& ppath_lmax,
    const Numeric& ppath_lraytrace,
//    const Tensor3& p_path_maxlength,
    const Vector& p_grid,
    const Tensor3& z_field,
    const ConstTensor3View& t_field,
    const ConstTensor4View& vmr_field,
    const Vector& refellipsoid,
    const Vector& f_grid,
    const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  out2
      << "  EstimatePressurePoints1D for absorption\n";
  out2 << "  ------------------------------------------------------------- \n";

  // Number of zenith angles.
  const Index N_za = za_grid.nelem();
  const Index N_p = p_grid.nelem();


  // If theta is between 90° and the limiting value, the intersection point
  // is exactly at the same level as the starting point (cp. AUG)
  Numeric theta_lim =
      180. - asin((refellipsoid[0] + z_field(0, 0, 0)) /
                  (refellipsoid[0] + z_field(N_p-1, 0, 0))) *
             RAD2DEG;


  //Get the highest zenith angle index , which is still a limb angle
  MaxLimbIndex= N_za-1;
  while (MaxLimbIndex>=0 && za_grid[MaxLimbIndex]>=theta_lim){
    MaxLimbIndex--;
  }

  if (za_grid[MaxLimbIndex]<=90){
    MaxLimbIndex=-1;
  }



  PressureArray.resize(N_p*N_za);
  TemperatureArray.resize(N_p*N_za);
  VmrArray.resize(N_p*N_za);
  InterpWeightsArray.resize(N_p*N_za);
  InterpWeightsZenithArray.resize(N_p*N_za);
  GposArray.resize(N_p*N_za);
  GposZenithArray.resize(N_p*N_za);
  LstepArray.resize(N_p*N_za);

//  const bool adaptive = (p_path_maxlength.npages() > 0);
  Numeric ppath_lmax_temp = ppath_lmax;
  Numeric ppath_lraytrace_temp = ppath_lraytrace;


  //Loop over all directions, defined by za_grid
  for (Index i_za = 0; i_za < N_za; i_za++) {



    // uplooking angles
    if (za_grid[i_za] <= 90.) {
      // Loop over all positions inside the cloud box defined by the
      // cloudbox_limits excluding the upper boundary. For uplooking
      // directions, we start from cloudbox_limits[1]-1 and go down
      // to cloudbox_limits[0] to do a sequential update of the
      // radiation field
      for (Index i_p = N_p - 2; i_p >= 0; i_p--) {
//        if (adaptive) {
//          ppath_lmax_temp = p_path_maxlength(i_p, 0, 0);
//          ppath_lraytrace_temp = p_path_maxlength(i_p, 0, 0);
//        }

        const Index idx = subscript2index(i_za, i_p, N_za);

        Vector Pressures;
        Vector Temperatures;
        Matrix Vmrs;
        ArrayOfGridPos cloud_gp_p;
        ArrayOfGridPos cloud_gp_za;
        Vector lstep;
        Matrix itw;
        Matrix itw_za;
        CloudPropagationPath1D(ws,
                             Pressures,
                             Temperatures,
                             Vmrs,
                             cloud_gp_p,
                             cloud_gp_za,
                             lstep,
                             itw,
                             itw_za,
                             i_p,
                             i_za,
                             za_grid,
                             cloudbox_limits,
                             ppath_step_agenda,
                             ppath_lmax_temp,
                             ppath_lraytrace_temp,
                             p_grid,
                             z_field,
                             t_field,
                             vmr_field,
                             refellipsoid,
                             f_grid,
                             verbosity);

        PressureArray[idx] = Pressures;
        TemperatureArray[idx] = Temperatures;
        VmrArray[idx] = Vmrs;
        InterpWeightsArray[idx] = itw;
        InterpWeightsZenithArray[idx] = itw_za;
        GposArray[idx] = cloud_gp_p;
        GposZenithArray[idx] = cloud_gp_za;
        LstepArray[idx] = lstep;
      }
    } else if (za_grid[i_za] >= theta_lim) {
      //
      // downlooking angles
      //
      for (Index i_p = 0 + 1; i_p <= N_p-1; i_p++) {
//        if (adaptive) {
//          ppath_lmax_temp = p_path_maxlength(i_p, 0, 0);
//          ppath_lraytrace_temp = p_path_maxlength(i_p, 0, 0);
//        }

        const Index idx = subscript2index(i_za, i_p, N_za);

        Vector Pressures;
        Vector Temperatures;
        Matrix Vmrs;
        ArrayOfGridPos cloud_gp_p;
        ArrayOfGridPos cloud_gp_za;
        Vector lstep;
        Matrix itw;
        Matrix itw_za;
        CloudPropagationPath1D(ws,
                               Pressures,
                               Temperatures,
                               Vmrs,
                               cloud_gp_p,
                               cloud_gp_za,
                               lstep,
                               itw,
                               itw_za,
                               i_p,
                               i_za,
                               za_grid,
                               cloudbox_limits,
                               ppath_step_agenda,
                               ppath_lmax_temp,
                               ppath_lraytrace_temp,
                               p_grid,
                               z_field,
                               t_field,
                               vmr_field,
                               refellipsoid,
                               f_grid,
                               verbosity);

        PressureArray[idx] = Pressures;
        TemperatureArray[idx] = Temperatures;
        VmrArray[idx] = Vmrs;
        InterpWeightsArray[idx] = itw;
        InterpWeightsZenithArray[idx] = itw_za;
        GposArray[idx] = cloud_gp_p;
        GposZenithArray[idx] = cloud_gp_za;
        LstepArray[idx] = lstep;

      }  // Close loop over p_grid (inside cloudbox).
    }    // end if downlooking.

      //
      // Limb looking:
      // We have to include a special case here, as we may miss the endpoints
      // when the intersection point is at the same level as the aactual point.
      // To be save we loop over the full cloudbox. Inside the function
      // cloud_ppath_update1D it is checked whether the intersection point is
      // inside the cloudbox or not.

    else {
      for (Index i_p = 0; i_p <= N_p - 1; i_p++) {
        // For this case the cloudbox goes down to the surface and we
        // look downwards. These cases are outside the cloudbox and
        // not needed. Switch is included here, as ppath_step_agenda
        // gives an error for such cases.

        if (i_p != 0) {
//          if (adaptive) {
//            ppath_lmax_temp = p_path_maxlength(i_p, 0, 0);
//            ppath_lraytrace_temp = p_path_maxlength(i_p, 0, 0);
//          }

          const Index idx = subscript2index(i_za, i_p, N_za);

          Vector Pressures;
          Vector Temperatures;
          Matrix Vmrs;
          ArrayOfGridPos cloud_gp_p;
          ArrayOfGridPos cloud_gp_za;
          Vector lstep;
          Matrix itw;
          Matrix itw_za;
          CloudPropagationPath1D(ws,
                                 Pressures,
                                 Temperatures,
                                 Vmrs,
                                 cloud_gp_p,
                                 cloud_gp_za,
                                 lstep,
                                 itw,
                                 itw_za,
                                 i_p,
                                 i_za,
                                 za_grid,
                                 cloudbox_limits,
                                 ppath_step_agenda,
                                 ppath_lmax_temp,
                                 ppath_lraytrace_temp,
                                 p_grid,
                                 z_field,
                                 t_field,
                                 vmr_field,
                                 refellipsoid,
                                 f_grid,
                                 verbosity);

          PressureArray[idx] = Pressures;
          TemperatureArray[idx] = Temperatures;
          VmrArray[idx] = Vmrs;
          InterpWeightsArray[idx] = itw;
          InterpWeightsZenithArray[idx] = itw_za;
          GposArray[idx] = cloud_gp_p;
          GposZenithArray[idx] = cloud_gp_za;
          LstepArray[idx] = lstep;

        }
      }
    }
  }
}

void CloudPropagationPath1D(Workspace& ws,
                          Vector& Pressures,
                          Vector& Temperatures,
                          Matrix& Vmrs,
                          ArrayOfGridPos& cloud_gp_p,
                          ArrayOfGridPos& cloud_gp_za,
                          Vector& lstep,
                          Matrix& itw,
                          Matrix& itw_za,
                          const Index& p_index,
                          const Index& za_index,
                          const ConstVectorView& za_grid,
                          const ArrayOfIndex& cloudbox_limits,
                          const Agenda& ppath_step_agenda,
                          const Numeric& ppath_lmax,
                          const Numeric& ppath_lraytrace,
                          const ConstVectorView& p_grid,
                          const ConstTensor3View& z_field,
                          const ConstTensor3View& t_field,
                          const ConstTensor4View& vmr_field,
                          const ConstVectorView& refellipsoid,
                          const ConstVectorView& f_grid,
                          const Verbosity& verbosity) {

  CREATE_OUT3;

  Ppath ppath_step;

  //Initialize ppath for 1D.
  ppath_init_structure(ppath_step, 1, 1);
  // See documentation of ppath_init_structure for understanding
  // the parameters.

  // Assign value to ppath.pos:
  ppath_step.pos(0, 0) = z_field(p_index, 0, 0);
  ppath_step.r[0] = refellipsoid[0] + z_field(p_index, 0, 0);

  // Define the direction:
  ppath_step.los(0, 0) = za_grid[za_index];

  // Define the grid positions:
  ppath_step.gp_p[0].idx = p_index+cloudbox_limits[0];
  ppath_step.gp_p[0].fd[0] = 0;
  ppath_step.gp_p[0].fd[1] = 1;

  // Call ppath_step_agenda:
  ppath_step_agendaExecute(ws,
                           ppath_step,
                           ppath_lmax,
                           ppath_lraytrace,
                           f_grid,
                           ppath_step_agenda);


  // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.

  if ((cloudbox_limits[0] <= ppath_step.gp_p[1].idx &&
       cloudbox_limits[1] > ppath_step.gp_p[1].idx) ||
      (cloudbox_limits[1] == ppath_step.gp_p[1].idx &&
       abs(ppath_step.gp_p[1].fd[0]) < 1e-6)) {

    // Ppath_step normally has 2 points, the starting
    // point and the intersection point.
    // But there can be points in between, when a maximum
    // lstep is given. We have to interpolate on all the
    // points in the ppath_step.
    cloud_gp_p = ppath_step.gp_p;

    //get the right indices for the cloudbox
    for (Index i = 0; i < ppath_step.np; i++)
      cloud_gp_p[i].idx -= cloudbox_limits[0];

    // Grid index for points at upper limit of cloudbox must be shifted
    const Index n1 = cloudbox_limits[1] - cloudbox_limits[0];
    gridpos_upperend_check(cloud_gp_p[0], n1);
    gridpos_upperend_check(cloud_gp_p[ppath_step.np - 1], n1);

    //Additionally to the gridpos, we need only lstep from ppath_step
    lstep = ppath_step.lstep;

    //Precalculate interpolation weights and gridpos for angular interpolation
    Vector los_grid = ppath_step.los(joker, 0);
    cloud_gp_za.resize(los_grid.nelem());
    gridpos(cloud_gp_za, za_grid, los_grid);

    itw_za.resize(cloud_gp_p.nelem(), 4);
    interpweights(itw_za, cloud_gp_p, cloud_gp_za);


    // interpolate pressure, temperature and vmr
    itw.resize(cloud_gp_p.nelem(), 2);
    interpweights(itw, cloud_gp_p);

    //pressure interpolation
    out3 << "Interpolate pressure for ppath\n";
    Pressures.resize(ppath_step.np);
    itw2p(Pressures, p_grid, cloud_gp_p, itw);

    //temperature interpolation
    out3 << "Interpolate temperature for ppath\n";
    Temperatures.resize(ppath_step.np);
    interp(Temperatures,itw, t_field(joker,0,0), cloud_gp_p);

    Vmrs.resize(vmr_field.nbooks(),ppath_step.np);
    for (Index i_v = 0; i_v < vmr_field.nbooks(); i_v++){
      interp(Vmrs(i_v,joker),itw, vmr_field(i_v,joker,0,0), cloud_gp_p);
    }




  }  //end if inside cloudbox
}

//------------------------------------------------------------------------------
// Auxilary functions-----------------------------------------------------------
//
//------------------------------------------------------------------------------


void FlattenMeshGrid(Matrix& flattened_meshgrid,
                     const Vector& grid1,
                     const Vector& grid2) {

  const Index N_grid1 = grid1.nelem();
  const Index N_grid2 = grid2.nelem();

  // preparing output
  flattened_meshgrid.resize(N_grid1 * N_grid2, 2);

  //Flatten grids
  Index Idir_idx = 0;
  for (Index i_aa = 0; i_aa < N_grid2; i_aa++) {
    for (Index i_za = 0; i_za < N_grid1; i_za++) {
      flattened_meshgrid(Idir_idx, 0) = grid1[i_za];
      flattened_meshgrid(Idir_idx, 1) = grid2[i_aa];

      Idir_idx += 1;
    }
  }
}

void FlattenMeshGridIndex(ArrayOfIndex& index_grid1,
                          ArrayOfIndex& index_grid2,
                     const Vector& grid1,
                     const Vector& grid2) {

  const Index N_grid1 = grid1.nelem();
  const Index N_grid2 = grid2.nelem();

  // preparing output
  index_grid1.resize(N_grid1*N_grid2);
  index_grid2.resize(N_grid1*N_grid2);

  //Flatten grids
  Index Idir_idx = 0;
  for (Index i_aa = 0; i_aa < N_grid2; i_aa++) {
    for (Index i_za = 0; i_za < N_grid1; i_za++) {
      index_grid1[Idir_idx] = i_za;
      index_grid2[Idir_idx] = i_aa;

      Idir_idx += 1;
    }
  }
}


// Subscript to index functions-------------------------------------------------

//2d array
Index subscript2index(Index i, Index j, Index Ni) { return j * Ni + i; }

//3d array
Index subscript2index(Index i, Index j, Index k, Index Ni, Index Nj) {
  return k * Nj * Ni + subscript2index(i, j, Ni);
}

//4d array
Index subscript2index(
    Index i, Index j, Index k, Index l, Index Ni, Index Nj, Index Nk) {
  return l * Nk * Nj * Ni + subscript2index(i, j, k, Ni, Nj);
}

//5d array
Index subscript2index(Index i,
                     Index j,
                     Index k,
                     Index l,
                     Index m,
                     Index Ni,
                     Index Nj,
                     Index Nk,
                     Index Nl) {
  return m * Nl * Nk * Nj * Ni + subscript2index(i, j, k, l, Ni, Nj, Nk);
}


// Index to subscript functions-------------------------------------------------
//2d array
Tuple2D index2subscript(Index idx, Index Ni) {
  const Index j = int(idx / Ni);

  const Index i = idx - j * Ni;

  return std::make_tuple(i, j);
}

//3d array
Tuple3D index2subscript(Index idx, Index Ni, Index Nj) {
  const Index k = int(idx / (Ni * Nj));

  Tuple2D t = index2subscript(idx - k * Ni * Nj, Ni);

  Index i = std::get<0>(t);
  Index j = std::get<1>(t);

  return std::make_tuple(i, j, k);
}

//4d array
Tuple4D index2subscript(Index idx, Index Ni, Index Nj, Index Nk) {
  const Index l = int(idx / (Ni * Nj * Nk));

  Tuple3D t = index2subscript(idx - l * Ni * Nj * Nk, Ni, Nj);

  Index i = std::get<0>(t);
  Index j = std::get<1>(t);
  Index k = std::get<2>(t);

  return std::make_tuple(i, j, k, l);
}

//5d array
Tuple5D index2subscript(Index idx, Index Ni, Index Nj, Index Nk, Index Nl) {
  const Index m = int(idx / (Ni * Nj * Nk * Nl));

  Tuple4D t = index2subscript(idx - m * Ni * Nj * Nk * Nl, Ni, Nj, Nk);

  Index i = std::get<0>(t);
  Index j = std::get<1>(t);
  Index k = std::get<2>(t);
  Index l = std::get<3>(t);

  return std::make_tuple(i, j, k, l, m);
}