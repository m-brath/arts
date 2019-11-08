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

#include "newdoit.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
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

void Initialize_doit_i_field(  //Output
    Tensor7& doit_i_field,
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
    doit_i_field.resize(Nf, Np_cloud, 1, 1, Nrza, 1, Ns);
  } else if (atmosphere_dim == 3) {
    const Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
    const Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
    const Index Nraa = rt_aa_grid.nelem();

    doit_i_field.resize(Nf, Np_cloud, Nlat_cloud, Nlon_cloud, Nrza, Nraa, Ns);
  }

  doit_i_field = NAN;
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

void SetClearsky_doit_i_field(Tensor7& doit_i_field,
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
    const Index nf = f_grid.nelem() ? doit_i_field.nlibraries() : 1;

    for (Index f_index = 0; f_index < nf; f_index++) {
      Index N_za = doit_i_field.npages();
      Index N_aa = doit_i_field.nrows();
      Index N_i = doit_i_field.ncols();

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
          doit_i_field(f_index, 0, joker, joker, joker, joker, joker);
      scat_i_p(1, joker, joker, joker, joker, joker) =
          doit_i_field(f_index,
                       doit_i_field.nvitrines() - 1,
                       joker,
                       joker,
                       joker,
                       joker,
                       joker);

      for (Index za_index = 0; za_index < N_za; ++za_index) {
        for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
          for (Index i = 0; i < N_i; ++i) {
            VectorView target_field = doit_i_field(
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
          "Error in doit_i_fieldSetClearsky: For 3D "
          "all_frequencies option is not implemented \n");

    for (Index f_index = 0; f_index < doit_i_field.nvitrines(); f_index++) {
      Index N_p = doit_i_field.nvitrines();
      Index N_lat = doit_i_field.nshelves();
      Index N_lon = doit_i_field.nbooks();
      Index N_za = doit_i_field.npages();
      Index N_aa = doit_i_field.nrows();
      Index N_i = doit_i_field.ncols();

      Tensor6 scat_i_p(2, N_lat, N_lon, N_za, N_aa, N_i);
      scat_i_p(0, joker, joker, joker, joker, joker) =
          doit_i_field(f_index, 0, joker, joker, joker, joker, joker);
      scat_i_p(1, joker, joker, joker, joker, joker) =
          doit_i_field(f_index, N_p - 1, joker, joker, joker, joker, joker);

      Tensor6 scat_i_lat(N_p, 2, N_lon, N_za, N_aa, N_i);
      scat_i_lat(joker, 0, joker, joker, joker, joker) =
          doit_i_field(f_index, joker, 0, joker, joker, joker, joker);
      scat_i_lat(joker, 1, joker, joker, joker, joker) =
          doit_i_field(f_index, joker, N_lat - 1, joker, joker, joker, joker);

      Tensor6 scat_i_lon(N_p, N_lat, 2, N_za, N_aa, N_i);
      scat_i_lon(joker, joker, 0, joker, joker, joker) =
          doit_i_field(f_index, joker, joker, 0, joker, joker, joker);
      scat_i_lon(joker, joker, 1, joker, joker, joker) =
          doit_i_field(f_index, joker, joker, N_lon - 1, joker, joker, joker);

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
                VectorView target_field = doit_i_field(f_index,
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
                VectorView target_field = doit_i_field(f_index,
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
                VectorView target_field = doit_i_field(f_index,
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
                          Tensor7& doit_i_field,
                          const Agenda& iy_main_agenda,
                          const Index& atmosphere_dim,
                          const Vector& lat_grid,
                          const Vector& lon_grid,
                          const Tensor3& z_field,
                          const Tensor4& nlte_field,
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

    //--- Get doit_i_field at lower and upper boundary
    //    (boundary=0: lower, boundary=1: upper)
    for (Index boundary = 0; boundary <= 1; boundary++) {
      const Index boundary_index = boundary ? doit_i_field.nvitrines() - 1 : 0;
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
      doit_i_field(joker, boundary_index, 0, 0, 0, 0, joker) = iy;

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

        doit_i_field(joker, boundary_index, 0, 0, za_index, 0, joker) = iy;

        for (Index fi = 0; fi < Nf; fi++) {
          if (doit_i_field(fi, boundary_index, 0, 0, za_index - 1, 0, 0) /
                      doit_i_field(fi, boundary_index, 0, 0, za_index, 0, 0) >
                  maxratio ||
              doit_i_field(fi, boundary_index, 0, 0, za_index - 1, 0, 0) /
                      doit_i_field(fi, boundary_index, 0, 0, za_index, 0, 0) <
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
                 << doit_i_field(fi, boundary_index, 0, 0, za_index - 1, 0, 0)
                 << " and "
                 << doit_i_field(fi, boundary_index, 0, 0, za_index, 0, 0)
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

    //--- Get doit_i_field at lower and upper boundary
    //    (boundary=0: lower, boundary=1: upper)
    for (Index boundary = 0; boundary <= 1; boundary++) {
      const Index boundary_index = boundary ? doit_i_field.nvitrines() - 1 : 0;
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

              doit_i_field(joker,
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
      const Index boundary_index = boundary ? doit_i_field.nshelves() - 1 : 0;
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

              doit_i_field(joker,
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
      const Index boundary_index = boundary ? doit_i_field.nbooks() - 1 : 0;
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

              doit_i_field(joker,
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
    ;
  }
}

void NewDoitMonoCalc(Workspace& ws,
    //Input and Output:
                     Tensor6& doit_i_field_mono,
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
                     const Index& f_index,
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     const Index& t_interp_order,
                     const String& iy_unit,
                     const Vector& refellipsoid,
                     const Vector& epsilon,
                     const Index& max_num_iterations,
                     const Numeric& tau_max,
                     const Index& accelerated,
                     const Numeric& ppath_lmax,
                     const Numeric& ppath_lraytrace,
                     const Verbosity& verbosity)

{

  CREATE_OUT0;

  //calculate gas extinction
  CalcGasExtinction(ws,
                    gas_extinction,
                    propmat_clearsky_agenda,
                    t_field,
                    vmr_field,
                    p_grid,
                    lat_grid,
                    lon_grid,
                    f_mono);

  ostringstream os, os1, os2;
  os << "gas absorption calculated \n";
  out0 << os.str();
  os.clear();

  //calculate par_optpropCalc_doit
  CalcParticleOpticalProperties(extinction_matrix,
                                absorption_vector,
                                scat_data,
                                za_grid,
                                f_index,
                                pnd_field,
                                t_field,
                                stokes_dim);

  os1 << "particle optical properties calculated \n";
  out0 << os1.str();
  os1.clear();

  //calculate sca_optpropCalc_doit
  CalcScatteringProperties(scattering_matrix,
                           t_field,
                           f_index,
                           scat_data,
                           pnd_field,
                           stokes_dim,
                           atmosphere_dim,
                           za_grid,
                           aa_grid,
                           scat_za_grid,
                           scat_aa_grid,
                           t_interp_order,
                           verbosity);

  os2 << "particle optical properties calculated \n";
  out0 << os2.str();
  os2.clear();


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

  Tensor3 p_path_maxlength;
  if (tau_max > 0) {
    //calculate local ppath_lmax
    CalcPropagationPathMaxLength(
        p_path_maxlength,
        extinction_matrix,  //(Np,Nlat,Nlon,ndir,nst,nst)
        p_grid,
        lat_grid,
        lon_grid,
        za_grid,
        tau_max);
  }

  //run new doit
  RunNewDoit(ws,
      doit_i_field_mono,
      convergence_flag,
      iteration_counter,
      gas_extinction,
      extinction_matrix,
      absorption_vector,
      scattering_matrix,
      cloudbox_limits,
      propmat_clearsky_agenda,
      surface_rtprop_agenda,
      ppath_step_agenda,
      atmosphere_dim,
      stokes_dim,
      t_field,
      z_field,
      z_surface,
      p_grid,
      lat_grid,
      lon_grid,
      za_grid,
      aa_grid,
      scat_za_grid,
      scat_aa_grid,
      f_mono,
      f_index,
      iy_unit,
      ppath_lmax,
      ppath_lraytrace,
      p_path_maxlength,
      refellipsoid,
      epsilon,
      max_num_iterations,
      accelerated,
      verbosity);

}

void CalcGasExtinction(Workspace& ws,
                     Tensor3& gas_extinct,
                     const Agenda& propmat_clearsky_agenda,
                     const ConstTensor3View& t_field,
                     const ConstTensor4View& vmr_field,
                     const ConstVectorView& p_grid,
                     const ConstVectorView& lat_grid,
                     const ConstVectorView& lon_grid,
                     const ConstVectorView& f_mono) {
  // Initialization
  gas_extinct = 0.;

  const Index Np = p_grid.nelem();
  const Index Nlat = lat_grid.nelem() > 0 ? lat_grid.nelem() : 1;
  const Index Nlon = lon_grid.nelem() > 0 ? lon_grid.nelem() : 1;

  // Local variables to be used in agendas

  ArrayOfPropagationMatrix propmat_clearsky_local;

  const Vector rtp_temperature_nlte_local_dummy(0);

  // Calculate layer averaged gaseous extinction
  for (Index ip = 0; ip < Np; ip++) {
    for (Index ilat = 0; ilat < Nlat; ilat++) {
      for (Index ilon = 0; ilon < Nlon; ilon++) {
        const Vector rtp_mag_dummy(3, 0);
        const Vector ppath_los_dummy;

        ArrayOfStokesVector nlte_dummy;
        ArrayOfPropagationMatrix partial_dummy;
        ArrayOfStokesVector partial_source_dummy, partial_nlte_dummy;
        propmat_clearsky_agendaExecute(ws,
                                       propmat_clearsky_local,
                                       nlte_dummy,
                                       partial_dummy,
                                       partial_source_dummy,
                                       partial_nlte_dummy,
                                       ArrayOfRetrievalQuantity(0),
                                       (Vector) f_mono,  // monochromatic calculation
                                       rtp_mag_dummy,
                                       ppath_los_dummy,
                                       p_grid[ip],
                                       t_field(ip, ilat, ilon),
                                       rtp_temperature_nlte_local_dummy,
                                       vmr_field(joker, ip, ilat, ilon),
                                       propmat_clearsky_agenda);

        //Assuming non-polarized light and only one frequency
        //TODO: Check if polarization is needed for absorption
        if (propmat_clearsky_local.nelem()) {
          for (Index j = 0; j < propmat_clearsky_local.nelem(); j++) {
            gas_extinct(ip, ilat, ilon) += propmat_clearsky_local[j].Kjj()[0];
          }
        }
      }
    }
  }
}

void CalcParticleOpticalProperties(Tensor6& extinction_matrix,//(Np,Nlat,Nlon,ndir,nst,nst)
                    Tensor5& absorption_vector,//(Np,Nlat,Nlon,ndir,nst)
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    const Vector& scat_za_grid,
                    const Index& f_index,
                    const ConstTensor4View& pnd_field,
                    const ConstTensor3View& t_field,
                    const Index& stokes_dim) {
  // Initialization
  extinction_matrix = 0.;
  absorption_vector = 0.;

  const Index Np = pnd_field.npages();
  const Index Nlat = pnd_field.nrows();
  const Index Nlon = pnd_field.ncols();


  assert(absorption_vector.nshelves() == Np);
  assert(extinction_matrix.nvitrines() == Np);

  // preparing input data
  Matrix dir_array(scat_za_grid.nelem(), 2, 0.);
  dir_array(joker, 0) = scat_za_grid;

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
                          f_index);

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

      extinction_matrix(joker, ilat, ilon, joker, joker, joker) =
          ext_mat_bulk_ii(0, joker, joker, joker, joker);

      absorption_vector(joker, ilat, ilon, joker, joker) =
          abs_vec_bulk_ii(0, joker, joker, joker);
    }
  }
}

void CalcScatteringProperties(  //Output
    Tensor7& scattering_matrix,
    const Tensor3& t_field,
    const Index& f_index,
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
  idir_array.resize(scat_za_grid.nelem() * scat_aa_grid.nelem(), 2);

  //Flatten incoming directions
  Index Idir_idx = 0;
  for (Index i_aa = 0; i_aa < scat_aa_grid.nelem(); i_aa++) {
    for (Index i_za = 0; i_za < scat_za_grid.nelem(); i_za++) {
      idir_array(Idir_idx, 0) = scat_za_grid[i_za];
      idir_array(Idir_idx, 1) = scat_aa_grid[i_aa];

      Idir_idx += 1;
    }
  }

  Matrix pdir_array;
  if (atmosphere_dim == 1) {
    pdir_array.resize(za_grid.nelem(), 2);
    pdir_array(joker, 0) = za_grid;
    pdir_array(joker, 1) = 0.;

  } else {
    pdir_array.resize(scat_za_grid.nelem() * scat_aa_grid.nelem(), 2);

    //Flatten outgoing/propagation directions
    Index Pdir_idx = 0;
    for (Index i_aa = 0; i_aa < aa_grid.nelem(); i_aa++) {
      for (Index i_za = 0; i_za < za_grid.nelem(); i_za++) {
        pdir_array(Idir_idx, 0) = za_grid[i_za];
        pdir_array(Idir_idx, 1) = aa_grid[i_aa];

        Pdir_idx += 1;
      }
    }
  }

  // making output containers
  ArrayOfArrayOfTensor6 sca_mat_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor6 sca_mat_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor6 scat_mat_bulk_ii;
  Index ptype_bulk_ii;

  if (atmosphere_dim == 1){
    scattering_matrix.resize(Np,
                             Nlat,
                             Nlon,
                             pdir_array.nrows(),
                             scat_za_grid.nelem(),
                             stokes_dim,
                             stokes_dim);
  } else {
    scattering_matrix.resize(Np,
                             Nlat,
                             Nlon,
                             pdir_array.nrows(),
                             idir_array.nrows(),
                             stokes_dim,
                             stokes_dim);
  }
  scattering_matrix=0.;


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
          f_index,
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


        os << "As we have a 1D atmosphere, we can integrate over \n"
           << "incoming azimuth. \n";
        out2 << os.str();
        os.clear();


        Tensor6 scat_mat_bulk_ii_temp;
        Index idx_i = 0;
        Index idx_i1 = 0;
        Numeric delta_phi = 360. *DEG2RAD / (scat_aa_grid.nelem() - 1);

        scat_mat_bulk_ii_temp.resize(1,
                                     Np,
                                     za_grid.nelem(),
                                     scat_za_grid.nelem(),
                                     stokes_dim,
                                     stokes_dim);
        scat_mat_bulk_ii_temp=0.;

        //integrate over incoming azimuth
        for (Index i_za = 0; i_za < scat_za_grid.nelem(); i_za++) {
          for (Index i_aa = 0; i_aa < (scat_aa_grid.nelem() - 1); i_aa++) {
            //flattened indices
            idx_i = i_aa * scat_za_grid.nelem() + i_za;
            idx_i1 = (i_aa + 1) * scat_za_grid.nelem() + i_za;

            scat_mat_bulk_ii_temp(joker, joker, joker, i_za, joker, joker) +=
                scat_mat_bulk_ii(joker, joker, joker, idx_i, joker, joker);

            scat_mat_bulk_ii_temp(joker, joker, joker, i_za, joker, joker) +=
                scat_mat_bulk_ii(joker, joker, joker, idx_i1, joker, joker);

            scat_mat_bulk_ii_temp(joker, joker, joker, i_za, joker, joker) /=
                2.;

            scat_mat_bulk_ii_temp(joker, joker, joker, i_za, joker, joker) *=
                delta_phi;
          }
        }

        scattering_matrix(joker, ilat, ilon, joker, joker, joker, joker) =
            scat_mat_bulk_ii_temp(0, joker, joker, joker, joker, joker);

      } else {
        scattering_matrix(joker, ilat, ilon, joker, joker, joker, joker) =
            scat_mat_bulk_ii(0, joker, joker, joker, joker, joker);
      }
    }
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
                      const ConstVectorView f_grid,
                      const ConstVectorView za_grid,
                      const ConstVectorView aa_grid,
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
            //for 3d atmosphere sensor_los has to angles.
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
    Tensor6& extinction_matrix,  //(Np,Nlat,Nlon,ndir,nst,nst)
    const ConstVectorView& p_grid,
    const ConstVectorView& lat_grid,
    const ConstVectorView& lon_grid,
    const Vector& scat_za_grid,
    const Numeric& tau_max) {


  const Index Np = p_grid.nelem();
  const Index Nlat = lat_grid.nelem() > 0 ? lat_grid.nelem() : 1;
  const Index Nlon = lon_grid.nelem() > 0 ? lon_grid.nelem() : 1;
  const Index Ndir = extinction_matrix.npages();
  Numeric ext_mat_elem;

  //prepare output container
  p_path_maxlength.resize(Np, Nlat, Nlon);
  p_path_maxlength = 0.;

  for (Index ip = 0; ip < Np; ip++) {
    for (Index ilat = 0; ilat < Nlat; ilat++) {
      for (Index ilon = 0; ilon < Nlon; ilon++) {
        ext_mat_elem = 0.;

        // calculate average over zenith direction,
        // for azimuthal random orientation, this is not fully correct,
        // but since this is more an rough estimation we assume, that the
        // optical thickness is only a function of the position.
        for (Index idir = 0; idir < Ndir; idir++) {
          ext_mat_elem += (extinction_matrix(ip, ilat, ilon, idir, 0, 0) *
                           sin(scat_za_grid[idir] * DEG2RAD));
        }
        ext_mat_elem /= Numeric(Ndir);

        // maximum propagation path length for each grid point.
        p_path_maxlength(ip, ilat, ilon) = tau_max / ext_mat_elem;
      }
    }
  }
}

void RunNewDoit(Workspace& ws,
                //Input and Output:
                Tensor6& doit_i_field_mono,
                Index& convergence_flag,
                Index& iteration_counter,
                const ConstTensor3View& gas_extinction,
                const ConstTensor6View& extinction_matrix,
                const ConstTensor5View& absorption_vector,
                const ConstTensor7View& scattering_matrix,
                const ArrayOfIndex& cloudbox_limits,
                const Agenda& propmat_clearsky_agenda,
                const Agenda& surface_rtprop_agenda,
                const Agenda& ppath_step_agenda,
                const Index& atmosphere_dim,
                const Index& stokes_dim,
                const Tensor3& t_field,
                const Tensor3& z_field,
                const Matrix& z_surface,
                const Vector& p_grid,
                const Vector& lat_grid,
                const Vector& lon_grid,
                const Vector& za_grid,
                const Vector& aa_grid,
                const Vector& scat_za_grid,
                const Vector& scat_aa_grid,
                const Numeric& f_mono,
                const Index& f_index,
                const String& iy_unit,
                const Numeric& ppath_lmax,
                const Numeric& ppath_lraytrace,
                const Tensor3& p_path_maxlength,
                const Vector& refellipsoid,
                const Vector& epsilon,
                const Index& max_num_iterations,
                const Index& accelerated,
                const Verbosity& verbosity) {
  CREATE_OUT2;

  for (Index v = 0; v < doit_i_field_mono.nvitrines(); v++)
    for (Index s = 0; s < doit_i_field_mono.nshelves(); s++)
      for (Index b = 0; b < doit_i_field_mono.nbooks(); b++)
        for (Index p = 0; p < doit_i_field_mono.npages(); p++)
          for (Index r = 0; r < doit_i_field_mono.nrows(); r++)
            for (Index c = 0; c < doit_i_field_mono.ncols(); c++)
              if (std::isnan(doit_i_field_mono(v, s, b, p, r, c)))
                throw std::runtime_error(
                    "*doit_i_field_mono* contains at least one NaN value.\n"
                    "This indicates an improper initialization of *doit_i_field*.");

  //doit_i_field_mono can not be further checked here, because there is no way
  //to find out the size without including a lot more interface
  //variables
  //-----------End of checks--------------------------------------

  Tensor6 doit_i_field_mono_old;

  // Resize and initialize doit_scat_field,
  // which  has the same dimensions as doit_i_field
  Tensor6 doit_scat_field(doit_i_field_mono.nvitrines(),
                          doit_i_field_mono.nshelves(),
                          doit_i_field_mono.nbooks(),
                          doit_i_field_mono.npages(),
                          doit_i_field_mono.nrows(),
                          doit_i_field_mono.ncols(),
                          0.);

  convergence_flag = 0;
  iteration_counter = 0;
  // Array to save the last iteration steps
  ArrayOfTensor6 acceleration_input;
  if (accelerated) {
    acceleration_input.resize(4);
  }
  while (convergence_flag == 0) {
    // 1. Copy doit_i_field to doit_i_field_old.
    doit_i_field_mono_old = doit_i_field_mono;

    // 2.Calculate scattered field vector for all points in the cloudbox.

    // Calculate the scattered field.
    out2 << "  Calculate scattering field. \n";
    CalcScatteredField(ws,
                       doit_scat_field,
                       doit_i_field_mono,
                       scattering_matrix,
                       atmosphere_dim,
                       cloudbox_limits,
                       za_grid,
                       aa_grid,
                       scat_za_grid,
                       scat_aa_grid,
                       verbosity);

    // Update doit_i_field.
    out2 << "  Execute doit_rte_agenda. \n";
    Vector f_grid(1, f_mono);
    UpdateSpectralRadianceField(ws,
                                doit_i_field_mono,
                                doit_scat_field,
                                gas_extinction,
                                extinction_matrix,
                                absorption_vector,
                                cloudbox_limits,
                                za_grid,
                                aa_grid,
                                atmosphere_dim,
                                ppath_step_agenda,
                                ppath_lmax,
                                ppath_lraytrace,
                                p_path_maxlength,
                                p_grid,
                                z_field,
                                refellipsoid,
                                t_field,
                                f_grid,
                                f_index,
                                verbosity);

    //Convergence test.
    ChackConvergence(convergence_flag,
                     iteration_counter,
                     doit_i_field_mono,
                     doit_i_field_mono_old,
                     epsilon,
                     max_num_iterations,
                     iy_unit,
                     verbosity);

    // Convergence Acceleration, if wished.
    if (accelerated > 0 && convergence_flag == 0) {
      acceleration_input[(iteration_counter - 1) % 4] = doit_i_field_mono;
      // NG - Acceleration
      if (iteration_counter % 4 == 0) {
        doit_i_field_ngAcceleration(
            doit_i_field_mono, acceleration_input, accelerated, verbosity);
      }
    }
  }  //end of while loop, convergence is reached.
}

void CalcScatteredField(Workspace& ws,
                        // WS Output and Input
                        Tensor6& doit_scat_field,
                        //WS Input:
                        const Tensor6& doit_i_field_mono,
                        const Tensor7& scattering_matrix,
                        const Index& atmosphere_dim,
                        const ArrayOfIndex& cloudbox_limits,
                        const Vector& za_grid,
                        const Vector& aa_grid,
                        const Vector& scat_za_grid,
                        const Vector& scat_aa_grid,
                        const Verbosity& verbosity) {




}

void UpdateSpectralRadianceField(Workspace& ws,
                                 // WS Input and Output:
                                 Tensor6& doit_i_field_mono,
                                 Tensor6& doit_scat_field,
                                 // WS Input:
                                 const ConstTensor3View& gas_extinction,
                                 const ConstTensor6View& extinction_matrix,
                                 const ConstTensor5View& absorption_vector,
                                 const ArrayOfIndex& cloudbox_limits,
                                 const Vector& za_grid,
                                 const Vector& aa_grid,
                                 const Index& atmosphere_dim,
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
                                 const Index& f_index,
                                 const Verbosity& verbosity) {
  if (atmosphere_dim == 1) {
    UpdateSpectralRadianceField1D(ws,
                                  doit_i_field_mono,
                                  doit_scat_field,
                                  gas_extinction,
                                  extinction_matrix,
                                  absorption_vector,
                                  cloudbox_limits,
                                  za_grid,
                                  aa_grid,
                                  ppath_step_agenda,
                                  ppath_lmax,
                                  ppath_lraytrace,
                                  p_path_maxlength,
                                  p_grid,
                                  z_field,
                                  refellipsoid,
                                  t_field,
                                  f_grid,
                                  f_index,
                                  verbosity);
  } else if (atmosphere_dim == 3) {
    UpdateSpectralRadianceField3D(ws,
                                  doit_i_field_mono,
                                  doit_scat_field,
                                  gas_extinction,
                                  extinction_matrix,
                                  absorption_vector,
                                  cloudbox_limits,
                                  za_grid,
                                  aa_grid,
                                  ppath_step_agenda,
                                  ppath_lmax,
                                  ppath_lraytrace,
                                  p_path_maxlength,
                                  p_grid,
                                  z_field,
                                  refellipsoid,
                                  t_field,
                                  f_grid,
                                  f_index,
                                  verbosity);
  }
}

void UpdateSpectralRadianceField1D(Workspace& ws,
                                   // WS Input and Output:
                                   Tensor6& doit_i_field_mono,
                                   Tensor6& doit_scat_field,
                                   // WS Input:
                                   const ConstTensor3View& gas_extinction,
                                   const ConstTensor6View& extinction_matrix,
                                   const ConstTensor5View& absorption_vector,
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
                                   const Index& f_index,
                                   const Verbosity& verbosity) {

}

void UpdateSpectralRadianceField3D(Workspace& ws,
                                   // WS Input and Output:
                                   Tensor6& doit_i_field_mono,
                                   Tensor6& doit_scat_field,
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
                                   const Index& f_index,
                                   const Verbosity& verbosity) {

}

void ChackConvergence(  //WS Input and Output:
    Index& convergence_flag,
    Index& iteration_counter,
    Tensor6& doit_i_field_mono,
    // WS Input:
    const Tensor6& doit_i_field_mono_old,
    // Keyword:
    const Vector& epsilon,
    const Index& max_iterations,
    const String& iy_unit,
    const Verbosity& verbosity) {

}