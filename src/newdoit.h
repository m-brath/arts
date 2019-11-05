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
#include "agenda_class.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "doit.h"
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

//TODO:Add doxygen doc
void Initialize_doit_i_field(
    //Output
    Tensor7& doit_i_field,
    //Input
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Vector& f_grid,
    const Vector& rt_za_grid,
    const Vector& rt_aa_grid,
    const ArrayOfIndex& cloudbox_limits);

//TODO:Add doxygen doc
void SetAngularGrids(
    //Output
    Vector& za_grid,
    Vector& aa_grid,
    Vector& scat_za_grid,
    Vector& scat_aa_grid,
    //Input
    const Index& N_za_grid,
    const Index& N_aa_grid,
    const Index& N_scat_za_grid,
    const Index& N_scat_aa_grid,
    const String& za_grid_type);

//TODO:Add doxygen doc
void SetClearsky_doit_i_field(
    //Output
    Tensor7& doit_i_field,
    //Input
    const Vector& f_grid,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const ArrayOfIndex& cloudbox_limits,
    const Index& atmosphere_dim,
    const Verbosity& verbosity);

//TODO:Add doxygen doc
void GetIncomingRadiation(Workspace& ws,
                          //Output
                          Tensor7& doit_i_field,
                          //Input
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
                          const Verbosity& verbosity);

//TODO:Add doxygen doc
void LimitInputGridsAndFieldsToCloudbox(
                                        //Output
                                        Vector& p_grid_cldbx,
                                        Vector& lat_grid_cldbx,
                                        Vector& lon_grid_cldbx,
                                        Tensor3& t_field_cldbx,
                                        Tensor3& z_field_cldbx,
                                        Tensor4& vmr_field_cldbx,
                                        //Input
                                        ConstVectorView p_grid,
                                        ConstVectorView lat_grid,
                                        ConstVectorView lon_grid,
                                        ConstTensor3View t_field,
                                        ConstTensor3View z_field,
                                        ConstTensor4View vmr_field,
                                        const ArrayOfIndex& cloudbox_limits,
                                        const Verbosity& verbosity);

//TODO:Add doxygen doc
void NewDoitMonoCalc(Workspace& ws,
                     //Output
                     Tensor6& doit_i_field_mono,
                     Tensor3& gas_extinction,
                     Tensor6& extinction_matrix,
                     Tensor5& absorption_vector,
                     Tensor7& scattering_matrix,
                     //Input
                     const ArrayOfIndex& cloudbox_limits,
                     const Agenda& propmat_clearsky_agenda,
                     const Agenda& surface_rtprop_agenda,
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
                     const Index& max_lvl_optimize,
                     const Numeric& tau_scat_max,
                     const Numeric& sgl_alb_max,
                     const Verbosity& verbosity);

//TODO:Add doxygen doc
void CalcGasExtinction(Workspace& ws,
                       //Output
                       Tensor3& gas_extinct,
                       //Input
                       const Agenda& propmat_clearsky_agenda,
                       const ConstTensor3View& t_field,
                       const ConstTensor4View& vmr_field,
                       const ConstVectorView& p_grid,
                       const ConstVectorView& lat_grid,
                       const ConstVectorView& lon_grid,
                       const ConstVectorView& f_mono);

//TODO:Add doxygen doc
void CalcParticleOpticalProperties(
    //Output
    Tensor6& extinction_matrix,
    Tensor5& absorption_vector,
    //input
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Vector& scat_za_grid,
    const Index& f_index,
    const ConstTensor4View& pnd_field,
    const ConstTensor3View& t_field,
    const Index& stokes_dim);

//TODO:Add doxygen doc
void CalcScatteringProperties(
    //Output
    Tensor7& scattering_matrix,
    //Input
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
    const Verbosity& verbosity);

//TODO:Add doxygen doc
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
                           const Matrix& surface_field);