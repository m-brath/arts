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
                     Index& convergence_flag,
                     Index& iteration_counter,
                     //Input
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
                     const Numeric& ppath_lmax,
                     const Numeric& ppath_lraytrace,
                     const Verbosity& verbosity);

//TODO:Add doxygen doc
void CalcGasExtinction(Workspace& ws,
                       //Output
                       Tensor3& gas_extinct,
                       Vector& p_grid_abs,
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
    const Vector& za_grid,
    const ConstTensor4View& pnd_field,
    const ConstTensor3View& t_field,
    const Index& stokes_dim);

//TODO:Add doxygen doc
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
                           const ConstVectorView& f_grid,
                           const ConstVectorView& za_grid,
                           const ConstVectorView& aa_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           const Index& atmosphere_dim,
                           const Index& stokes_dim,
                           const Matrix& surface_field);

//TODO:Add doxygen doc
void CalcPropagationPathMaxLength(
    Tensor3& p_path_maxlength,
    const Tensor6View& extinction_matrix,  //(Np,Nlat,Nlon,ndir,nst,nst)
    const ConstVectorView& p_grid,
    const ConstTensor3View& gas_extinct,
    const ConstVectorView& p_grid_abs,
    const ConstVectorView& lat_grid,
    const ConstVectorView& lon_grid,
    const Vector& scat_za_grid,
    const Numeric& tau_max);

//TODO:Add doxygen doc
void RunNewDoit(Workspace& ws,
                Tensor6& doit_i_field_mono,
                Index& convergence_flag,
                Index& iteration_counter,
                const ConstTensor3View& gas_extinction,
                const ConstTensor6View& extinction_matrix,
                const ConstTensor5View& absorption_vector,
                const ConstTensor7View& scattering_matrix,
                const ConstTensor6View& surface_reflection_matrix,
                const ConstTensor5View& surface_emission,
                const ArrayOfIndex& cloudbox_limits,
                const Agenda& ppath_step_agenda,
                const Index& atmosphere_dim,
                const Tensor3& t_field,
                const Tensor3& z_field,
                const Vector& p_grid,
                const Vector& p_grid_abs,
                const Vector& za_grid,
                const Vector& aa_grid,
                const Vector& scat_za_grid,
                const Vector& scat_aa_grid,
                ArrayOfIndex& idir_idx0,
                ArrayOfIndex& idir_idx1,
                ArrayOfIndex& pdir_idx0,
                ArrayOfIndex& pdir_idx1,
                ArrayOfGridPos& gp_za_i,
                ArrayOfGridPos& gp_aa_i,
                Tensor3& itw,
                const Numeric& f_mono,
                const String& iy_unit,
                const Numeric& ppath_lmax,
                const Numeric& ppath_lraytrace,
                const Tensor3& p_path_maxlength,
                const Vector& refellipsoid,
                const Vector& epsilon,
                const Index& max_num_iterations,
                const Index& accelerated,
                const Verbosity& verbosity);


//TODO:Add doxygen doc
void CalcScatteredField(  // WS Output and Input
    Tensor6& doit_scat_field,
    //WS Input:
    const Tensor6& doit_i_field_mono,
    const Tensor7& scattering_matrix,
    const Index& atmosphere_dim,
    const Vector& scat_za_grid,
    const Vector& scat_aa_grid,
    ArrayOfIndex& idir_idx0,
    ArrayOfIndex& idir_idx1,
    ArrayOfIndex& pdir_idx0,
    ArrayOfIndex& pdir_idx1,
    ArrayOfGridPos& gp_za_i,
    ArrayOfGridPos& gp_aa_i,
    Tensor3& itw,
    const Verbosity& verbosity);

//TODO:Add doxygen doc
void CalcScatteredField1D(
    // Output
    Tensor6& doit_scat_field,
    // Input:
    const Tensor6& doit_i_field_mono,
    const Tensor7& scattering_matrix,
    const Vector& iza_grid,  // incoming direction
    ArrayOfIndex& pdir_idx0,//index array of propagation direction
    ArrayOfGridPos& gp_za_i, // grid pos for zenith angle interpolation
    Tensor3& itw, //interpolation weights
    const Verbosity& verbosity);

//TODO:Add doxygen doc
void CalcScatteredField3D(
    // WS Output and Input
    Tensor6& doit_scat_field,
    //WS Input:
    const Tensor6& doit_i_field_mono,
    const Tensor7& scattering_matrix,
    const Vector& iza_grid,  // incoming direction
    const Vector& iaa_grid,  // incoming direction
    const ArrayOfIndex& idir_idx0, //index array of flattened inc. direction
    const ArrayOfIndex& idir_idx1, //index array of flattened inc. direction
    const ArrayOfIndex& pdir_idx0, //index array of flattened propagation direction
    const ArrayOfIndex& pdir_idx1, //index array of flattened propagation direction
    const ArrayOfGridPos& gp_za_i, // grid pos for zenith angle interpolation
    const ArrayOfGridPos& gp_aa_i, // grid pos for azimuth angle interpolation
    const Tensor3& itw, //interpolation weights
    const Verbosity& verbosity);


//TODO:Add doxygen doc
void UpdateSpectralRadianceField(
    Workspace& ws,
    // WS Input and Output:
    Tensor6& doit_i_field_mono,
    Tensor6& doit_scat_field,
    // WS Input:
    const ConstTensor3View& gas_extinction,
    const ConstTensor6View& extinction_matrix,
    const ConstTensor5View& absorption_vector,
    const ConstTensor6View& surface_reflection_matrix,
    const ConstTensor5View& surface_emission,
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
    const Vector& p_grid_abs,
    const Tensor3& z_field,
    const Vector& refellipsoid,
    // Calculate thermal emission:
    const Tensor3& t_field,
    const Vector& f_grid,
    const Verbosity& verbosity);

//TODO:Add doxygen doc
void UpdateSpectralRadianceField1D(
    Workspace& ws,
    // WS Input and Output:
    Tensor6& doit_i_field_mono,
    Tensor6& doit_scat_field,
    // WS Input:
    const ConstTensor3View& gas_extinction,
    const ConstTensor6View& extinction_matrix,  //(Np,Nlat,Nlon,ndir,nst,nst)
    const ConstTensor5View& absorption_vector,  //(Np,Nlat,Nlon,ndir,nst)
    const ConstTensor6View& surface_reflection_matrix,
    const ConstTensor5View& surface_emission,
    const ArrayOfIndex& cloudbox_limits,
    const Vector& za_grid,
    // Propagation path calculation:
    const Agenda& ppath_step_agenda,
    const Numeric& ppath_lmax,
    const Numeric& ppath_lraytrace,
    const Tensor3& p_path_maxlength,
    const Vector& p_grid,
    const Vector& p_grid_abs,
    const Tensor3& z_field,
    const Vector& refellipsoid,
    // Calculate thermal emission:
    const Tensor3& t_field,
    const Vector& f_grid,
    const Verbosity& verbosity);

//TODO:Add doxygen doc
void UpdateSpectralRadianceField3D(
    Workspace& ws,
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
    const Verbosity& verbosity);

//TODO:Add doxygen doc
void UpdateCloudPropagationPath1D(
    Workspace& ws,
    Tensor6View& doit_i_field_mono,
    const Index& p_index,
    const Index& za_index,
    const ConstVectorView& za_grid,
    const ArrayOfIndex& cloudbox_limits,
    const ConstTensor6View& doit_scat_field,
    const Agenda& ppath_step_agenda,
    const Numeric& ppath_lmax,
    const Numeric& ppath_lraytrace,
    const ConstVectorView& p_grid,
    const ConstVectorView& p_grid_abs,
    const ConstTensor3View& z_field,
    const ConstVectorView& refellipsoid,
    const ConstTensor3View& t_field,
    const ConstVectorView& f_grid,
    const ConstTensor5View& ext_mat_field,
    const ConstTensor4View& abs_vec_field,
    const ConstTensor3View& gas_extinction,
    const ConstTensor6View& surface_reflection_matrix,
    const ConstTensor5View& surface_emission,
    const Verbosity& verbosity);

//TODO:Add doxygen doc
void InterpolateOnPropagation1D(  //Output
    VectorView&  gas_abs_int,
    Tensor3View& ext_mat_int,
    MatrixView& abs_vec_int,
    MatrixView& sca_vec_int,
    MatrixView& doit_i_field_mono_int,
    VectorView& t_int,
    const ConstTensor3View& gas_extinction,
    const ConstTensor5View& ext_mat_field,
    const ConstTensor4View& abs_vec_field,
    const ConstTensor6View& doit_scat_field,
    const ConstTensor6View& doit_i_field_mono,
    const ConstVectorView& p_grid,
    const ConstVectorView& p_grid_abs,
    const ConstTensor3View& t_field,
    const Ppath& ppath_step,
    const ArrayOfIndex& cloudbox_limits,
    const ConstVectorView& za_grid,
    const Verbosity& verbosity);

//TODO:Add doxygen doc
void RTStepInCloudNoBackground(Tensor6View doit_i_field_mono,
                               const Ppath& ppath_step,
                               const ConstVectorView& t_int,
                               const VectorView& gas_abs_int,
                               const ConstTensor3View& ext_mat_int,
                               const ConstMatrixView& abs_vec_int,
                               const ConstMatrixView& sca_vec_int,
                               const ConstMatrixView& doit_i_field_mono_int,
                               const ArrayOfIndex& cloudbox_limits,
                               const ConstVectorView& f_grid,
                               const Index& p_index,
                               const Index& lat_index,
                               const Index& lon_index,
                               const Index& za_index,
                               const Index& aa_index,
                               const Verbosity& verbosity);

//TODO:Add doxygen doc
void RadiativeTransferStep(  //Output and Input:
    VectorView stokes_vec,
    MatrixView trans_mat,
    //Input
    const PropagationMatrix& ext_mat_av,
    const StokesVector& abs_vec_av,
    const ConstVectorView& sca_vec_av,
    const Numeric& lstep,
    const Numeric& rtp_planck_value,
    const bool& trans_is_precalc = false);

//TODO:Add doxygen doc
void CheckConvergence(  //WS Input and Output:
    Index& convergence_flag,
    Index& iteration_counter,
    Tensor6&  doit_i_field_mono,
    const Tensor6& doit_i_field_mono_old,
    const Numeric& f_mono,
    const Vector& epsilon,
    const Index& max_iterations,
    const String& iy_unit,
    const Verbosity& verbosity);

//TODO:Add doxygen doc
void FlattenMeshGrid(Matrix& flattened_meshgrid,
                     const Vector& grid1,
                     const Vector& grid2);

//TODO:Add doxygen doc
void FlattenMeshGridIndex(ArrayOfIndex& index_grid1,
                          ArrayOfIndex& index_grid2,
                          const Vector& grid1,
                          const Vector& grid2);
