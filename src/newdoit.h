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
  \file   newdoit.h
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

class RTDomain {
 public:
  //TODO: Add doxygen doc and remove unneeded methods and constructors.

  RTDomain() = default;

  // constructor for 3d
//  RTDomain(const Vector& p_grid,
//           const Vector& lat_grid,
//           const Vector& lon_grid,
//           const Vector& za_grid,
//           const Vector& aa_grid,
//           const ArrayOfVector& PressureArray,
//           const ArrayOfVector& TemperatureArray,
//           const ArrayOfVector& GasExtinctionArray,
//           const ArrayOfMatrix& InterpWeightsArray,
//           const ArrayOfMatrix& InterpWeightsAngleArray,
//           const ArrayOfArrayOfGridPos& GposPArray,
//           const ArrayOfArrayOfGridPos& GposLatArray,
//           const ArrayOfArrayOfGridPos& GposLonArray,
//           const ArrayOfArrayOfGridPos& GposZenithArray,
//           const ArrayOfArrayOfGridPos& GposAzimuthArray,
//           const ArrayOfVector& LStepArray,
//           const Index& MaxLimbIndex,
//           const Index& AtmosphereDim)
//      : mp_grid(p_grid),
//        mlat_grid(lat_grid),
//        mlon_grid(lon_grid),
//        mza_grid(za_grid),
//        maa_grid(aa_grid),
//        mPressureArray(PressureArray),
//        mTemperatureArray(TemperatureArray),
//        mGasExtinctionArray(GasExtinctionArray),
//        mInterpWeightsArray(InterpWeightsArray),
//        mInterpWeightsAngleArray(InterpWeightsAngleArray),
//        mGposPArray(GposPArray),
//        mGposLatArray(GposLatArray),
//        mGposLonArray(GposLonArray),
//        mGposZenithArray(GposZenithArray),
//        mGposAzimuthArray(GposAzimuthArray),
//        mLStepArray(LStepArray),
//        mMaxLimbIndex(MaxLimbIndex),
//        mAtmosphereDim(AtmosphereDim) {}

  // constructor for 1d
//  RTDomain(const Vector& p_grid,
//           const Vector& lat_grid,
//           const Vector& lon_grid,
//           const Vector& za_grid,
//           const Vector& aa_grid,
//           const ArrayOfVector& PressureArray,
//           const ArrayOfVector& TemperatureArray,
//           const ArrayOfVector& GasExtinctionArray,
//           const ArrayOfMatrix& InterpWeightsArray,
//           const ArrayOfMatrix& InterpWeightsAngleArray,
//           const ArrayOfArrayOfGridPos& GposPArray,
//           const ArrayOfArrayOfGridPos& GposZenithArray,
//           const ArrayOfVector& LStepArray,
//           const Index& MaxLimbIndex,
//           const Index& AtmosphereDim)
//      : mp_grid(p_grid),
//        mlat_grid(lat_grid),
//        mlon_grid(lon_grid),
//        mza_grid(za_grid),
//        maa_grid(aa_grid),
//        mPressureArray(PressureArray),
//        mTemperatureArray(TemperatureArray),
//        mGasExtinctionArray(GasExtinctionArray),
//        mInterpWeightsArray(InterpWeightsArray),
//        mInterpWeightsAngleArray(InterpWeightsAngleArray),
//        mGposPArray(GposPArray),
//        mGposZenithArray(GposZenithArray),
//        mLStepArray(LStepArray),
//        mMaxLimbIndex(MaxLimbIndex),
//        mAtmosphereDim(AtmosphereDim) {}

  // constructor (only) grids
  RTDomain(const Vector& p_grid,
           const Vector& lat_grid,
           const Vector& lon_grid,
           const Vector& za_grid,
           const Vector& aa_grid,
           const Index& AtmosphereDim)
      : mp_grid(p_grid),
        mlat_grid(lat_grid),
        mlon_grid(lon_grid),
        mza_grid(za_grid),
        maa_grid(aa_grid),
        mAtmosphereDim(AtmosphereDim) {}

  // set and get methods of variables
  void set_p_grid(const Vector& p_grid) { mp_grid = p_grid; }

  void set_lat_grid(const Vector& lat_grid) { mlat_grid = lat_grid; }

  void set_lon_grid(const Vector& lon_grid) { mlon_grid = lon_grid; }

  void set_za_grid(const Vector& za_grid) { mza_grid = za_grid; }

  void set_aa_grid(const Vector& aa_grid) { maa_grid = aa_grid; }

  void set_PressureArray(const ArrayOfVector& PressureArray) {
    mPressureArray = PressureArray;
  }

  void set_TemperatureArray(const ArrayOfVector& TemperatureArray) {
    mTemperatureArray = TemperatureArray;
  }

  void set_GasExtinctionArray(const ArrayOfVector& GasExtinctionArray) {
    mGasExtinctionArray = GasExtinctionArray;
  }

  void set_InterpWeightsArray(const ArrayOfMatrix& InterpWeightsArray) {
    mInterpWeightsArray = InterpWeightsArray;
  }

  void set_InterpWeightsAngleArray(
      const ArrayOfMatrix& InterpWeightsAngleArray) {
    mInterpWeightsAngleArray = InterpWeightsAngleArray;
  }

  void set_GposPArray(const ArrayOfArrayOfGridPos& GposPArray) {
    mGposPArray = GposPArray;
  }

  void set_GposLatArray(const ArrayOfArrayOfGridPos& GposLatArray) {
    mGposLatArray = GposLatArray;
  }

  void set_GposLonArray(const ArrayOfArrayOfGridPos& GposLonArray) {
    mGposLonArray = GposLonArray;
  }

  void set_GposZenithArray(const ArrayOfArrayOfGridPos& GposZenithArray) {
    mGposZenithArray = GposZenithArray;
  }

  void set_GposAzimuthArray(const ArrayOfArrayOfGridPos& GposAzimuthArray) {
    mGposAzimuthArray = GposAzimuthArray;
  }

  void set_LStepArray(const ArrayOfVector& LStepArray) {
    mLStepArray = LStepArray;
  }

  void set_MaxLimbIndex(const Index& MaxLimbIndex) {
    mMaxLimbIndex = MaxLimbIndex;
  }

  void set_AtmosphereDim(const Index& AtmosphereDim) {
    mAtmosphereDim = AtmosphereDim;
  }

  void set_CloudboxField(const Tensor6& CloudboxField) {
    mCloudboxField = CloudboxField;
  }

  void set_CloudboxScatteringField(const Tensor6& CloudboxScatteringField) {
    mCloudboxScatteringField = CloudboxScatteringField;
  }

  //get-methods of variables----------------------------------------------------
  const Vector& get_p_grid() const { return mp_grid; }
  Vector& get_p_grid() { return mp_grid; }

  const Vector& get_lat_grid() const { return mlat_grid; }
  Vector& get_lat_grid() { return mlat_grid; }

  const Vector& get_lon_grid() const { return mlon_grid; }
  Vector& get_lon_grid() { return mlon_grid; }

  const Vector& get_za_grid() const { return mza_grid; }
  Vector& get_za_grid() { return mza_grid; }

  const Vector& get_aa_grid() const { return maa_grid; }
  Vector& get_aa_grid() { return maa_grid; }

  const ArrayOfVector& get_PressureArray() const { return mPressureArray; }
  ArrayOfVector& get_PressureArray() { return mPressureArray; }

  const ArrayOfVector& get_TemperatureArray() const {
    return mTemperatureArray;
  }
  ArrayOfVector& get_TemperatureArray() { return mTemperatureArray; }

  const ArrayOfVector& get_GasExtinctionArray() const {
    return mGasExtinctionArray;
  }
  ArrayOfVector& get_GasExtinctionArray() { return mGasExtinctionArray; }

  const ArrayOfMatrix& get_InterpWeightsArray() const {
    return mInterpWeightsArray;
  }
  ArrayOfMatrix& get_InterpWeightsArray() { return mInterpWeightsArray; }

  const ArrayOfMatrix& get_InterpWeightsAngleArray() const {
    return mInterpWeightsAngleArray;
  }
  ArrayOfMatrix& get_InterpWeightsAngleArray() {
    return mInterpWeightsAngleArray;
  }

  const ArrayOfArrayOfGridPos& get_GposPArray() const { return mGposPArray; }
  ArrayOfArrayOfGridPos& get_GposPArray() { return mGposPArray; }

  const ArrayOfArrayOfGridPos& get_GposLatArray() const {
    return mGposLatArray;
  }
  ArrayOfArrayOfGridPos& get_GposLatArray() { return mGposLatArray; }

  const ArrayOfArrayOfGridPos& get_GposLonArray() const {
    return mGposLonArray;
  }
  ArrayOfArrayOfGridPos& get_GposLonArray() { return mGposLonArray; }

  const ArrayOfArrayOfGridPos& get_GposZenithArray() const {
    return mGposZenithArray;
  }
  ArrayOfArrayOfGridPos& get_GposZenithArray() { return mGposZenithArray; }

  const ArrayOfArrayOfGridPos& get_GposAzimuthArray() const {
    return mGposAzimuthArray;
  }
  ArrayOfArrayOfGridPos& get_GposAzimuthArray() { return mGposAzimuthArray; }

  const ArrayOfVector& get_LStepArray() const { return mLStepArray; }
  ArrayOfVector& get_LStepArray() { return mLStepArray; }

  const Index& get_MaxLimbIndex() const { return mMaxLimbIndex; }
  Index& get_MaxLimbIndex() { return mMaxLimbIndex; }

  const Index& get_AtmosphereDim() const { return mAtmosphereDim; }
  Index& get_AtmosphereDim() { return mAtmosphereDim; }

  const Tensor7& get_ExtinctionMatrix() const { return mExtinctionMatrix; }
  Tensor7& get_ExtinctionMatrix() { return mExtinctionMatrix; }

  const Tensor6& get_AbsorptionVector() const { return mAbsorptionVector; }
  Tensor6& get_AbsorptionVector() { return mAbsorptionVector; }

  const Tensor6& get_CloudboxField() const { return mCloudboxField; }
  Tensor6& get_CloudboxField() { return mCloudboxField; }

  const Tensor6& get_CloudboxScatteringField() const {
    return mCloudboxScatteringField;
  }
  Tensor6& get_CloudboxScatteringField() { return mCloudboxScatteringField; }

 protected:
  Vector mp_grid;
  Vector mlat_grid;
  Vector mlon_grid;
  Vector mza_grid;
  Vector maa_grid;
  ArrayOfVector mPressureArray;
  ArrayOfVector mTemperatureArray;
  ArrayOfVector mGasExtinctionArray;
  ArrayOfMatrix mInterpWeightsArray;
  ArrayOfMatrix mInterpWeightsAngleArray;
  ArrayOfArrayOfGridPos mGposPArray;
  ArrayOfArrayOfGridPos mGposLatArray;
  ArrayOfArrayOfGridPos mGposLonArray;
  ArrayOfArrayOfGridPos mGposZenithArray;
  ArrayOfArrayOfGridPos mGposAzimuthArray;
  ArrayOfVector mLStepArray;
  Index mMaxLimbIndex;
  Index mAtmosphereDim;
  Tensor7 mExtinctionMatrix;
  Tensor6 mAbsorptionVector;
  Tensor6 mCloudboxField;
  Tensor6 mCloudboxScatteringField;
};

typedef Array<RTDomain> ArrayOfRTDomain;

class RTDomainScatteringProperties {
 public:
  //TODO: Add doxygen doc and remove unneeded methods and constructors.

  RTDomainScatteringProperties() = default;

  /**Constructor for DomainScatteringProperties
 *
  * @param[in] ScatteringMatrix scattering matrix (Np,Nlat,Nlon,npdir,nidir,nst,nst)
 * @param[in] IdirIdx0 index of flattened incidence zenith angle meshgrid
 * @param[in] IdirIdx1 index of flattened incidence azimuth angle meshgrid
 * @param[in] PdirIdx0 index of flattened propagation zenith angle meshgrid
 * @param[in] PdirIdx1 index of flattened propagation azimuth angle meshgrid
 * @param[in] GpZaI interpolation gridpoints zenith incidence angle
 * @param[in] GpAaI interpolation gridpoints azimuth incidence angle
 * @param[in] Itw interpolation weight for incidence angles
 * @param[in] ScatZaGrid Zenith angle grid for the scattering calculation.
 * @param[in] ScatAaGrid Azimuth angle grid for the scattering calculation.
 */
  RTDomainScatteringProperties(
                               const Tensor7& ScatteringMatrix,
                               const ArrayOfIndex& IdirIdx0,
                               const ArrayOfIndex& IdirIdx1,
                               const ArrayOfIndex& PdirIdx0,
                               const ArrayOfIndex& PdirIdx1,
                               const ArrayOfGridPos& GpZaI,
                               const ArrayOfGridPos& GpAaI,
                               const Tensor3& Itw,
                               const Vector& ScatZaGrid,
                               const Vector& ScatAaGrid)
      : mScatteringMatrix(ScatteringMatrix),
        mIdirIdx0(IdirIdx0),
        mIdirIdx1(IdirIdx1),
        mPdirIdx0(PdirIdx0),
        mPdirIdx1(PdirIdx1),
        mGpZaI(GpZaI),
        mGpAaI(GpAaI),
        mItw(Itw),
        mScatZaGrid(ScatZaGrid),
        mScatAaGrid(ScatAaGrid) {}

  //Set methods
  void set_ScatZaGrid(const Vector& ScatZaGrid) { mScatZaGrid = ScatZaGrid; }

  void set_ScatAaGrid(const Vector& ScatAaGrid) { mScatAaGrid = ScatAaGrid; }

  //Get methods
  const Tensor7& get_ScatteringMatrix() const { return mScatteringMatrix; }
  Tensor7& get_ScatteringMatrix() { return mScatteringMatrix; }

  const ArrayOfIndex& get_IdirIdx0() const { return mIdirIdx0; }
  ArrayOfIndex& get_IdirIdx0() { return mIdirIdx0; }

  const ArrayOfIndex& get_IdirIdx1() const { return mIdirIdx1; }
  ArrayOfIndex& get_IdirIdx1() { return mIdirIdx1; }

  const ArrayOfIndex& get_PdirIdx0() const { return mPdirIdx0; }
  ArrayOfIndex& get_PdirIdx0() { return mPdirIdx0; }

  const ArrayOfIndex& get_PdirIdx1() const { return mPdirIdx1; }
  ArrayOfIndex& get_PdirIdx1() { return mPdirIdx1; }

  const ArrayOfGridPos& get_GpZaI() const { return mGpZaI; }
  ArrayOfGridPos& get_GpZaI() { return mGpZaI; }

  const ArrayOfGridPos& get_GpAaI() const { return mGpAaI; }
  ArrayOfGridPos& get_GpAaI() { return mGpAaI; }

  const Tensor3& get_Itw() const { return mItw; }
  Tensor3& get_Itw() { return mItw; }

  const Vector& get_ScatZaGrid() const { return mScatZaGrid; }
  Vector& get_ScatZaGrid() { return mScatZaGrid; }

  const Vector& get_ScatAaGrid() const { return mScatAaGrid; }
  Vector& get_ScatAaGrid() { return mScatAaGrid; }

  /** Calculates the scattered field
  *
  * @param[out] cloudbox_scat_field Monochromatic scattered radiation field inside the
  *              cloudbox at iteration step.
  * @param[in] cloudbox_field_mono Monochromatic radiation field inside the
  *              cloudbox.
  * @param[in] atmosphere_dim The atmospheric dimensionality (1-3).
  * @param[in] verbosity Verbosity setting.
  */
  void CalcScatteredField(Tensor6& cloudbox_scat_field,
                          //WS Input:
                          const Tensor6& cloudbox_field_mono,
                          const Index& atmosphere_dim,
                          const Verbosity& verbosity);

 private:
  /** Calculates the scattered field for 1D atmosphere
  *
  * @param[out] cloudbox_scat_field Monochromatic scattered radiation field inside the
  *              cloudbox at iteration step.
  * @param[in] cloudbox_field_mono Monochromatic radiation field inside the
  *              cloudbox.
  * @param[in] verbosity Verbosity setting.
  */
  void CalcScatteredField1D(
      // Output
      Tensor6& cloudbox_scat_field,
      // Input:
      const ConstTensor6View& cloudbox_field_mono,
      const Verbosity& verbosity);

  /** Calculates the scattered field for 3D atmosphere
  *
  * @param[out] cloudbox_scat_field Monochromatic scattered radiation field inside the
  *              cloudbox at iteration step.
  * @param[in] cloudbox_field_mono Monochromatic radiation field inside the
  *              cloudbox.
  * @param[in] verbosity Verbosity setting.
  */
  void CalcScatteredField3D(
      // Output
      Tensor6& cloudbox_scat_field,
      // Input:
      const ConstTensor6View& cloudbox_field_mono,
      const Verbosity& verbosity);

 protected:
  Tensor7 mScatteringMatrix;
  ArrayOfIndex mIdirIdx0;
  ArrayOfIndex mIdirIdx1;
  ArrayOfIndex mPdirIdx0;
  ArrayOfIndex mPdirIdx1;
  ArrayOfGridPos mGpZaI;
  ArrayOfGridPos mGpAaI;
  Tensor3 mItw;
  Vector mScatZaGrid;
  Vector mScatAaGrid;
};

typedef Array<RTDomainScatteringProperties> ArrayOfRTDomainScatteringProperties;

/** Initialises variables for DOIT scattering calculations.
 *
 * Note that cloudbox_field is Nan-initialzed
 *
 * @param[out] cloudbox_field Spectral radiance field inside the cloudbox
 * @param[in] stokes_dim The dimensionality of the Stokes vector (1-4).
 * @param[in] atmosphere_dim The atmospheric dimensionality (1-3).
 * @param[in] f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 * @param[in] rt_za_grid Zenith angle grid.
 * @param[in] rt_aa_grid Azimuth angle grid.
 * @param[in] cloudbox_limits The limits of the cloud box.
 */
void Initialize_cloudbox_field(
    //Output
    Tensor7& cloudbox_field,
    //Input
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Vector& f_grid,
    const Vector& rt_za_grid,
    const Vector& rt_aa_grid,
    const ArrayOfIndex& cloudbox_limits);

/** Set the angular grids for DOIT scattering calculations.
 *
 * @param[out] za_grid Zenith angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[out] aa_grid Azimuth angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[out] scat_za_grid Zenith angle grid for the scattering calculation.
 * @param[out] scat_aa_grid Azimuth angle grid for the scattering calculation.
 * @param[in] N_za_grid Number of zenith angle grid points.
 * @param[in] N_aa_grid Number of azimuth angle grid points.
 * @param[in] N_scat_za_grid Number of zenith angle grid points for the
 *              scattering calculation.
 * @param[in] N_scat_aa_grid Number of azimuth angle grid points for the
 *              scattering calculation.
 * @param[in] za_grid_type String with the type of the grid for za_grid.
 */
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

/** Interpolate clearsky field on all gridpoints in cloudbox.
 *
 * @param[out] cloudbox_field Spectral radiance field inside the cloudbox.
 * @param[in] f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 * @param[in] p_grid Pressure grid.
 * @param[in] lat_grid Latitude grid.
 * @param[in] lon_grid Longitude grid.
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] atmosphere_dim The atmospheric dimensionality (1-3).
 * @param[in] verbosity Verbosity setting.
 */
void SetClearsky_cloudbox(
    //Output
    Tensor7& cloudbox_field,
    //Input
    const Vector& f_grid,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const ArrayOfIndex& cloudbox_limits,
    const Index& atmosphere_dim,
    const Verbosity& verbosity);

/** Calculates incoming radiation field of the cloudbox.
 *
 * The method performs monochromatic pencil beam calculations for
 * all grid positions on the cloudbox boundary, and all directions
 * given by angle grids (*za/aa_grid*). Found radiances
 * are stored in cloudbox_field which can be used as boundary
 * conditions when scattering inside the cloud box is solved.
 *
 * @param[in,out] ws Current workspace.
 * @param[out] cloudbox_field Spectral radiance field inside the cloudbox.
 * @param[in] iy_main_agenda Agenda calculating the single monochromatic pencil
 *              beam spectrum.
 * @param[in] atmosphere_dim The atmospheric dimensionality (1-3).
 * @param[in] lat_grid Latitude grid.
 * @param[in] lon_grid Longitude grid.
 * @param[in] z_field Field of geometrical altitudes.
 * @param[in] nlte_field The field of NLTE temperatures and/or ratios.
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 * @param[in] za_grid Zenith angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[in] aa_grid Azimuth angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[in] verbosity Verbosity setting.
 */
void GetIncomingRadiation(Workspace& ws,
                          //Output
                          Tensor7& cloudbox_field,
                          //Input
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
                          const Verbosity& verbosity);

/** Limits the atmospheric input fields to the cloudbox.
 *
 * @param[out] p_grid_cldbx Pressure grid within the cloudbox.
 * @param[out] lat_grid_cldbx Latitude grid within the cloudbox.
 * @param[out] lon_grid_cldbx Longitude grid within the cloudbox.
 * @param[out] t_field_cldbx Temperature field within the cloudbox.
 * @param[out] z_field_cldbx Field of geometrical altitudes within the cloudbox.
 * @param[out] vmr_field_cldbx Field of volume mixing ratios within the cloudbox.
 * @param[in] p_grid Pressure grid.
 * @param[in] lat_grid Latitude grid.
 * @param[in] lon_grid Longitude grid.
 * @param[in] t_field Temperature field
 * @param[in] z_field Field of geometrical altitudes.
 * @param[in] vmr_field Field of volume mixing ratios.
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] verbosity Verbosity setting.
 */
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

/** Main function of the DOIT solver.
 *
 * Iterative solution of the VRTE (DOIT method). A solution for the RTE
 * with scattering is found using the DOIT method:
 * 1. Calculate scattering integral.
 * 2. Calculate RT with fixed scattered field using
 * 3. Convergence check
 *
 * Note: The atmospheric dimensionality atmosphere_dim can be
 * either 1 or 3. To these dimensions the method adapts
 * automatically. 2D scattering calculations are not
 * supported.
 *
 * @param[in,out] ws Current workspace.
 * @param[out] cloudbox_field_mono Monochromatic radiation field inside the
 *              cloudbox.
 * @param[out] gas_extinction Field with the gas extinction, this is not used
 *              for the actual RT calculation. It is just used for the adaptive
 *              ppath length.
 * @param[out] extinction_matrix Extinction matrix field
 *              (Np,Nlat,Nlon,ndir,nst,nst).
 * @param[out] absorption_vector absorption vector field (Np,Nlat,Nlon,ndir,nst)
 * @param[out]scattering_matrix Scattering matrix field
 *              (Np,Nlat,Nlon,pdir,idir,nst,nst).
 * @param[out] convergence_flag Flag for the convergence test.
 * @param[out] iteration_counter  Counter for number of iterations.
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] propmat_clearsky_agenda Agenda to calculate the absorption
 *              coefficient matrix.
 * @param[in] surface_rtprop_agenda Agenda to calculate radiative properties of
 *              the surface.
 * @param[in] ppath_step_agenda Agenda to calculate a propagation path step.
 * @param[in] atmosphere_dim The atmospheric dimensionality (1 or 3).
 * @param[in] stokes_dim The dimensionality of the Stokes vector (1-4).
 * @param[in] pnd_field Particle number density field.
 * @param[in] t_field Temperature field.
 * @param[in] z_field Field of geometrical altitudes.
 * @param[in] vmr_field Field of volume mixing ratios.
 * @param[in] z_surface Surface altitude field.
 * @param[in] p_grid Pressure grid.
 * @param[in] lat_grid Latitude grid.
 * @param[in] lon_grid Longitude grid.
 * @param[in] za_grid Zenith angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[in] aa_grid Azimuth angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[in] scat_za_grid Zenith angle grid for the scattering calculation.
 * @param[in] scat_aa_grid Azimuth angle grid for the scattering calculation.
 * @param[in] f_mono Monochromatic frequency.
 * @param[in] scat_data Array of single scattering data.
 * @param[in] t_interp_order Temperature interpolation order.
 * @param[in] iy_unit Unit in which the convergence is checked.
 * @param[in] refellipsoid Reference ellipsoid.
 * @param[in] epsilon Limits for convergence. A vector with length matching
 *              stokes_dim with unit of iy_unit.
 * @param[in] max_num_iterations Maximum number of iterations.
 * @param[in] tau_min Maximum optical thickness per propagation step.
 * @param[in] accelerated Index wether to accelerate only the intensity (1) or
 *              the whole Stokes Vector(>1).
 * @param[in] ForwardCorrectionFlag Index wether to use the forward scattering
 *              corretion (1) or not (0).
 * @param[in] ppath_lmax Maximum length between points describing propagation
 *              paths.
 * @param[in]ppath_lraytrace Maximum length of ray tracing steps when
 *              determining propagation paths.
 * @param[in] verbosity Verbosity setting.
 */
void NewDoitMonoCalc(Workspace& ws,
                     //Output
                     Tensor6& cloudbox_field_mono,
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
                     const Numeric& refinement_criterion,
                     const Numeric& tau_min,
                     const Index& N_decade,
                     const Index& accelerated,
                     const Index& ForwardCorrectionFlag,
                     const Numeric& ppath_lmax,
                     const Numeric& ppath_lraytrace,
                     const Verbosity& verbosity);

/** Calculates the gas extinction for an atmospheric column
 *
 *
 * @param[in,out] ws Current workspace.
 * @param[out] gas_extinct Column with the gas extinction.
 * @param[in] p_grid Pressure grid.
 * @param[in] t_vector Temperature field
 * @param[in] vmr_matrix Column of n volume mixing ratios.
 * @param[in] propmat_clearsky_agenda Agenda to calculate the absorption
 *              coefficient matrix.
 * @param[in] f_mono Monochromatic frequency.
 */
void CalcGasExtinction(Workspace& ws,
                       Vector& gas_extinct,
                       const ConstVectorView& p_grid,
                       const ConstVectorView& t_vector,
                       const ConstMatrixView& vmr_matrix,
                       const Agenda& propmat_clearsky_agenda,
                       const ConstVectorView& f_mono);

/** Calculates gas extinction for a field.
 *
 * @param[in,out] ws Current workspace.
 * @param[out] gas_extinct Field with the gas extinction, this is not used
 *              for the actual RT calculation. It is just used for the adaptive
 *              ppath length.
 * @param[in] p_grid Pressure grid.
 * @param[in] lat_grid Latitude grid.
 * @param[in] lon_grid Longitude grid.
 * @param[in] t_field Temperature field.
 * @param[in] vmr_field Field of volume mixing ratios.
 * @param[in] propmat_clearsky_agenda Agenda to calculate the absorption
 *              coefficient matrix.
 * @param[in] f_mono Monochromatic frequency.
 */
void CalcGasExtinctionField(Workspace& ws,
                            Tensor3& gas_extinct,
                            const ConstVectorView& p_grid,
                            const ConstVectorView& lat_grid,
                            const ConstVectorView& lon_grid,
                            const ConstTensor3View& t_field,
                            const ConstTensor4View& vmr_field,
                            const Agenda& propmat_clearsky_agenda,
                            const ConstVectorView& f_mono);

//TODO: Add doxygen doc
/**
 * @param DomainScatteringProperties
 * @param t_field
 * @param scat_data
 * @param pnd_field
 * @param stokes_dim
 * @param atmosphere_dim
 * @param za_grid
 * @param aa_grid
 * @param scat_za_grid
 * @param scat_aa_grid
 * @param t_interp_order
 * @param ForwardCorrectionFlag
 * @param verbosity
 */
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
    const Verbosity& verbosity);

/** Calculates bulk optical properties from per-scat-species bulk properties.
 *
 * @param[out] extinction_matrix Bulk extinction matrix (Np,Nlat,Nlon,Nza,Naa,nst,nst)
 * @param[out] absorption_vector absorption vector (Np,Nlat,Nlon,Nza,Naa,nst)
 * @param[in] scat_data Array of single scattering data.
 * @param[in] za_grid Zenith angle grid of Spectral radiance field inside the
 *             cloudbox.
 * @param[in] pnd_field Particle number density field.
 * @param[in] t_field Temperature field.
 * @param[in] stokes_dim The dimensionality of the Stokes vector (1-4).
 */
void CalcParticleOpticalProperties(
    //Output
    Tensor7& extinction_matrix,
    Tensor6& absorption_vector,
    //input
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Vector& za_grid,
    const ConstTensor4View& pnd_field,
    const ConstTensor3View& t_field,
    const Index& stokes_dim);

/**Calculates bulk scattering properties from per-scat-species bulk properties.
 *
 * @param[out] scattering_matrix scattering matrix (Np,Nlat,Nlon,npdir,nidir,nst,nst)
 * @param[out] idir_idx0 index of flattened incidence zenith angle meshgrid
 * @param[out] idir_idx1 index of flattened incidence azimuth angle meshgrid
 * @param[out] pdir_idx0 index of flattened propagation zenith angle meshgrid
 * @param[out] pdir_idx1 index of flattened propagation azimuth angle meshgrid
 * @param[out] gp_za_i interpolation gridpoints zenith incidence angle
 * @param[out] gp_aa_i interpolation gridpoints azimuth incidence angle
 * @param[out] itw interpolation weight for incidence angles
 * @param[in] t_field Temperature field.
 * @param[in] scat_data Array of single scattering data.
 * @param[in] pnd_field Particle number density field.
 * @param[in] stokes_dim The dimensionality of the Stokes vector (1-4).
 * @param[in] atmosphere_dim The atmospheric dimensionality (1-3).
 * @param[in] za_grid Zenith angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[in] aa_grid Azimuth angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[in] scat_za_grid Zenith angle grid for the scattering calculation.
 * @param[in] scat_aa_grid Azimuth angle grid for the scattering calculation.
 * @param[in] t_interp_order Temperature interpolation order.
 * @param[in] verbosity Verbosity setting.
*/
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


/** Applies forward scattering correction
 *
 * @param scattering_matrix
 * @param extinction_matrix
 * @param absorption_vector
 * @param idir_idx0
 * @param idir_idx1
 * @param pdir_idx0
 * @param pdir_idx1
 * @param atmosphere_dim
 * @param za_grid
 * @param aa_grid
 * @param scat_za_grid
 * @param scat_aa_grid
 * @param f_mono
 */
void ForwardScatteringCorrection(  //Output
    Tensor7& scattering_matrix,
    Tensor7& extinction_matrix,        //(Np,Nlat,Nlon,Nza,Naa,nst,nst)
    Tensor6& absorption_vector,  //(Np,Nlat,Nlon,Nza,Naa,nst)
    const ArrayOfIndex& idir_idx0,
    const ArrayOfIndex& idir_idx1,
    const ArrayOfIndex& pdir_idx0,
    const ArrayOfIndex& pdir_idx1,
    const Index& atmosphere_dim,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Vector& scat_za_grid,
    const Vector& scat_aa_grid);


/**Calculates surface properties
 *
 * @param ws[in,out] Current workspace.
 * @param[out] surface_skin_t Surface skin temperature
 * @param[out] surface_los Downwelling radiation directions to consider in surface
 *              reflection.
 * @param[out] surface_reflection_matrix The reflection coefficients for the directions
 *              given by surface_los to the direction of interest.
 * @param[out] surface_emission The emission from the surface.
 * @param[in] surface_rtprop_agenda Agenda providing radiative properties of the surface.
 * @param[in] f_grid f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 * @param[in] za_grid Zenith angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[in] aa_grid Azimuth angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[in] lat_grid Latitude grid.
 * @param[in] lon_grid Longitude grid.
 * @param[in] atmosphere_dim The atmospheric dimensionality (1-3).
 * @param[in] stokes_dim The dimensionality of the Stokes vector (1-4).
 * @param[in] surface_field Surface altitude field.
 */
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

/**Calculates maximum propagation path length according to tau_max
 *
 * @param[out] p_path_maxlength maximum length of propagation path
 * @param[in] extinction_matrix extinction_matrix Bulk extinction matrix
 *              (Np,Nlat,Nlon,ndir,nst,nst)
 * @param[in] gas_extinct Field with the gas extinction, this is not used
 *              for the actual RT calculation. It is just used for the adaptive
 *              ppath length.
 * @param[in] p_grid Pressure grid.
 * @param[in] lat_grid Latitude grid.
 * @param[in] lon_grid Longitude grid.
 * @param[in] scat_za_grid Zenith angle grid for the scattering calculation.
 * @param[in] ppath_lmax Maximum length between points describing propagation
 *              paths outside of cloudbox
 * @param[in] ppath_lraytrace Maximum length of ray tracing steps when determining
 *              propagation paths otside of cloudbox
 * @param[in] tau_max Maximum optical thickness per propagation step.
 */
void CalcPropagationPathMaxLength(
    Tensor3& p_path_maxlength,
    const Tensor7View& extinction_matrix,  //(Np,Nlat,Nlon,ndir,nst,nst)
    const ConstTensor3View& gas_extinct,
    const ConstVectorView& p_grid,
    const ConstVectorView& lat_grid,
    const ConstVectorView& lon_grid,
    const Vector& scat_za_grid,
    const Numeric& ppath_lmax,
    const Numeric& ppath_lraytrace,
    const Numeric& tau_max);

/**Checks where refinement is needed.
 *
 * @param[out] refine Flag which indicate which grid CELL needs refinement.
 * @param[out] DlogK Difference in log10 of the mean extinction in 1 to 3 directions.
 * @param[out] tau_p mean optical thickness in in 1 to 3 directions.
 * @param[in] extinction_matrix extinction_matrix extinction_matrix Bulk extinction matrix
 *              (Np,Nlat,Nlon,ndir,nst,nst)
 * @param[in] gas_extinct Field with the gas extinction, this is not used
 *              for the actual RT calculation. It is just used for the adaptive
 *              ppath length.
 * @param[in] p_grid Pressure grid.
 * @param[in] lat_grid Latitude grid.
 * @param[in] lon_grid Longitude grid.
 * @param[in] scat_za_grid Zenith angle grid for the scattering calculation.
 * @param[in] z_field Field of geometrical altitudes.
 * @param[in] refellipsoid Reference ellipsoid.
 * @param[in] atmosphere_dim The atmospheric dimensionality (1 or 3).
 * @param[in] refine_crit Refinement criterium in terms of Log10 of ratio of two adjacent
 *            points of the total extinction.
 * @param[in] tau_crit Minimum optical thickness for refinement. Cells with smaller optical
 *            thickness will not be refined.
 */
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
    const Numeric& tau_crit);

/**Does the refinement and prepares the subdomains.
 *
 * @param[in,out] ws Current workspace.
 * @param[out] Subdomains Array of sub radiative transfer domain
 * @param[out] SubdomainsScatteringProperties Array of subdomain scattering properties
 * @param[in] MainDomain Radiative transfer domain
 * @param[in] MainDomainScatteringProperties Main domain scattering properties.
 * @param[in] t_field Temperature field.
 * @param[in] z_field Altitude field.
 * @param[in] vmr_field Field of volume mixing ratios.
 * @param[in] pnd_field Particle number density field.
 * @param[in] scat_data Array of single scattering data.
 * @param[in] t_interp_order Temperature interpolation order.
 * @param[in] ForwardCorrectionFlag Index wether to use the forward scattering
 *              corretion (1) or not (0).
 * @param[in] refine Flag which indicate which grid CELL will be refinement.
 * @param[in] DlogK Difference in log10 of the mean extinction in 1 to 3 directions.
 * @param[in] N_decade Number of points per decade of change of extinction to be inserted.
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] propmat_clearsky_agenda Agenda to calculate the absorption
 *              coefficient matrix.
 * @param[in] ppath_step_agenda Agenda to calculate a propagation path step.
 * @param[in] ppath_lmax Maximum length between points describing propagation
 *              paths.
 * @param[in]ppath_lraytrace Maximum length of ray tracing steps when
 *              determining propagation paths.
 * @param[in] z_surface Surface altitude.
 * @param[in] refellipsoid Reference ellipsoid.
 * @param[in] f_grid f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 */
void PrepareSubdomains(
    Workspace& ws,
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
    const ArrayOfIndex& cloudbox_limits,
    const Agenda& propmat_clearsky_agenda,
    const Agenda& ppath_step_agenda,
    const Numeric& ppath_lmax,
    const Numeric& ppath_lraytrace,
    const Matrix& z_surface,
    const Vector& refellipsoid,
    const Vector& f_grid,
    const Verbosity& verbosity);

/**Does the refinement and prepares the subdomains for 1D atmospheres.
 *
 * @param[in,out] ws Current workspace.
 * @param[out] Subdomains Array of sub radiative transfer domain
 * @param[out] SubdomainsScatteringProperties Array of subdomain scattering properties
 * @param[in] MainDomain Radiative transfer domain
 * @param[in] MainDomainScatteringProperties Main domain scattering properties.
 * @param[in] t_field Temperature field.
 * @param[in] z_field Altitude field.
 * @param[in] vmr_field Field of volume mixing ratios.
 * @param[in] pnd_field Particle number density field.
 * @param[in] scat_data Array of single scattering data.
 * @param[in] t_interp_order Temperature interpolation order.
 * @param[in] ForwardCorrectionFlag Index wether to use the forward scattering
 *              corretion (1) or not (0).
 * @param[in] refine Flag which indicate which grid CELL will be refinement.
 * @param[in] DlogK Difference in log10 of the mean extinction in 1 to 3 directions.
 * @param[in] N_decade Number of points per decade of change of extinction to be inserted.
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] propmat_clearsky_agenda Agenda to calculate the absorption
 *              coefficient matrix.
 * @param[in] ppath_step_agenda Agenda to calculate a propagation path step.
 * @param[in] ppath_lmax Maximum length between points describing propagation
 *              paths.
 * @param[in]ppath_lraytrace Maximum length of ray tracing steps when
 *              determining propagation paths.
 * @param[in] z_surface Surface altitude.
 * @param[in] refellipsoid Reference ellipsoid.
 * @param[in] f_grid f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 */
void PrepareSubdomains1D(
    Workspace& ws,
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
    const ArrayOfIndex& cloudbox_limits,
    const Agenda& propmat_clearsky_agenda,
    const Agenda& ppath_step_agenda,
    const Numeric& ppath_lmax,
    const Numeric& ppath_lraytrace,
    const Matrix& z_surface,
    const Vector& refellipsoid,
    const Vector& f_grid,
    const Verbosity& verbosity);

/** The actual DOIT scattering solver
 *
 * @param[out] cloudbox_field_mono Monochromatic radiation field inside the
 *              cloudbox.
 * @param[out] convergence_flag Flag for the convergence test.
 * @param[out] iteration_counter  Counter for number of iterations.
 * @param[in] extinction_matrix Extinction matrix field
 *              (Np,Nlat,Nlon,ndir,nst,nst).
 * @param[in] absorption_vector absorption vector field (Np,Nlat,Nlon,ndir,nst)
 * @param[in] scattering_matrix Scattering matrix field
 *              (Np,Nlat,Nlon,pdir,idir,nst,nst).
 * @param[in] surface_reflection_matrix The reflection coefficients for the directions
 *              given by surface_los to the direction of interest.
 * @param[in] surface_emission The emission from the surface.
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] atmosphere_dim The atmospheric dimensionality (1-3).
 * @param[in] scat_za_grid Zenith angle grid for the scattering calculation.
 * @param[in] scat_aa_grid Azimuth angle grid for the scattering calculation.
 * @param[in] PressureArray Array with the pressures of each propagation step
 * @param[in] TemperatureArray Array with the temperature at each propagation step
 * @param[in] GasExtinctionArray Array with the gas extinction at each propagation step
 * @param[in] InterpWeightsArray Array with the pressure interpolation weights
 *              at each propagation step
 * @param[in] InterpWeightsZenithArray Array with the zenith angle interpolation
 *              weights at each propagation step
 * @param[in] GposArray Array with the pressure gridpoints at each propagation step
 * @param[in] GposZenithArray Array with the zenith angle gridpoints at each
 *              propagation step
 * @param[in] LstepArray Array with the length of each propagation step
 * @param[in] MaxLimbIndex Index with index of the highest limb angle (-1 if no
 *              limb angle)
 * @param[in] idir_idx0 index of flattened incidence zenith angle meshgrid
 * @param[in] idir_idx1 index of flattened incidence azimuth angle meshgrid
 * @param[in] pdir_idx0 index of flattened propagation zenith angle meshgrid
 * @param[in] pdir_idx1 index of flattened propagation azimuth angle meshgrid
 * @param[in] gp_za_i interpolation gridpoints zenith incidence angle
 * @param[in] gp_aa_i interpolation gridpoints azimuth incidence angle
 * @param[in] itw interpolation weight for incidence angles
 * @param[in] f_mono Monochromatic frequency.
 * @param[in] iy_unit Unit in which the convergence is checked.
 * @param[in] epsilon Limits for convergence. A vector with length matching
 *              stokes_dim with unit of iy_unit.
 * @param[in] max_num_iterations Maximum number of iterations.
 * @param[in] tau_max Maximum optical thickness per propagation step.
 * @param[in] accelerated Index wether to accelerate only the intensity (1) or
 *              the whole Stokes Vector(>1).
 * @param[in] verbosity Verbosity setting.
 */
void RunNewDoit(//Input and Output:
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
    const Verbosity& verbosity);

/** RT calculation at iteration step
 *
 * @param[in,out] Domain Radiative transfer domain
 * @param[in] surface_reflection_matrix The reflection coefficients for the directions
 *              given by surface_los to the direction of interest.
 * @param[in] surface_emission The emission from the surface.
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] atmosphere_dim The atmospheric dimensionality (1-3).
 * @param[in] f_grid f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 * @param[in] verbosity Verbosity setting.
 */
void UpdateSpectralRadianceField(//Input and Output:
    RTDomain& Domain,
    //Input:
    const ConstTensor6View& surface_reflection_matrix,
    const ConstTensor5View& surface_emission,
    const ArrayOfIndex& cloudbox_limits,
    const Index& atmosphere_dim,
    const Vector& f_grid,
    const Verbosity& verbosity);

/** RT calculation for 1D atmosphere at iteration step
 *
 * @param[in,out] Domain Radiative transfer domain
 * @param[in] surface_reflection_matrix The reflection coefficients for the directions
 *              given by surface_los to the direction of interest.
 * @param[in] surface_emission The emission from the surface.
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] f_grid f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 * @param[in] verbosity Verbosity setting.
 */
void UpdateSpectralRadianceField1D(
    //Input and Output:
    RTDomain& Domain,
    const ConstTensor6View& surface_reflection_matrix,
    const ConstTensor5View& surface_emission,
    const ArrayOfIndex& cloudbox_limits,
    const Vector& f_grid,
    const Verbosity& verbosity);

//TODO:Add doxygen doc
void UpdateSpectralRadianceField3D(
    Workspace& ws,
    // WS Input and Output:
    Tensor6& cloudbox_field_mono,
    Tensor6& cloudbox_scat_field,
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

/** Calculates the radiation at propagation path for a grid point and direction in 1D cloudbox
 *
 * @param[in] cloudbox_field_mono Monochromatic radiation field inside the
 *              cloudbox at iteration step.
 * @param[in] p_index Index of pressure gridpoint
 * @param[in] za_index index of zenith angle gridpoint
 * @param[in] idx Index of ppath
 * @param[in] f_grid f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 * @param[in] surface_reflection_matrix The reflection coefficients for the directions
 *              given by surface_los to the direction of interest.
 * @param[in] surface_emission The emission from the surface.
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] verbosity Verbosity setting.
 */
void UpdateCloudPropagationPath1D(
    RTDomain& Domain,
    const Index& p_index,
    const Index& za_index,
    const Index& idx,
    const ConstVectorView& f_grid,
    const ConstTensor6View& surface_reflection_matrix,
    const ConstTensor5View& surface_emission,
    const ArrayOfIndex& cloudbox_limits,
    const Verbosity& verbosity);

/** Interpolates fields on propagation path
 *
 * @param[out] ext_mat_int Extinction matrix field on propagation path
 * @param[out] abs_vec_int Absorption vector field on propagation path
 * @param[out] sca_vec_int Scattered field on propagation path
 * @param[out] cloudbox_field_mono_int Monochromatic radiation field on propagation
 *              path
 * @param[in] ext_mat_field Extinction matrix field for given direction
 *              (Np,Nlat,Nlon,nst,nst).
 * @param[in] abs_vec_field Absorption vector field for given direction
 *              (Np,Nlat,Nlon,ndir,nst)
 * @param[in] cloudbox_scat_field Monochromatic scattered radiation field
 *              inside the cloudbox at iteration step.
 * @param[in] cloudbox_field_mono Monochromatic radiation field inside the
 *              cloudbox at iteration step.
 * @param[in] cloud_gp_p grid pos for pressure interpolation on propagation path
 * @param[in] cloud_gp_za grid pos for zenith angle interpolation on propagation path
 * @param[in] itw pressure interpolation weights
 * @param[in] itw_za zenith angle interpolation weights
 * @param[in] verbosity verbosity Verbosity setting.
 */
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
    const Verbosity& verbosity);

/** Calculates the radiation transport along propagation path
 *
 * @param[out] cloudbox_field_mono Monochromatic radiation field inside the
 *              cloudbox at iteration step.
 * @param[in] lstep_ppath Length of ppath steps
 * @param[in] temperature_ppath Temperature at ppath steps
 * @param[in] pressure_ppath Pressure at ppath steps
 * @param[in] gas_extinction_ppath Gas extinction at ppath steps
 * @param[in] ext_mat_int Extinction matrix field on propagation path
 * @param[in] abs_vec_int Absorption vector field on propagation path
 * @param[in] sca_vec_int Scattered field on propagation path
 * @param[in] cloudbox_field_mono_int Monochromatic radiation field on propagation
 *              path
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] f_grid f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 * @param[in] p_index Index of pressure gridpoint
 * @param[in] lat_index Index of latitude gridpoint
 * @param[in] lon_index Index of longitude gridpoint
 * @param[in] za_index index of zenith angle gridpoint
 * @param[in] aa_index index of azimuth angle gridpoint
 * @param[in] verbosity verbosity Verbosity setting.
 */
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
                               const Verbosity& verbosity);

/** Calculates the radiation transport for one propagation path
 *
 * The function can be used for cloudbox calculations.
 *
 * The function is best explained by considering a homogenous layer. That is,
 * the physical conditions inside the layer are constant. In reality they
 * are not constant, so in practical all coefficients have to be averaged
 * before calling this function.
 * Total extinction and absorption inside the layer are described by
 * *ext_mat_av* and *abs_vec_av* respectively,
 * the blackbdody radiation of the layer is given by *rte_planck_value*
 * and the propagation path length through the layer is *lstep*.
 *
 * There is an additional scattering source term in the
 * VRTE, the scattering integral term. For this function a constant
 * scattering term is assumed. The radiative transfer step is only a part
 * the iterative solution of the scattering problem, for more
 * information consider AUG. In the clearsky case this variable has to be
 * set to 0.
 *
 * When calling the function, the vector *stokes_vec* shall contain the
 * Stokes vector for the incoming radiation. The function returns this
 * vector, then containing the outgoing radiation on the other side of the
 * layer.
 *
 * The function performs the calculations differently depending on the
 * conditions to improve the speed. There are three cases: <br>
 *  1. Scalar absorption (stokes_dim = 1). <br>
 *  2. The matrix ext_mat_gas is diagonal (unpolarised absorption). <br>
 *  3. The total general case.
 *
 * @param[in,out] stokes_vec Stokes vector of propagation step.
 * @param[in,out] trans_mat Transmission matrix of propagation step.
 * @param[in] ext_mat_av Averaged extinction matrix of propagation step.
 * @param[in] abs_vec_av Averaged absorption vector.
 * @param[in] sca_vec_av averaged scattering vector.
 * @param[in] lstep The length of propagation step.
 * @param[in] rtp_planck_value Blackbody radiation.
 * @param[in] trans_is_precalc FIXMEDOC
 *
 * \author Richard Larsson, Manfred Brath
    \date   2017-08-14, 2020-02-10

 */
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

//TODO add doxygen doc
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
    const Verbosity& verbosity);

//TODO add doxygen doc
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
    const Verbosity& verbosity);

//TODO add doxygen doc
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
    const Verbosity& verbosity);

//TODO add doxygen doc
void CalcSourceForRTStep(  //Output:
    VectorView source,
    VectorView j, //emission density
    //Input
    const PropagationMatrix& ext_mat,
    const StokesVector& abs_vec,
    const ConstVectorView& sca_vec,
    const Numeric& rtp_planck_value);

//TODO add doxygen doc
void RTSolver2ndOrder(  //in and out:
    VectorView stokes_vec,
    ConstVectorView& source_0,
    ConstVectorView& source_1,
    ConstMatrixView& trans_mat_01,
    const Numeric& tau_01,
    ConstVectorView& Qmax);

//TODO add doxygen doc
void RTSolver3rdOrder(//in and out:
    VectorView stokes_vec,
    ConstVectorView& source_0,
    ConstVectorView& source_1,
    ConstVectorView& source_2,
    ConstMatrixView& trans_mat_01,
    const Numeric& tau_01,
    const Numeric& tau_12,
    ConstVectorView& Qmax);

/** Checks the convergence of the DOIT iteration
 *
 * @param[out] convergence_flag Flag for the convergence test.
 * @param[out] iteration_counter  Counter for number of iterations.
 * @param[in] cloudbox_field_mono Monochromatic radiation field inside the
 *              cloudbox of the current iteration step.
 * @param[in] cloudbox_field_mono_old Monochromatic radiation field inside the
 *              cloudbox of the previous iteration step.
 * @param[in] f_mono Monochromatic frequency.
 * @param[in] epsilon Limits for convergence. A vector with length matching
 *              stokes_dim with unit of iy_unit.
 * @param[in] max_iterations Maximum number of iterations.
 * @param[in] iy_unit Unit in which the convergence is checked.
 * @param[in] verbosity Verbosity setting.
 */
void CheckConvergence(  //WS Input and Output:
    Index& convergence_flag,
    Index& iteration_counter,
    Tensor6&  cloudbox_field_mono,
    const Tensor6& cloudbox_field_mono_old,
    const Numeric& f_mono,
    const Vector& epsilon,
    const Index& max_iterations,
    const String& iy_unit,
    const Verbosity& verbosity);

/** Precalculates the propagation path steps and quantities needed for interpolation
 *
 * @param[in,out] ws Current workspace.
 * @param[out] PressureArray Array with the pressures of each propagation step
 * @param[out] TemperatureArray Array with the temperature at each propagation step
 * @param[out] VmrArray Array with the VMR at each propagation step
 * @param[out] InterpWeightsArray Array with the pressure interpolation weights
 *              at each propagation step
 * @param[out] InterpWeightsZenithArray Array with the zenith angle interpolation
 *              weights at each propagation step
 * @param[out] GposArray Array with the pressure interpolation gridpoints at
 *              each propagation step
 * @param[out] GposZenithArray Array with the zenith angle interpolatrion
 *              gridpoints at each propagation step
 * @param[out] LstepArray Array with the length of each propagation step
 * @param[out] MaxLimbIndex Index with index of the highest limb angle (-1 if no
 *              limb angle)
 * @param[in] cloudbox_limits The limits of the cloud box.
 * @param[in] p_index Index of pressure gridpoint
 * @param[in] za_grid Zenith angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[in] ppath_step_agenda Agenda to calculate a propagation path step.
 * @param[in] ppath_lmax Maximum length between points describing propagation
 *              paths.
 * @param[in]ppath_lraytrace Maximum length of ray tracing steps when
 *              determining propagation paths.
 * @param[in] p_grid Pressure grid.
 * @param[in] lat_grid Latitude grid.
 * @param[in] lon_grid Longitude grid.
 * @param[in] z_field Field of geometrical altitudes.
 * @param[in] t_field Temperature field.
 * @param[in] vmr_field Field of volume mixing ratios.
 * @param[in] z_surface Surface altitude.
 * @param[in] refellipsoid Reference ellipsoid.
 * @param[in] f_grid f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 * @param[in] verbosity Verbosity setting.
 */
void EstimatePPathElements1D(Workspace& ws,
                             ArrayOfVector& PressureArray,
                             ArrayOfVector& TemperatureArray,
                             ArrayOfMatrix& VmrArray,
                             ArrayOfMatrix& InterpWeightsArray,
                             ArrayOfMatrix& InterpWeightsZenithArray,
                             ArrayOfArrayOfGridPos& GposArray,
                             ArrayOfArrayOfGridPos& GposZenithArray,
                             ArrayOfVector& LstepArray,
                             Index& MaxLimbIndex,
                             const Vector& za_grid,
                             const Agenda& ppath_step_agenda,
                             const Numeric& ppath_lmax,
                             const Numeric& ppath_lraytrace,
                             const Vector& p_grid,
                             const Vector& lat_grid,
                             const Vector& lon_grid,
                             const Tensor3& z_field,
                             const ConstTensor3View& t_field,
                             const ConstTensor4View& vmr_field,
                             const ConstMatrixView& z_surface,
                             const Vector& refellipsoid,
                             const Vector& f_grid,
                             const Verbosity& verbosity);

/** Calculates a propagation path for 1D
 *
 * @param ws[in,out] Current workspace.
 * @param[out] Pressures Pressure points along propagation path
 * @param[out] Temperatures Temperature along propagation path
 * @param[out] Vmrs Volume mixing ration along propagation path
 * @param[out] cloud_gp_p Pressure interpolation gridpoints along propagation path
 * @param[out] cloud_gp_za Zenith angle interpolation gridpoints along
 *              propagation path
 * @param[out] lstep Length of each propagation step along propagation path
 * @param[out] itw Pressure interpolation weights along propagation path
 * @param[out] itw_za Zenith angle interpolation weights along propagation path
 * @param[in] p_index Index of pressure gridpoint
 * @param[in] za_index index of zenith angle gridpoint
 * @param[in] za_grid Zenith angle grid of Spectral radiance field inside the
 *              cloudbox.
 * @param[in] ppath_step_agenda Agenda to calculate a propagation path step.
 * @param[in] ppath_lmax Maximum length between points describing propagation
 *              paths.
 * @param[in]ppath_lraytrace Maximum length of ray tracing steps when
 *              determining propagation paths.
 * @param[in] p_grid Pressure grid.
 * @param[in] lat_grid Latitude grid.
 * @param[in] lon_grid Longitude grid.
 * @param[in] z_field Field of geometrical altitudes.
 * @param[in] t_field Temperature field.
 * @param[in] vmr_field Field of volume mixing ratios.
 * @param[in] z_surface Surface altitude.
 * @param[in] refellipsoid Reference ellipsoid.
 * @param[in] f_grid f_grid The frequency grid for monochromatic pencil beam
 *              calculations.
 * @param[in] verbosity Verbosity setting.
 */
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
                            const Agenda& ppath_step_agenda,
                            const Numeric& ppath_lmax,
                            const Numeric& ppath_lraytrace,
                            const ConstVectorView& p_grid,
                            const ConstVectorView& lat_grid,
                            const ConstVectorView& lon_grid,
                            const ConstTensor3View& z_field,
                            const ConstTensor3View& t_field,
                            const ConstTensor4View& vmr_field,
                            const ConstMatrixView& z_surface,
                            const ConstVectorView& refellipsoid,
                            const ConstVectorView& f_grid,
                            const Verbosity& verbosity);

//------------------------------------------------------------------------------
// Auxilary functions-----------------------------------------------------------
//
//------------------------------------------------------------------------------


//TODO:Add doxygen doc
void FlattenMeshGrid(Matrix& flattened_meshgrid,
                     const Vector& grid1,
                     const Vector& grid2);

//TODO:Add doxygen doc
void FlattenMeshGridIndex(ArrayOfIndex& index_grid1,
                          ArrayOfIndex& index_grid2,
                          const Vector& grid1,
                          const Vector& grid2);

// Subscript to index functions-------------------------------------------------
//TODO:Add doxygen doc
//2d
Index subscript2index(Index i, Index j, Index Ni);

//3d
Index subscript2index(Index i, Index j, Index k, Index Ni, Index Nj);

//4d
Index subscript2index(
    Index i, Index j, Index k, Index l, Index Ni, Index Nj, Index Nk);

//5d
Index subscript2index(Index i,
                     Index j,
                     Index k,
                     Index l,
                     Index m,
                     Index Ni,
                     Index Nj,
                     Index Nk,
                     Index Nl);


// Index to subscript functions-------------------------------------------------
//TODO:Add doxygen doc
//2d array
typedef std::tuple<int, int> Tuple2D;
Tuple2D index2subscript(Index idx, Index Ni);

//3d array
typedef std::tuple<int, int, int> Tuple3D;
Tuple3D index2subscript(Index idx, Index Ni, Index Nj);

//4d array
typedef std::tuple<int, int, int, int> Tuple4D;
Tuple4D index2subscript(Index idx, Index Ni, Index Nj, Index Nk);

//5d array
typedef std::tuple<int, int, int, int, int> Tuple5D;
Tuple5D index2subscript(Index idx, Index Ni, Index Nj, Index Nk, Index Nl);