/* Copyright (C) 2021
   Jon Petersen <jon.petersen@studium.uni-hamburg.de>
   Manfred Brath  <manfred.brath@uni-hamburg.de>

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
   USA. */

/*===========================================================================
  ===  File description
  ===========================================================================*/

#include "agenda_class.h"
#include "arts.h"
#include "rte.h"

/*!
  \file   m_gas_scattering.cc
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>,
          Manfred Brath  <manfred.brath@.uni-hamburg.de>
  \date   2021-02-08

  \brief  Workspace functions related to gas scattering.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric BOLTZMAN_CONST;

/*===========================================================================
  === The functions
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scatteringOff(Index& gas_scattering_do,
                       Agenda& gas_scattering_agenda,
                       const Verbosity&) {
  // set flag to False (default)
  gas_scattering_do = 0;

  gas_scattering_agenda = Agenda();
  gas_scattering_agenda.set_name("gas_scattering_agenda");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scatteringCoefXsecConst(PropagationMatrix& sca_coef,
                                 const Vector& f_grid,
                                 const Numeric& rtp_pressure,
                                 const Numeric& rtp_temperature,
                                 const Index& stokes_dim,
                                 const Numeric& ConstXsec,
                                 const Verbosity&) {
  // Some basic sizes
  const Index nf = f_grid.nelem();

  // Number density
  Numeric N;
  N = rtp_pressure / rtp_temperature / BOLTZMAN_CONST;

  //Vector of constant cross sections
  Vector Xsec(nf, ConstXsec);

  // set coefficients
  PropagationMatrix sca_coef_temp(nf, stokes_dim);
  sca_coef_temp.SetZero();
  sca_coef_temp.Kjj() += Xsec;
  sca_coef_temp.Kjj() *= N;

  sca_coef=sca_coef_temp;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scatteringMatrixIsotropic(TransmissionMatrix& sca_mat,
                                   const Vector& f_grid,
                                   const Index& stokes_dim,
                                   const Verbosity&) {

  //TODO add in_los,out_los as input for the machanism that if in_los or out_los
  // is empty then sca_mat is empty.

  TransmissionMatrix sca_mat_temp(f_grid.nelem(), stokes_dim);
  sca_mat_temp.setIdentity();
  sca_mat_temp *= 1 / (4 * PI);

  sca_mat=sca_mat_temp;
}