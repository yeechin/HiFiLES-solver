/*!
 * \file flux.h
 * \author - Original code: SD++ developed by Patrice Castonguay, Antony Jameson,
 *                          Peter Vincent, David Williams (alphabetical by surname).
 *         - Current development: Aerospace Computing Laboratory (ACL)
 *                                Aero/Astro Department. Stanford University.
 * \version 0.1.0
 *
 * High Fidelity Large Eddy Simulation (HiFiLES) Code.
 * Copyright (C) 2014 Aerospace Computing Laboratory (ACL).
 *
 * HiFiLES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HiFiLES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HiFiLES.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "arrayt.h"

/*! calculate inviscid flux in 2D */
void calc_invf_2d(arrayt<double>& in_u, arrayt<double>& out_f);

/*! calculate inviscid flux in 3D */
void calc_invf_3d(arrayt<double>& in_u, arrayt<double>& out_f);

/*! calculate viscous flux in 2D */
void calc_visf_2d(arrayt<double>& in_u, arrayt<double>& in_grad_u, arrayt<double>& out_f);

/*! calculate viscous flux in 3D */
void calc_visf_3d(arrayt<double>& in_u, arrayt<double>& in_grad_u, arrayt<double>& out_f);

/*!
 * \brief calculate & add addtional ALE flux term in 2D
 * \param[in] in_u - Solution vector
 * \param[in] in_v - Grid velocity
 * \param[in,out] out_f - Modified flux vector
 */
void calc_alef_2d(arrayt<double>& in_u, arrayt<double>& in_v, arrayt<double>& out_f);

/*!
 * \brief calculate & add addtional ALE flux term in 3D
 * \param[in] in_u - Solution vector
 * \param[in] in_v - Grid velocity
 * \param[in,out] out_f - Modified flux vector
 */
void calc_alef_3d(arrayt<double>& in_u, arrayt<double>& in_v, arrayt<double>& out_f);
