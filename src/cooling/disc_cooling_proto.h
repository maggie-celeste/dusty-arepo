/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/cooling/disc_cooling_proto.h
 * \date        10/2023
 * \brief       Header for cooling functions.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 26.09.2024 Prepared file for public release -- Maggie Celeste
 */

#ifndef INLINE_FUNC
#define INLINE_FUNC
#endif /* #ifndef INLINE_FUNC */

double star_mass;
double star_x;
double star_y;
double star_z;


#ifdef BETA_COOLING
double DoCooling_Beta(double u_old, double rho, double dt, double *ne_guess, double omega, double Radius);
double GetCoolingTime_Beta(double u_old, double rho, double *ne_guess, double omega, double Radius);
void cooling_only_Beta(void);
void cool_cell_Beta(int i);
#endif

#ifdef MOD_LOMBARDI_COOLING
void MOD_LOMBARDI_COOLING(void);
#endif