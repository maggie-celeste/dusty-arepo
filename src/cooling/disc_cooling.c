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
 * \file        src/cooling/disc_cooling.c
 * \date        10/2023
 * \brief       Module for gas radiative cooling
 * \details     contains functions:
 *                double DoCooling_Beta(double u_old, double rho, double dt, 
 *                                      double *ne_guess, double omega, double Radius);
 *                double GetCoolingTime_Beta(double u_old, double rho, double *ne_guess, 
 *                                      double omega, double Radius);
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 01.09.2024 Prepared file for public release -- Maggie Celeste
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef DISC_COOLING

#ifdef BETA_COOLING
  static double BETA = BETA_COOLING;
#ifdef BETA_COOLING_POWER
  static double BETA_POWER = BETA_COOLING_POWER;
#else
  // default to fixed-beta, with the option to have beta vary as a function of R
  static double BETA_POWER = 0;
#endif // #ifdef BETA_COOLING_POWER

#ifdef BETA_COOLING_R0
  static double R_0 = BETA_COOLING_R0;
#else
  static double R_0 = 1.496e11; 
  // default to R_0 for beta(R) to 1AU
#endif // #ifdef BETA_COOLING_R0
#endif // #ifdef BETA_COOLING

#ifdef MOD_LOMBARDI_COOLING
  static double (*opacities) = NULL; 
  
#ifdef LOMBARDI_TEQ
  static double T_EQ = LOMBARDI_TEQ;
#else
  static double T_EQ = 5.0; // default background temperature
#endif //#ifdef LOMBARDI_TEQ
#endif //#ifdef MOD_LOMBARDI_COOLING



#ifdef BETA_COOLING
/*! \brief Computes the new internal energy per unit mass.
 *
 *  Computes cooling rate as in Lodato & Rice 2004
 *
 *  \param[in] u_old the initial (before cooling is applied) internal energy
 *             per unit mass of the gas cell.
 *  \param[in] rho   the proper density of the gas cell.
 *  \param[in] dt    the duration of the time step.
 *  \param[in] ne_guess electron number density relative to hydrogen number
 *             density (for molecular weight computation).
 *  \param[in] omega angular frequency (got from star mass + position)
 *  \return The new internal energy per unit mass of the gas cell.
 */
double DoCooling_Beta(double u_old, double rho, double dt, double *ne_guess, double omega, double Radius)
{
  double u ;

  DoCool.u_old_input    = u_old;
  DoCool.rho_input      = rho;
  DoCool.dt_input       = dt;
  DoCool.ne_guess_input = *ne_guess;

  if(!gsl_finite(u_old))
    terminate("invalid input: u_old=%g\n", u_old);

  /* bracketing */
  double t_cool = GetCoolingTime_Beta(u_old, rho, ne_guess, omega, Radius) ;
  
  u = u_old / (1 + dt/t_cool) ;

  return u;
}

/*! \brief Returns the cooling time.
 *
 *  If we actually have heating, a cooling time of 0 is returned.
 *
 *  \param[in] u_old The initial (before cooling is applied) internal energy
 *             per unit mass of the gas cell.
 *  \param[in] rho The proper density of the gas cell.
 *  \param[in] ne_guess Electron number density relative to hydrogen number
 *             density (for molecular weight computation).
 *
 *  \return Cooling time; 0 if heating.
 */
double GetCoolingTime_Beta(double u_old, double rho, double *ne_guess, double omega, double Radius)
{
  // If we set beta to vary with radius st Beta = (R/R_0)^1/2 * Beta_0
  // => Beta(100AU) = (100/5)**1/2 * Beta_0 = 4.5 * Beta_0
  double beta_mod = pow((Radius/R_0), BETA_POWER) * BETA;
  return beta_mod / omega ;
}

/*! \brief Apply the beta cooling to all the active gas cells.
 *
 *  \return void
 */
void cooling_only_Beta(void) /* normal cooling routine when star formation is disabled */
{
  int idx, i;

  CPU_Step[CPU_MISC] += measure_time();
  
  star_mass = PartSpecialListGlobal[0].mass;
  star_x = PartSpecialListGlobal[0].pos[0];
  star_y = PartSpecialListGlobal[0].pos[1];
  star_z = PartSpecialListGlobal[0].pos[2];
  
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue; /* skip cells that have been swallowed or eliminated */

          cool_cell_Beta(i);
        }
    }
  CPU_Step[CPU_COOLINGSFR] += measure_time();
}

/*! \brief Apply the isochoric cooling to a given gas cell.
 *
 *  This function applies the normal isochoric cooling to a single gas cell.
 *  Once the cooling has been applied according to one of the cooling models
 *  implemented, the internal energy per unit mass, the total energy and the
 *  pressure of the cell are updated.
 *
 *  \param[in] i Index of the gas cell to which cooling is applied.
 *
 *  \return void
 */
void cool_cell_Beta(int i)
{
  double dt, dtime, ne = 1;
  double unew, dens, dtcool, dx, dy, dz;

  dens = SphP[i].Density;
  #ifdef DUST_INCLUDE
    dens += SphP[i].DustDensity;
  #endif
  
  dx = P[i].Pos[0] - star_x;
  dy = P[i].Pos[1] - star_y;
  dz = P[i].Pos[2] - star_z;
  
  double radius = sqrt( dx*dx + dy*dy + dz*dz );
  
  double omega = sqrt( All.G * star_mass / pow(radius, 3));

  dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

  dtime = All.cf_atime * dt / All.cf_time_hubble_a;

  dtcool = dtime;

  ne         = SphP[i].Ne; /* electron abundance (gives ionization state and mean molecular weight) */
  unew       = DoCooling_Beta(dmax(All.MinEgySpec, SphP[i].Utherm), dens * All.cf_a3inv, dtcool, &ne, omega, radius);
  SphP[i].Ne = ne;

  if(unew < 0)
    terminate("invalid temperature: Thistask=%d i=%d unew=%g\n", ThisTask, i, unew);

  double du = unew - SphP[i].Utherm;

  if(unew < All.MinEgySpec)
    du = All.MinEgySpec - SphP[i].Utherm;

  SphP[i].Utherm += du;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

  #ifdef DUST_INCLUDE
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].DustMass;
  #endif

#ifdef OUTPUT_COOLHEAT
  if(dtime > 0)
    SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif /* #ifdef OUTPUT_COOLHEAT */

  set_pressure_of_cell(i);
}
#endif // #ifdef BETA_COOLING



#ifdef MOD_LOMBARDI_COOLING
/*! \brief Apply the Lombardi cooling to all the active gas cells, as described in Lombardi et al 2015.
 *
 *  \return void
 */
void MOD_LOMBARDI_COOLING(void) /* lombardi cooling routine for protoplanetary discs */
{
  double Teq = T_EQ;
  if (opacities == NULL){
    // read in opacity table
    FILE *fdcool;
    if(!(fdcool = fopen(All.LombardiOpacityFile, "r")))
      terminate("COOLING: cannot read opacity table in file `%s'\n", All.LombardiOpacityFile);
    //double opacities[100][120]; 
    opacities = mymalloc_movable(&opacities, "opacities", 120 * 100 * sizeof(double));
    int row, col;
    for (row=0; row<100; row++) {
      for (col=0; col<120; col++) {
        fscanf(fdcool, "%lf", &opacities[table_index(row, col, 120)]);
      }
    }
   fclose(fdcool);

   int indx;
   for (indx=0; indx< 100 * 120; indx++)
   {
    opacities[indx] *= 1 / ( All.UnitLength_in_cm * All.UnitLength_in_cm / All.UnitMass_in_g );
   }
  }

  if (pseudo_opacities == NULL){
    // read in opacity table
    FILE *fdcool2;
    if(!(fdcool2 = fopen(All.LombardiPseudoOpacityFile, "r")))
      terminate("COOLING: cannot read opacity table in file `%s'\n", All.LombardiPseudoOpacityFile);
    //double opacities[100][120]; 
    pseudo_opacities = mymalloc_movable(&pseudo_opacities, "pseudo_opacities", 120 * 100 * sizeof(double));
    int row, col;
    for (row=0; row<100; row++) {
      for (col=0; col<120; col++) {
        fscanf(fdcool2, "%lf", &pseudo_opacities[table_index(row, col, 120)]);
      }
    }
   fclose(fdcool2);

   int indx;
   for (indx=0; indx< 100 * 120; indx++)
   {
    pseudo_opacities[indx] *= 1 / (All.UnitLength_in_cm * All.UnitLength_in_cm / All.UnitMass_in_g);
   }
  }
  //(double *)mymalloc("opacities", sizeof(double[100][120]));

  // Apply cooling to each active cell
  double STEFAN_BOLTZMANN = 5.6704e-5 / All.UnitEnergy_in_cgs * pow(All.UnitLength_in_cm,2);
  double meanweight = 2.3; //4.0 / (1 + 3 * HYDROGEN_MASSFRAC); // note: assuming NEUTRAL GAS
  int idx, i;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue; // skip cells that have been swallowed or eliminated 

          // Start by computing the pseudo-mean column density, rho_column
          double star_mass = PartSpecialListGlobal[0].mass;
          
          double Posx, Posy, Radius, Omega;
          Posx = P[i].Pos[0] - PartSpecialListGlobal[star_idx].pos[0];
          Posy = P[i].Pos[1] - PartSpecialListGlobal[star_idx].pos[1];
          Radius = sqrt( pow(Posx,2) + pow(Posy,2));
          Omega = sqrt( All.G * star_mass / pow(Radius, 3)); 

          double c_s = sqrt(SphP[i].Utherm * (GAMMA) * (GAMMA-1));
          double Temp = SphP[i].Utherm  * meanweight * PROTONMASS / All.UnitMass_in_g  * (GAMMA-1) / BOLTZMANN * All.UnitEnergy_in_cgs;

          if (Temp < Teq){ 
            continue;
          }

          #ifdef TWODIMS
          double H = c_s / Omega;
          double rho_column = SphP[i].Density;
          #else
          double H_cs = c_s / Omega;
          double Q_3D = Omega * Omega / (All.G * 4 * 3.14 * SphP[i].Density);

          double H_disc = H_cs * 1.253 / (1 + 1 / sqrt(1.253 * Q_3D));    // MODIFYING H FOR LOMBARDI
          double H_pres = SphP[i].scaleheight;                            // ORIGINAL LOMBARDI H -- saved in voronoi_gradients_lsf.c, BEFORE slope limiters are applied.
          double H = 1.0 / sqrt( 1 / (H_disc * H_disc) + 1 / (H_pres * H_pres)); // COMBINED IN HARMONIC QUADRATURE

          if (H > FLT_MAX){
            H = FLT_MAX;
          }
          double rho_column = SphP[i].Density * 1.014 * H ;
          #endif

          // Interpolate pseudo-mean opacity from the opacities table
          // k1, j1 are the index for the opacity of Temp, Rho below the actual Temp, Rho;
          // k2, j2 are the index for the opacity of Temp, Rho above the actual Temp, Rho;
          // Use bilinear interpolation to estimate the opacity at i,j

          double k = log10(Temp)/0.04;

          #ifdef TWODIMS
          double j = log10( SphP[i].Density/H * All.UnitMass_in_g)/0.2 + 100;
          #else
          double j = log10(SphP[i].Density * All.UnitMass_in_g) / 0.2 + 100;
          #endif

          if(k > 99 ){k=99.;}
          if(k < 0 ){k=0.0;}
          int k_1 = floor(k);
          if (k_1 > 98){ k_1=98;}
          int k_2 = k_1 + 1;

          if(j > 119){j=119.0;}
          if(j < 0){j=0;}
          int j_1 = floor(j);
          if (j_1 > 118){j_1 = 118;}
          int j_2 = j_1 + 1;          

          if (k!= k){
            printf("NAN IN k i=%i!\n", i);
          }

          if (j!= j){
            printf("NAN IN j i=%i!\n", i);
          }
          
          double op11 = opacities[table_index(k_1, j_1, 120)];
          double op12 = opacities[table_index(k_1, j_2, 120)];
          double op21 = opacities[table_index(k_2, j_1, 120)];
          double op22 = opacities[table_index(k_2, j_2, 120)];

          double planck_mean_opacity = (j_2 - j) * ( (k_2 - k) * op11  +  (k - k_1) * op21 )
                                     + (j - j_1) * ( (k_2 - k) * op12  +  (k - k_1) * op22 );

          double pseudo_op11 = pseudo_opacities[table_index(k_1, j_1, 120)];
          double pseudo_op12 = pseudo_opacities[table_index(k_1, j_2, 120)];
          double pseudo_op21 = pseudo_opacities[table_index(k_2, j_1, 120)];
          double pseudo_op22 = pseudo_opacities[table_index(k_2, j_2, 120)];

          double pseudo_mean_opacity = (j_2 - j) * ( (k_2 - k) * pseudo_op11  +  (k - k_1) * pseudo_op21 )
                                     + (j - j_1) * ( (k_2 - k) * pseudo_op12  +  (k - k_1) * pseudo_op22 );

          // Now we can calculate the radiative cooling rate:
          double dt_cell = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval ; // the timestep of the cell 

          double du_dt_rad = 4 * STEFAN_BOLTZMANN * (pow(Teq, 4) - pow(Temp, 4)) / 
                             ( pow(rho_column, 2) * pseudo_mean_opacity + 1/(planck_mean_opacity));

          // then use the radiative cooling rate to compute the new utherm as in eq 21
          double Utherm_eq = Teq * (All.UnitMass_in_g / (meanweight * PROTONMASS )) * (1 / (GAMMA-1)) * BOLTZMANN / All.UnitEnergy_in_cgs;

          double t_therm_alt = ( Utherm_eq - SphP[i].OldUtherm ) / du_dt_rad;

          double C_v = (All.UnitMass_in_g / (meanweight * PROTONMASS )) * (1 / (GAMMA-1)) * BOLTZMANN / All.UnitEnergy_in_cgs;
          double temp_term = 1 / (4 * STEFAN_BOLTZMANN * (Teq + Temp) * (Teq * Teq + Temp * Temp));
          double f_tau = rho_column * rho_column * pseudo_mean_opacity + 1/(planck_mean_opacity); 

          double t_therm = C_v * temp_term * f_tau; 

          double du_hydro = SphP[i].Utherm - SphP[i].OldUtherm;

          // Then update SphP values... Thermal energy, total energy, and Pressure.
          double pre_cool_U = SphP[i].Utherm;
          SphP[i].Utherm = SphP[i].OldUtherm * exp(-dt_cell / t_therm) 
                           + Utherm_eq * (1 - exp(-dt_cell / t_therm)) 
                           + du_hydro / (1 + dt_cell / t_therm); 
          
          double old_therm_val = SphP[i].OldUtherm * exp(-dt_cell / t_therm) ;
          double eq_val = Utherm_eq * (1 - exp(-dt_cell / t_therm)) ;
          double du_hydro_val = du_hydro / (1 + dt_cell / t_therm); ;

          
          //printf("\n COOLING: Utherm %g, OLD Utherm %g, EQ Utherm %g, du_hydro %g, Radius %g, t_therm %g, t_therm_alt %g, t_therm/t_therm_alt %g",
          //       SphP[i].Utherm, SphP[i].OldUtherm, Utherm_eq, du_hydro, Radius, t_therm, t_therm_alt, t_therm / t_therm_alt);

          if (SphP[i].Utherm < 0){
            printf("Cooling Utherm %g Old Utherm %g Eq_utherm %g du_hydro %g Radius %g i=%i t_therm %g \n", 
                    SphP[i].Utherm, SphP[i].OldUtherm, Utherm_eq, du_hydro, Radius, i, t_therm);

            printf("old therm val = %g eq val = %g du_hydro val = %g t_therm %g dt_cell %g \n", 
                    old_therm_val, eq_val, du_hydro_val, t_therm, dt_cell);
            fflush(stdout);
              SphP[i].Utherm = pre_cool_U;
              printf("COOLING FAILED \n");
          }


          SphP[i].Energy += All.cf_atime * All.cf_atime * (SphP[i].Utherm - pre_cool_U) * P[i].Mass;
#ifdef DUST_INCLUDE
          //SphP[i].Energy += All.cf_atime * All.cf_atime * (SphP[i].Utherm - pre_cool_U) * SphP[i].DustMass;
#endif

          set_pressure_of_cell(i);
          
        }
    }
}
#endif

/*! \brief Apply the isochoric cooling to a given gas cell.
 *
 *  This function applies the normal isochoric cooling to a single gas cell.
 *  Once the cooling has been applied according to one of the cooling models
 *  implemented, the internal energy per unit mass, the total energy and the
 *  pressure of the cell are updated.
 *
 *  \param[in] i Index of the gas cell to which cooling is applied.
 *
 *  \return void
 */
void cool_cell(int i)
{
  double dt, dtime, ne = 1;
  double unew, dens, dtcool, dx, dy, dz;

  dens = SphP[i].Density;
  #ifdef DUST_INCLUDE
    dens += SphP[i].DustDensity;
  #endif
  
  dx = P[i].Pos[0] - star_x;
  dy = P[i].Pos[1] - star_y;
  dz = P[i].Pos[2] - star_z;
  
  double radius = sqrt( dx*dx + dy*dy + dz*dz );
  
  double omega = sqrt( All.G * star_mass / pow(radius, 3));

  dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

  dtime = All.cf_atime * dt / All.cf_time_hubble_a;

  dtcool = dtime;

  ne         = SphP[i].Ne; /* electron abundance (gives ionization state and mean molecular weight) */
  unew       = DoCooling_Beta(dmax(All.MinEgySpec, SphP[i].Utherm), dens * All.cf_a3inv, dtcool, &ne, omega, radius);
  SphP[i].Ne = ne;

  if(unew < 0)
    terminate("invalid temperature: Thistask=%d i=%d unew=%g\n", ThisTask, i, unew);

  double du = unew - SphP[i].Utherm;

  if(unew < All.MinEgySpec)
    du = All.MinEgySpec - SphP[i].Utherm;

  SphP[i].Utherm += du;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

  #ifdef DUST_INCLUDE
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].DustMass;
  #endif

#ifdef OUTPUT_COOLHEAT
  if(dtime > 0)
    SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif /* #ifdef OUTPUT_COOLHEAT */

  set_pressure_of_cell(i);
}

#endif /* #ifdef DISC_COOLING */
