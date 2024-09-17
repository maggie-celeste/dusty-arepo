#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief This routine returns K*rho given a particular Stokes number.
 * Assumes Keplerian orbit about a star of specified mass and position.
 * stopping time t_s defined as: del (v_d - v_g)/del t = (v_d-v_g) / t_s
 * we make use of the drag coefficient, K, where K is defined as follows...
 * rate of change of dust momentum per unit volume: del p_d/del t = K * (v_d-v_g) * rho_g * rho_d
 * and rate of change of gas momentum per unit vol: del p_g/del t = - FB* K * (v_d-v_g) * rho_g * rho_d
 * where FB is a switch for dust backreaction onto gas.
 * thus, in this code: K =  1 / ( t_s * (rho_g + FB* rho_d))
 *
 *  \param[in] star_mass Mass of star
 *  \param[in] pID      ID Number of particle (for pulling up particle positions)
 *  \param[in] rho_d    Dust density per unit volume
 *  \param[in] rho_g    Gas density per unit volume
 *  \param[in] FB       feedback of dust onto gas.
 *  \param[in] star_x   X position of star
 *  \param[in] star_y   Y position of star
 *  \param[in] star_z   Z position of star
 * 
 *  \return K
 */

//TODO: EPSTEIN DRAG FUNCTION


#ifdef DUST_STOKES
double K_from_stokes(double Radius, 
                     double rho_d, double rho_g, int FB)
{     
      double Omega, t_stop, K;

      double stokes = DUST_STOKES;
      // Dynamical timescale
      double v0 = 1.0;
      double R0 = 1.0;
      double eps_d = 17/3.;
      double vmax = v0 / sqrt(1 - exp(-eps_d*R0));

      
      
      if (Radius < 1e-3)
      {
            //printf("grav_external.c: Correcting a division by zero...\n R=%f, dx=%f, dy=%f", Radius, Posx, Posy);
            Radius= eps_d * 1e-3;
      }
      double V_gal = vmax * sqrt( 1 - exp(-eps_d*Radius));

      Omega = V_gal/Radius;
      //Omega=1.0;
      
      t_stop = stokes / Omega;
      
      //SUSPECT IT'S DRIVEN BY LARGE OR SMALL VALUES OF K.
      K = 1.0 / (t_stop * ( rho_g  + FB*rho_d));
      return K;
};
#endif


#ifdef DUST_SIZE
double K_from_size(double c_s, double rho_d, double rho_g, int FB)
{     
      double t_stop, K;
#ifdef DUST_RHO_GRAIN
      double rho_grain = DUST_RHO_GRAIN *  (All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm / All.UnitMass_in_g);
#else
      double rho_grain = 2 * (All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm / All.UnitMass_in_g); //about 2 gcm^-3
#endif
      double size_grain = DUST_SIZE;
      
      // t_stop from Epstein drag law:
      // TODO: directly compute K instead
      t_stop = rho_grain * size_grain / (rho_g * c_s);
      K = 1.0 / (t_stop * ( rho_g  + FB*rho_d));
      
      return K;
};
#endif

void check_all_particles_for_NAN(void)
{
      int idx, i;
      int is_nan=0;
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;
          if (SphP[i].DustMomentum[0]!=SphP[i].DustMomentum[0] ||
              SphP[i].DustMomentum[1]!=SphP[i].DustMomentum[1] ||
              SphP[i].DustMomentum[2]!=SphP[i].DustMomentum[2] ||
              SphP[i].Momentum[0] != SphP[i].Momentum[0] ||
              SphP[i].Momentum[1] != SphP[i].Momentum[1] ||
              SphP[i].Momentum[2] != SphP[i].Momentum[2]
            )
            {is_nan=is_nan+1;
             printf("        NAN DETECTED FOR PARTICLE:\n Numgas=%i, pID=%i,\n velx=%10f, vely=%10f, velz=%10f, \n dustvelx=%10f, dustvely=%10f, dustvelz=%10f, \n density=%10f, dustdensity=%10f  \n posx=%8f, posy=%8f, posz=%8f",
                    NumGas, i, P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], SphP[i].DustVel[0], SphP[i].DustVel[1], SphP[i].DustVel[2],
                    SphP[i].Density, SphP[i].DustDensity, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);}  

          if (P[i].Vel[0]!=P[i].Vel[0] ||
            P[i].Vel[1]!=P[i].Vel[1] ||
            P[i].Vel[2]!=P[i].Vel[2] ||
            SphP[i].DustVel[0] != SphP[i].DustVel[0] ||
            SphP[i].DustVel[1] != SphP[i].DustVel[1] ||
            SphP[i].DustVel[2] != SphP[i].DustVel[2]
            )
            {is_nan=is_nan+1;
             printf("        NAN DETECTED FOR PARTICLE:\n Numgas=%i, pID=%i,\n velx=%10f, vely=%10f, velz=%10f, \n dustvelx=%10f, dustvely=%10f, dustvelz=%10f, \n density=%10f, dustdensity=%10f  \n posx=%8f, posy=%8f, posz=%8f",
                    NumGas, i, P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], SphP[i].DustVel[0], SphP[i].DustVel[1], SphP[i].DustVel[2],
                    SphP[i].Density, SphP[i].DustDensity, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);}      

        if (SphP[i].Density!= SphP[i].Density)
            {is_nan=is_nan+1;
             printf("        NAN DETECTED FOR PARTICLE:\n Numgas=%i, pID=%i,\n velx=%10f, vely=%10f, velz=%10f, \n dustvelx=%10f, dustvely=%10f, dustvelz=%10f, \n density=%10f, dustdensity=%10f  \n posx=%8f, posy=%8f, posz=%8f",
                    NumGas, i, P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], SphP[i].DustVel[0], SphP[i].DustVel[1], SphP[i].DustVel[2],
                    SphP[i].Density, SphP[i].DustDensity, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);}  

        if (SphP[i].DustDensity!= SphP[i].DustDensity)
            //print out all the properties prior to terminating. 
            //TODO: tally up nans, terminate after
            {is_nan=is_nan+1;
             printf("        NAN DETECTED FOR PARTICLE:\n Numgas=%i, pID=%i,\n velx=%10f, vely=%10f, velz=%10f, \n dustvelx=%10f, dustvely=%10f, dustvelz=%10f, \n density=%10f, dustdensity=%10f  \n posx=%8f, posy=%8f, posz=%8f",
                    NumGas, i, P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], SphP[i].DustVel[0], SphP[i].DustVel[1], SphP[i].DustVel[2],
                    SphP[i].Density, SphP[i].DustDensity, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);}  
        
        
        if (SphP[i].Density < 0)
            //print out all the properties prior to terminating. 
            //TODO: tally up nans, terminate after
            {is_nan=is_nan+1;
             printf("        NEGATIVE GAS DENSITY DETECTED FOR PARTICLE:\n Numgas=%i, pID=%i,\n velx=%10f, vely=%10f, velz=%10f, \ndustvelx=%10f, dustvely=%10f, dustvelz=%10f, \n density=%10f, dustdensity=%10f  \n posx=%8f, posy=%8f, posz=%8f",
                    NumGas, i, P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], SphP[i].DustVel[0], SphP[i].DustVel[1], SphP[i].DustVel[2],
                    SphP[i].Density, SphP[i].DustDensity, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);}  

        if (SphP[i].DustDensity < 0)
            //print out all the properties prior to terminating. 
            //TODO: tally up nans, terminate after
            {is_nan=is_nan+1;
             printf("        NEGATIVE DUST DENSITY DETECTED FOR PARTICLE:\n Numgas=%i, pID=%i,\n velx=%10f, vely=%10f, velz=%10f, \n dustvelx=%10f, dustvely=%10f, dustvelz=%10f, \n density=%10f, dustdensity=%10f  \n posx=%8f, posy=%8f, posz=%8f",
                    NumGas, i, P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], SphP[i].DustVel[0], SphP[i].DustVel[1], SphP[i].DustVel[2],
                    SphP[i].Density, SphP[i].DustDensity, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);}  
        }

        if (is_nan != 0)
          { fflush(stdout);
            terminate("NANs detected \n");}
};


void check_state_for_NAN(struct state *st)
{
      int is_nan=0;

      if ((st->velx != st->velx) ||
          (st->vely != st->vely) ||
          (st->velz != st->velz) )
          {
            is_nan = is_nan+1;
            printf("    NAN DETECTED IN VELX   \n");
            printf("rho=%g velx=%g vely=%g velz=%g press=%g\n", st->rho, st->velx, st->vely, st->velz, st->press);
            printf("DUST rho=%g velx=%g vely=%g velz=%g\n", st->rhodust, st->velx_dust, st->vely_dust, st->velz_dust);
          }


      
      if ((st->velDust[0] != st->velDust[0]) ||
          (st->velDust[1] != st->velDust[1]) ||
          (st->velDust[2] != st->velDust[2]) )
          {
            is_nan = is_nan+1;
            printf("    NAN DETECTED in VELDUST   \n");
            printf("rho=%g velx=%g vely=%g velz=%g press=%g\n", st->rho, st->velx, st->vely, st->velz, st->press);
            printf(" vel_Gas=%g, %g, %g, vel_dust = %g, %g, %g \n", st->velGas[0],  st->velGas[1], st->velGas[2], st->velDust[0],  st->velDust[1], st->velDust[2]);
            printf("DUST rho=%g velx=%g vely=%g velz=%g\n", st->rhodust, st->velx_dust, st->vely_dust, st->velz_dust);
          }

        if (is_nan != 0)
          { fflush(stdout);
            terminate("NANs detected \n");}
};