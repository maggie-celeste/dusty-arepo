/*!
 * 
 * \file        src/hydro/dust.h
 * \date        10/2022
 * \brief       Header for dust stokes functions.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef INLINE_FUNC
#define INLINE_FUNC
#endif /* #ifndef INLINE_FUNC */

#ifdef DUST_STOKES
double K_from_stokes(double Radius, 
                     double rho_d, double rho_g, int FB);

#endif

#ifdef DUST_SIZE
double K_from_size(double c_s, double rho_d, double rho_g, int FB);
#endif

void check_all_particles_for_NAN(void);

void check_state_for_NAN(struct state *st);