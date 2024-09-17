
/* 
  Maggie Celeste 2024

  The functions opac and opac_ice are adapted from Zhaohuan Zhu et al 2012, available here: https://www.physics.unlv.edu/~zhzhu/Opacity.html
  opac() and opac_ice() generate opacities from power-law fits.

  The rest of the code is for turning opacities generated in opac or opac_ice into opacity tables or pseudo-mean-opacity tables,
  for use with the modified Lombardi cooling algorithm.
  
  You don't need to re-generate these tables unless you want to change them.
*/


#include <stdio.h>
#include <math.h>

// Function by Zhu et al. 2012 to generate opacity for a given density and temperature
double opac(double rho, double T)
{
	double xlop, xlp, xlt, rosstabd,pre;
	pre=rho*T*8.314472*1.e7/2.4;
	if (pre < 0. || T < 0.) {
		fprintf(stderr, "error: pre or T negative\n");
		fprintf(stderr, "pre: %g T: %g\n", pre, T);
		return(10.);
	}
	if (pre == 0 || T == 0)
		return (1.);
	xlp = log10(pre);
	xlt = log10(T);


        if (xlt<3.+0.03*(xlp+4.)){
          xlop=-1.27692+0.73846*xlt;
	}
        else if (xlt<3.08+0.028084*(xlp+4)){
          xlop=129.88071-42.98075*xlt+(142.996475-129.88071)*0.1*(xlp+4);
	}
        else if (xlt<3.28+xlp/4.*0.12){
          xlop=-15.0125+4.0625*xlt;
	}
        else if (xlt<3.41+0.03328*xlp/4.){
          xlop=58.9294-18.4808*xlt+(61.6346-58.9294)*xlp/4.;
	}
        else if (xlt<3.76+(xlp-4)/2.*0.03){
          xlop=-12.002+2.90477*xlt+(xlp-4)/4.*(13.9953-12.002);
	}
        else if (xlt<4.07+(xlp-4)/2.*0.08){
          xlop=-39.4077+10.1935*xlt+(xlp-4)/2.*(40.1719-39.4077);
	}
        else if (xlt<5.3715+(xlp-6)/2.*0.5594){
          xlop=17.5935-3.3647*xlt+(xlp-6)/2.*(17.5935-15.7376);
	}
        else{
        xlop=-0.48;
        }
        if (xlop<3.586*xlt-16.85&&xlt<4.){xlop=3.586*xlt-16.85;}
        if (xlt<2.9){xlop=-1.27692+0.73846*xlt;}
        rosstabd=pow(10.,xlop);
        return(rosstabd);

}

// Function by Zhu et al. 2012 to generate opacity for a given density and temperature
// -- this time with ice
double opac_ices(double rho,double T) 
{
	double xlop, xlp, xlt, rosstabd, pre;
	pre=rho*T*8.314472*1.e7/2.4;
	if (pre < 0. || T < 0.) {
		fprintf(stderr, "error: pre or T negative\n");
		fprintf(stderr, "pre: %g T: %g\n", pre, T);
	}
	if (pre == 0 || T == 0)
		return (1.);
	xlp = log10(pre);
	xlt = log10(T);


	if(xlt<2.23567+0.01899*(xlp-5.)){
	  xlop=1.5*(xlt-1.16331)-0.736364;
	}
	else if (xlt<2.30713+0.01899*(xlp-5.)){
	  xlop=-3.53154212*xlt+8.767726-(7.24786-8.767726)*(xlp-5.)/16.;
	}
	else if (xlt<2.79055){
          xlop=1.5*(xlt-2.30713)+0.62 ;
	}
	else if  (xlt<2.96931){
	   xlop=-5.832*xlt+17.7;
	}
	else if (xlt<3.29105+(3.29105-3.07651)*(xlp-5.)/8.){
	   xlop=2.129*xlt-5.9398;
	}
        else if (xlt<3.08+0.028084*(xlp+4)){
          xlop=129.88071-42.98075*xlt+(142.996475-129.88071)*0.1*(xlp+4);
	}
        else if (xlt<3.28+xlp/4.*0.12){
          xlop=-15.0125+4.0625*xlt;
	}
        else if (xlt<3.41+0.03328*xlp/4.){
          xlop=58.9294-18.4808*xlt+(61.6346-58.9294)*xlp/4.;
	}
        else if (xlt<3.76+(xlp-4)/2.*0.03){
          xlop=-12.002+2.90477*xlt+(xlp-4)/4.*(13.9953-12.002);
	}
        else if (xlt<4.07+(xlp-4)/2.*0.08){
          xlop=-39.4077+10.1935*xlt+(xlp-4)/2.*(40.1719-39.4077);
	}
        else if (xlt<5.3715+(xlp-6)/2.*0.5594){
          xlop=17.5935-3.3647*xlt+(xlp-6)/2.*(17.5935-15.7376);
	}
        else{
        xlop=-0.48;
        }
        if (xlop<3.586*xlt-16.85&&xlt<4.){xlop=3.586*xlt-16.85;}
        rosstabd=pow(10.,xlop);
        return(rosstabd);
}



// These functions are used in calculating the pseudo-mean opacity (see: Lombardi et al 2015)
void theta_f(double eps_val, double n, double delta, double& theta, double& theta_prime)
{
  double phi = 0;
  theta = 1;
  double eps = 0;
  while (eps < eps_val)
  {
    eps += delta;
    phi += pow(theta, n) * pow(eps, 2) * delta;
    theta += - phi / pow(eps, 2) * delta;
  }
  theta_prime = - phi / pow(eps, 2);
}

double integrand(double rho, double T, double n, double eps1, double eps2, double delta)
{
  double theta_eps1, theta_eps2, theta_eps1_prime, theta_eps2_prime;
  theta_f(eps1, n, delta, theta_eps1, theta_eps1_prime);
  theta_f(eps2, n, delta, theta_eps2, theta_eps2_prime);

  double k = opac(rho * pow(theta_eps2 / theta_eps1, n), T * theta_eps2 / theta_eps1);
  double integrand = ((theta_eps1_prime / theta_eps1) * pow(eps1, 2))
                    * (k * pow(theta_eps2, n));
  return integrand;
}

double pseudo_mean_opacity(double rho, double T)
{
  // Now we need to also generate our pseudo-mean opacity table for a given polytropic index, n.
  double n = 1.5;         // polytropic index -- for protoplanetary discs, somewhere between 1 and 1.5
  double Zeta_prime = 1.014; //NB: THIS CHANGES WITH n, see table 1. of Lombardi et al 2015
  double delta = 1e-2;      // step for numerical integration / ODE solving

  // First we need to find eps_B for our chosen value of polytropic index n, by finding theta(eps) = 0 when eps=eps_B
  // note -- we solve the Lane-Emden function as a pair of ODEs; the wikipedia page is pretty helpful...double eps = 0;
  double phi = 0;
  double theta = 1;
  double eps = 0;
  double delta_phi, delta_theta;

  while (theta >= 0)
  {
    eps += delta;
    delta_phi = pow(theta, n) * pow(eps, 2) * delta;
    phi += delta_phi;
    delta_theta = - phi / pow(eps, 2) * delta;
    theta += delta_theta;
  }
  double eps_B = eps - delta; //go back to the last value before 0, since current value is below 0
  double theta_prime_B = - (phi - delta_phi) / pow(eps_B, 2);
  
  // we're going to need the derivative of theta at eps_B too:
  double multiplied_part = (n + 1) / (Zeta_prime * pow(eps_B, 2) * theta_prime_B);
  double theta_prime = 0;   // differential of theta wrt eps
  double f_a, f_b, g_a, g_b;  // for trapezoidal rule
  double integral=0;
  double first_integral, second_integral;
  double eps_prime, theta_eps_prime, phi_eps_prime;  // note: theta_eps_prime = theta(eps_prime)
  double delta_theta_eps_prime, delta_phi_eps_prime;

  // IC of integral:
  phi = 0;
  theta = 1;
  eps = 0;
  f_a = 0;

  while (eps < eps_B)
  { // integrate first integral from eps=0 to eps=eps_B
    //printf("integration is %g %% done:\n", eps/eps_B * 100);
    eps += delta;
    delta_phi = pow(theta, n) * pow(eps, 2) * delta;
    phi += delta_phi;
    delta_theta = - phi / pow(eps, 2) * delta;
    theta += delta_theta;
    theta_prime = - phi / pow(eps,2);

    f_b = (theta_prime / theta) * pow(eps, 2);
    first_integral = delta * (f_a + f_b) / 2.0;
    f_a = f_b;
    // now we hold eps constant while we integrate over eps_prime:
    // IC:
    second_integral=0;
    eps_prime = eps;
    theta_eps_prime = theta;
    phi_eps_prime = phi;
    g_a = opac(rho, T) * pow(theta_eps_prime, n);   // since rho_adjusted = rho when theta_eps_prime = theta
    while (eps_prime < eps_B)
    { // integrate second integral from eps to eps_B
      eps_prime += delta;
      delta_phi_eps_prime = pow(theta_eps_prime, n) * pow(eps_prime, 2) * delta;
      phi_eps_prime += delta_phi_eps_prime;
      delta_theta_eps_prime = - phi_eps_prime / pow(eps_prime, 2) * delta;
      theta_eps_prime += delta_theta_eps_prime;
      //printf("adjusted opacity = %g vs %g\n", opac(rho * pow(theta_eps_prime / theta, n), T * (theta_eps_prime / theta)), opac(rho, T));
      g_b = opac(rho * pow(theta_eps_prime / theta, n), T * (theta_eps_prime / theta)) * pow(theta_eps_prime, n);
      second_integral += delta * (g_a + g_b) / 2.0;
      g_a = g_b;
    }
    integral += (first_integral * second_integral);
  }
  double pseudo_mean_opac = multiplied_part * integral;
  return (pseudo_mean_opac);
}


double pseudo_mean_opacity_ice(double rho, double T)
{
  // Now we need to also generate our pseudo-mean opacity table for a given polytropic index, n.
  double n = 1.5;         // polytropic index -- for protoplanetary discs, somewhere between 1 and 1.5
  double Zeta_prime = 1.014; //NB: THIS CHANGES WITH n, see table 1. of Lombardi et al 2015
  double delta = 1e-2;      // step for numerical integration / ODE solving

  // First we need to find eps_B for our chosen value of polytropic index n, by finding theta(eps) = 0 when eps=eps_B
  // note -- we solve the Lane-Emden function as a pair of ODEs; the wikipedia page is pretty helpful...double eps = 0;
  double phi = 0;
  double theta = 1;
  double eps = 0;
  double delta_phi, delta_theta;

  while (theta >= 0)
  {
    eps += delta;
    delta_phi = pow(theta, n) * pow(eps, 2) * delta;
    phi += delta_phi;
    delta_theta = - phi / pow(eps, 2) * delta;
    theta += delta_theta;
  }
  double eps_B = eps - delta; //go back to the last value before 0, since current value is below 0
  double theta_prime_B = - (phi - delta_phi) / pow(eps_B, 2);
  
  // we're going to need the derivative of theta at eps_B too:
  double multiplied_part = (n + 1) / (Zeta_prime * pow(eps_B, 2) * theta_prime_B);
  double theta_prime = 0;   // differential of theta wrt eps
  double f_a, f_b, g_a, g_b;  // for trapezoidal rule
  double integral=0;
  double first_integral, second_integral;
  double eps_prime, theta_eps_prime, phi_eps_prime;  // note: theta_eps_prime = theta(eps_prime)
  double delta_theta_eps_prime, delta_phi_eps_prime;

  // IC of integral:
  phi = 0;
  theta = 1;
  eps = 0;
  f_a = 0;

  while (eps < eps_B)
  { // integrate first integral from eps=0 to eps=eps_B
    //printf("integration is %g %% done:\n", eps/eps_B * 100);
    eps += delta;
    delta_phi = pow(theta, n) * pow(eps, 2) * delta;
    phi += delta_phi;
    delta_theta = - phi / pow(eps, 2) * delta;
    theta += delta_theta;
    theta_prime = - phi / pow(eps,2);

    f_b = (theta_prime / theta) * pow(eps, 2);
    first_integral = delta * (f_a + f_b) / 2.0;
    f_a = f_b;
    // now we hold eps constant while we integrate over eps_prime:
    // IC:
    second_integral=0;
    eps_prime = eps;
    theta_eps_prime = theta;
    phi_eps_prime = phi;
    g_a = opac_ices(rho, T) * pow(theta_eps_prime, n);   // since rho_adjusted = rho when theta_eps_prime = theta
    while (eps_prime < eps_B)
    { // integrate second integral from eps to eps_B
      eps_prime += delta;
      delta_phi_eps_prime = pow(theta_eps_prime, n) * pow(eps_prime, 2) * delta;
      phi_eps_prime += delta_phi_eps_prime;
      delta_theta_eps_prime = - phi_eps_prime / pow(eps_prime, 2) * delta;
      theta_eps_prime += delta_theta_eps_prime;
      //printf("adjusted opacity = %g vs %g\n", opac(rho * pow(theta_eps_prime / theta, n), T * (theta_eps_prime / theta)), opac(rho, T));
      g_b = opac_ices(rho * pow(theta_eps_prime / theta, n), T * (theta_eps_prime / theta)) * pow(theta_eps_prime, n);
      second_integral += delta * (g_a + g_b) / 2.0;
      g_a = g_b;
    }
    integral += (first_integral * second_integral);
  }
  double pseudo_mean_opac = multiplied_part * integral;
  return (pseudo_mean_opac);
}



int main()
{
FILE *fp = fopen("opacity_table.txt", "w"); //a not w!

for (int i=0; i < 100; i++) 
{
    double T = pow(10.0, i/25.0);

    for (int j=-100; j < 20; j++){
        double rho = pow(10, j/5.0);
        
        double opacity = opac(rho, T);

        fprintf(fp, "%lf ", opacity );
    }

    fprintf(fp, "\n");
}
fclose(fp);


FILE *fp_ice = fopen("opacity_table_ice.txt", "w"); //a not w!

for (int i=0; i < 100; i++) 
{
    double T = pow(10.0, i/25.0);

    for (int j=-100; j < 20; j++){
        double rho = pow(10, j/5.0);
        
        double opacity = opac_ices(rho, T);

        fprintf(fp_ice, "%lf ", opacity );
    }

    fprintf(fp_ice, "\n");
}
fclose(fp_ice);


FILE *fp2 = fopen("pseudo_mean_opacity_table_ice.txt", "w"); 

for (int i=0; i < 100; i++) 
{
    double T = pow(10.0, i/25.0);

    for (int j=-100; j < 20; j++){
        printf("%i i, %i j \n", i, j);
        double rho = pow(10, j/5.0);
        
        double pseudo_mean_opac = pseudo_mean_opacity_ice(rho, T);
        fprintf(fp2, "%lf ", pseudo_mean_opac );
    }

    fprintf(fp2, "\n");
}
fclose(fp2);


/* The following is an example of how to look these tables up later on.

double pseudo_mean_opac = pseudo_mean_opacity_long(100, 10, 1.5, 1e-2);
double opacity = opac(100, 10);
double pseudo_mean = pseudo_mean_opacity(100, 10);
printf("OPAC=%g PSEUDO_OPAC=%g or %g \n", opacity, pseudo_mean_opac, pseudo_mean);

double opacities[100][120];
int row, col;
char dummy;

FILE *fp = fopen("opacity_table.txt", "r");

for (row=0; row<100; row++) {
    for (col=0; col<120; col++) {
        fscanf(fp, "%lf", &opacities[row][col]);
    }
    printf("row %d op %lf dummy %c ", row, opacities[row][119], dummy);
  }
*/

/*
i = log_10(T)/0.04
i_1 = int(i)
i_2 = i_1 + 1
j = log_10(rho)/0.2 + 10
j_1 = int(j)
j_2 = j_1 + 1


op11 = opacity(i_1, j_1)
op12 = opacity(i_1, j_2)
op21 = opacity(i_2, j_1)
op22 = opacity(i_2, j_2)

//interpolate along i, between op11 and op21, and between op12 and op22
half_interp_op1 = (i_2 - i) / (i_2 - i_1) * op11  +  (i - i_1) / (i_2 - i_1) * op21
half_interp_op2 = (i_2 - i) / (i_2 - i_1) * op21  +  (i - i_1) / (i_2 - i_1) * op22

//interpolate the i-interpolated values along j 
interp_op = (j_2 - j) / (j_2 - j_1) * interp_op1  +  (j - j_1) / (j_2 - j_1) * interp_op2
//since we've already normalised to an evenly spaced integer grid, j_2 - j_1 and i_2 - i_1 will both be 1
*/
return(1);
}