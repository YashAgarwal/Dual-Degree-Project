#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef S_SPLINT_S
#include <complex.h>
#endif

#define TIME_STEP_INITIAL 50
#define MAX_FRAC_PPT 0.7

#include "../headers/functions.h"

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

void evolve(int n_x, int n_y, double delta_x, double delta_y, double kappa_c,
            double kappa_eta, double P1, double P2, double P3, double P4,
            double P5, double P6, double P7, double P8, double P9, double P10,
            double P11, double P12, double P13, double P14, double P15,
            double P16, double P17, double P18, double P19, double P20,
            double P21, double P22, double P23, double P24, double M0,
            double L0, double A, double B, double P, double c_zero,
            double r_zero, double noise_str, double gamma_i, double gamma_a,
            double beta11, double beta22, double beta12, double coeff,
            double delta_t, int time_steps, double delta_t_M0,
            int time_steps_M0, int file_timer, int theta_int, int theta_start,
            int theta_steps, int cont, int INDEX_start, char * output) {

  fftw_complex *comp, *eta; /* FFTW variables : composition, order-parameter(eta) */
  double ave_comp;
  gsl_rng *ran_num;
  const gsl_rng_type *Taus;

  FILE *fpr, *fpw, *fpw1, *fpw2;

  int INDEX = 0, t_start_updated;    /* time indexes  */
  int i1, i2, i, i_new;              /* space indexes  */
  int half_nx, half_ny;              /* Fourier transformation variables  */
  double kx, ky, delta_kx, delta_ky; /* Fourier transformation variables  */
  double kx2, ky2, kx4, ky4, kx6, ky6, k2, k3, k4; /* Fourier transformation variables  */
  double denom_eta, denom_c, factor; /* Denominator in semi-implicit fourier spectral technique */
  char NAME[50];

  char output_folder_name[50];
  strcpy(output_folder_name, output);

  double *c, *eta_r, *mod, *mod_e, *M, *L; /* variables to reproduce extended
                                              Cahn-Hilliard model (Abi-Haider,
                                              PhilMag 2001) eta in real space */
  double r;                                /* radius of the precipitates */
  double *fm_c, *fp_c; /* bulk free eng density at matrix and precipitates */
  double *W_eta, *W_eta1, *wang_f; /* w__eta */
  double pr_x1, pr_y1, pr_x2, pr_y2, pr_xy1, pr_yx1, pr_xy2, pr_yx2, *pr_t,
         *pr_ta, *aspect_ratio; /* precipitates radius */
  int r_counter, cnt_x, cnt_y;
  double r0_new, i_ini, *mx1, *mx2, *mx1_final, maxM, mod_cut,
         mod_limit; /* initial radius */
  double *alpha, *zeta;
  unsigned FLAG;
  int *ic, *jc;

  double P25, P26, P27, P28;

  P25 = 2. * (P20 - P24);
  P26 = 3. * P16 - P17;
  P27 = P16 + 1. / 9. * (P21 - P17);
  P28 = 3. * P16 - 2. / 3. * P21 - 1. / 3. * P17;

  fftw_complex *g_c, *g_eta, *M_mod, *h, *hx, *hy, *qx, *qy, *nx, *ny, *nex,
               *ney, *nxx, *nyy, *d2c; /* FFTW variables for partial derivative of f0 wrt c & eta */
  fftw_plan plan1, plan2;

  g_c = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  g_eta = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  comp = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  eta = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  d2c = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  h = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  nx = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  ny = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  nex = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  ney = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  nxx = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  nyy = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  M_mod = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  hx = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  hy = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  qx = fftw_malloc(n_x * n_y * sizeof(fftw_complex));
  qy = fftw_malloc(n_x * n_y * sizeof(fftw_complex));

  FLAG = FFTW_ESTIMATE;
  // FLAG = FFTW_EXHAUSTIVE;

  fftw_plan_with_nthreads(8);

  plan1 = fftw_plan_dft_2d(n_x, n_y, comp, comp, FFTW_FORWARD, FLAG);
  plan2 = fftw_plan_dft_2d(n_x, n_y, comp, comp, FFTW_BACKWARD, FLAG);

  (void)gsl_rng_env_setup();
  Taus = gsl_rng_taus;
  ran_num = gsl_rng_alloc(Taus);

  c = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  fm_c = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  fp_c = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  wang_f = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  W_eta = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  W_eta1 = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  eta_r = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  mod = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  mod_e = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  M = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  L = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  mx1 = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  mx2 = (double *)malloc((size_t)n_x * n_y * sizeof(double));
  mx1_final = (double *)malloc((size_t)n_x * n_y * sizeof(double));

  pr_t = (double *)malloc((size_t)(time_steps + 1) * sizeof(double));
  pr_ta = (double *)malloc((size_t)(time_steps + 1) * sizeof(double));
  aspect_ratio = (double *)malloc((size_t)(time_steps + 1) * sizeof(double));
  alpha = (double *)malloc((size_t)(time_steps + 1) * sizeof(double));
  zeta = (double *)malloc((size_t)(time_steps + 1) * sizeof(double));

  ic = (int *)malloc((size_t)n_x * n_y * sizeof(int));
  jc = (int *)malloc((size_t)n_x * n_y * sizeof(int));

  printf("continuation %d\n", cont);

  //Creating the initial profile
  if (cont == 0) {
    ave_comp = 0.0;

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {

        i = i2 + n_y * i1;

        r = sqrt((1. * n_x * delta_x / 2.0 - i1 * delta_x) *
                 (1. * n_x * delta_x / 2.0 - i1 * delta_x) +
                 (1. * n_y * delta_y / 2.0 - i2 * delta_y) *
                 (1. * n_y * delta_y / 2.0 - i2 * delta_y));

        if (r <= r_zero) {
          __real__ comp[i] = 1.0 + noise_str *(1 - 2 * gsl_rng_uniform(ran_num));
          __imag__ comp[i] = 0.0;
          __real__ eta[i] = 1.0 + noise_str * (1 - 2 * gsl_rng_uniform(ran_num));
          __imag__ eta[i] = 0.0;
        } else {
          __real__ comp[i] = c_zero + noise_str * (1 - 2 * gsl_rng_uniform(ran_num));
          __imag__ comp[i] = 0.0;
          __real__ eta[i] = 0.0 + noise_str * (1 - 2 * gsl_rng_uniform(ran_num));
          __imag__ eta[i] = 0.0;
        }
        ave_comp = ave_comp + creal(comp[i]);
      }
    }
    ave_comp = 1. * ave_comp / (n_x * n_y);

    printf("Initial Average composition = %lf", ave_comp);

    /* Calculation of initial radius r0 */
    i_ini = ((r_zero + 1.0 * delta_x * n_x / 2) / delta_x);
    r0_new = i_ini * delta_x +
             ((i_ini + 1) * delta_x -
             i_ini * delta_x) * (0.5 - 1.0) / (c_zero - 1.0) -
             delta_x * (int)n_x / 2;

    printf("r_zero = %le || i_ini = %le || r0_new = %le\n\n",
           r_zero, i_ini, r0_new);

    /*
    Copying the complex type vectors comp and eta into double type vectors c and eta
    */
    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;
        c[i] = __real__ comp[i];
        eta_r[i] = __real__ eta[i];
      }
    }

    /*
    Save the initial c and eta profile
    */
    sprintf(NAME, "%s/initial.dat", output_folder_name);
    fpw = fopen(NAME, "w");
    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;
        fprintf(fpw, "%d %d %le %le\n", i1, i2, c[i], eta_r[i]);
      }
    }
    fclose(fpw);

    half_nx = (int)n_x / 2;
    half_ny = (int)n_y / 2;

    delta_kx = (2.0 * M_PI) / (n_x * delta_x);
    delta_ky = (2.0 * M_PI) / (n_y * delta_y);

    // TIME LOOP
    sprintf(NAME, "%s/r2_t-constM.dat", output_folder_name);
    fpw1 = fopen(NAME, "w");
    sprintf(NAME, "%s/mod_r2_t-constM.dat", output_folder_name);
    fpw2 = fopen(NAME, "w");

    for (INDEX = 1; INDEX < time_steps_M0 + 1; ++INDEX){

      for (i1 = 0; i1 < n_x; ++i1) {
        for (i2 = 0; i2 < n_y; ++i2) {
          i = i2 + n_y * i1;

          c[i] = __real__ comp[i];
          eta_r[i] = __real__ eta[i];

          fm_c[i] = A * c[i] * c[i];
          fp_c[i] = B * (1 - c[i]) * (1 - c[i]);

          if (eta_r[i] < 0.0) {
            W_eta[i] = 0.0;
            W_eta1[i] = 0.0;
          }

          else if (eta_r[i] > 1.0) {
            W_eta[i] = 1.0;
            W_eta1[i] = 0.0;
          }

          else if (eta_r[i] >= 0.0 && eta_r[i] <= 1.0) {
            W_eta[i] = eta_r[i] * eta_r[i] * eta_r[i] *
                       (10 - 15 * eta_r[i] + 6 * eta_r[i] * eta_r[i]);

            //DOUBT: Is this W_eta calc wrong, 10 ?
            W_eta1[i] = 3 * eta_r[i] * eta_r[i] *
                        (10 - 15 * eta_r[i] + 6 * eta_r[i] * eta_r[i]) +
                        eta_r[i] * eta_r[i] * eta_r[i] * (12 * eta_r[i] - 15);
          }
          //DOUBT: Is this formula wrong,
          g_eta[i] = fm_c[i] * (-W_eta1[i]) + fp_c[i] * W_eta1[i] +
                     2.0 * P * eta_r[i] * (1.0 - eta_r[i]) * (1.0 - 2.0 * eta_r[i]) +
                     0.0 * _Complex_I;
          g_c[i] = 2.0 * A * c[i] * (1 - W_eta[i]) -
                   2.0 * B * (1 - c[i]) * W_eta[i] + 0.0 * _Complex_I;
          //g_c[i] = 2.0*A*comp[i]*(1.0-comp[i])*(1.0-2.0*comp[i]) +
          //0.0*_Complex_I;
        }
      }

      fftw_execute_dft(plan1, comp, comp);
      fftw_execute_dft(plan1, eta, eta);
      fftw_execute_dft(plan1, g_c, g_c);
      fftw_execute_dft(plan1, g_eta, g_eta);

      for (i1 = 0; i1 < n_x; ++i1) {
        if (i1 < half_nx)
          kx = i1 * delta_kx;
        else
          kx = (i1 - n_x) * delta_kx;
        kx2 = kx * kx;

        for (i2 = 0; i2 < n_y; ++i2) {

          if (i2 < half_ny)
              ky = i2 * delta_ky;
          else
              ky = (i2 - n_y) * delta_ky;

          ky2 = ky * ky;

          k2 = kx2 + ky2;
          k4 = k2 * k2;

          i = i2 + n_y * i1;

          denom_eta = 1.0 + 2.0 * delta_t_M0 * kappa_eta * k2;
          eta[i] = 1.0 * (eta[i] - delta_t_M0 * g_eta[i]) / denom_eta;

          denom_c = 1.0 + 2.0 * delta_t_M0 * k4 * kappa_c; // Note : We put M to be 1
          comp[i] = 1.0 * (comp[i] - k2 * delta_t_M0 * g_c[i]) /
                    denom_c; // Note : We put M to be 1
        }
      }

      fftw_execute_dft(plan2, comp, comp);
      fftw_execute_dft(plan2, eta, eta);

      for (i1 = 0; i1 < n_x; ++i1) {
        for (i2 = 0; i2 < n_y; ++i2) {
          i = i2 + n_y * i1;

          eta[i] = 1. * __real__ eta[i] / (n_x * n_y);
          __imag__ eta[i] = 0.0;

          comp[i] = 1. * __real__ comp[i] / (n_x * n_y);
          __imag__ comp[i] = 0.0;

          c[i] = __real__ comp[i];
          eta_r[i] = __real__ eta[i];
        }
      }

      //Calculate the Rx and Ry of the central phase
      double *R = findRadius(c, n_x, n_y, delta_x, delta_y);

      //Calculate (Rx^2 + Ry^2) / 2
      pr_t[INDEX] = (R[0] * R[0] / 4.0 + R[1] * R[1] / 4.0) / 2.0;

      /*
      fprintf(fpw1, "%le %le %le %le %le %le %le\n", INDEX * delta_t_M0,
                                            (R[2] - R[3])/2,
                                            (R[4] - R[5])/2,
                                            (pr_t[INDEX] - r0_new * r0_new),
                                            (pr_t[INDEX] - r0_new * r0_new) / (delta_x * delta_x),
                                            pr_t[INDEX]);
      */
      fprintf(fpw1, "%le %le %le %le %le %le\n", INDEX * delta_t_M0,
                                            pr_t[INDEX],
                                            R[2], R[3], R[4], R[5]);

      //Find number of Blocks inside the precipitate
      r_counter = 0;
      for (i1 = 0; i1 < n_x; ++i1) {
        for (i2 = 0; i2 < n_y; ++i2) {
          i = i2 + n_y * i1;

          if (c[i] >= 0.5) {
            r_counter = r_counter + 1;
          }
        }
      }

      //Calculate R^2 from the number blocks inside the precipitate
      pr_ta[INDEX] = 1. * r_counter * (delta_x * delta_y) / M_PI;

      /*
      fprintf(fpw2, "%le %le %le\n", INDEX * delta_t_M0,
              (pr_ta[INDEX] - r0_new * r0_new),
              1. * (pr_ta[INDEX] - r0_new * r0_new) / (delta_x * delta_x));
      */
      fprintf(fpw2, "%le %le\n", INDEX * delta_t_M0, pr_ta[INDEX]);

      //Output the state data for future analysis after every TIME_STEP_INITIAL
      if (INDEX % TIME_STEP_INITIAL == 0) {

        // copy composition in c and Print composition in file
        sprintf(NAME, "%s/data/c-time_%09d.dat", output_folder_name, INDEX);
        fpw = fopen(NAME, "w");
        for (i1 = 0; i1 < n_x; ++i1) {
          for (i2 = 0; i2 < n_y; ++i2) {
            i = i2 + n_y * i1;
            c[i] = __real__ comp[i];

            fprintf(fpw, "%d %d %.10le\n", i1, i2, c[i]);
          }
          fprintf(fpw, "\n");
        }
        fclose(fpw);
        //Fin

        //Copy eta in eta_r and Print eta in a file
        sprintf(NAME, "%s/data/eta-time_%09d.dat",output_folder_name, INDEX);
        fpw = fopen(NAME, "w");
        for (i1 = 0; i1 < n_x; ++i1) {
          for (i2 = 0; i2 < n_y; ++i2) {
            i = i2 + n_y * i1;
            eta_r[i] = __real__ eta[i];

            fprintf(fpw, "%d %d %.10le\n", i1, i2, eta_r[i]);
          }
          fprintf(fpw, "\n");
        }
        fclose(fpw);
        //Fin

        //printf("%d %.10le %.10le %.10le\n", INDEX, pr_t[INDEX], r0_new * r0_new,
        //       pr_t[INDEX] - r0_new * r0_new);

        printf("%d %le\n", INDEX, pr_ta[INDEX]);

      }

    }
    fclose(fpw1);
    fclose(fpw2);

    t_start_updated = time_steps_M0 + 1;
  }
  //If continuation was ON
  else if (cont == 1) {
    sprintf(NAME, "input/data/c-time_%09d.dat", INDEX_start);
    if ((fpr = fopen(NAME, "r")) == NULL) {
      printf("Unable to open the data file to read.\n");
      printf(
          "Continuation is not possible. Exiting from the following file %s\n",
          NAME);
      exit(0);
    } else {
      fpr = fopen(NAME, "r");
    }

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;

        fscanf(fpr, "%d %d %le\n", &ic[i], &jc[i], &c[i]);
      }
      fscanf(fpr, "\n");
    }
    fclose(fpr);

    sprintf(NAME, "input/data/eta-time_%09d.dat", INDEX_start);
    if ((fpr = fopen(NAME, "r")) == NULL) {
      printf("Unable to open the data file to read.\n");
      printf(
          "Continuation is not possible. Exiting from the following file %s\n",
          NAME);
      exit(0);
    } else {
      fpr = fopen(NAME, "r");
    }
    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;

        fscanf(fpr, "%d %d %le\n", &ic[i], &jc[i], &eta_r[i]);
      }
      fscanf(fpr, "\n");
    }
    fclose(fpr);

    ave_comp = 0.0;

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {

        i = i2 + n_y * i1;

        __real__ comp[i] = c[i];
        __imag__ comp[i] = 0.0;

        __real__ eta[i] = eta_r[i];
        __imag__ eta[i] = 0.0;

        ave_comp = ave_comp + creal(comp[i]);
      }
    }
    ave_comp = 1. * ave_comp / (n_x * n_y);
    if (INDEX % 50 == 0) {
      printf("Initial Average composition = %lf", ave_comp);
    }

    sprintf(NAME, "%s/initial.dat", output_folder_name);
    fpw = fopen(NAME, "w");
    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;
        fprintf(fpw, "%d %d %le %le\n", i1, i2, c[i], eta_r[i]);
      }
    }
    fclose(fpw);

    t_start_updated = INDEX_start + 1;
  }


  /************************************************************/
  /*               Variable Mobility : MAIN TIME LOOP         */
  /************************************************************/

  sprintf(NAME, "%s/r2_t-varM.dat", output_folder_name);
  fpw1 = fopen(NAME, "w");
  sprintf(NAME, "%s/mod_r2_t-varM.dat", output_folder_name);
  fpw2 = fopen(NAME, "w");

  for (INDEX = t_start_updated; INDEX < time_steps + 1; ++INDEX) {

    ave_comp = 0.0;

    /*********       g(c) Calculation       **********/
    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;

        comp[i] = comp[i] + 0.0 * _Complex_I;
        eta[i] = eta[i] + 0.0 * _Complex_I;

        c[i] = __real__ comp[i];
        eta_r[i] = __real__ eta[i];

        fm_c[i] = A * c[i] * c[i];
        fp_c[i] = B * (1 - c[i]) * (1 - c[i]);

        if (eta_r[i] < 0.0) {
          W_eta[i] = 0.0;
          W_eta1[i] = 0.0;
        }

        else if (eta_r[i] > 1.0) {
          W_eta[i] = 1.0;
          W_eta1[i] = 0.0;
        }

        else if (eta_r[i] >= 0.0 && eta_r[i] <= 1.0) {
          W_eta[i] = eta_r[i] * eta_r[i] * eta_r[i] *
                     (10 - 15 * eta_r[i] + 6 * eta_r[i] * eta_r[i]);

          W_eta1[i] = 3 * eta_r[i] * eta_r[i] *
                          (10 - 15 * eta_r[i] + 6 * eta_r[i] * eta_r[i]) +
                      eta_r[i] * eta_r[i] * eta_r[i] * (12 * eta_r[i] - 15);
        }

        __real__ g_eta[i] =
            fm_c[i] * (-W_eta1[i]) + fp_c[i] * W_eta1[i] +
            2.0 * P * eta_r[i] * (1.0 - eta_r[i]) * (1.0 - 2.0 * eta_r[i]) +
            0.0 * _Complex_I;
        __imag__ g_eta[i] = 0.0;
        __real__ g_c[i] = 2.0 * A * c[i] * (1 - W_eta[i]) -
                          2.0 * B * (1 - c[i]) * W_eta[i] + 0.0 * _Complex_I;
        __imag__ g_c[i] = 0.0;
        //  g_c[i] = 2.0*A*c[i]*(1.0-c[i])*(1.0-2.0*c[i]) + 0.0*_Complex_I;
        ave_comp = ave_comp + creal(comp[i]);
      }
    }

    ave_comp = 1. * ave_comp / (n_x * n_y);

    //Values gone OUT OF RANGE
    //TERMINATE PROGRAM
    if(ave_comp >= 1){
      printf("Oops Explosion!!");
      exit(0);
    }

    if (INDEX % 50 == 0) {
      printf("%d  Ave_comp  %lf\n", INDEX, ave_comp);
    }

    /*********       d2c(2nd gradient of Composition) Calculation        **********/

    fftw_execute_dft(plan1, comp, comp);
    fftw_execute_dft(plan1, eta, eta);

    for (i1 = 0; i1 < n_x; ++i1) {
      if (i1 < half_nx)
        kx = i1 * delta_kx;
      else
        kx = (i1 - n_x) * delta_kx;
      kx2 = kx * kx;

      for (i2 = 0; i2 < n_y; ++i2) {
        if (i2 < half_ny)
          ky = i2 * delta_ky;
        else
          ky = (i2 - n_y) * delta_ky;
        ky2 = ky * ky;

        k2 = kx2 + ky2;
        k3 = k2 * (kx + ky);
        k4 = k2 * k2;

        i = i2 + n_y * i1;

        d2c[i] = -2. * comp[i] * (kappa_c * k2);
      }
    }

    /*********  Backward FFTW of d2c and mu (h here) Calculation    **********/
    fftw_execute_dft(plan2, d2c, d2c);

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {

        i = i2 + n_y * i1;

        __real__ comp[i] = c[i];
        __imag__ comp[i] = 0.0;

        __real__ eta[i] = eta_r[i];
        __imag__ eta[i] = 0.0;
      }
    }

    /*
    sprintf(NAME, "%s/initial.dat", output_folder_name);
    fpw = fopen(NAME, "w");
    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;
        fprintf(fpw, "%d %d %le %le\n", i1, i2, c[i], eta_r[i]);
      }
    }
    fclose(fpw);
    */

    t_start_updated = INDEX_start + 1;
  }

  /************************************************************/
  /*                  Variable L : time loop                  */
  /************************************************************/

  sprintf(NAME, "%s/r2_t-varM.dat", output_folder_name);
  fpw1 = fopen(NAME, "w");
  sprintf(NAME, "%s/mod_r2_t-varM.dat", output_folder_name);
  fpw2 = fopen(NAME, "w");

  for (INDEX = t_start_updated; INDEX < time_steps + 1; ++INDEX) {

    ave_comp = 0.0;

    /*********       g(c) Calculation       **********/
    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;

        comp[i] = comp[i] + 0.0 * _Complex_I;
        eta[i] = eta[i] + 0.0 * _Complex_I;

        c[i] = __real__ comp[i];
        eta_r[i] = __real__ eta[i];

        fm_c[i] = A * c[i] * c[i];
        fp_c[i] = B * (1 - c[i]) * (1 - c[i]);

        if (eta_r[i] < 0.0) {
          W_eta[i] = 0.0;
          W_eta1[i] = 0.0;
        }

        else if (eta_r[i] > 1.0) {
          W_eta[i] = 1.0;
          W_eta1[i] = 0.0;
        }

        else if (eta_r[i] >= 0.0 && eta_r[i] <= 1.0) {
          W_eta[i] = eta_r[i] * eta_r[i] * eta_r[i] *
                     (10 - 15 * eta_r[i] + 6 * eta_r[i] * eta_r[i]);

          W_eta1[i] = 3 * eta_r[i] * eta_r[i] *
                          (10 - 15 * eta_r[i] + 6 * eta_r[i] * eta_r[i]) +
                      eta_r[i] * eta_r[i] * eta_r[i] * (12 * eta_r[i] - 15);
        }

        __real__ g_eta[i] =
            fm_c[i] * (-W_eta1[i]) + fp_c[i] * W_eta1[i] +
            2.0 * P * eta_r[i] * (1.0 - eta_r[i]) * (1.0 - 2.0 * eta_r[i]) +
            0.0 * _Complex_I;
        __imag__ g_eta[i] = 0.0;
        __real__ g_c[i] = 2.0 * A * c[i] * (1 - W_eta[i]) -
                          2.0 * B * (1 - c[i]) * W_eta[i] + 0.0 * _Complex_I;
        __imag__ g_c[i] = 0.0;
        //  g_c[i] = 2.0*A*c[i]*(1.0-c[i])*(1.0-2.0*c[i]) + 0.0*_Complex_I;
        ave_comp = ave_comp + creal(comp[i]);
      }
    }

    ave_comp = 1. * ave_comp / (n_x * n_y);

    if (INDEX % 50 == 0) {
      printf("%d  Ave_comp  %lf\n", INDEX, ave_comp);
    }

    /*********       d2c Calculation        **********/

    fftw_execute_dft(plan1, comp, comp);
    fftw_execute_dft(plan1, eta, eta);

    for (i1 = 0; i1 < n_x; ++i1) {
      if (i1 < half_nx)
        kx = i1 * delta_kx;
      else
        kx = (i1 - n_x) * delta_kx;
      kx2 = kx * kx;

      for (i2 = 0; i2 < n_y; ++i2) {
        if (i2 < half_ny)
          ky = i2 * delta_ky;
        else
          ky = (i2 - n_y) * delta_ky;
        ky2 = ky * ky;

        k2 = kx2 + ky2;
        k3 = k2 * (kx + ky);
        k4 = k2 * k2;

        i = i2 + n_y * i1;

        d2c[i] = -2. * comp[i] * (kappa_c * k2);
      }
    }

    /*********  Backward FFTW of d2c and mu (h here) Calculation    **********/
    fftw_execute_dft(plan2, d2c, d2c);

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {

        i = i2 + n_y * i1;

        d2c[i] = 1. * __real__ d2c[i] / (n_x * n_y);
        __imag__ d2c[i] = 0.0;

        h[i] = (g_c[i] - d2c[i]) + _Complex_I * 0.0; // mu calculation
      }
    }

    /*********  del_mu Calculation (it would be a vector hx and hy) **********/

    fftw_execute_dft(plan1, g_c, g_c);

    for (i1 = 0; i1 < n_x; ++i1) {
      if (i1 < half_nx)
        kx = i1 * delta_kx;
      else
        kx = (i1 - n_x) * delta_kx;
      kx2 = kx * kx;
      kx4 = kx2 * kx * kx;
      kx6 = kx4 * kx * kx;

      for (i2 = 0; i2 < n_y; ++i2) {
        if (i2 < half_ny)
          ky = i2 * delta_ky;
        else
          ky = (i2 - n_y) * delta_ky;
        ky2 = ky * ky;
        ky4 = ky2 * ky * ky;
        ky6 = ky4 * ky * ky;

        k2 = kx2 + ky2;
        k3 = k2 * (kx + ky);
        k4 = k2 * k2;

        i = i2 + n_y * i1;

        hx[i] = _Complex_I * kx *
                (g_c[i] +
                 comp[i] * (2. * kappa_c * k2

                            + 2. * P16 * kx6 + P26 * ky2 * kx4

                            + 2. * P27 * ky6 + P28 * kx2 * ky4

                            + 2. * P21 * ky2 * kx4 + P28 * kx2 * ky4

                            + 2. * P17 * kx2 * ky4 + P26 * ky2 * kx4));

        hy[i] = _Complex_I * ky *
                (g_c[i] +
                 comp[i] * (2. * kappa_c * k2

                            + 2. * P16 * kx6 + P26 * ky2 * kx4

                            + 2. * P27 * ky6 + P28 * kx2 * ky4

                            + 2. * P21 * ky2 * kx4 + P28 * kx2 * ky4

                            + 2. * P17 * kx2 * ky4 + P26 * ky2 * kx4));
      }
    }

    fftw_execute_dft(plan2, hx, hx);
    fftw_execute_dft(plan2, hy, hy);

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {

        i = i2 + n_y * i1;

        hx[i] = 1. * __real__ hx[i] / (n_x * n_y);
        __imag__ hx[i] = 0.0;

        hy[i] = 1. * __real__ hy[i] / (n_x * n_y);
        __imag__ hy[i] = 0.0;
      }
    }

    /*********          nx and ny Calculation           **********/

    for (i1 = 0; i1 < n_x; ++i1) {
      if (i1 < half_nx)
        kx = i1 * delta_kx;
      else
        kx = (i1 - n_x) * delta_kx;
      kx2 = kx * kx;

      for (i2 = 0; i2 < n_y; ++i2) {
        if (i2 < half_ny)
          ky = i2 * delta_ky;
        else
          ky = (i2 - n_y) * delta_ky;
        ky2 = ky * ky;

        k2 = kx2 + ky2;
        k4 = k2 * k2;

        i = i2 + n_y * i1;

        nx[i] = _Complex_I * kx * comp[i];
        ny[i] = _Complex_I * ky * comp[i];

        nex[i] = _Complex_I * kx * eta[i];
        ney[i] = _Complex_I * ky * eta[i];

        nxx[i] = -kx2 * comp[i];
        nyy[i] = -ky2 * comp[i];
      }
    }

    fftw_execute_dft(plan2, nx, nx);
    fftw_execute_dft(plan2, ny, ny);
    fftw_execute_dft(plan2, nex, nex);
    fftw_execute_dft(plan2, ney, ney);
    fftw_execute_dft(plan2, nxx, nxx);
    fftw_execute_dft(plan2, nyy, nyy);

    fftw_execute_dft(plan2, comp, comp);
    fftw_execute_dft(plan2, eta, eta);

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {

        i = i2 + n_y * i1;

        nx[i] = 1. * __real__ nx[i] / (n_x * n_y);
        ny[i] = 1. * __real__ ny[i] / (n_x * n_y);

        __imag__ nx[i] = 0.0;
        __imag__ ny[i] = 0.0;

        nex[i] = 1. * __real__ nex[i] / (n_x * n_y);
        ney[i] = 1. * __real__ ney[i] / (n_x * n_y);

        __imag__ nex[i] = 0.0;
        __imag__ ney[i] = 0.0;

        nxx[i] = 1. * __real__ nxx[i] / (n_x * n_y);
        nyy[i] = 1. * __real__ nyy[i] / (n_x * n_y);

        __imag__ nxx[i] = 0.0;
        __imag__ nyy[i] = 0.0;

        comp[i] = 1. * __real__ comp[i] / (n_x * n_y);
        __imag__ comp[i] = 0.0;

        eta[i] = 1. * __real__ eta[i] / (n_x * n_y);
        __imag__ eta[i] = 0.0;
      }
    }

    /*********          Mobiity (M) and q = M*del_mu Calculation **********/

    mod_cut = 1.0e-2;
    mod_limit = 1.0e-4;

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;

        mod[i] = (nx[i] * nx[i] + ny[i] * ny[i]);
        mod_e[i] = (nex[i] * nex[i] + ney[i] * ney[i]);

        // M[i] = M0;
        L[i] = L0;

        if (mod[i] > mod_cut) {
          M[i] = M0 + beta11 +
                 1. * (beta12 * (nx[i] * nx[i]) *
                       (nx[i] * nx[i] - 3. * ny[i] * ny[i]) *
                       (nx[i] * nx[i] - 3. * ny[i] * ny[i]) /
                       (mod[i] * mod[i] * mod[i]));
        }

        else if (mod[i] <= mod_cut) {
          if (mod[i] > mod_limit && mod[i] <= mod_cut) {

            M[i] = M0 + beta11 +
                   1. * (beta12 * (nx[i] * nx[i]) *
                         (nx[i] * nx[i] - 3. * ny[i] * ny[i]) *
                         (nx[i] * nx[i] - 3. * ny[i] * ny[i]) /
                         (mod[i] * mod[i] * mod[i]));

            maxM = 1. * mod[i] / (mod_cut);

            if (maxM < 0) {
              wang_f[i] = 0.0;
            }

            else if (maxM > 1.0) {
              wang_f[i] = 1.0;
            }

            else if (maxM >= 0 && maxM <= 1.0) {
              wang_f[i] =
                  maxM * maxM * maxM * (10 - 15 * maxM + 6 * maxM * maxM);
            }
            M[i] = M0 + wang_f[i] * (M[i] - M0);
          }

          else if (mod[i] <= mod_limit || __real__ comp[i] > 1.0) {
            M[i] = M0;
          }
        }

        if (mod_e[i] > mod_cut) {

          L[i] = L0 + beta11 +
                 1. * (beta12 * (nx[i] * nx[i]) *
                       (nx[i] * nx[i] - 3. * ny[i] * ny[i]) *
                       (nx[i] * nx[i] - 3. * ny[i] * ny[i]) /
                       (mod[i] * mod[i] * mod[i]));
        }

        qx[i] = M[i] + _Complex_I * 0.0;
        qy[i] = M[i] + _Complex_I * 0.0;

        M_mod[i] = M[i] + _Complex_I * 0.0;
      }
    }

    fftw_execute_dft(plan1, qx, qx);
    fftw_execute_dft(plan1, qy, qy);

    for (i1 = 0; i1 < n_x; ++i1) {
      if (i1 < half_nx)
        kx = i1 * delta_kx;
      else
        kx = (i1 - n_x) * delta_kx;
      kx2 = kx * kx;

      for (i2 = 0; i2 < n_y; ++i2) {
        if (i2 < half_ny)
          ky = i2 * delta_ky;
        else
          ky = (i2 - n_y) * delta_ky;
        ky2 = ky * ky;

        k2 = kx2 + ky2;
        k4 = k2 * k2;

        i = i2 + n_y * i1;

        qx[i] = _Complex_I * kx * qx[i];
        qy[i] = _Complex_I * ky * qy[i];
      }
    }

    fftw_execute_dft(plan2, qx, qx);
    fftw_execute_dft(plan2, qy, qy);

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {

        i = i2 + n_y * i1;

        qx[i] = 1. * __real__ qx[i] / (n_x * n_y);
        __imag__ qx[i] = 0.0;

        qy[i] = 1. * __real__ qy[i] / (n_x * n_y);
        __imag__ qy[i] = 0.0;

        qx[i] = hx[i] * qx[i];

        i = i2 + n_y * i1;

        d2c[i] = 1. * __real__ d2c[i] / (n_x * n_y);
        __imag__ d2c[i] = 0.0;

        h[i] = (g_c[i] - d2c[i]) + _Complex_I * 0.0; // mu calculation
      }
    }

    /*********  del_mu Calculation (it would be a vector hx and hy) **********/

    fftw_execute_dft(plan1, g_c, g_c);

    for (i1 = 0; i1 < n_x; ++i1) {
      if (i1 < half_nx)
        kx = i1 * delta_kx;
      else
        kx = (i1 - n_x) * delta_kx;
      kx2 = kx * kx;
      kx4 = kx2 * kx * kx;
      kx6 = kx4 * kx * kx;

      for (i2 = 0; i2 < n_y; ++i2) {
        if (i2 < half_ny)
          ky = i2 * delta_ky;
        else
          ky = (i2 - n_y) * delta_ky;
        ky2 = ky * ky;
        ky4 = ky2 * ky * ky;
        ky6 = ky4 * ky * ky;

        k2 = kx2 + ky2;
        k3 = k2 * (kx + ky);
        k4 = k2 * k2;

        i = i2 + n_y * i1;

        hx[i] = _Complex_I * kx *
                (g_c[i] +
                 comp[i] * (2. * kappa_c * k2

                            + 2. * P16 * kx6 + P26 * ky2 * kx4

                            + 2. * P27 * ky6 + P28 * kx2 * ky4

                            + 2. * P21 * ky2 * kx4 + P28 * kx2 * ky4

                            + 2. * P17 * kx2 * ky4 + P26 * ky2 * kx4));

        hy[i] = _Complex_I * ky *
                (g_c[i] +
                 comp[i] * (2. * kappa_c * k2

                            + 2. * P16 * kx6 + P26 * ky2 * kx4

                            + 2. * P27 * ky6 + P28 * kx2 * ky4

                            + 2. * P21 * ky2 * kx4 + P28 * kx2 * ky4

                            + 2. * P17 * kx2 * ky4 + P26 * ky2 * kx4));
      }
    }

    fftw_execute_dft(plan2, hx, hx);
    fftw_execute_dft(plan2, hy, hy);

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {

        i = i2 + n_y * i1;

        hx[i] = 1. * __real__ hx[i] / (n_x * n_y);
        __imag__ hx[i] = 0.0;

        hy[i] = 1. * __real__ hy[i] / (n_x * n_y);
        __imag__ hy[i] = 0.0;
      }
    }

    /*********          nx and ny Calculation           **********/

    for (i1 = 0; i1 < n_x; ++i1) {
      if (i1 < half_nx)
        kx = i1 * delta_kx;
      else
        kx = (i1 - n_x) * delta_kx;
      kx2 = kx * kx;

      for (i2 = 0; i2 < n_y; ++i2) {
        if (i2 < half_ny)
          ky = i2 * delta_ky;
        else
          ky = (i2 - n_y) * delta_ky;
        ky2 = ky * ky;

        k2 = kx2 + ky2;
        k4 = k2 * k2;

        i = i2 + n_y * i1;

        nx[i] = _Complex_I * kx * comp[i];
        ny[i] = _Complex_I * ky * comp[i];

        nex[i] = _Complex_I * kx * eta[i];
        ney[i] = _Complex_I * ky * eta[i];

        nxx[i] = -kx2 * comp[i];
        nyy[i] = -ky2 * comp[i];
      }
    }

    fftw_execute_dft(plan2, nx, nx);
    fftw_execute_dft(plan2, ny, ny);
    fftw_execute_dft(plan2, nex, nex);
    fftw_execute_dft(plan2, ney, ney);
    fftw_execute_dft(plan2, nxx, nxx);
    fftw_execute_dft(plan2, nyy, nyy);

    fftw_execute_dft(plan2, comp, comp);
    fftw_execute_dft(plan2, eta, eta);

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {

        i = i2 + n_y * i1;

        nx[i] = 1. * __real__ nx[i] / (n_x * n_y);
        ny[i] = 1. * __real__ ny[i] / (n_x * n_y);

        __imag__ nx[i] = 0.0;
        __imag__ ny[i] = 0.0;

        nex[i] = 1. * __real__ nex[i] / (n_x * n_y);
        ney[i] = 1. * __real__ ney[i] / (n_x * n_y);

        __imag__ nex[i] = 0.0;
        __imag__ ney[i] = 0.0;

        nxx[i] = 1. * __real__ nxx[i] / (n_x * n_y);
        nyy[i] = 1. * __real__ nyy[i] / (n_x * n_y);

        __imag__ nxx[i] = 0.0;
        __imag__ nyy[i] = 0.0;

        comp[i] = 1. * __real__ comp[i] / (n_x * n_y);
        __imag__ comp[i] = 0.0;

        eta[i] = 1. * __real__ eta[i] / (n_x * n_y);
        __imag__ eta[i] = 0.0;
      }
    }

    /*********          Mobiity (M) and q = M*del_mu Calculation **********/

    mod_cut = 1.0e-2;
    mod_limit = 1.0e-4;

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;

        mod[i] = (nx[i] * nx[i] + ny[i] * ny[i]);
        mod_e[i] = (nex[i] * nex[i] + ney[i] * ney[i]);

        // M[i] = M0;
        L[i] = L0;

        if (mod[i] > mod_cut) {
          M[i] = M0 + beta11 +
                 1. * (beta12 * (nx[i] * nx[i]) *
                       (nx[i] * nx[i] - 3. * ny[i] * ny[i]) *
                       (nx[i] * nx[i] - 3. * ny[i] * ny[i]) /
                       (mod[i] * mod[i] * mod[i]));
        }

        else if (mod[i] <= mod_cut) {
          if (mod[i] > mod_limit && mod[i] <= mod_cut) {

            M[i] = M0 + beta11 +
                   1. * (beta12 * (nx[i] * nx[i]) *
                         (nx[i] * nx[i] - 3. * ny[i] * ny[i]) *
                         (nx[i] * nx[i] - 3. * ny[i] * ny[i]) /
                         (mod[i] * mod[i] * mod[i]));

            maxM = 1. * mod[i] / (mod_cut);

            if (maxM < 0) {
              wang_f[i] = 0.0;
            }

            else if (maxM > 1.0) {
              wang_f[i] = 1.0;
            }

            else if (maxM >= 0 && maxM <= 1.0) {
              wang_f[i] =
                  maxM * maxM * maxM * (10 - 15 * maxM + 6 * maxM * maxM);
            }
            M[i] = M0 + wang_f[i] * (M[i] - M0);
          }

          else if (mod[i] <= mod_limit || __real__ comp[i] > 1.0) {
            M[i] = M0;
          }
        }

        if (mod_e[i] > mod_cut) {

          L[i] = L0 + beta11 +
                 1. * (beta12 * (nx[i] * nx[i]) *
                       (nx[i] * nx[i] - 3. * ny[i] * ny[i]) *
                       (nx[i] * nx[i] - 3. * ny[i] * ny[i]) /
                       (mod[i] * mod[i] * mod[i]));
        }

        qx[i] = M[i] + _Complex_I * 0.0;
        qy[i] = M[i] + _Complex_I * 0.0;

        M_mod[i] = M[i] + _Complex_I * 0.0;
      }
    }

    fftw_execute_dft(plan1, qx, qx);
    fftw_execute_dft(plan1, qy, qy);

    for (i1 = 0; i1 < n_x; ++i1) {
      if (i1 < half_nx)
        kx = i1 * delta_kx;
      else
        kx = (i1 - n_x) * delta_kx;
      kx2 = kx * kx;

      for (i2 = 0; i2 < n_y; ++i2) {
        if (i2 < half_ny)
          ky = i2 * delta_ky;
        else
          ky = (i2 - n_y) * delta_ky;
        ky2 = ky * ky;

        k2 = kx2 + ky2;
        k4 = k2 * k2;

        i = i2 + n_y * i1;

        qx[i] = _Complex_I * kx * qx[i];
        qy[i] = _Complex_I * ky * qy[i];
      }
    }

    fftw_execute_dft(plan2, qx, qx);
    fftw_execute_dft(plan2, qy, qy);

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {

        i = i2 + n_y * i1;

        qx[i] = 1. * __real__ qx[i] / (n_x * n_y);
        __imag__ qx[i] = 0.0;

        qy[i] = 1. * __real__ qy[i] / (n_x * n_y);
        __imag__ qy[i] = 0.0;

        qx[i] = hx[i] * qx[i] + _Complex_I * 0.0; // mu calculation
        qy[i] = hy[i] * qy[i] + _Complex_I * 0.0;
      }
    }

    /*************************************************************************/
    /*********      Evolution Equation          **********/
    /*************************************************************************/

    fftw_execute_dft(plan1, comp, comp);
    fftw_execute_dft(plan1, M_mod, M_mod);
    fftw_execute_dft(plan1, qx, qx);
    fftw_execute_dft(plan1, qy, qy);

    fftw_execute_dft(plan1, eta, eta);
    fftw_execute_dft(plan1, g_eta, g_eta);

    for (i1 = 0; i1 < n_x; ++i1) {
      if (i1 < half_nx)
        kx = i1 * delta_kx;
      else
        kx = (i1 - n_x) * delta_kx;
      kx2 = kx * kx;

      for (i2 = 0; i2 < n_y; ++i2) {
        if (i2 < half_ny)
          ky = i2 * delta_ky;
        else
          ky = (i2 - n_y) * delta_ky;
        ky2 = ky * ky;

        k2 = kx2 + ky2;
        k4 = k2 * k2;

        i = i2 + n_y * i1;

        denom_eta =
            1.0 + 2.0 * delta_t * 1.0 * kappa_eta * k2; // Instead of 1.0, L[i]
                                                        // will make variable
                                                        // relaxation parameter
        eta[i] =
            1.0 * (eta[i] - 1.0 * delta_t * g_eta[i]) /
            denom_eta; // Instead of 1.0, L[i] will make variable relaxation
                       // parameter

        factor = (kappa_c * k2);
        denom_c = 1.0 + 2.0 * delta_t * k2 * M[i] * factor;

        comp[i] = (comp[i] + delta_t * (qx[i] + qy[i] - k2 * M[i] * g_c[i])) /
                  denom_c;

        // comp[i] = comp[i] -_Complex_I*delta_t*(kx*qx[i] + ky*qy[i]);
      }
    }

    fftw_execute_dft(plan2, comp, comp);
    fftw_execute_dft(plan2, eta, eta);

    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;

        comp[i] = 1. * __real__ comp[i] / (n_x * n_y);
        __imag__ comp[i] = 0.0;

        eta[i] = 1. * __real__ eta[i] / (n_x * n_y);
        __imag__ eta[i] = 0.0;

        c[i] = __real__ comp[i];
        eta_r[i] = __real__ eta[i];
      }
    }

    /*********  x-radius Calculation   *************/

    for (i1 = (int)n_x / 2; i1 < n_x; ++i1) {

      i2 = (int)n_y / 2;

      if (__real__ comp[i2 + n_y * i1] < 0.5) {
        pr_x1 = (i1 - 1) * delta_x +
                1.0 * (i1 * delta_x - (i1 - 1) * delta_x) *
                    (0.5 - __real__ comp[i2 + n_y * (i1 - 1)]) /
                    (__real__ comp[i2 + n_y * i1] -
                     __real__ comp[i2 + n_y * (i1 - 1)]);

        break;
      }
    }

    r_counter = 0;
    for (i1 = (int)n_x / 2; i1 < n_x; ++i1) {

      i2 = (int)n_y / 2;

      r_counter = r_counter + 1;
      i_new = (int)n_x / 2 - r_counter;

      if (__real__ comp[i2 + n_y * (i_new)] < 0.5) {
        pr_x2 = ((i_new + 1) * delta_x +
                 1.0 * ((i_new)*delta_x - (i_new + 1) * delta_x) *
                     (0.5 - __real__ comp[i2 + n_y * (i_new + 1)]) /
                     (__real__ comp[i2 + n_y * (i_new)] -
                      __real__ comp[i2 + n_y * ((i_new) + 1)]));

        break;
      }
    }

    for (i1 = (int)n_x / 2; i1 < n_x; ++i1) {
      for (i2 = (int)n_y / 2; i2 < n_y; ++i2) {
        if (i1 == i2 && __real__ comp[i2 + n_y * i1] < 0.5) {

          pr_xy1 = (i1 - 1) * delta_x +
                   1.0 * (i1 * delta_x - (i1 - 1) * delta_x) *
                   (0.5 - __real__ comp[i2 + n_y * (i1 - 1)]) /
                   (__real__ comp[i2 + n_y * i1] -
                    __real__ comp[i2 + n_y * (i1 - 1)]);
          break;
        }
      }
    }

    r_counter = 0;
    for (i1 = (int)n_x / 2; i1 < n_x; ++i1) {
      for (i2 = (int)n_y / 2; i2 < 0; --i2) {

        r_counter = r_counter + 1;
        i_new = (int)n_x / 2 - r_counter;

        if (i_new == i2 && __real__ comp[i2 + n_y * (i_new)] < 0.5) {
          pr_xy2 = ((i_new + 1) * delta_x +
                    1.0 * ((i_new)*delta_x - (i_new + 1) * delta_x) *
                        (0.5 - __real__ comp[i2 + n_y * (i_new + 1)]) /
                        (__real__ comp[i2 + n_y * (i_new)] -
                         __real__ comp[i2 + n_y * ((i_new) + 1)]));

          break;
        }
      }
    }

    /*********  y-radius Calculation   *************/

    for (i2 = (int)n_y / 2; i2 < n_y; ++i2) {

      i1 = (int)n_x / 2;
      i = i2 + n_y * i1;

      if (__real__ comp[i] < 0.5) {
        pr_y1 = (i2 - 1) * delta_y +
                1.0 * (i2 * delta_y - (i2 - 1) * delta_y) *
                    (0.5 - __real__ comp[i - 1]) /
                    (__real__ comp[i] - __real__ comp[i - 1]);

        break;
      }
    }

    r_counter = 0;

    for (i2 = (int)n_y / 2; i2 < n_y; ++i2) {

      i1 = (int)n_x / 2;
      i = i2 + n_y * i1;

      r_counter = r_counter + 1;
      i_new = (int)n_y / 2 - r_counter;

      if (__real__ comp[i_new + n_y * i1] < 0.5) {
        pr_y2 = ((i_new + 1) * delta_y +
                 1.0 * (i_new * delta_y - (i_new + 1) * delta_y) *
                     (0.5 - __real__ comp[(i_new + 1) + n_y * i1]) /
                     (__real__ comp[i_new + n_y * i1] -
                      __real__ comp[(i_new + 1) + n_y * i1]));

        break;
      }
    }

    /*********      Data Generation         **********/

    if (INDEX % file_timer == 0) {

      sprintf(NAME, "%s/Mx-time_%09d.dat",output_folder_name, INDEX);
      fpw = fopen(NAME, "w");
      for (i1 = 0; i1 < n_x; ++i1) {
        for (i2 = 0; i2 < n_y; ++i2) {
          i = i2 + n_y * i1;

          if (i2 == (int)n_y / 2) {

            fprintf(fpw, "%d %le %le %le\n", i1, c[i], eta_r[i], M[i]);
          }
        }
      }
      fclose(fpw);

      sprintf(NAME, "%s/data/M-time_%09d.dat",output_folder_name, INDEX);
      fpw = fopen(NAME, "w");
      for (i1 = 0; i1 < n_x; ++i1) {
        for (i2 = 0; i2 < n_y; ++i2) {
          i = i2 + n_y * i1;
          fprintf(fpw, "%d %d %le\n", i1, i2, M[i]);
        }
        fprintf(fpw, "\n");
      }
      fclose(fpw);

      sprintf(NAME, "%s/My-time_%09d.dat",output_folder_name, INDEX);
      fpw = fopen(NAME, "w");
      for (i1 = 0; i1 < n_x; ++i1) {
        for (i2 = 0; i2 < n_y; ++i2) {
          i = i2 + n_y * i1;
          if (i1 == (int)n_x / 2) {
            fprintf(fpw, "%d %le %le %le\n", i2, c[i], eta_r[i], M[i]);
          }
        }
      }
      fclose(fpw);

      sprintf(NAME, "%s/data/c-time_%09d.dat",output_folder_name, INDEX);
      fpw = fopen(NAME, "w");
      for (i1 = 0; i1 < n_x; ++i1) {
        for (i2 = 0; i2 < n_y; ++i2) {
          i = i2 + n_y * i1;
          c[i] = __real__ comp[i];
          eta_r[i] = __real__ eta[i];

          fprintf(fpw, "%d %d %le %le\n", i1, i2, c[i], eta_r[i]);
        }
        fprintf(fpw, "\n");
      }
      fclose(fpw);

      sprintf(NAME, "%s/preci-rx_%09d.dat",output_folder_name, INDEX);
      fpw = fopen(NAME, "w");
      for (i1 = 0; i1 < n_x; ++i1) {
        for (i2 = 0; i2 < n_y; ++i2) {
          i = i2 + n_y * i1;

          /*********      Calculation of precipitates radius **********/

          if (i2 == (int)n_y / 2) {

            fprintf(fpw, "%le\n", __real__ comp[i]);
          }
        }
      }
      fclose(fpw);
    }

    pr_t[INDEX] = 1. * (1. * (pr_x1 - pr_x2) * (pr_x1 - pr_x2) / 4.0 +
                        1. * (pr_y1 - pr_y2) * (pr_y1 - pr_y2) / 4.0) /
                  2.0;
    aspect_ratio[INDEX] = sqrt(
        1. * ((pr_x1 - pr_x2) * (pr_x1 - pr_x2) / 4.0) /
        (sqrt(2.0) * (pr_xy1 - pr_xy2) * sqrt(2.0) * (pr_xy1 - pr_xy2) / 4.0));

    /*fprintf(fpw1, "%le %le %le %le %le\n", INDEX * delta_t,
            (pr_t[INDEX] - r0_new * r0_new),
            1. * (pr_t[INDEX] - r0_new * r0_new) / (delta_x * delta_x),
            pr_t[INDEX], aspect_ratio[INDEX]);
            */
    fprintf(fpw1, "%le %le %le %le %le %le %le %le\n", INDEX * delta_t,
                                          (pr_x1 - pr_x2)/2,
                                          (pr_y1 - pr_y2)/2,
                                          (double)(pr_y1 - pr_y2)/(pr_y1 - pr_y2),
                                          (pr_t[INDEX] - r0_new * r0_new),
                                          1. * (pr_t[INDEX] - r0_new * r0_new) / (delta_x * delta_x),
                                          pr_t[INDEX],
                                          aspect_ratio[INDEX]);

    r_counter = 0;
    for (i1 = 0; i1 < n_x; ++i1) {
      for (i2 = 0; i2 < n_y; ++i2) {
        i = i2 + n_y * i1;

        if (c[i] >= 0.5) {
          r_counter = r_counter + 1;
        }
      }
    }

    pr_ta[INDEX] = 1. * r_counter * (delta_x * delta_y) / M_PI;

    fprintf(fpw2, "%le %le %le\n", INDEX * delta_t,
            (pr_ta[INDEX] - r0_new * r0_new),
            1. * (pr_ta[INDEX] - r0_new * r0_new) / (delta_x * delta_x));
  }
  fclose(fpw1);
  fclose(fpw2);

  sprintf(NAME, "%s/alpha_t.dat", output_folder_name);
  fpw1 = fopen(NAME, "w");

  for (INDEX = 1; INDEX < time_steps + 1; ++INDEX) {

    if (INDEX > (int)(time_steps * delta_t / 20) &&
        INDEX < (time_steps - (int)(time_steps * delta_t / 20))) {
      alpha[INDEX] =
          sqrt(1. * (pr_t[INDEX + (int)(time_steps * delta_t / 20)] -
                     pr_t[INDEX - (int)(time_steps * delta_t / 20)]) /
               (2. * A * time_steps * delta_t / 20));
      fprintf(fpw1, "%le %le\n", INDEX * delta_t, alpha[INDEX]);
    }
  }
  fclose(fpw1);

  fftw_free(comp);
  fftw_free(eta);
  gsl_rng_free(ran_num);
  fftw_free(g_c);
  fftw_free(g_eta);
  fftw_free(nx);
  fftw_free(ny);
  fftw_free(nex);
  fftw_free(ney);
  fftw_free(nxx);
  fftw_free(nyy);

  fftw_free(M_mod);
  fftw_free(d2c);
  fftw_free(h);
  fftw_free(hx);
  fftw_free(qx);
  fftw_free(hy);
  fftw_free(qy);

  free(fm_c);
  free(fp_c);
  free(wang_f);
  free(W_eta);
  free(W_eta1);
  free(c);
  free(eta_r);
  free(mod);
  free(mod_e);
  free(M);
  free(L);
  free(mx1);
  free(mx2);
  free(mx1_final);

  free(pr_t);
  free(pr_ta);
  free(aspect_ratio);
  free(alpha);
  free(zeta);

  free(ic);
  free(jc);

  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
}
