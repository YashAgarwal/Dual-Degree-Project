#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include "../headers/functions.h"

int main (void){

    setbuf(stdout, NULL);
    FILE *fpr, *fpw;
    int n_x, n_y;
    double delta_x, delta_y;
    double kappa_c, kappa_eta;
    double A, B;
    double M0, L0;
    double P;
    int time_steps, time_steps_M0, file_timer;
    double delta_t, delta_t_M0;
    int theta_int, theta_start, theta_steps;
    double c_zero, r_zero;
    double noise_str;

    double gamma_i, gamma_a;
    double beta11, beta22, beta12;
    double coeff;

    double P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15,
    P16, P17, P18, P19, P20, P21, P22, P23, P24;

    int cont, INDEX_start;

    int multi_output = 1;

    /***************************
    Preparing the Output directory
    ***************************/
    /*Generating the output folder name*/
    char output_folder_name[50];

    if(multi_output == 1){
        time_t t = time(NULL);
        struct tm tm = *localtime(&t);

        sprintf(output_folder_name, "output/%d_%d_%d/%d_%d_%d",
                tm.tm_mday,
                tm.tm_mon + 1,
                tm.tm_year + 1900,
                tm.tm_hour,
                tm.tm_min,
                tm.tm_sec);

        printf("Creating output in %s\n", output_folder_name);
        fflush(stdout);
    }
    else{
        (void) system ("rm -rf output/*");
        sprintf(output_folder_name, "output");
    }

    char str[100];
    sprintf(str, "mkdir -p %s/data", output_folder_name);
    (void) system (str);

    /*
    The input values used for any particular run is written on to the file
    "README" in the directory named "output/".  I open it here.
    */
    sprintf(str, "%s/README", output_folder_name);
    fpw = fopen (str, "w");
    if (fpw == NULL){
        printf ("Unable to open output/README. Exiting");
        exit (0);
    }

    /***************************
    Reading the required constants
    ***************************/

    /*system size and grid size*/
    fpr = fopen ("input/system_size", "r");
    if (fpr == NULL){
        printf ("Unable to open input/system_size. Exiting");
        exit (0);
    }

    (void) fscanf (fpr, "%d%d%le%le", &n_x, &n_y, &delta_x, &delta_y);
    (void) fclose (fpr);
    fprintf (fpw, "n_x = %d\n", n_x);
    fprintf (fpw, "n_y = %d\n", n_y);
    fprintf (fpw, "delta_x = %le\n", delta_x);
    fprintf (fpw, "delta_y = %le\n", delta_y);

    /*
    Chemical free energy related constants
    */
    fpr = fopen ("input/constants", "r");
    if (fpr == NULL)
    {
      printf ("Unable to open input/constants. Exiting");
      exit (0);
    }

    (void) fscanf (fpr, "%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le",
     &kappa_c, &kappa_eta, &M0, &L0, &A, &B, &P, &c_zero, &r_zero,
     &noise_str, &gamma_i, &gamma_a, &beta11, &beta22, &beta12,
     &coeff);
    (void) fclose (fpr);
    fprintf (fpw, "kappa_c = %le\n", kappa_c);
    fprintf (fpw, "kappa_eta = %le\n", kappa_eta);
    fprintf (fpw, "M0 = %le\n", M0);
    fprintf (fpw, "L0 = %le\n", L0);
    fprintf (fpw, "A = %le\n", A);
    fprintf (fpw, "B = %le\n", B);
    fprintf (fpw, "P = %le\n", P);
    fprintf (fpw, "c_zero = %le\n", c_zero);
    fprintf (fpw, "r_zero = %le\n", r_zero);
    fprintf (fpw, "noise_str = %le\n", noise_str);
    fprintf (fpw, "gamma_i = %le\n", gamma_i);
    fprintf (fpw, "gamma_a = %le\n", gamma_a);
    fprintf (fpw, "beta11 = %le\n", beta11);
    fprintf (fpw, "beta22 = %le\n", beta22);
    fprintf (fpw, "beta12 = %le\n", beta12);
    fprintf (fpw, "coeff = %le\n", coeff);

    /*
    time step and number of time steps
    */
    fpr = fopen ("input/time_information", "r");
    if (fpr == NULL)
    {
      printf ("Unable to open input/time_information. Exiting");
      exit (0);
    }

    (void) fscanf (fpr, "%le%d%le%d%d%d%d%d", &delta_t, &time_steps,
     &delta_t_M0, &time_steps_M0, &file_timer, &theta_int,
     &theta_start, &theta_steps);
    (void) fclose (fpr);
    fprintf (fpw, "delta_t = %le\n", delta_t);
    fprintf (fpw, "time_steps = %d\n", time_steps);
    fprintf (fpw, "delta_t_M0 = %le\n", delta_t_M0);
    fprintf (fpw, "time_steps_M0 = %d\n", time_steps_M0);
    fprintf (fpw, "file_timer = %d\n", file_timer);
    fprintf (fpw, "theta_interval = %d\n", theta_int);
    fprintf (fpw, "theta_start = %d\n", theta_start);
    fprintf (fpw, "theta_steps = %d\n", theta_steps);

    /*
    Let us read the continuation variable and INDEX_start
    */
    fpr = fopen ("input/continuation", "r");
    if (fpr == NULL)
    {
      printf ("Unable to open input/continuation. Exiting");
      exit (0);
    }

    (void) fscanf (fpr, "%d%d", &cont, &INDEX_start);
    (void) fclose (fpr);
    fprintf (fpw, "cont = %d\n", cont);
    fprintf (fpw, "INDEX_start = %d\n", INDEX_start);

    /*
    Let us read the interfacial parameters
    */
    fpr = fopen ("input/intf_param", "r");
    if (fpr == NULL)
    {
      printf ("Unable to oepn input/intf_param. Exiting");
      exit (0);
    }

    (void) fscanf (fpr,
     "%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le",
     &P1, &P2, &P3, &P4, &P5, &P6, &P7, &P8, &P9, &P10, &P11,
     &P12, &P13, &P14, &P15, &P16, &P17, &P18, &P19, &P20, &P21,
     &P22, &P23, &P24);
    (void) fclose (fpr);

    fprintf (fpw, "P1 = %le\n", P1);
    fprintf (fpw, "P2 = %le\n", P2);
    fprintf (fpw, "P3 = %le\n", P3);
    fprintf (fpw, "P4 = %le\n", P4);
    fprintf (fpw, "P5 = %le\n", P5);
    fprintf (fpw, "P6 = %le\n", P6);
    fprintf (fpw, "P7 = %le\n", P7);
    fprintf (fpw, "P8 = %le\n", P8);
    fprintf (fpw, "P9 = %le\n", P9);
    fprintf (fpw, "P10 = %le\n", P10);
    fprintf (fpw, "P11 = %le\n", P11);
    fprintf (fpw, "P12 = %le\n", P12);
    fprintf (fpw, "P13 = %le\n", P13);
    fprintf (fpw, "P14 = %le\n", P14);
    fprintf (fpw, "P15 = %le\n", P15);
    fprintf (fpw, "P16 = %le\n", P16);
    fprintf (fpw, "P17 = %le\n", P17);
    fprintf (fpw, "P18 = %le\n", P18);
    fprintf (fpw, "P19 = %le\n", P19);
    fprintf (fpw, "P20 = %le\n", P20);
    fprintf (fpw, "P21 = %le\n", P21);
    fprintf (fpw, "P22 = %le\n", P22);
    fprintf (fpw, "P23 = %le\n", P23);
    fprintf (fpw, "P24 = %le\n", P24);
    (void) fclose(fpw);

    //Start Simulation
    evolve (n_x, n_y, delta_x, delta_y, kappa_c, kappa_eta, P1, P2, P3, P4, P5,
      P6, P7, P8, P9, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19,
      P20, P21, P22, P23, P24, M0, L0, A, B, P, c_zero, r_zero, noise_str,
      gamma_i, gamma_a, beta11, beta22, beta12, coeff, delta_t,
      time_steps, delta_t_M0, time_steps_M0, file_timer, theta_int,
      theta_start, theta_steps, cont, INDEX_start, output_folder_name);

    //Analysis
    //printf("Simulaiton completed, do you want to start the analysis (y/n)?");
    char ans = 'y';
    //getch(ans);
    /*
    if(ans == 'y' || ans == 'Y'){
      fpw = fopen("analysis.sh", "w");

      if(fpw == NULL){
        printf("Couldn't create file analysis.sh, Aborting....\n");
        exit(0);
      }

      fprintf(fpw, "mkdir %s/analysis\n", output_folder_name);

      //Generating images for the c profile
      fprintf(fpw, "echo \"outputPath='%s/analysis'; name='%s'; N1=%d; N2=%d; S=%d; load 'makeMovie.gnu';\" | gnuplot\n",
              output_folder_name,
              "c",
              time_steps,
              time_steps_M0,
              file_timer);

      //Making a movie out of the images generated above
      fprintf(fpw, "cd %s/analysis\n",output_folder_name);
      fprintf(fpw, "ffmpeg -framerate 5 -pattern_type glob -i '*.png' -c:v libx264 %s_animation.mp4\n","c");
      fprintf(fpw, "rm -rf *.png\n");
      (void) fclose(fpw);

      (void) system("./analysis.sh\n");
      (void) system("chmod +x analysis.sh\n");
      //(void) system("rm analysis.sh")
    }
    */
    return 0;
}
