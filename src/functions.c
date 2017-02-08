#include "../headers/functions.h"
#include <stdlib.h>
#define c(x,y) c[(y) + n_y * (x)]
#define interpolate(x0, y0, x2, y2, y1) (((y1)-(y0)) * (((double)(x2)-(x0))/((y2)-(y0)))) + (x0)

double* findRadius(double *c, int n_x, int n_y, double delta_x, double delta_y){

  double pr_x1, pr_x2, pr_y1, pr_y2;

  for (int i1 = n_x / 2; i1 < n_x; ++i1) {
    int i2 = n_y / 2;
    if (c(i1, i2) < 0.5) {
      //pr_x1 = (i1-1) + (c(i1-1,i2) - 0.5)/(c(i1-1,i2) - c(i1, i2));
      pr_x1 = interpolate(i1-1, c(i1-1,i2), i1, c(i1,i2), 0.5);
      pr_x1 *= delta_x;
      break;
    }
  }

  for (int i1 = n_x / 2; i1 > 0; --i1) {
    int i2 = n_y / 2;
    if (c(i1, i2) < 0.5) {
      pr_x2 = interpolate(i1, c(i1, i2), i1 + 1, c(i1+1, i2), 0.5);
      pr_x2 *= delta_x;
      break;
    }
  }

  for (int i2 = n_y / 2; i2 < n_y; ++i2) {
    int i1 = n_x / 2;
    if (c(i1, i2) < 0.5) {
      //pr_x1 = (i1-1) + (c(i1-1,i2) - 0.5)/(c(i1-1,i2) - c(i1, i2));
      pr_y1 = interpolate(i2-1, c(i1,i2-1), i2, c(i1,i2), 0.5);
      pr_y1 *= delta_y;
      break;
    }
  }

  for (int i2 = n_y / 2; i2 > 0; --i2) {
    int i1 = n_x / 2;
    if (c(i1, i2) < 0.5) {
      //pr_x1 = (i1-1) + (c(i1-1,i2) - 0.5)/(c(i1-1,i2) - c(i1, i2));
      pr_y2 = interpolate(i2+1, c(i1,i2+1), i2, c(i1,i2), 0.5);
      pr_y2 *= delta_y;
      break;
    }
  }

  double *R = (double*)malloc(sizeof(double)*6);
  R = (double[]){pr_x1 - pr_x2, pr_y1 - pr_y2, pr_x1, pr_x2, pr_y1, pr_y2};
  return R;
}
