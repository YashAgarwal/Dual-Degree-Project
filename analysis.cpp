#include <stdio.h>
#include <bits/stdc++.h>
#include <iostream>
#include <string>
#define c(x,y) c[(y) + n_y * (x)]
//#define c(x,y) c[x][y]
#define interpolate(x0, y0, x2, y2, y1) (((y1)-(y0)) * (((double)(x2)-(x0))/((y2)-(y0)))) + (x0)

using namespace std;

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
		//R = (double[]){pr_x1 - pr_x2, pr_y1 - pr_y2, pr_x1, pr_x2, pr_y1, pr_y2};
		//R = (double[]){pr_x1 - pr_x2, pr_y1 - pr_y2, pr_x1, pr_x2, pr_y1, pr_y2};
		R[0] = pr_x1 - pr_x2;
		R[1] = pr_y1 - pr_y2;
		R[2] = pr_x1;
		R[3] = pr_x2;
		R[4] = pr_y1;
		R[5] = pr_y2;

		return R;
}

int main(int argc, char ** argv){
		if(argc != 3){
			cout << argc << " inputs given\n";
			cout << "\nERROR: usage ./analysis <DATA_ID> time\n";
			return 0;
		}

		//TODO: check if the readme file is there in the folder
		//TODO: Read the parameters
		int n_x=1024, n_y=1024;
		double delta_y=0.2, delta_x=0.2;

		//Make the filename
		char filename[50];
		sprintf(filename, "./output/%s/data/c-time_%09d.dat", argv[1], stoi(argv[2]));

		//open the filename
		ifstream f;
		f.open(filename);

		double *c = new double[n_x * n_y];
		int t;
		double d;
		for(int i=0; i<n_x; i++){
			for(int j=0; j<n_y; j++){
				f >> t;
				f >> t;
				f >> c(i,j);
				f >> d;
			}
		}

		//calculate the radius
		double *R = findRadius(c, 1024, 1024, 0.2, 0.2);

		for(int i=0; i<6; i++){
			cout << R[i] << endl;
		}

		return 0;
}
