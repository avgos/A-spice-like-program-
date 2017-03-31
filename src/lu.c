
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_complex.h>

#include <gsl/gsl_linalg.h> //this requires gsl

#define COMPLEX(a,b) a+(b*_Complex_I)

int main (void)
{
  double complex a_data[] = { COMPLEX(0.18,11), 0.60, 0.57, 0.96,
                      0.41, 0.24, 0.99, COMPLEX(0.58,12),
                      0.14, COMPLEX(0.30,7.30), 0.97, 0.66,
                      0.51, 0.13, COMPLEX(0.19,12.3), 0.85 };

  double complex  b_data[] = { 1.0, 2.0, 3.0, 4.0 };




  gsl_complex b_;
  gsl_matrix_complex * a;
  //a = malloc(sizeof(gsl_complex)*16);
  a = gsl_matrix_complex_alloc (4,4);
  gsl_vector_complex * b;
  b = gsl_vector_complex_alloc(4);
  gsl_vector_complex * x;
  x = gsl_vector_complex_alloc(4);
  //for(int i=0;i<16;i++){
  for(int i=0;i<4;i++){
  	b_.dat[0] = 11.241-i; //creal(11)
  	b_.dat[1] = 3.213431+i; //cimag
  	gsl_vector_complex_set(b,i,b_);
  	for(int j=0;j<4;j++){
  		b_.dat[0] = creal(a_data[i*4+j]);
  		b_.dat[1] = cimag(a_data[i*4+j]);
  		gsl_matrix_complex_set(a, i, j, b_);
  	}
  }




//--------------------------------------------------
 //choleksy start
	gsl_linalg_complex_cholesky_decomp (a);
	gsl_linalg_complex_cholesky_solve (a,b,x);
//cholesky end
//--------------------------------------------------

  gsl_permutation * p;
  p = gsl_permutation_alloc(4);
  int skata;
  gsl_linalg_complex_LU_decomp(a,p,&skata);
  gsl_linalg_complex_LU_solve (a, p, b,x);



  for(int i =0;i<4;i++){

  	for(int j=0;j<4;j++){
  		gsl_complex temp =  gsl_matrix_complex_get(a, i, j);

  		printf("%.4f %.4fi\t",GSL_REAL(temp), GSL_IMAG(temp));

//  		printf("%.4f %.4fi\t", a[i*4+j].dat[0], a[i*4+j].dat[1]);
  	}
  	  gsl_complex b_i = gsl_vector_complex_get (b, i);
  	  gsl_complex x_i = gsl_vector_complex_get (x, i);

  	printf("\t| %6.4f %6.4fi = %6.4f %6.4fi",
  		GSL_REAL(b_i),GSL_IMAG(b_i),
  		GSL_REAL(x_i),GSL_IMAG(x_i)

  	 );
  	printf("\n");
  }









  int p_data[4];
  for(int i =0;i<4;i++) p_data[i] = i;
  double l_data[16];
  int vector_size = 4;
  int n = 4;
  for(int k=0; k<vector_size;k++){
  	//double x = fabs(a_data[k*vector_size+k]);
  	double x = sqrt(creal(a_data[k*vector_size+k])*creal(a_data[k*vector_size+k]) + cimag(a_data[k*vector_size+k])*cimag(a_data[k*vector_size+k]) ) ;
  	int m;
  	int changed;
  	changed = 0;
  	for(int i=k;i<vector_size;i++){
  		if(sqrt(creal(a_data[i*vector_size+k])*creal(a_data[i*vector_size+k]) + cimag(a_data[i*vector_size+k])*cimag(a_data[i*vector_size+k]) ) > x){
  	//	if(fabs(a_data[i*vector_size+k]) > x){
  			m = i;
  			changed = 1;
  		}
  	}
  	if(changed){
  		int t = p_data[k];
  		p_data[k] = p_data[m];
	  	p_data[m] = t;
  	}
	printf("%d %d %d\n",changed,k,m );
  	for(int i=0;i<vector_size && changed;i++){
  		double complex temp = a_data[m*vector_size+i];
  		a_data[m*vector_size+i] = a_data[k*vector_size+i];
  		a_data[k*vector_size+i] = temp;
  	}
  	for(int i=k+1;i<n;i++){
  		a_data[i*vector_size+k] = l_data[i*vector_size+k] = a_data[i*vector_size+k]/a_data[k*vector_size+k];
  	}
  	for(int i=k+1;i<n;i++){
  		for(int j=k+1;j<n;j++){
  			a_data[i*vector_size+j] -= l_data[i*vector_size+k]*a_data[k*vector_size+j];
  		}
  	}
  }
  for(int i =0;i<4;i++){
  	printf("%d\t", p_data[i]);
  	for(int j=0;j<4;j++){
  		printf("%.4f %.4fi\t", creal(a_data[i*vector_size+j]), cimag(a_data[i*vector_size+j]));
  	}
  	printf("\n");
  }


  return 0;
}



