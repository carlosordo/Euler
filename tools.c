//aqui se crea el codigo para el Jacobi 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int Jacobi(int n, double **matrizA, double *vectorb, int *xiter, double tol, int niter){
	
	int i, j, k = 0;//contadores usados para nuestro ciclo for.
	double *x=NULL;  //apuntador de nuestro vector solucion aproximado 
	x = (double *) malloc(n*sizeof(double));
	if (x==NULL){
		perror("ERROR. No hay sufuciente memoria");
		exit(EXIT_FAILURE);
	}
	//variable para analizar el error de la aproximación.
	double err;

	while (k<niter,k++){
		for (i=0;i<n;i++){
			//variables para guardar la suma o la diferencia de las dos aproximaciones guardadas en x
			double sum, diff;
			sum=vectorb[i];
			for (j=0;j<n;j++){
				sum=sum-matrizA[i][j]*xiter[j];  //operacion dada por el método
			}
			sum=sum+matrizA[i][i]*xiter[i];
			x[i]=sum/matrizA[i][i];
			
			diff=x[i]-xiter[i];
			err=err+diff*diff;
		}
		//se comprueba si cumple la condicion de convergencia 
		if (err<tol){
			break;
		//se prepara para realizar la siguiente iteración	
		} else{
			xiter[i] = x[i];
			
		}
	}
	return k;
}


