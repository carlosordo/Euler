//aqui se crea el codigo para el metodo de euler con el que se soluciona la edp 
//metodo explicito 
#include <stdio.h>
#include <math.h>
#include <stdlib.h> //para el malloc
#include "tools.h"

//llamamos los prototipos de las condiciones iniciales 
double T0(double x);
double q(double x, double t);
//escribimos el int main tanto para el metodo explicito e implicito
int main(){
	int i, j,k; //para los ciclos for 
	int n=100; //n=100 debido a que h=0.01
	int tmax=1, m=10000; //tiempo maximo y m se toma por delta t 
	
	double h = 0.01, dt=1e-4,lambda; //condiciones dadas por el problema y condición para que converga el método 	
	lambda=(dt)/(h*h);// como se define lambda segun la teoria
	//apuntadores para las matrices dinamicas
	double **matrizT,**matrizTit,**matrizA;
	
	//alloc de la matriz heat para el metodo explicito e implicito.
	matrizT = (double * *)malloc(n * sizeof(double));
	matrizTit = (double * *)malloc(n * sizeof(double));
	matrizA = (double * *)malloc(n * sizeof(double));
	for (i=0; i<n; i++){
		matrizT[i] = (double *) malloc(m*sizeof(double) );
		matrizTit[i] = (double *) malloc(m*sizeof(double) );
		matrizA[i] = (double *) malloc(n*sizeof(double) );
		//verificamos que existe ese espacio de memoria 
		if(matrizT[i] == NULL){
			perror("ERROR. No hay suficiente memoria ");
			exit(EXIT_FAILURE);
		}
		if(matrizTit[i] == NULL){
			perror("ERROR. No hay suficiente memoria ");
			exit(EXIT_FAILURE);
		}
		if(matrizA[i] == NULL){
			perror("ERROR. No hay suficiente memoria ");
			exit(EXIT_FAILURE);	
		}
	}
	
	//se escribe la matriz A
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			if (i==j){
				matrizA[i][j]=1-2*lambda;
			}
			if (j==i+1){
				matrizA[i][j]= -lambda;
			}
			if (j==i-1){
				matrizA[i][j]= lambda;
			}
			else{
				matrizA[i][j]=0;
			}
		}
	}	
		
	//se da la memoria para los vectores (usaremos para el metodo explicito e implicito pues no cambian)
	double *x = malloc(n*sizeof(double));
	double *t = malloc(m*sizeof(double));
	
	//Vector de posición 
	x[0] = 0;
	x[n] = 1;
	for( i = 1; i < n; i++ ) {
      	x[i] = x[i-1] + h ;
	}
	
	//Vector de tiempo 
	t[0] = 0;
	t[m] = tmax;
	for( i = 1; i < m; i++ ) {
      	t[i] = t[i-1] + dt ; 
	}
	
	//se usan las funciones para colocar nuestros valores iniciales dados por el problema
	//Colocamos la condición inicial t=0 en nuestra matriz T
	for (j=0; j < n; j++){
		matrizT[j][0] = T0(x[j]);
		matrizTit[j][0] = T0(x[j]);
		
	}
	//Colocamos las condiciones de frontera x=0 y x=1
	for (i=0;i<m;i++){
		matrizT[0][i] = 0; //no usamos una funcion, pues solamente es 0 
		matrizTit[0][i] = 0;
	}
	
	for (i=0;i<m;i++){
		matrizT[n-1][i] = 0;   //no se usa una funcion pues es 0
		matrizTit[n-1][i] = 0;    
	}
	//Se aplica el calculo del metodo de Euler explicito
	for (i = 1; i < n-1; i++){
		for (j = 1; j < m-1; j++){
			matrizT[i][j]=matrizT[i][j-1]+(lambda*(matrizT[i+1][j-1]-2*matrizT[i][j-1]+matrizT[i-1][j-1]))+dt*q(x[i],t[j]);
		}	
	}
	// Exportación de datos a archivo de texto plano, datos para usar en Gráfica.
	FILE *eulerdata=fopen("data2.txt", "w");
	for (i=0; i<n; i++){
		for (j =0 ; j < m; j++)
		{
			fprintf(eulerdata, "%f, %f, %f\n", x[i], t[j], matrizT[i][j]);
		}
	}
	fclose(eulerdata);
	
	
	//creando la nueva matriz para el metodo implicito
	while (k<20,k++){
		for (i = 1; i < n-1; i++){
			for (j = 1; j < m-1; j++){
				
				matrizTit[i][j]=(matrizT[i][j-1]+ lambda*matrizT[i+1][j]+ lambda*matrizT[i-1][j])/(2*lambda +1)+(dt*q(x[i],t[j]));
				matrizT[i][j]=matrizTit[i][j];
				
			}
				
		}			
	}
	free (matrizT);
	free (matrizTit);
	return 0;
}


//definimos las funciones para las condiciones iniciales

//se define la funcion para t(x,0) como e^x
double T0(double x){
	return exp(x);
}
//definimos q(x,t) como se nos da en el ejemplo 
double q(double x, double t){
	return  cos(M_PI*t) * sin(2*M_PI*x);
}





