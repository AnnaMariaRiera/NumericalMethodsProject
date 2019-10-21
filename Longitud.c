//Càlcul de la longitud L de l'elipse, aproximant la integral amb la fórmula de quadratura de Gauss-Txebishev
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//Definim el nombre pi amb 17 decimals
#define P 3.14159265359879323

int main()
{
  //Totes les variables necessàries al problema les fem en precisió doble
  double integral; //Aquest serà el càlcul de la integral aproximada
  int n=9; //Nombre de nodes
  double omega=P/n; //Les wi
  double x; //Nodes de Txebishev
  int i;
  
  //Càlcul aproximat de la integral
  integral=0;
  x=0;
  for(i=0;i<n;i++)
  {
      x = cos(((2*i-1)/((double)2*n))*P); //Obliguem al denominador a ser de tipus double, ja que si no ens divideix per un enter
      integral = omega*(sqrt(1-((2*x)*(2*x))/4)*sqrt(1+(((4*x*x)/64)*(1/(1-(4*x*x)/4))))) + integral; //Fem la suma de wi*f(2*xi)
  }
  
  integral = 4*integral;
  
  printf("El valor de la longitud de l'elipse es: %.16f\n", integral);
  
  return 0;
}