//Anna Maria Riera Escandell, NIU 1359293
//Càlcul del punt (x,y), amb y>0, de l'el·lipse (x/2)^2+(2y)^2=1 (x pertanyent a [-2,2]), tal que la proporció entre les longituds dels arcs mesurats des del punt (2,0) sigui 4:1.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//Definim el nombre pi amb 17 decimals
#define P 3.14159265359879323
#define ht 0.01
//Precisió del resultat de Gauss-Seidel (resolució del sistema per trobar el vector de moments dels splines)
#define T 0.000001
 

//Per a imprimir les matrius per pantalla
//Passem com a arguments la matriu i la seva dimensió
void PintaMatriu ( double** matriu, int n)
{
	int i, j;
	printf("\n");
	for (i=0 ; i<n ; i++)
	{ 
		for ( j=0 ; j<n ; j++)
		{
			printf("%f ", matriu[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

//Funció per a calcular la longitud total de l'el·lipse (integral de la corba aproximada mitjançant quadratura de Gauss-Txebyshev)
double LongTotal()
{
	double integral;
	int n=9; //Nombre de nodes
	double omega=P/n; //Les wi
	double x; //Nodes de Txebishev
	int i;
  
	//Càlcul aproximat de la integral
	integral=0;
	x=0;
	for(i=0;i<n;i++)
	{
		x = cos(((2*i-1)/(2.*n))*P); //Obliguem al denominador a ser de tipus double, ja que si no ens divideix per un enter
		integral = omega*(sqrt(1-((2*x)*(2*x))/4)*sqrt(1+(((4*x*x)/64)*(1/(1-(4*x*x)/4))))) + integral; //Fem la suma de wi*f(2*xi)
	}

	integral = 4*integral;

	return integral;
} 

//Resolució del sistema lineal A*x=b mitjançant el mètode de Gauss-Seidel. Torna el vector de la solució (x)
//Passem com a arguments la matriu (A), el vector de termes independents (b) i la dimensió de la matriu 
double* Resol_Sistema(double** a, double* b, int n)
{
	int k; //Compta iteracions
	int Kmax=50; //Nombre màxim d'iteracions
	int i, j;
	double e; //Criteri de parada 
	double *x; //Vector de la solució, de dimensió n (de 0 fins a n-1)
	double *xa; //Vector auxiliar de la solució
	double s; //Variable per anar fent els sumatoris
	
	x=(double*)malloc(n*sizeof(double));
	xa=(double*)malloc(n*sizeof(double));
	
	//Inicialitzem les coses
	for(i=0; i<n; i++)
	{
		x[i]=0;
		xa[i]=0;
	}
	k=0; 
	//Fem les iteracions
	do
	{
		k=k+1; 
		e=-1;
		for(i=0; i<n; i++) //Guardem totes les components del vector x[i]
		{
			s=0;
			for (j=0; j<=(i-1); j++)
			{
				s = s + a[i][j]*x[i];
			}
			for (j=(i+1); j<n; j++)
			{
				s = s + a[i][j]*x[i];
			}
			x[i] = (b[i]-s)/(a[i][i]);
		}
		for(i=0; i<(n-1); i++)
		{
			if (fabs(xa[i]-x[i])>e)
			{
				e=fabs(xa[i]-x[i]);  
			}
			xa[i]=x[i];
		}
		
	} while ((k<Kmax)&&(e>=T)); //Fa el procés mentre es compleixin les condicions. Si alguna deixa de complir-se, pararà
	
	
	free(xa);
	xa = NULL;
	
	return(x);
}


//Funció que troba la matriu dels coeficients de cada spline (alpha_i, beta_i, delta_i i gamma_i)
//Passem com a arguments el vector dels nodes t_i (t_i=x+ih) i el nombre n de components del vector sense comptar el 0 (de manera que hi ha n+1 nodes)
double** Splines (double *t, int n)
{
	int	i, j; //Variables auxiliars
	double *y; //Coordenades de suport
	double *h; //vector amb hi=t_i-t_{i-1}
	double *d, *lambda, *mu;
	double *M; //Vector dels moments (de 0 fins a n)
	double **matriu; //Guarda la matriu amb 2 a la diagonal, lambda i mu
	double *resultat; //Moments de 1 fins a n-1
	double **Coefs_Splines; //Matriu que guardarà els coeficients finals de cada spline
		
	//Reservem memòria
	y=(double*)malloc((n+1)*sizeof(double));
	h=(double*)malloc(n*sizeof(double)); //El vector en sí va de 0 fins n-1 (com sempre en C), però l'entenem com de 1 fins a n (guardem en les posicions 0 fins n-1 els valors de h per 1 fins n). El mateix amb d, lambda i mu, però amb una posició menys (van de 1 fins n-1 i els escrivim com 0 fins n-2)
	d=(double*)malloc((n-1)*sizeof(double));
	lambda=(double*)malloc((n-1)*sizeof(double));
	mu=(double*)malloc((n-1)*sizeof(double));
	M=(double*)malloc((n+1)*sizeof(double));
	resultat = (double*)malloc((n-1)*sizeof(double));
	matriu = (double**) malloc ((n-1)*sizeof(double*));
	for(i=0; i<(n-1); i++)
	{
		matriu[i] = (double*) malloc ((n-1)*sizeof(double));
	}
	Coefs_Splines = (double**) malloc (n*sizeof(double*));
	for(i=0; i<n; i++)
	{
		Coefs_Splines[i] = (double*) malloc (4*sizeof(double));
	}
	
	//Rellenem els vectors yi, hi, di, lambdai, mui
    for(i=0;i<=n;i++)
    {
	   y[i]=sqrt((64-15*t[i]*t[i])/(64-16*t[i]*t[i]));
    }
	for(i=0;i<=n-1;i++)
    {
	   h[i]=t[i+1]-t[i];
    }
	for(i=0;i<=n-2;i++)
    {
	   d[i]=(6/(h[i]+h[i+1]))*((y[i+2]-y[i+1])/(h[i+1])-(y[i+1]-y[i])/(h[i]));
	   lambda[i]=h[i+1]/(h[i]+h[i+1]);
	   mu[i]=h[i]/(h[i]+h[i+1]);
    }
	
	//Trobem els moments M: fem la matriu i resolem el sistema
	//Emplenem la matriu
	for (i=0 ; i<=n-2 ; i++)
	{ 
		for (j=0 ; j<=n-2 ; j++)
		{
			matriu[i][j]=0;
		}
		for (j=0 ; j<=n-2 ; j++)
		{
			if (j==i) {
				matriu[i][j] = 2; //A la diagonal tot 2
			}
	        if (j+1==i) {
				matriu[i][j] = mu[i]; //La diagonal de baix tot mu
			}
		    if (j-1==i) {
				matriu[i][j] = lambda[i]; //La diagonal de dalt tot lambda
			}
		}
	}
	//Resolem el sistema:
	resultat = Resol_Sistema(matriu, d, n-1);
	
	//Ara tenim el resultat guardat a "resultat", així que el posem a "M" i hi afegim els termes que falten per acabar d'emplenar el vector de moments:
	M[0]=0;
	M[n]=0;
	for(i=0;i<(n-1);i++)
    {
	   M[i+1]=resultat[i];
    }
	
	//Construïm els coeficients dels splines, ara que ja tenim els moments, i els emmagatzemem tots en una matriu
	for(i=0;i<=n-1;i++)
    {
		Coefs_Splines[i][0]=y[i];
		Coefs_Splines[i][1]=((y[i+1]-y[i])/(h[i]))-((2*M[i]+M[i+1])/(6.))*h[i];
		Coefs_Splines[i][2]=M[i]/(2.);
		Coefs_Splines[i][3]=(M[i+1]-M[i])/(6.*h[i]);
    }

	//Alliberem memòria de tot el que no tornarem a utilizar
	free(y);
	y = NULL;
	free(h);
	h = NULL;
	free(d);
	d = NULL;
	free(lambda);
	lambda = NULL;
	free(mu);
	mu = NULL;
	free(M);
	M = NULL;
	free(resultat);
	resultat = NULL;
	for (i=0; i<(n-1); i++){
		free(matriu[i]);
		matriu[i] = NULL;
	}
	free(matriu);
	matriu = NULL;

	return (Coefs_Splines);
}

//Funció que calcula F(x) i si es vol, imprimeix per pantalla x, F(x) i Lea
//Passem com a argument x (el valor en què avaluem F) i un enter que si val 1 farà que s'imprimeixin per pantalla els valors anteriors
double F_x (double x, int i)
{
	int j;
	double L;
	int n; //Nombre de nodes (que és n+1), que depèn de cada x
	double h0; //h0=(2-x)/n
	double *t; //Els nodes t_i
	double Fx; //Vector que guardarà F(xk)=L-5Lea
	double **Coefs_Splines; //Matriu on es guarden els coeficients de tots els splines per a cada x
	double Lea;
	
	L = LongTotal();
	
	n = floor((1.99-x)/(ht)); //La funció floor ens dóna la part entera
	h0 = (1.99-x)/n;
	
	//Fem un vector de les ti
	t=(double*)malloc((n+1)*sizeof(double)); 
	for(j=0;j<=n;j++)
	{
		t[j]=x+j*h0;
	}
	
	//Rervem memòria per la matriu dels Splines
	Coefs_Splines = (double**) malloc (n*sizeof(double*));
	for(j=0; j<=(n-1); j++)
	{
		Coefs_Splines[j] = (double*) malloc (4*sizeof(double));
	}
	
	//Ara cridem la funció Splines que ens empleni la matriu
	Coefs_Splines=Splines(t, n);
	
	//Càlcul de Lea (aproximació de la integral per una suma de 0 fins a n-1)
	Lea=0;
	for(j=0; j<=n-1; j++)
	{
		Lea = Lea + Coefs_Splines[j][0]*(t[j+1]-t[j]) + Coefs_Splines[j][1]*((t[j+1]-t[j])*(t[j+1]-t[j]))/(2.) + Coefs_Splines[j][2]*((t[j+1]-t[j])*(t[j+1]-t[j])*(t[j+1]-t[j]))/(3.) + Coefs_Splines[j][3]*((t[j+1]-t[j])*(t[j+1]-t[j])*(t[j+1]-t[j])*(t[j+1]-t[j]))/(4.);
	}
		
	Fx= L -5*Lea;
	
	if (i==1)
	{
		printf("x=%.16f\n", x);
		printf("F(%.16f)=%.16f\n", x, Fx);
		printf("Lea=%.16f\n", Lea);
	}
	
	
	//Alliberem memòria
	free(t);
	t = NULL;
	for (i=0; i<=(n-1); i++){
		free(Coefs_Splines[i]);
		Coefs_Splines[i] = NULL;
	}
	free(Coefs_Splines);
	Coefs_Splines = NULL;

	return(Fx);
}


int main()
{
	int i;
	int si=1; //Per imprimir els valors de la funció F_x
	int no=0; //Per no imprimir els valors de la funció F_x
	double L; 
	double x0=0.5; //Llavor inicial per fer Newton
	double x; //Els xk
	double Fx; //Vector que guardarà F(xk)=L-5Lea
	double xa; //xn-1 que utilitzem com a variable auxiliar al criteri de parada
	double E=0.000001; //serà l'épsilon del criteri de parada, és el que ens diu l'enunciat (10^-6)
	
	L = LongTotal();
	printf("\n\nEl valor de la longitud total de l'elipse es %.16f\n", L);
	
	//Inicialitzem
	x=x0;
	xa=x0;
	printf("Prenem com a llavor inicial x=%.1f\n\n", x);
	
	//Fem el mètode de Newton
	printf("Iteracio 0:\n"); //Imprimim els resultats de la iteració 0: x0, F(x0) i Lea
	Fx = F_x(x, si);
	printf("\n");
	for (i=0; i<=20; i++) //Posem un nombre màxim d'iterats, per exemple 50
	{
		printf("Iteracio %d:\n", i+1);
		
		Fx = F_x(x, no); //Fem F(x) per a poder calcular l'iterat
		
		//Calculem el valor de cada iterat x_{k+1}
		x = x - Fx/((5.)*sqrt((64-15*x*x)/(64-16*x*x))); //x_{k+1}=x_{k}-F(xk)/F'(xk)
		
		Fx = F_x(x, si); //Avaluem F en la nova x per a imprimir per pantalla els resultats
		printf("El valor absolut de la diferencia entre la x anterior i aquesta es %.16f\n\n", fabs(x-xa));
		
		if (fabs(x-xa)<=E) //Criteri de parada
		{
			printf("RESULTAT FINAL: despres de fer %d iteracions amb el metode de Newton, obtenim el valor de x*=%.16f.\n\n", i+1, x);
			return 0;
		}
		
		xa=x;
	}

	return 0;
}