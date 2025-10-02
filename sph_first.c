#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<malloc.h>
#define X 0
#define Y 1
#define _USE_MATH_DEFINES

typedef struct
{

  int id;
  double pos[2];
  double vel[2];
  double accel[2];
  double mass;
  double rho;
  double h;
  double p;
  double c;
  double du;
  double u;
  int *nn;
  int nNeighbors;
  double *dx;
  double *dy;
  double *r;
  double *W;
  double *dWx;
  double *dWy;
  int type;
  
}Particles;

Particles *part, *auxPart;
int nFluid, nPart;

/* Llamado a funciones secundarias */

void ics(int nx, int ny, double dx, double dy, double Lx,
	 double Ly); // Condiciones iniciales

double W(double r, double h);
double dW(double r, double dx, double h);
void testKernel(void);
void NN(int i);
void test_NN(void);
void density(void);
void eos(void);
void navierStokes(void);
void viscosity(double dx);
void boundaryInteraction(double dx);
void meanVelocity(void);
void acceleration(double dx);
void drift(double dt);
void kick(double dt);
void printState(char *outfile);

int main(int argc, char *argv[])
{

  int i, nx, ny, counter;
  double Lx, Ly, dx, dy;
  double dt = 5e-5;
  double t, tTotal = atoi(argv[1])*dt;
  char outfiles[500];
  //double t, tTotal = 4000*dt;
  
  printf("Corriendo durante %d pasos, un tiempo total de %lf s\n",
	 atoi(argv[1]),tTotal);

  nx = 40;
  ny = 40;
  Lx = 1e-3;
  Ly = 1e-3;
  dx = Lx/nx;
  dy = Ly/ny;

  nFluid = nx*ny;

  part = (Particles *)malloc((size_t)nFluid*sizeof(Particles));

  if( part==NULL )
    {
      printf("Error alocando part\n");
      exit(0);
    }

  // Create the initial conditions
  ics( nx, ny, dx, dy, Lx, Ly);

  // testing kernel function
  testKernel();
  counter = 0;
  t = 0;

  // printting system initial state
  sprintf(outfiles,"./output/state_%.4d",counter);
  printState(outfiles);

  // main loop

  while( t<=tTotal )
    {
      // searching near neighbors for all fuid particles
      for( i=0; i<nFluid; i++ )
	NN(i);

      // testing near neighbors searching
      if(counter==0)
	test_NN();

      // computing density
      density();

      // drift in leap−frog integration
      drift(dt);

      // computing acceleration
      acceleration(dx);

      // kick in leap−frog integration
      kick(dt);

      // drift in leap−frog integration
      drift(dt);
      t = t + dt;
      counter++;
      
      // printting system state
      sprintf(outfiles,"./output/state_%.4d",counter);
      printState(outfiles);
      printf("step = %d \n",counter);

    }
  
  free(part);
  return 0;
}

/* 2) */
double W(double r, double h){
  double R = r/h;
  double alpha = 15.0/(7.0*M_PI*h*h);
  if( (R >= 0.0) && (R < 1.0) )
    return alpha*((2.0/3.0) - R*R + 0.5*R*R*R);
  if( (R >= 1.0) && (R <= 2.0) )
    return alpha*((1.0/6.0)*(2.0-R)*(2.0-R)*(2.0-R));
  if( R>2.0)
    return 0.0;
  return 0.0;
}

double dW(double r, double dx, double h){
  double R = r/h;
  double alpha = 15.0/(7.0*M_PI*h*h);
  if( (R >= 0.0) && (R < 1.0) )
    return alpha*(-2.0 + 1.5*R)*dx/(h*h);
  if( (R >= 1.0) && (R <= 2.0) )
    return alpha*(-0.5*(2.0-R)*(2.0-R))*dx/(h*h*R);
  if( R>2.0)
    return 0.0;
  return 0.0;
}

void testKernel(void){
  
  double r, w, dw;
  FILE *fKernelTest;
  fKernelTest = fopen("kernel_test.output","w");

  for( r=-3.0; r<=3.0; r = r + 0.1){
    w = W( fabs(r), 1.0);
    dw = dW( fabs(r), r/sqrt(3.0), 1.0);
    fprintf(fKernelTest,"%16.10lf %16.10lf %16.10lf\n",r,w,dw);
    }
  
  fclose(fKernelTest);

}

// Searching the near neighbors
void NN(int i)
{
  
  double kappa = 2.0;
  double xij, yij, rij, hij;
  double *auxDouble;
  int j, *auxInt, nNeighbors;

  // liberamos la información anterior para cambiar el valor con el paso del tiempo
  free(part[i].nn);
  free(part[i].dx);
  free(part[i].dy);
  free(part[i].r);
  free(part[i].W);
  free(part[i].dWx);
  free(part[i].dWy);

  // apuntamos nuevamente a NULL para volver a llenar
  part[i].nn = NULL;
  part[i].dx = NULL;
  part[i].dy = NULL;
  part[i].r = NULL;
  part[i].W = NULL;
  part[i].dWx = NULL;
  part[i].dWy = NULL;
  
  // por cada ciclo 
  nNeighbors = 0;

  // Buscamos los vecinos para todas las partículas
  for( j=0; j<nPart; j++ )
    {

      // excepto para las de la frontera, sin embargo, estas partículas pueden ser vecinos de las de fluido
      if( i!=j )
	{
	  
	  xij = part[i].pos[X] - part[j].pos[X];
	  yij = part[i].pos[Y] - part[j].pos[Y];
	  rij = sqrt( xij*xij + yij*yij );
	  hij = 0.5*(part[i].h+part[j].h);

	  if( rij < kappa*hij + 1e-10 ) // La suma de 1e-10 hace que los problemas del borde se eliminen al colocar una distancia mucho
	    // mucho más pequeña adicional
	    {
	      // Si está dentro del radio entonces aumentamos el número de vecinos en 1
	      nNeighbors++;

	      
	      // add neighbor id
	      auxInt = NULL;
	      auxInt = (int *)realloc(part[i].nn,
				      (size_t)(nNeighbors)*
				      sizeof(int));
	      
	      part[i].nn = auxInt;
	      auxInt = NULL;

	      
	      // add neighbor dx
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dx,
					    (size_t)(nNeighbors)*
					    sizeof(double));

	      part[i].dx = auxDouble;
	      auxDouble = NULL;


	      // add neighbor dy
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dy,
					    (size_t)(nNeighbors)*
					    sizeof(double));

	      part[i].dy = auxDouble;
	      auxDouble = NULL;


	      // add neighbor r
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].r,
					    (size_t)(nNeighbors)*
					    sizeof(double));

	      part[i].r = auxDouble;
	      auxDouble = NULL;

	      
	      // add neighbor W
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].W,
					    (size_t)(nNeighbors)*
					    sizeof(double));

	      part[i].W = auxDouble;
	      auxDouble = NULL;


	      // add neighbor dWx
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dWx,
					    (size_t)(nNeighbors)*
					    sizeof(double));

	      part[i].dWx = auxDouble;
	      auxDouble = NULL;


	      // add neighbor dWy
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dWy,
					    (size_t)(nNeighbors)*
					    sizeof(double));

	      part[i].dWy = auxDouble;
	      auxDouble = NULL;

	      // El número de vecinos va hasta 0
	      part[i].nn[nNeighbors-1] = j;
	      part[i].dx[nNeighbors-1] = xij;
	      part[i].dy[nNeighbors-1] = yij;
	      part[i].r[nNeighbors-1] = rij;
	      part[i].W[nNeighbors-1] = W( rij, hij );
	      part[i].dWx[nNeighbors-1] = dW( rij, xij, hij);
	      part[i].dWy[nNeighbors-1] = dW( rij, yij, hij);

	    } // final if r_ij
	} // final if i!=j
    } // final for j
  
  part[i].nNeighbors = nNeighbors;
  
} // Final función del número de Vecino (NN)

/* Analizamos los vecinos de 20 partículas aleatorias para determinar
   si nuestro código anterior está correcto por lo que en necesario
   graficar y observar la estructura
*/

void test_NN(void)
{
  int i,j,k;
  FILE *fTestNN;
  fTestNN = fopen("NN_test.output","w");
  srand(time(NULL));
  for( k=0; k<20; k++)
    {
      i = rand() % nFluid;
      printf("testing for particle %d\n",i);
      printf("with %d neighbors\n",part[i].nNeighbors);
      fprintf(fTestNN,"%16d %16.10lf %16.10lf\n",	      
	      part[i].id,
	      part[i].pos[X],
	      part[i].pos[Y]);
      
      for( j=0; j<part[i].nNeighbors; j++ )
	fprintf(fTestNN,"%16d %16.10lf %16.10lf\n",
		part[i].nn[j],
		part[part[i].nn[j]].pos[X],
		part[part[i].nn[j]].pos[Y]);
      
      fprintf(fTestNN,"\n");
      /*Este salto de página permite hacer cortes en la información
	que entra a la gráfica de gnuplot*/
    }
  fclose(fTestNN);
}

void density(void)
{
  int i, j;
  double wii, norm;
  for( i=0; i<nFluid; i++ ){
      
    // self density
    wii = W( 0.0, part[i].h );
    
    // computing density
    part[i].rho = part[i].mass*wii;

    for( j=0; j<part[i].nNeighbors; j++ ){
      part[i].rho =
	part[i].rho + part[part[i].nn[j]].mass*part[i].W[j];
    }

    // normalizing the density
    norm = (part[i].mass/part[i].rho)*wii;

    for( j=0; j<part[i].nNeighbors; j++ ){
      norm =norm +
	(part[part[i].nn[j]].mass/part[part[i].nn[j]].rho)
	*part[i].W[j];
    } // final for j
    
    part[i].rho = part[i].rho/norm;
  } // final for i
  
  printf("density computed\n");
} // final función density

void eos(void){
  int i;
  for( i=0; i<nPart; i++ )
    {
      part[i].c = 0.01;
      part[i].p = part[i].c*part[i].c*part[i].rho;
    }
} // final función eos (Ecuación de Estado)

void navierStokes(void){

  int i, j, k;
  double pij, vdw;

  // computing sound speed and pression
  eos();

  // computing acceleration
  for( i=0; i<nFluid; i++ ){
    part[i].accel[X] = part[i].accel[Y] = 0.0; // Añadir las fuerzas gravitacionales acá
    part[i].du = 0.0;
    
    for( k=0; k<part[i].nNeighbors; k++ ){

      //Truco para simplificar los corchetes en la escritura
      //afecta tiempo de computo
      j = part[i].nn[k];  

      pij = ( part[i].p/(part[i].rho*part[i].rho) )
	+ ( part[j].p/(part[j].rho*part[j].rho) );
      
      part[i].accel[X] = part[i].accel[X]
	- part[j].mass*pij*part[i].dWx[k];
      
      part[i].accel[Y] = part[i].accel[Y]
	- part[j].mass*pij*part[i].dWy[k];
      
      vdw = (part[i].vel[X]-part[j].vel[X])*part[i].dWx[k]
	+ (part[i].vel[Y]-part[j].vel[Y])*part[i].dWy[k];
      part[i].du = part[i].du + 0.5*part[j].mass*pij*vdw;
      
    } //W = integral( f*dr ) = integral( f*dr*(dt/dt) ) = integral( f*v*dt )
}
  
printf("acceleration computed\n");
}

void viscosity(double dx){
  
  int i, j, k;
  double xij, yij, vxij, vyij, vijrij, vdw;
  double hij, cij, phiij, rhoij, Piij;
  double alphapi = 1.0;
  double betapi = 1.0;
  double eps = dx;
  double eps2 = 0.01*eps*eps;
  for( i=0; i<nFluid; i++ ){

    for( k=0; k<part[i].nNeighbors ; k++ ){

      j = part[i].nn[k];

      xij = part[i].pos[X] - part[j].pos[X];
      yij = part[i].pos[Y] - part[j].pos[Y];
      vxij = part[i].vel[X] - part[j].vel[X];
      vyij = part[i].vel[Y] - part[j].vel[Y];
      vijrij = vxij*xij + vyij*yij;

      if( vijrij < 0.0 ){
	
	hij = 0.5*(part[i].h+part[j].h);
	phiij = (hij*vijrij)/( xij*xij + yij*yij + eps2);
	cij = 0.5*(part[i].c+part[j].c);
	rhoij = 0.5*(part[i].rho+part[j].rho);
	
	Piij = ( -alphapi*cij*phiij + betapi*phiij*phiij )/( rhoij );

	part[i].accel[X] = part[i].accel[X]
	  - part[j].mass*Piij*part[i].dWx[k];
	
	part[i].accel[Y] = part[i].accel[Y]
	  - part[j].mass*Piij*part[i].dWy[k];

	vdw = (part[i].vel[X]-part[j].vel[X])*part[i].dWx[k]
	  + (part[i].vel[Y]-part[j].vel[Y])*part[i].dWy[k];

	part[i].du = part[i].du + 0.5*part[j].mass*Piij*vdw;
	
      } // final if
      
    } // final for k
    
  } // final for i
  
  printf("viscosity computed\n");
} // final función viscosity

void boundaryInteraction(double dx){

  int i, j;
  int n1 = 12, n2 = 4;
  double r0 = dx/2.0, D = 0.01;
  double xij, yij, rij, PBxij, PByij;

  for( i=0; i<nFluid; i++ ){
    
    for( j=0; j<part[i].nNeighbors; j++ ){

      if( part[part[i].nn[j]].type==-1 ){

	xij = part[i].pos[X] - part[part[i].nn[j]].pos[X];
	yij = part[i].pos[Y] - part[part[i].nn[j]].pos[Y];
	rij = sqrt( xij*xij + yij*yij );

	if( rij<r0 ){

	  PBxij = D*( pow((r0/rij),n1)
		      - pow((r0/rij),n2) )*(xij/(rij*rij));
	  
	  PByij = D*( pow((r0/rij),n1)
		      - pow((r0/rij),n2) )*(yij/(rij*rij));
	  
	  part[i].accel[X] = part[i].accel[X] + PBxij;
	  part[i].accel[Y] = part[i].accel[Y] + PByij;
	} // final if rij

      } // final if part 

    } // final for j
  } // final for i
  
  printf("interaction with boundary computed\n");
} // final 

void meanVelocity(void){
  
  int i, j;

  //Parámetro de corrección de discontinuidades XSPH
  double epsilon = 0.3; 

  double vxMean, vyMean;
  double vxij, vyij, rhoij;

  for( i=0; i<nFluid; i++ ){

    vxMean = 0.0;
    vyMean = 0.0;

    for( j=0; j<part[i].nNeighbors; j++ ){

      vxij = part[i].vel[X] - part[part[i].nn[j]].vel[X];
      vyij = part[i].vel[Y] - part[part[i].nn[j]].vel[Y];
      rhoij = 0.5*(part[i].rho + part[part[i].nn[j]].rho);

      vxMean = vxMean + (part[part[i].nn[j]].mass/rhoij)
	*vxij*part[i].W[j];

      vyMean = vyMean + (part[part[i].nn[j]].mass/rhoij)
	*vyij*part[i].W[j];
    } // final for j

    part[i].vel[X] = part[i].vel[X] - epsilon*vxMean;
    part[i].vel[Y] = part[i].vel[Y] - epsilon*vyMean;

  } // final for i
}

void acceleration(double dx){
  
  // computing acceleration and change of energy
  navierStokes();
  
  // computing viscosity contribution
  viscosity(dx);
  
  // computing interaction with boundary
  boundaryInteraction(dx);
  
  // correction to mean velocity
  meanVelocity(); //Modelo XSPH 
  printf("acceleration computed\n");
}

void drift(double dt)
{
  int i;
  for( i=0; i<nFluid; i++ )
    {
      part[i].pos[X] = part[i].pos[X] + 0.5*dt*part[i].vel[X];
      part[i].pos[Y] = part[i].pos[Y] + 0.5*dt*part[i].vel[Y];
      part[i].u = part[i].u + 0.5*dt*part[i].du;
    }
} // final función drift

void kick(double dt){
  int i;
  for( i=0; i<nFluid; i++ )
    {
      part[i].vel[X] = part[i].vel[X] + dt*part[i].accel[X];
      part[i].vel[Y] = part[i].vel[Y] + dt*part[i].accel[Y];
    }
} // final función kick

void printState(char *outfile){
  int i;
  FILE *fState;
  fState = fopen(outfile,"w");
 
  for( i=0; i<nPart; i++){
    fprintf(fState,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf%.10lf\n",
	    part[i].id,
	    part[i].pos[X],part[i].pos[Y],
	    part[i].vel[X],part[i].vel[Y],
	    part[i].accel[X],part[i].accel[Y],
	    part[i].rho,part[i].mass,
	    part[i].p,part[i].c,part[i].u);

}
  
fclose(fState);
}

/* Graficacion con gnuplot:
   
>plot "kernel_test.output" u 1:2 w lp
>set size square
>plot "kernel_test.output" u 1:2 w lp
>set grid
>replot
>plot "kernel_test.output" u 1:2 w lp, "" u 1:3 w lp  

// Las comillas "" reemplazan el nombre anterior entre comillas

>plot "fluid_ics.output" u 2:3 w p lc rgb "blue" not

>plot "fluid_ics.output" u 2:3 w p lc rgb "cyan" not,
 "NN_test.output" u 2:3 w p pt 7 lc rgb "blue" not,
 "" every ::0::0 u 2:3 w lp pt 7 lc rgb "red" not

*/
