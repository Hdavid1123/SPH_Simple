#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<malloc.h>

#define X 0
#define Y 1

/*
  Creamos el objeto celda
*/

typedef struct
{

  int id;
  double poscentro[2]; // Posición del centro de la celda
  int *idcelvec;         // Id de las cedas vecinas
  int ncelvec;         // Cantidad de celdas vecinas
  int *idmpart;     // Id de las partículas dentro de la celda
  int nmpart;       // Número de partículas dentro de la celda 
  
} Celdas;

/*
  Creamos el objeto partícula
*/

typedef struct
{
  
  int id;
  double pos[2];
  double vel[2];
  double accel[2];
  double mass, rho, h, p, c, du, u;
  int *nn;
  int nNeighbors;
  double *dx, *dy, *r, *W, *dWx, *dWy;
  int type;
  
} Particles;


// Variables globales
Celdas *celda;
Particles *part, *auxPart;
int nFluid, nPart, nCelda;

// Funciones secundarias

void iniconds(int nx, int ny, double dx, double dy,
	      double Lx, double Ly);

void Part_incel(int nGrid, int nFluid, double h);

double W(double r, double h);
double dW(double r, double dx, double h);
void testKernel(void);

void NN(int i, int ncel);
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



/* FUNCIÓN PRINCIPAL */


int main(int argc, char *argv[])
{

  int i, ncel, nx, ny, counter;
  double Lx, Ly, dx, dy;
  double dt = 5e-5; //Por defecto 5e-5
  double t, tTotal = atoi(argv[1])*dt;
  char outfiles[500];
  //double t, tTotal = 4000*dt;

  printf("Corriendo durante %d pasos, un tiempo total de %lf s\n",
         atoi(argv[1]),tTotal);

  nx = 40; // Número de partículas acomodadas
  ny = 40;
  Lx = 1e-3; // Tamaño del recuadro
  Ly = 1e-3;
  dx = Lx/nx;
  dy = Ly/ny;
  
  nFluid = nx*ny;

  // variables de las nuevas funciones
  int nGrid, partsincel;
  double h, cx, cy;

  h = dx;
  cx = Lx/(2*h);
  cy = Ly/(2*h);
  
  nGrid = cx * cy;
  
  // Alocando el espacio para las celdas
  celda = (Celdas *)malloc((size_t)nGrid*sizeof(Celdas));

  if( celda == NULL)
    {
      printf("Erro alocando las celdas\n");
      exit(0);
    }
  
  part = (Particles *)malloc((size_t)nFluid*sizeof(Particles));
  if( part==NULL )
    {
      printf("Error alocando part\n");
      exit(0);
    }

  /* Creando las condiciones iniciales */
  iniconds(nx, ny, dx, dy, Lx, Ly);

  /* Testeando la función de Kernell */
  testKernel();
  counter = 0;
  t = 0;

  /* Imprimiendo el estado inicial del sistema */
  sprintf(outfiles,"./ouput/state_%.4d",counter);
  printState(outfiles);

  
  /* Ciclo principal de la simulación */

  while( t<=tTotal)
    {
      
      Part_incel(nGrid, nFluid, h);

      // Buscamos los vecinos para todas las partículas de fluido
      for ( ncel=0; ncel<nGrid; ncel++)
	{
	  
	  partsincel = celda[i].nmpart;

	  for( i=0; i<partsincel; i++)
	    NN(i,ncel);
	    
	} // Final for i
    } // Final while
  
  
  // testing near neighbors searching
  if(counter==0)
    test_NN();
	      
  // computing density
  density();
  
  // drift in leap-frog integration
  drift(dt);
      
  // computing acceleration
  acceleration(dx);  
      
  // kick in leap-frog integration
  kick(dt);
      
  // drift in leap-frog integration
  drift(dt);
	
  t = t + dt;
  counter++;
  
  // printting system state
  sprintf(outfiles,"./output/state_%.4d",counter);
  printState(outfiles);
  
  printf("step = %d \n",counter);
  

  free(part);

  return 0;
  
} // FINAL FUNCIÓN PRINCIPAL


// ####################################################


/* FUNCIONES SECUNDARIAS */


// 1. Condiciones iniciales partículas

void iniconds(int nx, int ny, double dx, double dy,
	      double Lx, double Ly)
{

  int i, j, counter;
  double h, cx, cy;

  h = dx;
  cx = Lx/(2*h);
  cy = Ly/(2*h);

  FILE *fFluidIcs, *fbBorder, *frBorder, *ftBorder, *flBorder;
  fFluidIcs = fopen("fluid_ics.output","w");

  
  // 1.1 Para la malla de celdas

  // 1.1.1 Inicializando las variables
  
  for( i=0; i<cx; i++)
    {
      
      for( j=0; j<cy; j++)
	{
	  celda[counter].id = counter;
	  celda[counter].poscentro[X] = 2*i*h+h;
	  celda[counter].poscentro[Y] = 2*j*h+h;
	  celda[counter].idcelvec = NULL;
	  celda[counter].ncelvec = 0;
	  celda[counter].idmpart = NULL;
	  celda[counter].nmpart = 0;
	  counter++;
	} // Final for j
      
    } // Final for i

  
  // 1.1.2 Busqueda de celdas vecinas

  int Ncelvec, *auxInt;
  double distmax, kappa = 2.0;
  double xij, yij, rij;
  double eps = dx;
  double eps2 = 0.01*eps*eps;

  for( i=0; i<cx*cy; i++)
    {

      free(celda[i].idcelvec);
      celda[i].idcelvec = NULL;
      Ncelvec = 0;
      
      for( j=0; j<cx*cy; j++)
	{

	  // Calculamos la distancia entre los centros de la celda
	  xij = celda[i].poscentro[X] - celda[j].poscentro[X];
	  yij = celda[i].poscentro[Y] - celda[j].poscentro[Y];
	  rij = sqrt( xij*xij + yij*yij );

	  // Distancia máxima hallada con pitágoras + un eps
	  distmax = sqrt(2)*kappa*h + eps2;

	  if(rij <= distmax)
	    {

	      // Aumentamos en 1 el número de celdas vecinas
	      Ncelvec++;

	      auxInt = NULL;
	      auxInt = (int *)realloc(celda[i].idcelvec,
				      (size_t)(Ncelvec)*
				      sizeof(int));

	      celda[i].idcelvec = auxInt;
	      auxInt = NULL;

	      celda[i].idcelvec[Ncelvec-1] = j;
	    } // Final for if
	    
	} // Final for j

      celda[i].ncelvec = Ncelvec;
      
    } // Final for i
  
  
  // 1.2 Para las partículas de fluido

  counter = 0;
  for( j=0; j<ny; j++)
    {
      for( i=0; i<nx; i++)
	{	  
	  part[counter].id = counter;
	  part[counter].pos[X] = i*dx+dx/2.0;
	  part[counter].pos[Y] = j*dy+dy/2.0;
	  part[counter].vel[X] = 0.0;
	  part[counter].vel[Y] = 0.0;
	  part[counter].accel[X] = 0.0;
	  part[counter].accel[Y] = 0.0;
	  part[counter].rho = 1000;
	  part[counter].h = dx;
	  part[counter].mass = part[counter].rho*dx*dy;
	  part[counter].p = 0.0;
	  part[counter].c = 0.0;
	  part[counter].du = 0.0;
	  part[counter].u = 357.1;
	  part[counter].nn = NULL;
	  part[counter].nNeighbors = 0;
	  part[counter].dx = NULL;
	  part[counter].dy = NULL;
	  part[counter].r = NULL;
	  part[counter].W = NULL;
	  part[counter].dWx = NULL;
	  part[counter].dWy = NULL;
	  part[counter].type = 1;
	  counter++;
	  //El tipo fluido se representa como 1
	  
	} //Final for i
    } // Final for j
  
  
  // 1.2. Para las bartículas de los bordes.
  
  /* Velocidad que dan los bordes */
  double vBoundary = 1.5e-2;
  
  int npVirtI = 320;
  int npV = npVirtI/4;
  
  
  // Borde inferior, 81 puntos
  
  fbBorder = fopen("bottom_border.output","w");
  
  /* Las partículas de contorno se agregan al conjunto de
     partículas de estructura */
  
  nPart = nFluid;
  auxPart = NULL;
  
  auxPart = (Particles *)realloc(part, (size_t)(nPart+npV+1)*
				 sizeof(Particles));
  if (auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(0);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }

  // Contador = dimensiones de la caja, en este caso 40x40
  counter = nPart;
  printf("El contador en los bordes inicia en: %d", counter);

  
  for( i=0; i<=npV; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = i*dx/2.0;
      part[counter].pos[Y] = j*dy/2.0;
      part[counter].vel[X] = 0.0;
      //part[counter].vel[X] = -vBoundary;
      part[counter].vel[Y] = 0.0;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    } //Final for i
  
  for( i=nPart; i<nPart+npV+1; i++ )
    {
      fprintf(fbBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }

  fclose(fbBorder);

  // Borde derecho, 79 puntos

  frBorder = fopen("right_border.output","w");

  nPart = counter;
  auxPart = NULL;

  auxPart = (Particles *)realloc(part, (size_t)(nPart+npV-1)*
				 sizeof(Particles));

  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(0);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }

  for( i=0; i<npV-1; i++ )
    {
      part[counter].id = counter;
      part[counter].pos[X] = Lx;
      part[counter].pos[Y] = dy/2.0 + i*dy/2.0;
      part[counter].vel[X] = 0.0;
      part[counter].vel[Y] = 0.0;
      //part[counter].vel[Y] = -vBoundary;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }

  for( i=nPart; i<nPart+npV-1; i++ )
    {
      fprintf(frBorder, "%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }

  fclose(frBorder);

  // Borde superior, 81 puntos

  ftBorder = fopen("top_border.ouput","w");

  nPart = counter;
  auxPart = NULL;

  auxPart = (Particles *)realloc(part, (size_t)(nPart+npV+1)*
				 sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(0);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }

  for( i=0; i<=npV; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = i*dx/2.0;
      part[counter].pos[Y] = Ly;
      part[counter].vel[X] = vBoundary;
      part[counter].vel[Y] = 0.0;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }

  for( i=nPart; i<nPart+npV+1; i++)
    {
      fprintf(ftBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }

  fclose(ftBorder);

  // Borde izquierdo, 79 puntos

  flBorder = fopen("left_border.ouput","w");

  nPart = counter;
  auxPart = NULL;

  auxPart = (Particles *)realloc(part, (size_t)(nPart+npV-1)*
				 sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(0);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }

  for( i=0; i<npV-1; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = 0.0;
      part[counter].pos[Y] = dy/2.0 + i*dy/2.0;
      part[counter].vel[X] = 0.0;
      part[counter].vel[Y] = 0.0;
      //part[counter].vel[Y] = vBoundary;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }

  for( i=nPart; i<nPart+npV-1; i++)
    {
      fprintf(flBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }

  fclose(flBorder);

  // Imprimimos todas las partículas

  nPart = counter;

  for(i=0; i<nPart; i++)
    {
      fprintf(fFluidIcs,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    } // Final for i
  
  fclose(fFluidIcs);
  
} // Final func. cond. iniciales



// 2. Asignador de partículas para cada celda

void Part_incel(int nGrid, int nFluid, double h)
{

  int i,j, Nparts, *auxInt;
  double xc, yc, xp, yp, bool1, bool2;

  for( i=0; i<=nGrid; i++)
    {

      free(celda[i].idmpart);

      celda[i].idmpart = NULL;
      Nparts = 0;
      
      for( j=0; j<(nFluid+320); j++) // 320 No. part. de frontera
	{
	  
	  // Damos las posiciones de celdas y partículas
	  xc = celda[i].poscentro[X];
	  yc = celda[i].poscentro[Y];

	  xp = part[j].pos[X];
	  yp = part[j].pos[Y];

	  bool1 = fabs(xp-xc);
	  bool2 = fabs(yp-yc);

	  //if( (bool1 && bool2) < h + h*1e-5)
	  if( (bool1 < h+h/100000) && (bool2 < h+h/100000) )
	    {
	      
	      Nparts++;

	      auxInt = NULL;
	      auxInt = (int *)realloc(celda[i].idmpart,
				      (size_t)(Nparts)*
				      sizeof(int));
	      celda[i].idmpart = auxInt;
	      auxInt = NULL;

	      celda[i].idmpart[Nparts-1] = j;
	      
	    } // Final if
	 
	} // Final for j

      celda[i].nmpart = Nparts;
      
    } // Final for i
  
} // Final función asignadora de partículas



// 3. Función Kernell
double W(double r, double h)
{
  
  double R = r/h;
  double alpha = 15.0/(7.0*M_PI*h*h);

  if( (R >= 0.0) && (R < 1.0) )
    return alpha*((2.0/3.0) - R*R + 0.5*R*R*R);

  if( (R >= 1.0) && (R < 2.0) )
    return alpha*((1.0/6.0)*(2.0-R)*(2.0-R)*(2.0-R));

  if( R > 2.0 )
    return 0.0;

  return 0.0;
    
} // Final func. Kernell



// 4. Función derivada de Kernell
double dW(double r, double dx, double h)
{
  
  double R = r/h;
  double alpha = 15.0/(7.0*M_PI*h*h);

  if( (R >= 0.0) && (R < 1.0) )
    return alpha*(-2.0 + 1.5*R)*dx/(h*h);

  if( (R >= 1.0) && (R < 2.0) )
    return alpha*(-0.5*(2.0-R)*(2.0-R))*dx/(h*h*R);

  if( R > 2.0 )
    return 0.0;

  return 0.0;
  
} // Final func. der. Kernell



// 5. Test func. de Kernell
void testKernel(void)
{
  double r, w, dw;

  FILE *fKernelTest;
  fKernelTest = fopen("kernel_test.output","w");

  for( r=-3.0; r<3.0; r+=0.1)
    {
      w = W( fabs(r), 1.0);
      dw = dW( fabs(r), r/sqrt(3.0), 1.0);

      fprintf(fKernelTest,"%16.10lf %16.10lf %16.10lf\n",r,w,dw);
    }

  fclose(fKernelTest);
} // Final test func. Kernell



// 6. Buscador de vecinos. Usaremos linked-list

void NN(int i, int ncel)
{
  double kappa = 2.0;
  double xij, yij, rij, hij;
  double *auxDouble;
  int *auxInt, nNeighbors, partsin;
  int j, k, l, m;
  int allceldas, idvecina;
    
  k = celda[i].idmpart[ncel];
  free(part[k].nn);
  free(part[k].dx);
  free(part[k].dy);
  free(part[k].r);
  free(part[k].W);
  free(part[k].dWx);
  free(part[k].dWy);

  part[k].nn = NULL;
  part[k].dx = NULL;
  part[k].dy = NULL;
  part[k].r = NULL;
  part[k].W = NULL;
  part[k].dWx = NULL;
  part[k].dWy = NULL;
  nNeighbors = 0;

  // Determinamos el número de vecinas de la celda sobre la que
  // estamos
  allceldas = celda[i].ncelvec;

  for( l=0; l<allceldas; l++)
    {
      // Extraemos el id de cada una de las celdas vecinas
      idvecina = celda[i].idcelvec[l];

      // Y determinamos el número de partículas en esta
      partsin = celda[idvecina].nmpart;

      for( m=0; m<partsin; m++)
	{
	  
	  j = celda[idvecina].idmpart[m];
	  if (i!=j)
	    {
	      xij = part[k].pos[X] - part[j].pos[X];
	      yij = part[k].pos[Y] - part[j].pos[Y];
	      rij = sqrt( xij*xij + yij*yij );
	      hij = 0.5*(part[k].h + part[j].h);

	      if( rij < kappa*hij)
		{
		  nNeighbors++;

		  // add neighbor id
                  auxInt = NULL;
                  auxInt = (int *)realloc(part[k].nn,
					  (size_t)(nNeighbors)*
					  sizeof(int));
                  part[k].nn = auxInt;
                  auxInt = NULL;

		  
                  // add neighbor dx
                  auxDouble = NULL;
                  auxDouble = (double *)realloc(part[k].dx,
					       (size_t)(nNeighbors)*
					       sizeof(double));
                  part[k].dx = auxDouble;
                  auxDouble = NULL;

		  
                  // add neighbor dy
                  auxDouble = NULL;
                  auxDouble = (double *)realloc(part[k].dy,
					       (size_t)(nNeighbors)*
					       sizeof(double));
                  part[k].dy = auxDouble;
                  auxDouble = NULL;

		  
                  // add neighbor r
                  auxDouble = NULL;
                  auxDouble = (double *)realloc(part[k].r,
					       (size_t)(nNeighbors)*
					       sizeof(double));
                  part[k].r = auxDouble;
                  auxDouble = NULL;

		  
                  // add neighbor W
                  auxDouble = NULL;
                  auxDouble = (double *)realloc(part[k].W,
					       (size_t)(nNeighbors)*
					       sizeof(double));
                  part[k].W = auxDouble;
                  auxDouble = NULL;

		  
                  // add neighbor dWx
                  auxDouble = NULL;
                  auxDouble = (double *)realloc(part[k].dWx,
					       (size_t)(nNeighbors)*
					       sizeof(double));
                  part[k].dWx = auxDouble;
                  auxDouble = NULL;

		  
                  // add neighbor dWy
                  auxDouble = NULL;
                  auxDouble = (double *)realloc(part[k].dWy,
					       (size_t)(nNeighbors)*
					       sizeof(double));
                  part[k].dWy = auxDouble;
                  auxDouble = NULL;
		  
		  part[k].nn[nNeighbors-1] = j;
                  part[k].dx[nNeighbors-1] = xij;
                  part[k].dy[nNeighbors-1] = yij;
                  part[k].r[nNeighbors-1] = rij;
                  part[k].W[nNeighbors-1] = W( rij, hij ); 
                  part[k].dWx[nNeighbors-1] = dW( rij, xij, hij);
                  part[k].dWy[nNeighbors-1] = dW( rij, yij, hij);
		  
		} // Final if
	      
	    } // Final if

	} // Final for j

      part[k].nNeighbors = nNeighbors;
      
    } // Final for l
} // Final buscador de vecinos



// 7. Test del buscador de partículas

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
    }
  fclose(fTestNN);
    
} // Final test NN




// 8. Función de densidad

void density(void)
{
  int i, j;
  double wii, norm;
  
  for( i=0; i<nFluid; i++ )
    {
      // self density
      wii = W( 0.0, part[i].h );
      
      // computing density
      part[i].rho = part[i].mass*wii;
      for( j=0; j<part[i].nNeighbors; j++ )
	part[i].rho = part[i].rho + part[part[i].nn[j]].mass*
	  part[i].W[j];

      // normalizing the density
      norm = (part[i].mass/part[i].rho)*wii;
      for( j=0; j<part[i].nNeighbors; j++ )
	  norm = norm + (part[part[i].nn[j]].mass/
			 part[part[i].nn[j]].rho)*part[i].W[j];
      
      part[i].rho = part[i].rho/norm;
    }
  
  printf("density computed\n");
}



// 9. Función eos

void eos(void)
{
  int i;
  
  for( i=0; i<nPart; i++ )
    {
      part[i].c = 0.01;
      part[i].p = part[i].c*part[i].c*part[i].rho; 
    }
}



// 10. Función ecuaciones de Navier Stokes

void navierStokes(void)
{

  int i, j, k;
  double pij, vdw;
  // computing sound speed and pression
  eos();

  // computing acceleration
  for( i=0; i<nFluid; i++ )
    {
      
      part[i].accel[X] = part[i].accel[Y] = 0.0;
      part[i].du = 0.0;

      for( k=0; k<part[i].nNeighbors; k++ )
	{
	  j = part[i].nn[k];
	  
	  pij = ( part[i].p/(part[i].rho*part[i].rho) )
	    + ( part[j].p/(part[j].rho*part[j].rho) );

	  part[i].accel[X] = part[i].accel[X] - part[j].mass*
	    pij*part[i].dWx[k];

	  part[i].accel[Y] = part[i].accel[Y] - part[j].mass*
	    pij*part[i].dWy[k];
	  
	  vdw = (part[i].vel[X]-part[j].vel[X])*part[i].dWx[k]
	    + (part[i].vel[Y]-part[j].vel[Y])*part[i].dWy[k];

	  part[i].du = part[i].du + 0.5*part[j].mass*pij*vdw;
	  
	} 
    }
  
  printf("acceleration computed\n");
}



// 11. Función de viscosidad

void viscosity(double dx)
{
  
  int i, j, k;

  double xij, yij, vxij, vyij, vijrij, vdw;
  double hij, cij, phiij, rhoij, Piij;
  double alphapi = 1.0;
  double betapi = 1.0;
  double eps = dx;
  double eps2 = 0.01*eps*eps;

  for( i=0; i<nFluid; i++ )
    {
      for( k=0; k<part[i].nNeighbors ; k++ )
	{
	  
	  j = part[i].nn[k];
	  
	  xij = part[i].pos[X] - part[j].pos[X];
	  yij = part[i].pos[Y] - part[j].pos[Y];
	  vxij = part[i].vel[X] - part[j].vel[X];
	  vyij = part[i].vel[Y] - part[j].vel[Y];
	  vijrij = vxij*xij + vyij*yij;
	  
	  if( vijrij < 0.0 )
	    {
	      hij = 0.5*(part[i].h+part[j].h);
	      phiij = (hij*vijrij)/( xij*xij + yij*yij + eps2);
	      cij = 0.5*(part[i].c+part[j].c);
	      rhoij = 0.5*(part[i].rho+part[j].rho);
	      
	      Piij = ( -alphapi*cij*phiij + betapi*phiij*phiij )/
		( rhoij );

	      part[i].accel[X] = part[i].accel[X] - part[j].mass*
		Piij*part[i].dWx[k];

	      part[i].accel[Y] = part[i].accel[Y] - part[j].mass*
		Piij*part[i].dWy[k];

	      vdw = (part[i].vel[X]-part[j].vel[X])*part[i].dWx[k]
		+ (part[i].vel[Y]-part[j].vel[Y])*part[i].dWy[k];

	      part[i].du = part[i].du + 0.5*part[j].mass*Piij*vdw;
	      
	    }
	  
	}
    }
  printf("viscosity computed\n");
}



// 12. Función de interacción con los bordes

void boundaryInteraction(double dx)
{

  int i, j;
  int n1 = 12, n2 = 4;
  double r0 = dx/2.0, D = 0.01;
  double xij, yij, rij, PBxij, PByij;

  for( i=0; i<nFluid; i++ )
    {
      for( j=0; j<part[i].nNeighbors; j++ )
	{
	  if( part[part[i].nn[j]].type==-1 )
	    {
	      xij = part[i].pos[X] - part[part[i].nn[j]].pos[X];
	      yij = part[i].pos[Y] - part[part[i].nn[j]].pos[Y];
	      rij = sqrt( xij*xij + yij*yij );
	      
	      if( rij<r0 )
		{
		  PBxij = D*( pow((r0/rij),n1) - pow((r0/rij),n2))*
		    (xij/(rij*rij));

		  PByij = D*( pow((r0/rij),n1) - pow((r0/rij),n2))*
		    (yij/(rij*rij));
		  
		  part[i].accel[X] = part[i].accel[X] + PBxij;
		  part[i].accel[Y] = part[i].accel[Y] + PByij;
		  
		}
	    }
	}
    }
  printf("interaction with boundary computed\n");
}



// 13. Función de velocidad media

void meanVelocity(void)
{

  int i, j;
  double epsilon = 0.3;
  double vxMean, vyMean;
  double vxij, vyij, rhoij;

  for( i=0; i<nFluid; i++ )
    {
      
      vxMean = 0.0;
      vyMean = 0.0;
      
      for( j=0; j<part[i].nNeighbors; j++ )
	{
	  vxij = part[i].vel[X] - part[part[i].nn[j]].vel[X];
	  vyij = part[i].vel[Y] - part[part[i].nn[j]].vel[Y];
	  rhoij = 0.5*(part[i].rho+part[part[i].nn[j]].rho);

	  vxMean = vxMean + (part[part[i].nn[j]].mass/rhoij)*
	    vxij*part[i].W[j];

	  vyMean = vyMean + (part[part[i].nn[j]].mass/rhoij)*
	    vyij*part[i].W[j];
	}

      part[i].vel[X] = part[i].vel[X] - epsilon*vxMean;
      part[i].vel[Y] = part[i].vel[Y] - epsilon*vyMean;
   
    }
}



// 14. Función de aceleración

void acceleration(double dx)
{

  // computing acceleration and change of energy
  navierStokes();
  
  // computing viscosity contribution
  viscosity(dx);
  
  // computing interaction with boundary
  boundaryInteraction(dx);

  // correction to mean velocity
  meanVelocity();

  printf("acceleration computed\n");
  
}



// 15. Funciones drift, kick y printState

void drift(double dt)
{

  int i;
  
  for( i=0; i<nFluid; i++ )
    {
      part[i].pos[X] = part[i].pos[X] + 0.5*dt*part[i].vel[X];
      part[i].pos[Y] = part[i].pos[Y] + 0.5*dt*part[i].vel[Y];
      part[i].u = part[i].u + 0.5*dt*part[i].du;
    }
}

void kick(double dt)
{
  
  int i;
  
  for( i=0; i<nFluid; i++ )
    {
      part[i].vel[X] = part[i].vel[X] + dt*part[i].accel[X];
      part[i].vel[Y] = part[i].vel[Y] + dt*part[i].accel[Y];
    }
}

void printState(char *outfile)
{

  int i;

  FILE *fState;
  fState = fopen(outfile,"w");
 
  for( i=0; i<nPart; i++)
    {
      fprintf(fState,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  
  fclose(fState);

}
