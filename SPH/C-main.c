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