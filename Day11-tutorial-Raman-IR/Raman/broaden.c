
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0

#define MAX_STEPS (100000)

typedef int Boolean;

void
defaults ( double *, double *, double *, int *,
	   Boolean *, Boolean *, Boolean *, double * );
void arguments ( int , char *[], FILE **, FILE ** ,
		 double *, double *, double *, int *,
		 Boolean *, Boolean *, Boolean *,
		 double * );
Boolean get_next_data ( Boolean , Boolean ,
		double *, double *, FILE * );

#define MAX_LEN 133
char str[MAX_LEN];

int main ( int argc, char *argv[] )
{
  double Xmin, Xmax, width, cnorm;
  double xval, yval, Inte, rad, cnt, dx, minvalue, arg;
  double dist_counts[MAX_STEPS], xcrd[MAX_STEPS], weight;
  FILE *fp_in = stdin, *fp_out = stdout;
  char *filename;
  int i, j, count, ch, step, max_channel, begin_channel, nsteps;
  void usage(char *[]);
  Boolean integrate, bindata, useyval;
  
  defaults ( &Xmin, &Xmax, &width, &nsteps, &integrate, &bindata, &useyval,
	     &minvalue );
  
  arguments ( argc, argv, &fp_in, &fp_out, 
	      &Xmin, &Xmax, &width, &nsteps, 
	      &integrate, &bindata, &useyval,
	      &minvalue );
  
  for ( step = 0; step < nsteps; step++ )
    xcrd[step] = (double)step * ( Xmax - Xmin ) / ( (double)nsteps - 1.0 )
      + Xmin;
  
  /* So that 'width' = FWHM */
  width = width / ( 2.0 * sqrt ( (double)2.0 * log ( (double)2.0 ) ) );
  
  /*
   * Read in the data
   */
  
  for ( step = 0; step < nsteps; step++ )
    dist_counts[step] = 0;
  
  dx = ( Xmax - Xmin ) / ( (double)nsteps - 1.0 );
  
  cnorm = width * sqrt ( 2.0 * M_PI );
  
  count = 0;
  
  /*
  begin_channel = (int) ( (double)nsteps * 25.0 * sqrt ( (double)2.0 )
			  / (Xmax-Xmin) + 0.5 ) ;
  */
  
  begin_channel = (int)
    ( sqrt ( -2.0 * width*width * log ( minvalue / cnorm ) ) / dx );
  
  while ( get_next_data ( useyval, bindata, &xval, &yval, fp_in ) == TRUE ) {
    
    max_channel = (int) ( (double)nsteps * ( xval - Xmin )
			  / ( Xmax - Xmin ) + 0.5 ) ;
    
    if ( width > 0.0 ) { /* Gaussian */
      for ( ch = max_channel - begin_channel;
	    ch <= max_channel + begin_channel; ch++ ) {
	if ( ( ch >= 0 ) && ( ch < nsteps ) ) {
	  
	  cnt = (double)ch * ( Xmax - Xmin ) / ( (double)nsteps - 1.0 )
	    + Xmin - xval;
	  
	  arg = cnt / width;
	  weight = yval * exp ( - arg*arg / 2.0 ) / cnorm;
	  
	  dist_counts[ch] += weight;
	}
      }
      
    } else {
      weight = yval;
      ch = max_channel;
      
      if ( ( ch >= 0 ) && ( ch < nsteps ) )
	dist_counts[ch] += weight;
    } /* if ( width > 0 ) ... else ... */
    
    count++;
  }
  
  fprintf(stderr, "Number of values: %d\n", count);
  
  Inte = 0.0;
  
  for ( step = 0; step < nsteps; step++ ) {
    Inte += dist_counts[step];
    if ( integrate && ( step == 0 ) ) dist_counts[0] = 0.0;
    if ( integrate && ( step > 0 ) )
      dist_counts[step] += dist_counts[step-1];
    rad = ( Xmax - Xmin ) * (double)step / (double)(nsteps-1) + Xmin;
    fprintf ( fp_out, "  %14.7lf  %14.7lf\n", rad, dist_counts[step] );
  }
  
  fprintf ( stderr, "  Integral: %14.7lf\n", Inte * dx );
  
  if ( fp_in != stdin )
    fclose ( fp_in );
}

/*****************************************************************************/

void usage ( char *argv[] )
{
  
  fprintf ( stderr, "Usage: %s\n"
  " [ -int ]        > Integrate pair correlation function\n"
  " [ -bin ]        > Read in data in a binary format (2*double)\n"
  " [ -width val ]  > Width of Gaussian broadening\n"
  " [ -ch val ]     > Number of channels to collect distances (def 1000)\n"
  " [ -Xmin val ]   > Smallest energy to search for a eigs (def 10)\n"
  " [ -Xmax val ]   > Largest energy to search for a eigs (def 10)\n"
  " [ -f file ]     > Input filename  (def stdin)\n"
  " [ -o file ]     > Output filename (def stdout)\n"
  " [ -useyval ]    > The second column is used\n"
  " [ -minvalue ]   > Smallest value included\n"), argv[0];
  
  exit ( 99 );
}

/*****************************************************************************/

void arguments ( int argc, char *argv[], FILE **fp_in, FILE **fp_out,
  double *Xmin, double *Xmax, double *width, int *nsteps,
  Boolean *integrate, Boolean *bindata, Boolean *useyval,
  double *minvalue )
{
  char *filename;
  int i;
  
  for ( i = 1; i < argc; i++ ) {
    
    if ( ! strcmp( "-int", argv[i] ) )
      *integrate = TRUE;
    
    else if ( ! strcmp( "-width", argv[i] ) ) {
      if ( i == argc-1 ) usage ( argv );
      *width = atof(argv[i+1]);
      i++;
    }
    
    else if ( ! strcmp( "-ch", argv[i] ) ) {
      if ( i == argc-1 ) usage(argv);
      *nsteps = atoi(argv[i+1]);
      if ( *nsteps > MAX_STEPS ) {
	fprintf ( stderr, "Maximum number of steps %d\n", MAX_STEPS );
	exit ( 1 );
      }
      i++;
    }
    
    else if ( ! strcmp( "-Xmin", argv[i] ) ) {
      if ( i == argc-1 ) usage(argv);
      *Xmin = atof(argv[i+1]);
      i++;
    }
    
    else if ( ! strcmp( "-minvalue", argv[i] ) ) {
      if ( i == argc-1 ) usage(argv);
      *minvalue = atof(argv[i+1]);
      i++;
    }
    
    else if ( ! strcmp( "-bin", argv[i] ) ) {
      *bindata = TRUE;
    }
    
    else if ( ! strcmp( "-useyval", argv[i] ) ) {
      *useyval = TRUE;
    }
    
    else if ( ! strcmp( "-nouseyval", argv[i] ) ) {
      *useyval = FALSE;
    }
    
    else if ( ! strcmp( "-Xmax", argv[i] ) ) {
      if ( i == argc-1 ) usage(argv);
      *Xmax = atof(argv[i+1]);
      i++;
    }
    
    else if ( ! strcmp( "-f", argv[i] ) ) {
      if ( i == argc-1 ) usage(argv);
      filename = argv[i+1];
      if ( ( (*fp_in) = fopen(filename, "r" ) ) == NULL ) {
	fprintf(stderr, "Could not open file %s\n", filename);
	exit(1);
      }
      i++;
    }
    
    else if ( ! strcmp( "-o", argv[i] ) ) {
      if ( i == argc-1 ) usage(argv);
      filename = argv[i+1];
      if ( ( (*fp_out) = fopen(filename, "w" ) ) == NULL ) {
	fprintf(stderr,"Could not open file %s\n", filename);
	exit(1);
      }
      i++;
    }
    
    else usage ( argv );
  }
}

/*****************************************************************************/

void
defaults ( double *Xmin, double *Xmax, double *width, int *nsteps,
	   Boolean *integrate, Boolean *bindata, Boolean *useyval,
	   double *minvalue )
{
  *Xmin =  0.0;
  *Xmax = 10.0;
  *width = 0;
  *nsteps = 10000;
  *integrate = FALSE;
  *bindata = FALSE;
  *useyval = TRUE;
  *minvalue = 1.0e-8;
}

/*****************************************************************************/

Boolean
get_next_data ( Boolean useyval, Boolean bindata,
		double *xval, double *yval, FILE *fp_in )
{
  int expectedvalues;
  Boolean successfulreturn;
  
  expectedvalues = ( useyval == TRUE ) ? 2 : 1;
  
  successfulreturn = FALSE;
  
  if ( bindata ) {
    if ( useyval ) {
      if ( ( fread(xval, sizeof(double), 1, fp_in) == 1 ) &&
	   ( fread(yval, sizeof(double), 1, fp_in) == 1 ) )
	successfulreturn = TRUE;
      
    } else {
      *yval = 1.0;
      if ( fread ( xval, sizeof(double), 1, fp_in) == 1 )
	successfulreturn = TRUE;
    }
    
  } else {
    
    if ( fgets ( str, MAX_LEN, fp_in ) != NULL ) {
      
      if ( useyval ) {
	if ( sscanf ( str, "%lf %lf", xval, yval ) == 2 )
	  successfulreturn = TRUE;
	
      } else {
	if ( sscanf ( str, "%lf", xval ) == 1 )
	  successfulreturn = TRUE;
      }
    }
    
  }
  
  return successfulreturn;
}
