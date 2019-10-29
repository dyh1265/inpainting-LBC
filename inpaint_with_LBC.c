#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                    INPAINTING WITH LOCAL BINARY CONSTRAINTS              */
/*        This code implements the Algorithm for Level-n inpainting with    */
/*                           Local Binary Constraints                       */
/*                        (Andrei Sirazitdinov 10/2019)                     */
/*                                                                          */
/*--------------------------------------------------------------------------*/


void alloc_long_matrix

  (long   ***matrix,  /* matrix */
   long   n1,         /* size in direction 1 */
   long   n2)         /* size in direction 2 */

/*
   allocates memory for a matrix of size n1 * n2 in long format
*/

{
long i;    /* loop variable */

*matrix = (long **) malloc (n1 * sizeof(long *));

if (*matrix == NULL)
  {
  printf ("alloc_long_matrix: not enough memory available\n");
  exit (1);
  }

for (i=0; i<n1; i++)
  {
  (*matrix)[i] = (long *) malloc (n2 * sizeof(long));
  if ((*matrix)[i] == NULL)
    {
    printf ("alloc_long_matrix: not enough memory available\n");
    exit (1);
    }
  }

return;

}  /* alloc_long_matrix */

/*--------------------------------------------------------------------------*/

void alloc_double_matrix

  (double ***matrix,  /* matrix */
   long   n1,         /* size in direction 1 */
   long   n2)         /* size in direction 2 */

/*
   allocates memory for a matrix of size n1 * n2 in double format
*/

{
long i;    /* loop variable */

*matrix = (double **) malloc (n1 * sizeof(double *));

if (*matrix == NULL)
  {
  printf ("alloc_double_matrix: not enough memory available\n");
  exit (1);
  }

for (i=0; i<n1; i++)
  {
  (*matrix)[i] = (double *) malloc (n2 * sizeof(double));
  if ((*matrix)[i] == NULL)
    {
    printf ("alloc_double_matrix: not enough memory available\n");
    exit (1);
    }
  }

return;

}  /* alloc_double_matrix */


/*--------------------------------------------------------------------------*/

void alloc_long_cubix

     (long   ****cubix,  /* cubix */
      long   n1,         /* size in direction 1 */
      long   n2,         /* size in direction 2 */
      long   n3)         /* size in direction 3 */

/*
   allocates memory for a long format cubix of size n1 * n2 * n3
*/

{
long i, j;  /* loop variables */

*cubix = (long ***) malloc (n1 * sizeof(long **));

if (*cubix == NULL)
    {
    printf ("alloc_long_cubix: not enough memory available\n");
    exit (1);
    }

for (i=0; i<n1; i++)
    {
    (*cubix)[i] = (long **) malloc (n2 * sizeof(long *));
    if ((*cubix)[i] == NULL)
        {
    printf ("alloc_long_cubix: not enough memory available\n");
    exit (1);
    }
  for (j=0; j<n2; j++)
    {
    (*cubix)[i][j] = (long *) malloc (n3 * sizeof(long));
    if ((*cubix)[i][j] == NULL)
      {
      printf ("alloc_long_cubix: not enough memory available\n");
      exit (1);
      }
    }
  }

return;

}  /* alloc_long_cubix */

/*--------------------------------------------------------------------------*/

void alloc_double_cubix

  (double ****cubix,  /* cubix */
   long   n1,         /* size in direction 1 */
   long   n2,         /* size in direction 2 */
   long   n3)         /* size in direction 3 */

/*
   allocates memory for a double format cubix of size n1 * n2 * n3
*/

{
long i, j;  /* loop variables */

*cubix = (double ***) malloc (n1 * sizeof(double **));

if (*cubix == NULL)
  {
  printf ("alloc_double_cubix: not enough memory available\n");
  exit (1);
  }

for (i=0; i<n1; i++)
  {
  (*cubix)[i] = (double **) malloc (n2 * sizeof(double *));
  if ((*cubix)[i] == NULL)
    {
    printf ("alloc_double_cubix: not enough memory available\n");
    exit (1);
    }
  for (j=0; j<n2; j++)
    {
    (*cubix)[i][j] = (double *) malloc (n3 * sizeof(double));
    if ((*cubix)[i][j] == NULL)
      {
      printf ("alloc_double_cubix: not enough memory available\n");
      exit (1);
      }
    }
  }

return;

}  /* alloc_double_cubix */

/*--------------------------------------------------------------------------*/

void free_double_cubix

  (double ***cubix,   /* cubix */
   long   n1,         /* size in direction 1 */
   long   n2,         /* size in direction 2 */
   long   n3)         /* size in direction 3 */

/*
   frees memory for a double format cubix of size n1 * n2 * n3
*/

{
long i, j;   /* loop variables */

for (i=0; i<n1; i++)
  for (j=0; j<n2; j++)
    free (cubix[i][j]);

for (i=0; i<n1; i++)
  free (cubix[i]);

free (cubix);

return;

}  /* free_double_cubix */

/*--------------------------------------------------------------------------*/

void free_long_cubix

  (long   ***cubix,   /* cubix */
   long   n1,         /* size in direction 1 */
   long   n2,         /* size in direction 2 */
   long   n3)         /* size in direction 3 */

/*
   frees memory for a long format cubix of size n1 * n2 * n3
*/

{
long i, j;   /* loop variables */

for (i=0; i<n1; i++)
  for (j=0; j<n2; j++)
    free (cubix[i][j]);

for (i=0; i<n1; i++)
  free (cubix[i]);

free (cubix);

return;

}  /* free_long_cubix */

/*--------------------------------------------------------------------------*/
void alloc_double_vector

  (double **vector,   /* vector */
   long   n1)         /* size */

/*
   allocates memory for a double format vector of size n1
*/

{
*vector = (double *) malloc (n1 * sizeof(double));

if (*vector == NULL)
{
  printf ("alloc_double_vector: not enough memory available\n");
  exit (1);
}

return;

}  /* alloc_double_vector */

/*--------------------------------------------------------------------------*/

void alloc_long_vector

  (long   **vector,   /* vector */
   long   n1)         /* size */

/*
   allocates memory for a long format vector of size n1
*/

{
*vector = (long *) malloc (n1 * sizeof(long));

if (*vector == NULL)
{
  printf ("alloc_long_vector: not enough memory available\n");
  exit (1);
}

return;

}  /* alloc_long_vector */

/*--------------------------------------------------------------------------*/
void free_double_vector

  (double  *vector,    /* vector */
   long    n1)         /* size */

/*
   frees memory for a double format vector of size n1
*/

{

free (vector);

return;

}  /* free_double_vector */

/*--------------------------------------------------------------------------*/
void free_long_vector

  (long  *vector,    /* vector */
   long    n1)       /* size */

/*
   frees memory for a double format vector of size n1
*/

{

free (vector);

return;

}  /* free_double_vector */

/*--------------------------------------------------------------------------*/
void free_long_matrix

  (long    **matrix,   /* matrix */
   long    n1,         /* size in direction 1 */
   long    n2)         /* size in direction 2 */

/*
   frees memory for a matrix of size n1 * n2 in long format
*/

{
long i;   /* loop variable */

for (i=0; i<n1; i++)
  free (matrix[i]);

free (matrix);

return;

}  /* free_long_matrix */

/*--------------------------------------------------------------------------*/

void free_double_matrix

  (double  **matrix,   /* matrix */
   long    n1,         /* size in direction 1 */
   long    n2)         /* size in direction 2 */

/*
   frees memory for a matrix of size n1 * n2 in double format
*/

{
long i;   /* loop variable */

for (i=0; i<n1; i++)
  free (matrix[i]);

free (matrix);

return;

}  /* free_double_matrix */

/*--------------------------------------------------------------------------*/

void read_string

  (char *v)         /* string to be read */

/*
   reads a string v
*/

{
fgets (v, 80, stdin);

if (v[strlen (v)-1] == '\n')
  v[strlen (v)-1] = 0;

return;

}  /* read_string */

/*--------------------------------------------------------------------------*/

void read_long

  (long *v)         /* value to be read */

/*
   reads a long value v
*/

{
char   row[80];    /* string for reading data */

fgets (row, 80, stdin);
if (row[strlen (row)-1] == '\n')
  row[strlen (row)-1] = 0;
sscanf (row, "%ld", &*v);

return;

}  /* read_long */

/*--------------------------------------------------------------------------*/

void read_double

  (double *v)         /* value to be read */

/*
   reads a double value v
*/

{
char   row[80];    /* string for reading data */

fgets (row, 80, stdin);

if (row[strlen (row)-1] == '\n')
  row[strlen (row)-1] = 0;
sscanf (row, "%lf", &*v);

return;

}  /* read_double */

/*--------------------------------------------------------------------------*/

void skip_white_space_and_comments

  (FILE *inimage)  /* input file */

/*
   skips over white space and comments while reading the file
*/

{

int   ch = 0;   /* holds a character */
char  row[80];  /* for reading data */

/* skip spaces */
while (((ch = fgetc (inimage)) != EOF) && isspace (ch));

/* skip comments */
if (ch == '#')
{
  if (fgets (row, sizeof(row), inimage))
    skip_white_space_and_comments (inimage);
  else
    {
    printf ("skip_white_space_and_comments: cannot read file\n");
    exit (1);
    }
}
else
  fseek (inimage, -1, SEEK_CUR);

return;

} /* skip_white_space_and_comments */

/*--------------------------------------------------------------------------*/

void read_pgm_to_long

  (const char  *file_name,    /* name of pgm file */
   long        *nx,           /* image size in x direction, output */
   long        *ny,           /* image size in y direction, output */
   long        ***u)          /* image, output */

/*
   reads a greyscale image that has been encoded in pgm format P5 to
   an image u in long format;
   allocates memory for the image u;
   adds boundary layers of size 1 such that
   - the relevant image pixels in x direction use the indices 1,...,nx
   - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
char  row[80];      /* for reading data */
long  i, j;         /* image indices */
long  max_value;    /* maximum color value */
FILE  *inimage;     /* input file */

/* open file */
inimage = fopen (file_name, "rb");
if (inimage == NULL)
  {
  printf ("read_pgm_to_long: cannot open file '%s'\n", file_name);
  exit (1);
  }

/* read header */
if (fgets(row, 80, inimage) == NULL)
  {
  printf ("read_pgm_to_long: cannot read file\n");
  exit (1);
  }

/* image type: P5 */
if ((row[0] == 'P') && (row[1] == '5'))
  {
  /* P5: grey scale image */
  }
else
  {
  printf ("read_pgm_to_long: unknown image format\n");
  exit (1);
  }

/* read image size in x direction */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", nx))
  {
  printf ("read_pgm_to_long: cannot read image size nx\n");
  exit (1);
  }

/* read image size in y direction */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", ny))
  {
  printf ("read_pgm_to_long: cannot read image size ny\n");
  exit (1);
  }

/* read maximum grey value */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", &max_value))
  {
  printf ("read_pgm_to_long: cannot read maximal value\n");
  exit (1);
  }
fgetc(inimage);

/* allocate memory */
alloc_long_matrix (u, (*nx)+2, (*ny)+2);

/* read image data row by row */
for (j=1; j<=(*ny); j++)
  for (i=1; i<=(*nx); i++)
    (*u)[i][j] = (long) getc (inimage);

/* close file */
fclose (inimage);

return;

}  /* read_pgm_to_long */

/*--------------------------------------------------------------------------*/

void read_pgm_to_double

  (const char  *file_name,    /* name of pgm file */
   long        *nx,           /* image size in x direction, output */
   long        *ny,           /* image size in y direction, output */
   double      ***u)          /* image, output */

/*
   reads a greyscale image that has been encoded in pgm format P5 to
   an image u in double format;
   allocates memory for the image u;
   adds boundary layers of size 1 such that
   - the relevant image pixels in x direction use the indices 1,...,nx
   - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
char  row[80];      /* for reading data */
long  i, j;         /* image indices */
long  max_value;    /* maximum color value */
FILE  *inimage;     /* input file */

/* open file */
inimage = fopen (file_name, "rb");
if (inimage == NULL)
  {
  printf ("read_pgm_to_double: cannot open file '%s'\n", file_name);
  exit (1);
  }

/* read header */
if (fgets (row, 80, inimage) == NULL)
  {
  printf ("read_pgm_to_double: cannot read file\n");
  exit (1);
  }

/* image type: P5 */
if ((row[0] == 'P') && (row[1] == '5'))
  {
  /* P5: grey scale image */
  }
else
  {
  printf ("read_pgm_to_double: unknown image format\n");
  exit (1);
  }

/* read image size in x direction */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", nx))
  {
  printf ("read_pgm_to_double: cannot read image size nx\n");
  exit (1);
  }

/* read image size in x direction */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", ny))
  {
  printf ("read_pgm_to_double: cannot read image size ny\n");
  exit (1);
  }

/* read maximum grey value */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", &max_value))
   {
  printf ("read_pgm_to_double: cannot read maximal value\n");
  exit (1);
   }
fgetc(inimage);

/* allocate memory */
alloc_double_matrix (u, (*nx)+2, (*ny)+2);

/* read image data row by row */
for (j=1; j<=(*ny); j++)
  for (i=1; i<=(*nx); i++)
    (*u)[i][j] = (double) getc (inimage);

/* close file */
fclose (inimage);

return;

}  /* read_pgm_to_double */

/*--------------------------------------------------------------------------*/

void comment_line

  (char* comment,       /* comment string (output) */
   char* lineformat,    /* format string for comment line */
   ...)                 /* optional arguments */

/*
   Adds a line to the comment string comment. The string line can contain
   plain text and format characters that are compatible with sprintf.
   Example call:
   print_comment_line(comment, "Text %lf %ld", double_var, long_var).
   If no line break is supplied at the end of the input string, it is
   added automatically.
*/

{
char     line[80];
va_list  arguments;

/* get list of optional function arguments */
va_start (arguments, lineformat);

/* convert format string and arguments to plain text line string */
vsprintf (line, lineformat, arguments);

/* add line to total commentary string */
strncat (comment, line, 80);

/* add line break if input string does not end with one */
if (line[strlen (line)-1] != '\n')
  sprintf (comment, "%s\n", comment);

/* close argument list */
va_end (arguments);

return;

}  /* comment_line */

/*--------------------------------------------------------------------------*/

void write_double_to_pgm

  (double  **u,          /* image, unchanged */
   long    nx,           /* image size in x direction */
   long    ny,           /* image size in y direction */
   char    *file_name,   /* name of pgm file */
   char    *comments)    /* comment string (set 0 for no comments) */

/*
   writes a greyscale image in double format into a pgm P5 file
*/

{
FILE           *outimage;  /* output file */
long           i, j;       /* loop variables */
double         aux;        /* auxiliary variable */
unsigned char  byte;       /* for data conversion */

/* open file */
outimage = fopen (file_name, "wb");
if (NULL == outimage)
{
  printf ("Could not open file '%s' for writing, aborting\n", file_name);
  exit (1);
}

/* write header */
fprintf (outimage, "P5\n");                  /* format */
if (comments != 0)
  fprintf (outimage, "%s", comments);        /* comments */
fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");                 /* maximal value */

/* write image data */
for (j=1; j<=ny; j++)
  for (i=1; i<=nx; i++)
    {
    aux = u[i][j] + 0.499999;    /* for correct rounding */
    if (aux < 0.0)
      byte = (unsigned char)(0.0);
    else if (aux > 255.0)
      byte = (unsigned char)(255.0);
    else
      byte = (unsigned char)(aux);

    fwrite (&byte, sizeof(unsigned char), 1, outimage);
    }

/* close file */
fclose (outimage);

return;

}  /* write_double_to_pgm */

/*--------------------------------------------------------------------------*/

void write_long_to_pgm

  (long    **u,          /* image, unchanged */
   long    nx,           /* image size in x direction */
   long    ny,           /* image size in y direction */
   char    *file_name,   /* name of pgm file */
   char    *comments)    /* comment string (set 0 for no comments) */

/*
   writes a greyscale image in long format into a pgm P5 file;
*/

{
FILE           *outimage;  /* output file */
long           i, j;       /* loop variables */
unsigned char  byte;       /* for data conversion */

/* open file */
outimage = fopen (file_name, "wb");
if (NULL == outimage)
  {
  printf ("could not open file '%s' for writing, aborting\n", file_name);
  exit (1);
  }

/* write header */
fprintf (outimage, "P5\n");                  /* format */
if (comments != 0)
  fprintf (outimage, "%s", comments);        /* comments */
fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");                 /* maximal value */

/* write image data */
for (j=0; j<ny; j++)
  for (i=0; i<nx; i++)
    {
    if (u[i][j] < 0)
      byte = (unsigned char)(0);
    else if (u[i][j] > 255)
      byte = (unsigned char)(255);
    else
      byte = (unsigned char)(u[i][j]);
    fwrite (&byte, sizeof(unsigned char), 1, outimage);
    }

/* close file */
fclose (outimage);

return;

}  /* write_long_to_pgm */

/*--------------------------------------------------------------------------*/

void set_matrix_values_double

  (double  **u,           /* in+out :matrix (changed) */
   long    nx,            /* pixel number in x direction */
   long    ny,            /* pixel number in y direction */
   double  a)             /* new value */

/*
  set all entries of a 2-D array u to a given value a
*/

{
long  i,j;       /* loop variables */

for (i=1; i<=nx; i++)
  for (j=1; j<=ny; j++)
    u[i][j] = a;

return;

} /* set_matrix_values_double */

/*--------------------------------------------------------------------------*/

void set_matrix_values_long

  (long **u,           /* in+out :matrix (changed) */
   long nx,            /* pixel number in x direction */
   long ny,            /* pixel number in y direction */
   long a)             /* new value */

/*
  set all entries of a 2-D array u to a given value a
*/

{
long  i,j;       /* loop variables */

for (i=1; i<=nx; i++)
  for (j=1; j<=ny; j++)
      u[i][j] = a;

return;

} /* set_matrix_values_long */

/*--------------------------------------------------------------------------*/
void dummies_long

  (long **u,          /* image matrix */
   long   nx,         /* size in x direction */
   long   ny)         /* size in y direction */


/*
  creates dummy boundaries for image u in long format by mirroring
*/


{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
  {
  u[i][0]    = u[i][1];
  u[i][ny+1] = u[i][ny];
  }

for (j=0; j<=ny+1; j++)
  {
  u[0][j]    = u[1][j];
  u[nx+1][j] = u[nx][j];
  }

return;

}  /* dummies_long */

/* ------------------------------------------------------------------------- */
void dummies_double

  (double **u,        /* image matrix */
   long   nx,         /* size in x direction */
   long   ny)         /* size in y direction */

/*
  creates dummy boundaries for image u in double format by mirroring
*/

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
  {
  u[i][0]    = u[i][1];
  u[i][ny+1] = u[i][ny];
  }

for (j=0; j<=ny+1; j++)
  {
  u[0][j]    = u[1][j];
  u[nx+1][j] = u[nx][j];
  }

return;

}  /* dummies_double */

/* ------------------------------------------------------------------------- */

void set_bounds_long

  (long  **u,       /* image matrix */
   long   nx,       /* size in x direction */
   long   ny,       /* size in y direction */
   long   a)        /* boundary value */

/*
  set boundary entries of a 2-D array u to a given value a
*/

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
  {
  u[i][0]    = a;
  u[i][ny+1] = a;
  }

for (j=0; j<=ny+1; j++)
  {
  u[0][j]    = a;
  u[nx+1][j] = a;
  }

return;

}  /* set_bounds_long */

/* ------------------------------------------------------------------------*/

void set_bounds_double

  (double **u,          /* image matrix */
   long    nx,          /* size in x direction */
   long    ny,          /* size in y direction */
   double  a)           /* boundary value */

/*
 set boundary entries of a 2-D array u to a given value a
*/

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
  {
  u[i][0]    = a;
  u[i][ny+1] = a;
  }

for (j=0; j<=ny+1; j++)
  {
  u[0][j]    = a;
  u[nx+1][j] = a;
  }

return;

}  /* set_bounds_double */

/*--------------------------------------------------------------------------*/

void analyse_grey_double

  (double  **u,         /* image, unchanged */
   long    nx,          /* image size in x direction */
   long    ny,          /* image size in y direction */
   double  *min,        /* minimum, output */
   double  *max,        /* maximum, output */
   double  *mean,       /* mean, output */
   double  *std)        /* standard deviation, output */

/*
computes minimum, maximum, mean, and standard deviation of a greyscale
image u in double format
*/

{
long    i, j;       /* loop variables */
double  help1;      /* auxiliary variable */
double  help2;      /* auxiliary variable */

/* compute maximum, minimum, and mean */
*min  = u[1][1];
*max  = u[1][1];
help1 = 0.0;
for (i=1; i<=nx; i++)
  for (j=1; j<=ny; j++)
    {
    if (u[i][j] < *min) *min = u[i][j];
    if (u[i][j] > *max) *max = u[i][j];
    help1 = help1 + u[i][j];
    }
*mean = help1 / (nx * ny);

/* compute standard deviation */
*std = 0.0;
for (i=1; i<=nx; i++)
  for (j=1; j<=ny; j++)
    {
    help2  = u[i][j] - *mean;
    *std = *std + help2 * help2;
    }
*std = sqrt (*std / (nx * ny));

return;

}  /* analyse_grey_double */

/*------------------------------------------------------------------------*/

double MSE

     (long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      double   **u,       /* reconstructed image */
      double   **f)       /* original image */

/*
  mean squared error between the images u and f
*/

{
long     i, j, k;    /* loop variables */
double   aux;        /* time saver */
double   sum;        /* mean squared error */

sum = 0.0;
 for (i=1; i<=nx; i++)
  for (j=1; j<=ny; j++)
      {
      aux = u[i][j] - f[i][j];
      sum = sum + aux * aux;
      }
sum = sum / (double)(nx * ny);

return (sum);

} /* MSE */

/*--------------------------------------------------------------------------*/

long count_mask_points

  (long    nx,         /* image size in x direction */
   long    ny,         /* image size in y direction */
   long    **a)        /* binary inpainting mask, unchanged */

/* 
returns number of points in inpainting mask
*/

{
long    i, j;         /* loop variables */
long    counter;      /* counter for mask points */

/* count mask points */
counter = 0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
   if (a[i][j] == 1)
     counter++;

return (counter);

} /* count_mask_points */

/*--------------------------------------------------------------------------*/

void average_image

  (long   nx,       /* image size in x direction */
   long   ny,       /* image size in y direction */
   double **f,      /* input image */
   double **f_out)  /* output averaged image */

/*
The function averages pixels in a 3 by 3 area. The input is an initial
image with zero boundary. The output is an averaged image rounded
to the nearest integer  value.
*/

{
long           i, j, k, l; /* loop variables */
double         m;          /* counter */
double         help;       /* help variable */
double         aux;        /* auxiliary variable */
unsigned char  byte;       /* for data conversion */

for (i=1; i<=nx; i++)
  for (j=1; j<=ny; j++)
    {
    m = 0;
    help = 0;

    for (k=-1; k<=1; k++)
      for (l=-1; l<=1; l++)
        {
        m++;
        help += f[i+k][j+l];
        }

    f_out[i][j] = help / m;

    /* for correct rounding */
    aux = f_out[i][j] + 0.499999;
    if (aux < 0.0)
      byte = (unsigned char)(0.0);
    else if (aux > 255.0)
      byte = (unsigned char)(255.0);
    else
      byte = (unsigned char)(aux);
    f_out[i][j] = (double) byte;
    }

return;

} /* average image */

/*--------------------------------------------------------------------------*/

void reflect

  (long    nx,         /* image size in x direction */
   long    ny,         /* image size in y direction */
   long    i,          /* outer cycle iterator in x direction */
   long    j,          /* outer cycle iterator in y direction */
   long    m,          /* inner cycle iterator in x direction */
   long    n,          /* inner cycle iterator in y direction */
   long    d,          /* distance coefficient */
   long    *pi,        /* reflected x position */
   long    *pj)        /* reflected y position */

/*
computes position of reflected boundary
pixels by redirecting from mirrored boundary domain
to corresponding pixel in the image domain
Example: 1 2 3 | 3 2 1  - image and mirrored boundary.
*/

{
long hi, hj; /* help variables */

hi = i + m * d;
hj = j + n * d;

/* compute the new coordinates by redirecting */
*pi = hi <= 0 ? (1 - hi) : hi;
*pj = hj <= 0 ? (1 - hj) : hj;
*pi = hi >= (nx + 1) ? (1 - hi + 2 * nx) : *pi;
*pj = hj >= (ny + 1) ? (1 - hj + 2 * ny) : *pj;

return;

} /* compute_reflected_position */

/*------------------------------------------------------------------------*/

double err_compute_constraints

  (long   nx,       /* image size in x direction */
   long   ny,       /* image size in y direction */
   long   d,        /* distance to a patch centre */
   long   **bm,     /* binary constraints matrix */
   double **u,      /* inpainted image */
   double **f,      /* original image  */
   long   **a)      /* binary inpainting mask */

/*
computes the Constraint Satisfaction Rate (CSR)
*/

{
long   i, j, p, k, c; /* loop variables */
long   l;             /* counts number of mask points */
long   errn;          /* counts number of unsatisfied constraints */
double v;             /* help variable */
long   pi, pj;        /* mirrored positions */
double centre;        /* central pixel */
double neigh;         /* neighbouring pixel */
long   size;          /* length of the binary pattern */

/* to follow the clockwise direction */
long x_ids[8] = {-1,  0,  1,  1,  1,  0, -1, -1};
long y_ids[8] = {-1, -1, -1,  0,  1,  1,  1,  0};

/* initialise */
size = 8;
pi   = 0;
pj   = 0;
l    = 0;
errn = 0;

/* iterate over the reconstructed image */
for (i=1; i<=nx; i++)
  for(j=1; j<=ny; j++)
    if (a[i][j] != 0)
      {
      centre = u[i][j];

      /* visit neighbours of a mask pixel */
      for(c=0; c<size; c++)
        {
        p = x_ids[c];
        k = y_ids[c];

        /* dummy boundaries */
        reflect (nx, ny, i, j, p, k, d, &pi, &pj);

        /* find the neighbour */
        neigh = u[pi][pj];

        /* compare mask pixel value with the neighbour */
        if (neigh > centre)
          v = 1;
        else
          v = 0;

        /* compute the number of unsatisfied constraints */
        if (v - bm[l][c] != 0)
          errn++;
        }
      l++;
      }

return  (1 - errn / (double)(size * l));

} /* compute_constraints */

/*--------------------------------------------------------------------------*/

void encode_image

  (long    nx,        /* image size in x direction */
   long    ny,        /* image size in y direction */
   long    **a,       /* binary inpainting mask */
   double  **f,       /* initial image */
   long    levels,    /* number of levels */
   char    *fname)    /* file name */

/*
saves the image values and local binary patterns into the file
*/

{
long    size;           /* length of the binary pattern */
long    level;          /* loop variable */
long    n;              /* number of mask points */
double  ***F;           /* stores averaged images */
long    d;              /* distance to a patch centre */
long    i, j, k, p, c;  /* loop variables */
long    pi, pj;         /* mirrored positions */
long    v;              /* help variable */
double  centre;         /* central pixel */
double  neigh;          /* neighbouring pixel */
long    count;          /* counter */
FILE    *file;          /* file object */

/* to follow the clockwise direction */
long x_ids[8] = {-1,  0,  1,  1,  1,  0, -1, -1};
long y_ids[8] = {-1, -1, -1,  0,  1,  1,  1,  0};

/* open the file*/
file = fopen (fname, "w");

/* initialise */
size = 8;

/* allocate memory */
alloc_double_cubix (&F, levels, nx+2, ny+2);

/* count mask points */
n=count_mask_points (nx, ny, a);

/* set first array to original image */
for (i=1; i<=nx; i++)
  for(j=1; j<=ny; j++)
    F[0][i][j] = f[i][j];

/* mirror boundaries */
dummies_double (F[0], nx, ny);

/* write files */
for (level=0;level<levels;level++)
  {
  /* average images starting from the second level */
  if (level > 0)
    {
    average_image (nx, ny, F[level-1], F[level]);
    dummies_double (F[level], nx, ny);
    }

  /* distance to the neighbours */
  d = pow (3, level);

  /* initialisations */
  pi = 0;
  pj = 0;
  count=0;

  /* iterate */
  for (i=1; i<=nx; i++)
    for(j=1; j<=ny; j++)
      /* if we have a mask pixel */
      if (a[i][j] != 0)
        {
        /* save mask value */
        fprintf (file, "%u\n", (unsigned int)(F[level][i][j]));
        centre = F[level][i][j];

        /* visit centres of neighbouring regions */
        for(c=0; c<size; c++)
          {
          p = x_ids[c];
          k = y_ids[c];

          /* dummy boundaries */
          reflect (nx, ny, i, j, p, k, d, &pi, &pj);
          neigh = F[level][pi][pj];

          /* compare centres with its neighbours */
          if (neigh > centre)
            v=1;
          else
            v=0;

          /* save the result */
          fprintf (file, "%ld",v);
          }
        fprintf (file, " \n");
        count++;
        }
  }

/* close the file */
fclose (file);

/* free memory */
free_double_cubix (F, levels, nx+2, ny+2);

return;

} /* encode_image */

/*--------------------------------------------------------------------------*/

void read_file

    (long    nx,         /* image size in x direction */
     long    ny,         /* image size in y direction */
     long    levels,     /* number of levels */
     long    **a,        /* binary inpainting mask */
     long    ***bm,      /* binary constraints matrix */
     double  ***U,       /* mask image */
     char    *fname)     /* file name */

/*
read the file with information about constraints and known image values
*/

{
FILE *file;     /* file object */
long i, j, m;   /* iterators */
long val;       /* value */
long level;     /* loop variable */
long k;         /* counter */
long size;      /* length of the binary pattern */

size = 8;
char buff[size]; /* to read binary string */

/* open the file */
file = fopen (fname, "r");

/* for each level */
for (level=0; level<levels; level++)
  {
  k=0;

  /* loop over the image mask */
  for(i=1; i<=nx; i++)
    for(j=1; j<=ny; j++)
      if (a[i][j] != 0)
        {
        /* read initial image value from the file */
        fscanf (file, "%ld", &val);
        U[level][i][j] = (double) val;

        /* read pattern */
        fscanf (file, "%s", buff);
        for (m=0; m<size; m++)

          /* subtract 48 to convert char into long */
          bm[level][k][m] = (long)buff[m] - 48;

        k++;
        }
      else
        U[level][i][j] = 0.;
  }

/* close the file */
fclose (file);

return;

} /* read file */

/*--------------------------------------------------------------------------*/

/* data structure for saving a single matrix element */
typedef struct
{
long   row; /* row */
long   col; /* column */
double val; /* stored value */
} Elem;

/*--------------------------------------------------------------------------*/

/* data structure for saving a single sparse matrix */
typedef struct
{
long size;  /* number of non-zero elements */
long col;   /* vertical dimension */
long row;   /* horizontal dimension */
Elem *arr;  /* array of sparse matrix elements */
} Sparse;

/*--------------------------------------------------------------------------*/

/* data structure for saving the matrix of sparse matrices */
typedef struct
{
long   ncols;  /* number of columns */
long   nrows;  /* number of rows */
long   size;   /* number of matrices */
Sparse *arr;   /* sparse matrix element */
} SparseSparse;

/*--------------------------------------------------------------------------*/

double *sum_col

  (double *v,  /* input vector */
   Sparse s)   /* input sparse matrix */

/*
computes column sums of the input sparse matrix
*/

{
long i; /* iterator */

for (i=0; i<s.size; i++)
  v[s.arr[i].col] += fabs (s.arr[i].val);

return (v);

} /* sum_col */

/*--------------------------------------------------------------------------*/

double *sum_row

  (double *v,  /* input vector */
   Sparse s)   /* input sparse matrix */

/*
computes row sums of the input sparse matrix
*/

{
long i; /* iterator */

for (i=0; i<s.size; i++)
  v[s.arr[i].row] += fabs (s.arr[i].val);

return (v);

} /* sum_row */

/*--------------------------------------------------------------------------*/

double *sum_col_vec

  (double *v,     /* output vector */
   double *vec,   /* input vector */
   Sparse s)      /* input sparse matrix */

/*
computes transpose vector multiplication
*/

{
long i; /* iterator */

for (i=0; i<s.size; i++)
  v[s.arr[i].col] += s.arr[i].val * vec[s.arr[i].row];

return (v);

} /* sum_col_vec*/

/*--------------------------------------------------------------------------*/

double *sum_row_vec

  (double *v,     /* output vector */
   double *vec,   /* input vector */
   Sparse s)      /* input sparse matrix */

/*
computes matrix vector multiplication
*/

{
long i; /* iterator */

for (i=0; i<s.size; i++)
  v[s.arr[i].row] += s.arr[i].val * vec[s.arr[i].col];

return (v);

} /* sum_row_vec */

/*--------------------------------------------------------------------------*/

void create_constraint_matrix

   (long   nx,       /* image size in x direction */
    long   ny,       /* image size in y direction */
    long   level,    /* level of constraint */
    long   **a,      /* binary inpainting mask */
    double *b,       /* right part of the system Kx=b */
    long   **bm,     /* binary constraints matrix */
    double **u,      /* evolving image */
    Elem   *K)       /* storage for the constraint */

/*
fills in the constraint matrix A^(i)
*/

{
long   i, j, k, l, p, r, h ,c;  /* loop variables */
long   pi, pj;                  /* mirrored positions */
long   d;                       /* distance to a patch centre */
long   size;                    /* length of the binary pattern */
double epsilon;                 /* small constant for numerical stability */

/* to follow the clockwise direction */
long x_ids[8] = {-1,  0,  1,  1,  1,  0, -1, -1};
long y_ids[8] = {-1, -1, -1,  0,  1,  1,  1,  0};

/* initialise */
size    = 8;
l       = 0;
r       = 0;
d       = pow (3, level);
epsilon = 1e-6;
h       = 0;

/* loop over the image */
for(i=1; i<=nx; i++)
  for(j=1; j<=ny; j++)
    /* if we have a mask pixel */
    if (a[i][j] == 1)
      {
      /* visit centres of neighbouring regions */
      for(c=0; c<size; c++)
        {
        p = x_ids[c];
        k = y_ids[c];

        /* dummy boundaries */
        reflect (nx, ny, i, j, p, k, d, &pi, &pj);

        /* in case LBP = 1 */
        if (bm[h][c] == 1)
          {
          /* neigh >= centre */
          /* -neigh <= -centre */
          b[r] = -u[i][j];
          K[r].col = (pi-1)*ny + (pj-1);

          /* -neigh */
          K[r].val = -1.;
          K[r].row = r;
          r++;
          }

        /* in case LBP = 0 */
        if (bm[h][c] == 0)
          {
          /* neigh < centre */
          /* neigh <= centre - epsilon */
          b[r] = u[i][j] - epsilon;
          K[r].col = (pi-1)*ny + (pj-1);

          /*  neigh */
          K[r].val = 1.;
          K[r].row = r;
          r++;
          }

        }

        h++;
      }

return;

} /* create_constraint_matrix */

/*--------------------------------------------------------------------------*/

void create_averaging_constraint_matrix

   (long nx,       /* image size in x direction */
    long ny,       /* image size in y direction */
    long level,    /* the model level */
    long **a,      /* binary inpainting mask */
    long **bm,     /* binary constraints matrix */
    Elem *K)       /* averaging constraint matrix */

/*
fills in the averaging constraint matrix -B^(i)
*/

{
long i, j, k, l, m, n, p, r, h, c;  /* loop variables */
long pi, pj;                        /* mirrored positions */
long mi, mj;                        /* mirrored positions */
long d;                             /* distance to a patch centre */
long da;                            /* distance to a patch boundary */
long size;                          /* length of the binary pattern */

/* to follow the clockwise direction */
long x_ids[8] = {-1,  0,  1,  1,  1,  0, -1, -1};
long y_ids[8] = {-1, -1, -1,  0,  1,  1,  1,  0};

/* initialise */
size = 8;
l    = 0;
r    = 0;
d    = pow (3, level);
da   = 1;
h    = 0;

/* loop over the image */
for(i=1; i<=nx; i++)
  for(j=1; j<=ny; j++)
    /* if we have a mask pixel */
    if (a[i][j] == 1)
      {
      /* visit centres of neighbouring regions */
      for(c=0; c<size; c++)
        {
        p = x_ids[c];
        k = y_ids[c];

        /* dummy boundaries */
        reflect (nx, ny, i, j, p, k, d, &pi, &pj);

        /* loop around the centre of the neighbouring region */
        for (m=-da; m<=da; m++)
          for (n=-da; n<=da; n++)
            {
            /* reflect boundary */
            reflect (nx, ny, pi, pj, m, n, 1, &mi, &mj);

            /* in case LBP = 1 */
            if (bm[h][c] == 1)
              {
              K[r].col = (mi-1)*ny + (mj-1);
              K[r].val = 1. / 9.;
              K[r].row = l;
              r++;
              }

            /* in case LBP = 0 */
            if (bm[h][c] == 0)
              {
              K[r].col = (mi-1)*ny + (mj-1);
              K[r].val = -1. / 9.;
              K[r].row = l;
              r++;
              }
            }
          l++;
        }
      h++;
      }

return;

} /* create_averaging_constraint_matrix */

/*--------------------------------------------------------------------------*/

void inpaint_image_full

   (long   nx,     /* image size in x direction */
    long   ny,     /* image size in y direction */
    long   levels, /* number of levels */
    double eps,    /* stopping criterion */
    long   hx,     /* step size in x direction */
    long   hy,     /* step size in y direction */
    long   **a,    /* binary inpainting mask */
    double **f,    /* initial image */
    double ***U,   /* array to store evolving  images */
    long   ***bm)  /* binary constraints matrix */

/*
implements the algorithm for inpainting with LBC
*/

{
double  xp, xm, yp, ym; /* neighbourhood weights */
double  hx_2, hy_2;     /* time saver variables */
long    size;           /* length of the binary pattern */
long    i, j, k, l;     /* loop variables */
long    z;              /* current iteration number */
long    level;          /* loop variable */
long    Np;             /* number of mask points */
long    ncols;          /* number of columns in a matrix A or B*/
long    nrows;          /* number of rows in a matrix A or B*/
long    cols;           /* number of columns in a sparse matrix of matrices */
long    rows;           /* number of rows in a sparse matrix of matrices */
long    d;              /* distance to a patch centre */
long    a_elem;         /* number of elements in a sparse matrix A */
long    b_elem;         /* number of elements in a sparse matrix B */
long    num_mat_K1;     /* number of matrices in matrix K1 */
long    num_mat_K2;     /* number of matrices in matrix K2 */
double  *gradEu;        /* to store the Laplacian value */
double  *uold;          /* to store u on the step k */
double  *err;           /* to save CSR for each level */
double  perr;           /* dual variable error */
double  uerr;           /* primal variable error */
double  pprev;          /* previous error of the dual variable */
double  mse;            /* mean squared error */
double  ***F;           /* stores averaged images */
double  **Tau;          /* time step tau column sums */
double  **Sigma1;       /* time step sigma1 for K1 matrix*/
double  **Sigma2;       /* time step sigma2 for K2 matrix */
double  **P1;           /* to store p1 */
double  **P2;           /* to store p2 */
double  **K1_T_P1;      /* to store K^Tp1 */
double  **K2_T_P2;      /* to store K^Tp2 */
double  **U_bar;        /* to store u^bar */
double  **K1_U_bar;     /* to store K1*u^bar */
double  **K2_U_bar;     /* to store K2*u^bar */
double  **b;            /* right-hand side */
Sparse  *K1;            /* to store K1 matrices */
Sparse  *K2;            /* to store K2 matrices */
Elem    **A;            /* to store constraint elements */
Elem    **B;            /* to store averaging elements  */

/* initialise */

/* count number of mask points */
Np = count_mask_points (nx, ny, a);

/* length of the binary pattern */
size = 8;

/* number of non-zero elements in sparse matrices A and B */
a_elem = Np * size;
b_elem = Np * size * 9;

/* dimension of matrices A and B */
ncols = nx * ny;
nrows = size * Np;

/* number of sparse matrices in matrices K1 and K2 */
num_mat_K1 = levels;
num_mat_K2 = 2 * levels - 2;

/* dimension of a sparse matrix of sparse matrices K */
cols = levels;
rows = 2 * levels - 1;


/* ---- allocate memory ---- */

alloc_double_cubix  (&F, levels, nx+2, ny+2);
alloc_double_matrix (&P1, levels, nrows);
alloc_double_matrix (&P2, levels-1, nrows);
alloc_double_matrix (&U_bar, levels, nx*ny);
alloc_double_matrix (&K1_T_P1, cols, ncols);
alloc_double_matrix (&K2_T_P2, cols, ncols);
alloc_double_matrix (&K1_U_bar, levels, nrows);
alloc_double_matrix (&K2_U_bar, levels-1, nrows);
alloc_double_vector (&err, levels);
alloc_double_vector (&uold, levels);
alloc_double_vector (&gradEu, levels);
alloc_double_matrix (&Tau, cols, ncols);
alloc_double_matrix (&Sigma1, levels, nrows);
alloc_double_matrix (&Sigma2, levels-1, nrows);
alloc_double_matrix (&b, levels, nrows);

/* to store constraint matrices */
K1 = malloc (num_mat_K1 * sizeof(Sparse));
K2 = malloc (num_mat_K2 * sizeof(Sparse));

/* to combine K1 and K2 into one matrix K */
SparseSparse *K = malloc (2 * sizeof(Sparse));

/* to store inequality constraints */
A =  malloc (levels * sizeof(Elem *));
for (level=0; level<levels; level++)
  A[level] = malloc (a_elem * sizeof(Elem));

/* to store averaging constraints */
B = malloc ((levels-1) * sizeof(Elem *));
for (level=0; level<levels-1; level++)
  B[level] = malloc (b_elem * sizeof(Elem));


/* ---- initialisations ---- */

/*set the first array to the original image. We use it to compute the
   constraint satisfaction rate (CSR) */
for (i=1; i<=nx; i++)
  for(j=1; j<=ny; j++)
    F[0][i][j] = f[i][j];

/* average the initial image if needed*/
dummies_double (F[0], nx, ny);
for (level=0; level<levels; level++)
  {
  if (level>0)
    {
    average_image (nx, ny, F[level-1], F[level]);
    dummies_double (F[level], nx, ny);
    }
  }

/* mirror boundary values */
for (level=0; level<levels; level++)
  dummies_double (U[level], nx, ny);

/* fill in the  constraint matrix K1 */

/*            Matrix K1     */
/*     |A[0]  0  ...   0  | */
/* K1 =| 0  A[1] ...   0  | */
/*     |... ...  ...  ... | */
/*     | 0   0  ...  A[n] | */

for (level=0; level<levels; level++)
  {
  create_constraint_matrix (nx, ny, level, a,
      b[level], bm[level], U[level], A[level]);

  K1[level].arr  = A[level];
  K1[level].size = a_elem;
  K1[level].col  = level;
  K1[level].row  = level;
  }

/* fill in the constraint matrix K2 */

/*                  Matrix K2            */
/*      |-B[0]  A[1]    0     ...   0  | */
/* K2 = | 0    -B[1]  A[2]   ...    0  | */
/*      |...     ...   ...    ...  ... | */
/*      | 0   0  ...    0  -B[n-1] A[n]| */

/* fill in the averaging matrix */
l=0;
for (level=0; level<levels-1; level++)
  {
  create_averaging_constraint_matrix (nx, ny, level+1, a,
      bm[level+1], B[level]);

  K2[l].arr  = B[level];
  K2[l].size = b_elem;
  K2[l].col  = level;
  K2[l].row  = level;

  l++;

  K2[l].arr  = A[level+1];
  K2[l].size = a_elem;
  K2[l].col  = level+1;
  K2[l].row  = level;

  l++;
  }

/* create sparse matrix of sparse matrices */
K[0].arr   = K1;
K[0].size  = levels;
K[0].ncols = levels;
K[0].nrows = levels;

K[1].arr   = K2;
K[1].size  = 2*levels-2;
K[1].ncols = levels;
K[1].nrows = levels-1;

/* time step tau row sums */

/* initialise */
for (i=0; i<cols; i++)
  for (j=0; j<ncols; j++)
    Tau[i][j] = 0.;

/* compute row sums */
for (j=0; j<2; j++)
  for (i=0; i<K[j].size; i++)
    Tau[K[j].arr[i].col] = sum_col (Tau[K[j].arr[i].col], K[j].arr[i]);

/* get tau values */
for (i=0; i<cols; i++)
  for (j=0; j<ncols; j++)
    Tau[i][j]  =  1.0f / (4.0f+Tau[i][j]);

/* time step sigma column sums */

/* initialise */
for (i=0; i<levels; i++)
  for (j=0; j<nrows; j++)
    Sigma1[i][j] = 0;

for (i=0; i<levels-1; i++)
  for (j=0; j<nrows; j++)
    Sigma2[i][j] = 0;

/* compute row sums for matrix K1 */
for (i=0; i<num_mat_K1; i++)
  Sigma1[K1[i].row] = sum_row (Sigma1[K1[i].row], K1[i]);

/* compute row sums for matrix K2 */
for (i=0; i<num_mat_K2; i++)
  Sigma2[K2[i].row] = sum_row (Sigma2[K2[i].row], K2[i]);

/* get Sigma1 values */
for (i=0; i<levels; i++)
  for (j=0; j<nrows; j++)
    Sigma1[i][j]=1. / Sigma1[i][j];

/* get Sigma2 values */
for (i=0; i<levels-1; i++)
  for (j=0; j<nrows; j++)
    Sigma2[i][j]=1. / Sigma2[i][j];

/* initialise dual variables */
for (i=0; i<levels; i++)
  for (j=0; j<nrows; j++)
    P1[i][j] = 0.;

for (i=0; i<levels-1; i++)
  for (j=0; j<nrows; j++)
    P2[i][j] = 0.;

/* initialise u^bar = u^0 */
for(l=0; l<levels; l++)
  {
  k=0;
  for (i=1; i<=nx; i++)
    for (j=1; j<=ny; j++)
      {
      U_bar[l][k] = U[l][i][j];
      k++;
      }
  }

hx_2  = 1.0/(hx*hx);
hy_2  = 1.0/(hy*hy);
uerr  = 1;
perr  = 1;
pprev = 0;


/* ---- iterations ---- */

k=0;

/* while stopping criterion is not reached do iterations */
while (uerr>eps ||  perr>eps)
  {
  uerr = 0;
  perr = 0;

  /* initialise matrix vector products */
  for (i=0; i<levels; i++)
    for(j=0; j<nrows; j++)
      K1_U_bar[i][j] = 0;

  for (i=0; i<levels-1; i++)
    for(j=0; j<nrows; j++)
      K2_U_bar[i][j] = 0;

  /* compute K1*u^bar product */
  for (i=0; i<num_mat_K1; i++)
    K1_U_bar[K1[i].row] = sum_row_vec (K1_U_bar[K1[i].row],
        U_bar[K1[i].col], K1[i]);

  /* compute K2*u^bar product */
  for (i=0; i<num_mat_K2; i++)
    K2_U_bar[K2[i].row] = sum_row_vec (K2_U_bar[K2[i].row],
        U_bar[K2[i].col], K2[i]);

  /* update the dual variables p1 and p2 */

  /* in case of inequality constraints update p1 */
  /* Add 1e-6 for stability reasons */
  for (i=0; i<levels; i++)
    for (j=0; j<nrows; j++)
      {
      pprev = P1[i][j];
      P1[i][j] = fmax (0., pprev +
          Sigma1[i][j] * (K1_U_bar[i][j] - b[i][j]) + 1e-6);
      perr += (P1[i][j] - pprev) * (P1[i][j] - pprev);
      }

  /* in case of equality constraints update p2 */
  for (i=0; i<levels-1; i++)
    for (j=0; j<nrows; j++)
      {
      pprev = P2[i][j];
      P2[i][j] = pprev + Sigma2[i][j] * (K2_U_bar[i][j] + 1e-6);
      perr += (P2[i][j] - pprev) * (P2[i][j] - pprev);
      }

  /* initialise vector products with zeros */
  for (i=0; i<cols; i++)
    for (j=0; j<ncols; j++)
      {
      K1_T_P1[i][j] = 0.;
      K2_T_P2[i][j] = 0.;
      }

  /* compute K1^T*p1 product */
  for (i=0; i<num_mat_K1; i++)
    K1_T_P1[K1[i].col] = sum_col_vec (K1_T_P1[K1[i].col],
        P1[K1[i].row], K1[i]);

  /* compute K2^T*p2 product */
  for (i=0; i<num_mat_K2; i++)
    K2_T_P2[K2[i].col] = sum_col_vec (K2_T_P2[K2[i].col],
        P2[K2[i].row], K2[i]);

  /* update the primal variable u */
  z=0;
  for (i=1; i<=nx; i++)
    for (j=1; j<=ny; j++)
      {
      if (a[i][j] == 1)
        {
        for(l=0; l<levels; l++)
          U_bar[l][z] = U[l][i][j];
        }
      else
        {
        /* compute weights */
        xp =  (i<nx) * hx_2;
        xm =  (i>1)  * hx_2;
        yp =  (j<ny) * hy_2;
        ym =  (j>1)  * hy_2;

        for(l=0; l<levels; l++)
          {
          uold[l] = U[l][i][j];

          /* compute the Laplacian */
          gradEu[l] = (xm * U[l][i-1][j  ]+
                       ym * U[l][i  ][j-1]+
                       yp * U[l][i  ][j+1]+
                       xp * U[l][i+1][j  ]-
                       xm * U[l][i  ][j  ]-
                       xp * U[l][i  ][j  ]-
                       ym * U[l][i  ][j  ]-
                       yp * U[l][i  ][j  ]);

          U[l][i][j] = U[l][i][j] - Tau[l][z] *
          (-gradEu[l] + K1_T_P1[l][z] + K2_T_P2[l][z] + 1e-6) + 1e-6;

          /* update u^bar variable */
          U_bar[l][z]  = 2. * U[l][i][j] - uold[l];
          }
        uerr += (U[0][i][j] - uold[0]) * (U[0][i][j] - uold[0]);
        }
      z++;
      }

  /* compute MSE */
  mse = MSE (nx, ny, f, U[0]);

  /* compute CSR for each level */
  for (level=0; level<levels; level++)
    {
    d = pow (3, level);
    err[level] = err_compute_constraints (nx, ny,
        d, bm[level], U[level], F[level], a) ;
    }

  /* print CSR for each level */
  if (levels == 1)
    {
    printf ("u_err = %f, p_err = %f, lev_1_CSR = %f, ",
        uerr, perr, err[0]);
    printf ("MSE = %f, iter = %ld \n", mse, k);
    }
  else if (levels == 2)
    {
    printf ("u_err = %f, p_err = %f, lev_1_CSR = %f, ",
        uerr, perr, err[0]);
    printf ("lev_2_CSR = %f, ", err[1]);
    printf ("MSE = %f, iter = %ld \n", mse, k);
    }

  else
    {
    printf ("u_err = %f, p_err = %f, lev_1_CSR = %f, ",
        uerr, perr, err[0]);
    printf ("lev_2_CSR = %f, ", err[1]);
    printf ("lev_3_CSR = %f, ", err[2]);
    printf ("MSE = %f, iter = %ld \n", mse, k);
    }

  k++;
  }

/* clip to [0, 255]  in case outliers */
for (i=1; i<=nx; i++)
  for (j=1; j<=ny; j++)
	  {
	  if (U[0][i][j] > 255)
		  U[0][i][j] = 255;

	  if (U[0][i][j] < 0)
		  U[0][i][j] = 0;
	  }


/* ---- free memory ---- */

free(K1);
free(K2);
free(K);
for(level=0;level<levels;level++)
  free(A[level]);
free(A);
for(level=0;level<levels-1;level++)
  free(B[level]);
free(B);
free_double_matrix(b, levels, nrows);
free_double_matrix(Tau, cols,ncols);
free_double_matrix(Sigma1, levels, nrows);
free_double_matrix(Sigma2, levels-1,  nrows);
free_double_matrix(K1_T_P1, cols, ncols);
free_double_matrix(K2_T_P2, cols, ncols);
free_double_vector(err,levels);
free_double_vector(gradEu,levels);
free_double_vector(uold,levels);
free_double_matrix(P1, levels, nrows);
free_double_matrix(P2, levels-1, nrows);
free_double_matrix(U_bar, cols, ncols);
free_double_matrix(K1_U_bar, levels, nrows);
free_double_matrix(K2_U_bar, levels-1, nrows);
free_double_cubix(F, levels, nx+2, ny+2);

return;

} /* inpaint_image_full */

/*--------------------------------------------------------------------------*/
int main ()
{
char    in1[80];          /* filename of input image */
char    in2[80];          /* filename of input mask */
char    out1[80];         /* filename of output image */
char    fname[80];        /* filename for saving LBC */
double  **u;              /* evolving image */
double  **f;              /* original image */
long    **a;              /* inpainting mask */
long    i, j;             /* loop variables */
long    nx, ny;           /* image size in x, y direction */
double  hx = 1.0;         /* step size in x direction */
double  hy = 1.0;         /* step size in y direction */
double  max, min;         /* largest, smallest grey value */
double  mean;             /* average grey value */
double  std;              /* standard deviation */
long    n_mask_points;    /* number of mask points */
double  eps;              /* stopping criterion */
double  mse;              /* mean squared error */
double  ***U;             /* to store evolving  images */
long    ***bm;            /* to store binary constraints */
long    levels;           /* number of inequality constraints */
long    size;             /* length of the binary pattern */


printf ("\n");
printf ("INPAINTING WITH LOCAL BINARY CONSTRAINTS\n");
printf ("***************************************************\n\n");
printf ("    Copyright 2019 by Andrei Sirazitdinov          \n");
printf ("    Faculty of Mathematics and Computer Science    \n");
printf ("    Saarland University, Germany                   \n\n");
printf ("    All rights reserved. Unauthorized usage,       \n");
printf ("    copying, hiring, and selling prohibited.       \n\n");
printf ("    Send bug reports to                            \n");
printf ("    andreisirazitdinov@gmail.com                   \n\n");
printf ("***************************************************\n\n");


/* ---- read input image (pgm format p5) ---- */

printf ("input image (pgm):                             ");
read_string (in1);
read_pgm_to_double (in1, &nx, &ny, &f);

/* dummy boundary values */
dummies_double (f, nx, ny);


/* ---- read binary inpainting mask (pgm format p5) ---- */

printf ("inpainting mask (pgm):                         ");
read_string (in2);
read_pgm_to_long (in2, &nx, &ny, &a);

/* binarise inpainting mask */
n_mask_points = 0;
for (i=1; i<=nx; i++)
  for (j=1; j<=ny; j++)
    if (a[i][j] < 128)
		a[i][j] = 0;
    else
		a[i][j] = 1;


/* ---- read other parameters ---- */

printf ("number of levels (1, 2, or 3):                 ");
read_long (&levels);

printf ("stopping criterion for PDHG in ]0,1[:          ");
read_double (&eps);


/* ---- output image ---- */

printf ("output image (pgm):                            ");
read_string (out1);


/* ---- initialisations ---- */

/*size of a single constraint */
size = 8;

/* unit grid size*/
hx = hy = 1.;

/* count number of mask points */
n_mask_points = count_mask_points (nx, ny, a);


/* ---- allocate memory ---- */

alloc_double_matrix (&u, nx+2, ny+2);
alloc_long_cubix (&bm, levels, n_mask_points, size);
alloc_double_cubix (&U, levels, nx+2, ny+2);


/* ---- inpainting ---- */

/* save binary patterns in the file using the number of layers as
   a file name */
sprintf (fname, "%ld.txt", levels);
encode_image (nx, ny, a, f, levels, fname);

/* read the file containing LBC and  mask values */
read_file (nx, ny, levels, a, bm, U, fname);

/* inpainting with local binary constraints */
inpaint_image_full (nx, ny, levels, eps, hx, hy, a, f, U, bm);

/* extract the output image */
for (i=1; i<=nx; i++)
  for (j=1; j<=ny; j++)
	  u[i][j] = U[0][i][j];

/* check results */
mse = MSE(nx, ny, f, u);
analyse_grey_double (u, nx, ny, &min, &max, &mean, &std);
printf ("inpainted image\n");
printf ("minimum:         %10.2lf \n", min);
printf ("maximum:         %10.2lf \n", max);
printf ("mean:            %10.2lf \n", mean);
printf ("std. dev.:       %10.2lf \n\n", std);
printf ("mse:             %10.2lf \n\n", mse);


/* ---- write output image (pgm format P5) ---- */

write_double_to_pgm (u, nx, ny, out1, 0);
printf ("output image %s successfully written\n\n", out1);


/* ---- free memory ---- */

free_double_matrix (u, nx+2, ny+2);
free_double_matrix (f, nx+2, ny+2);
free_long_matrix (a, nx+2, ny+2);
free_long_cubix (bm, levels, n_mask_points, size);
free_double_cubix (U, levels, nx+2, ny+2);

return (0);

} /* main */

