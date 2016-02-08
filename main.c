/*
 * main.c
 *
 *  Created on: Jan 29, 2016
 *      Author: vova
 */

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include "ccarray.h"

#define UNUSED(x)               ((void)(x))
#define MAX_INPUT_LINE_LENGTH   4096
#define PI                      M_PI

/** known compression types */
typedef
enum compression_t {
    compression_unknown = -1,
    compression_none,
    compression_bzip,
    compression_gzip
} compression_t;



typedef
struct stars{
		double x;
		double f1;
		double f2;
		int count;
	}stars;


static FILE * open_file( const char * fname, compression_t * compression )
{
  char cmd[PATH_MAX] = {0};
  FILE * input = NULL;

  if ( *compression == compression_unknown )
  {
	const char * suffix;

	if ( (suffix = strstr(fname, ".bz2")) && *(suffix + 4) == 0 ) {
	  *compression = compression_bzip;
	}
	else if ( (suffix = strstr(fname, ".bz")) && *(suffix + 3) == 0 ) {
	  *compression = compression_bzip;
	}
	else if ( (suffix = strstr(fname, ".gz")) && *(suffix + 3) == 0 ) {
	  *compression = compression_gzip;
	}
	else {
	  *compression = compression_none;
	}
  }

  switch ( *compression )
  {
  case compression_bzip:
	snprintf(cmd, sizeof(cmd) - 1, "bzip2 -dc '%s'", fname);
	if ( !(input = popen(cmd, "r"))) {
	  fprintf(stderr, "popen('%s') fails: %s\n", cmd, strerror(errno));
	}
	break;

  case compression_gzip:
	snprintf(cmd, sizeof(cmd) - 1, "gzip -dc '%s'", fname);
	if ( !(input = popen(cmd, "r"))) {
	  fprintf(stderr, "popen('%s') fails: %s\n", cmd, strerror(errno));
	}
	break;

  case compression_none:
	if ( !(input = fopen(fname, "r")) ) {
	  fprintf(stderr, "fopen('%s') fails: %s\n", fname, strerror(errno));
	}
	break;

  default:
	fprintf(stderr, "BUG IN CODE: invalid compression tag=%d\n", *compression);
	break;
  }

  return input;
}
static int partition (double * m, int a, int b)
	{
	  int i = a;
	  int j;
	  for (j = a; j <= b; j++ )    // просматриваем с a по b
	   {
		 if (m[j] <= m[b])            // если элемент m[j] не превосходит m[b],
		  {
			//swap(m[i],m[j]);         // меняем местами m[j] и m[a], m[a+1], m[a+2] и так далее...
			double temp=m[i];						  // то есть переносим элементы меньшие m[b] в начало,
			m[i]=m[j];						  // а затем и сам m[b] «сверху»
			m[j]=temp;
			i++;                      // таким образом последний обмен: m[b] и m[i], после чего i++
		  }
	   }
	  return i-1;                     // в индексе i хранится <новая позиция элемента m[b]> + 1
	}

static void quicksort (double * m, int a, int b) // a - начало подмножества, b - конец
	{                                        // для первого вызова: a = 0, b = <элементов в массиве> - 1
	 if (a >= b) return;
	 int c = partition (m, a, b);
	 //printf("a=%d\tb=%d\tc=%d\n",a,b,c);

	 quicksort (m, a, c-1);
	 quicksort (m, c+1, b);

	}

double median_arry (double * m, int n ) // нахождение среднего значения массива размером n
{
	 if (n < 1) {
		 return 0;
	 }
	 quicksort(m,0,n-1);
	 if(n%2==0)
	   {
		int nn = n/2;
		return (m[nn]+m[nn-1])/2;
	   }
		else
		{
		  int nn = (int)((n-1)/2);
		  return m[nn];
		}
}

stars smedian_arry (stars * m, int n ) // нахождение среднего значения массива размером n
{
	stars median;
	median.x  = 0;
	median.f1 = 0;
	median.f2 = 0;
	 if (n < 1) {
		 return median;
	 }
	 double x[n];
	 double f1[n];
	 double f2[n];
	 int i;
	 for(i=0;i<n;++i)
	 {
		 x[i] = m[i].x;
		 f1[i]= m[i].f1;
		 f2[i]= m[i].f2;
	 }

	 median.x  = median_arry(x,n-1);
	 median.f1 = median_arry(f1,n-1);
	 median.f2 = median_arry(f2,n-1);
	 median.count = n;

	 return median;
}
stars sigma_arry(stars * m, stars msig, stars mean, int n ) // нахождение сигмы массива размером n
{

	if (n < 1) {
			 return msig;
		 }
	stars sig;
	sig.x	= 0;
	sig.f1	= 0;
	sig.f2	= 0;

	int count=0;
	 int i;
	 for( i=0; i< n; ++i)
	 {
		 if(fabs(mean.f1-m[i].f1)<3*msig.f1 && fabs(mean.f2-m[i].f2)<3*msig.f2)
		 {
			 sig.x	+= (m[i].x-mean.x)*(m[i].x-mean.x);
			 sig.f1	+= (m[i].f1-mean.f1)*(m[i].f1-mean.f1);
			 sig.f2	+= (m[i].f2-mean.f2)*(m[i].f2-mean.f2);
			 ++count;
		 }
	 }
	 sig.x	= sqrt(sig.x/count);
	 sig.f1	= sqrt(sig.f1/count);
	 sig.f2	= sqrt(sig.f2/count);
return sig;
}

stars mean_arry(stars * m, int n ) // нахождение среднего значения массива размером n
{
	stars mean;
	mean.x=0;
	mean.f1=0;
	mean.f2=0;

	 if (n < 1) {
		 return mean;
	 }
	 int ii=0;
	 int i;
	 for( i=0; i< n; ++i)
	 	 {
		 ++ii;
		 mean.x  += (m[i].x-mean.x)/ii;
		 mean.f1 += (m[i].f1-mean.f1)/ii;
		 mean.f2 += (m[i].f2-mean.f2)/ii;

	 	 }
	 mean.count = ii;
return mean;
}

stars mean_arry_sig(stars * m, stars sig, stars meand, int n ) // нахождение среднего значения массива размером n c заданной сигмой
{
	stars mean;
	mean.x=0;
	mean.f1=0;
	mean.f2=0;

	 if (n < 1) {
		 return mean;
	 }

	 int ii=0;
	 int i;

	 for( i=0; i< n; ++i)
	 	 {
		  if(fabs(meand.f1-m[i].f1)<3*sig.f1 && fabs(meand.f2-m[i].f2)<3*sig.f2)
			{
			  ++ii;
			  mean.x+=(m[i].x-mean.x)/ii;
			  mean.f1+=(m[i].f1-mean.f1)/ii;
			  mean.f2+=(m[i].f2-mean.f2)/ii;

			}
	 	 }

	 mean.count = ii;

return mean;
}

static int load_objects(FILE * input, int x, int f1, int f2, ccarray_t * objects)
{
  char line[MAX_INPUT_LINE_LENGTH];
  int ic;
  char * pc;
  size_t size;
  size_t capacity;

  stars * obj;

  size = ccarray_size(objects);
  capacity = ccarray_capacity(objects);

  while ( size < capacity && !feof(input) )
  {

	obj = ccarray_peek(objects, size);

    line[MAX_INPUT_LINE_LENGTH-1] = 0;
    if ( !fgets(line, MAX_INPUT_LINE_LENGTH, input) ) {
      break;
    }

    if ( line[MAX_INPUT_LINE_LENGTH - 1] != 0 ) {
      fprintf(stderr,"too long input line in this file\n");
      return -1;
    }

    /* remove trailing new line */
    line[strlen(line) - 1] = 0;


    for ( pc = line, ic = 1; ic < x; ++ic ) {
      if ( !(pc = strchr(pc + 1, '\t')) ) {
        break;
      }
    }
    if ( ic != x || sscanf(pc, " %lf", &obj->x) != 1 ) {
      continue;
    }
    for ( pc = line, ic = 1; ic < f1; ++ic ) {
	  if ( !(pc = strchr(pc + 1, '\t')) ) {
		break;
	  }
	}
	if ( ic != f1 || sscanf(pc, " %lf", &obj->f1) != 1 ) {
	  continue;
	}
	if(f2>0)
	{
	for ( pc = line, ic = 1; ic < f2; ++ic ) {
		  if ( !(pc = strchr(pc + 1, '\t')) ) {
			break;
		  }
		}
		if ( ic != f2 || sscanf(pc, " %lf", &obj->f2) != 1 ) {
		  continue;
		}
	}else{obj->f2=0;}

	ccarray_set_size(objects, ++size );

  }

  return size;
}


static void close_file( FILE * input, compression_t compression )
{
  if ( input && input != stdin ) {
	if ( compression > compression_none ) {
	  pclose(input);
	}
	else {
	  fclose(input);
	}
  }
}

static void show_usage( FILE * output, int argc, char * argv[] )
{
	  fprintf(output, "Sort txt file UTILITY\n");
	  fprintf(output, "USAGE:\n");
	  fprintf(output, "  %s OPTIONS FILE\n", basename(argv[0]));
	  fprintf(output, "OPTIONS:\n");
	  fprintf(output, "  x	=	<int>    one-based colum number for statistics\n");
	  fprintf(output, "  f1	=	<int>    one-based colum number for statistics\n");
	  fprintf(output, "  f2	=	<int>    one-based colum number for statistics\n");
	  fprintf(output, "  -sig              used mean in 3 the smallest sigma\n");
	  fprintf(output, "  -med              used median\n");
	  fprintf(output, "  -v                 Print some diagnostic messages to stderr (verbose mode)\n");
	  fprintf(output, "  \n");

	  UNUSED(argc);
}


int main(int argc, char *argv[])
{
	/* TODO: Fix this bug */

	//fprintf(stderr,"Hello, its me!!!\n");

	int beverbose = 0;
	int usedfiles = 0;

	compression_t compression =
		    { compression_unknown};

	const stars * obj;
	ccarray_t * list;
	size_t size,pos;
	size_t capacity = { 50000000 };

	const char * filename = 0;
	FILE * fp = {NULL};

	int x 	= -1;
	int f1 	= -1;
	int f2 	= 0;
	int sig = 0;
	int med = 0;

	int i;

		for (i = 1; i < argc; ++i )
		  {
			 if ( strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-help") == 0 ) {
			      show_usage(stdout, argc, argv);
			      return 0;
			    }

			/* read used columns from file1*/
			if ( strncmp(argv[i], "x=", 2) == 0 )
		        {
		          if ( sscanf(argv[i] + 2, "%d", &x) != 1 || x < 1 )
		          {
		            fprintf(stderr, "Invalid value of %s\n", argv[i]);
		            return 1;
		          }
		        }
			else if ( strncmp(argv[i], "f1=", 3) == 0 )
	        {
	          if ( sscanf(argv[i] + 3, "%d", &f1) != 1 || f1 < 1 )
	          {
	            fprintf(stderr, "Invalid value of %s\n", argv[i]);
	            return 1;
	          }
	        }
			else if ( strncmp(argv[i], "f2=", 3) == 0 )
	        {
	          if ( sscanf(argv[i] + 3, "%d", &f2) != 1 || f2 < 1 )
	          {
	            fprintf(stderr, "Invalid value of %s\n", argv[i]);
	            return 1;
	          }
	        }
		    else if ( strcmp(argv[i], "-v") == 0 ) {
		    	beverbose = 1;
		    	}
		    else if ( !filename ) {
		    	filename = argv[i];
		          usedfiles = 1;
		        }
		    else if ( strcmp(argv[i], "-sig") == 0 ) {
		    	sig = 1;
		        }
		    else if ( strcmp(argv[i], "-med") == 0 ) {
				med = 1;
				}

		    else
		        {
		          fprintf(stderr, "Invalid argument %s. Try %s --help\n", argv[i], argv[0]);
		          return 1;
		        }
		  }


			/* check command line inputs */
			if(usedfiles)
			{
				if ( !filename ) {
				    fprintf(stderr,"Input file name expected\n");
				    show_usage(stderr, argc, argv);
				    return -1;
				  }
				  if ( x < 1 ) {
				      fprintf(stderr,"x argument is mandatory\n");
				      show_usage(stderr, argc, argv);
				      return -1;
				    }
				  if ( f1 < 1 ) {
					  fprintf(stderr,"f1 argument is mandatory\n");
					  show_usage(stderr, argc, argv);
					  return -1;
					}

				  /* check if input files are readable */
				  if ( access(filename, R_OK) != 0 ) {
					  fprintf(stderr, "Can't read %s: %s\n", filename, strerror(errno));
					  return -1;
					}else{
						//FILE * fp = fopen(filename,"r");
					}
			}

				if ( beverbose ) {
					fprintf(stderr,"Used file '%s' and print satatistic colum number  %d, %d, %d \n",filename,x,f1,f2);

					if ( sig==1 )
					{
						fprintf(stderr,"Print satatistic colum number  %d, %d, %d \n",x,f1,f2);

					}
					if ( med==1 )
					{
						fprintf(stderr,"Print median for colum number  %d, %d, %d \n",x,f1,f2);

					}

				}

		/* allocate memory storage */
		  if ( !(list = ccarray_create(capacity, sizeof(stars))) ) {
			  fprintf(stderr, "ccarray_create(capacity=%zu) fails: %s\n", capacity, strerror(errno));
			  return -1;
			}


		/* load input files */
		  if ( beverbose ) {
			fprintf(stderr,"loading %s....\n", filename);
		  }

		  if ( !(fp = open_file(filename, &compression)) ) {
			fprintf(stderr, "Can't read '%s': %s\n", filename, strerror(errno));
			return -1;
		  }

		  int nstars = load_objects(fp, x,f1,f2,list);

		  if ( nstars <1){
			fprintf(stderr, "Can't load data from %s file\n", filename);
			return -1;
		  }

		  close_file( fp, compression );

		  if ( beverbose ) {
			fprintf(stderr,"load: %d objects\n", nstars);
		  }

		  //sort stars


		  /* found max and min*/
		size = ccarray_size(list);
		stars s_array[size];
		int n=0;
		for ( pos = 0; pos < size; ++pos )
		  {
			obj = ccarray_peek(list, pos);
			s_array[pos].x 	= obj->x;
			s_array[pos].f1 	= obj->f1;
			s_array[pos].f2 	= obj->f2;

			++n;

		  }

		stars mean;
		stars temp_mean;

		stars median;

		stars sigma;
		stars temp_sigma;

		temp_sigma.f1 = 1000;
		temp_sigma.f2 = 1000;

		if(n>0)
		{
			if(sig==1)
			{
				//printf("mean_x\tmean_f1\tmean_f2\tsigma_x\tsigma_f1\tsigma_f2\n");
				temp_mean=mean_arry(s_array,n);

				int k=0;
				do
				{
				sigma.f1 = temp_sigma.f1;
				sigma.f2 = temp_sigma.f2;
				mean=mean_arry_sig(s_array,sigma,temp_mean,n);

				temp_sigma=sigma_arry(s_array,sigma,mean,n);

				++k;

				if(beverbose)
				{
					fprintf(stderr,"count iteration = %d\n",k);
				}

				}while(temp_sigma.f1<sigma.f1 && temp_sigma.f2<sigma.f2);

				sigma.x	 = temp_sigma.x;
				sigma.f1 = temp_sigma.f1;
				sigma.f2 = temp_sigma.f2;

				if(k==2)
				{
					mean=mean_arry(s_array,n);
					sigma.x=0;
					sigma.f1=0;
					sigma.f2=0;

				}

				printf("%16.10f\t%16.10f\t%16.10f\t%16.10f\t%16.10f\t%16.10f\t%d\n",mean.x,mean.f1,mean.f2,sigma.x,sigma.f1,sigma.f2,mean.count);

			}
			else
			{
				mean=mean_arry(s_array,n);
				sigma=sigma_arry(s_array,temp_sigma,mean,n);
				printf("%16.10f\t%16.10f\t%16.10f\t%16.10f\t%16.10f\t%16.10f\t%d\n",mean.x,mean.f1,mean.f2,sigma.x,sigma.f1,sigma.f2,mean.count);
			}
			if(med==1)
			{
				//printf("median_x\tmedian_f1\tmedian_f2\n");
				median=smedian_arry(s_array,n-1);
				printf("%16.10f\t%16.10f\t%16.10f\t%d\n",median.x,median.f1,median.f2,median.count);
			}
			//print_res(s_array,n,decrease,head);
		}

	return 0;
}
