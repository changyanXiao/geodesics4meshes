/*
  qsmodel.c     headerfile for quasistatic probability model

  (c) Michael Schindler
  1997, 1998
  http://www.compressconsult.com/ or http://eiunix.tuwien.ac.at/~michael
  michael@compressconsult.com        michael@eiunix.tuwien.ac.at

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.  It may be that this
  program violates local patents in your country, however it is
  belived (NO WARRANTY!) to be patent-free here in Austria.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.

  Qsmodel is a quasistatic probability model that periodically
  (at chooseable intervals) updates probabilities of symbols;
  it also allows to initialize probabilities. Updating is done more
  frequent in the beginning, so it adapts very fast even without
  initialisation.

  it provides function for creation, deletion, query for probabilities
  and symbols and model updating.

  for usage see example.c
*/

#include "qsmodel.h"
#include <stdio.h>
#include <stdlib.h>

/* default tablesize 1<<TBLSHIFT */
#define TBLSHIFT 7

using namespace std;


// rescale frequency counts
static void dorescale( qsmodel *m)
{   int i, cf, missing;
    if (m->nextleft)  /* we have some more before actual rescaling */
    {   m->incr++;
        m->left = m->nextleft;
        m->nextleft = 0;
        return;
    }
    if (m->rescale < m->targetrescale)  /* double rescale interval if needed */
    {   m->rescale <<= 1;
        if (m->rescale > m->targetrescale)
            m->rescale = m->targetrescale;
    }
    cf = missing = m->cf[m->n];  /* do actual rescaling */
    for(i=m->n-1; i; i--)
    {   int tmp = m->newf[i];
        cf -= tmp;
        m->cf[i] = cf;
        tmp = tmp>>1 | 1;
        missing -= tmp;
        m->newf[i] = tmp;
    }
    if (cf!=m->newf[0])
    {   fprintf(stderr,"BUG: rescaling left %d total frequency\n",cf);
        deleteqsmodel(m);
        exit(1);
    }
    m->newf[0] = m->newf[0]>>1 | 1;
    missing -= m->newf[0];
    m->incr = missing / m->rescale;
    m->nextleft = missing % m->rescale;
    m->left = m->rescale - m->nextleft;
    if (m->search != NULL)
    {   i=m->n;
        while (i)
        {   int start, end;
            end = (m->cf[i]-1) >> m->searchshift;
            i--;
            start = m->cf[i] >> m->searchshift;
            while (start<=end)
            {   m->search[start] = i;
                start++;
            }
        }
    }
}


/* initialisation of qsmodel                           */
/* m   qsmodel to be initialized                       */
/* n   number of symbols in that model                 */
/* lg_totf  base2 log of total frequency count         */
/* rescale  desired rescaling interval, should be < 1<<(lg_totf+1) */
/* init  array of int's to be used for initialisation (NULL ok) */
/* compress  set to 1 on compression, 0 on decompression */
void initqsmodel( qsmodel *m, int n, int lg_totf, int rescale, int *init, int compress )
{   
    m->n = n;
    m->targetrescale = rescale;
    m->searchshift = lg_totf - TBLSHIFT;
    if (m->searchshift < 0)
        m->searchshift = 0;
    m->cf = (unsigned short*) malloc((n+1)*sizeof(uint2));
    m->newf = (unsigned short*) malloc((n+1)*sizeof(uint2));
    m->cf[n] = 1<<lg_totf;
    m->cf[0] = 0;
    if (compress)
        m->search = NULL;
    else
    {   m->search = (unsigned short*) malloc(((1<<TBLSHIFT)+1)*sizeof(uint2));
        m->search[1<<TBLSHIFT] = n-1;
    }
    resetqsmodel(m, init);
}


/* reinitialisation of qsmodel                         */
/* m   qsmodel to be initialized                       */
/* init  array of int's to be used for initialisation (NULL ok) */
void resetqsmodel( qsmodel *m, int *init)
{   int i, end, initval;
    m->rescale = m->n>>4 | 2;
    m->nextleft = 0;
    if (init == NULL)
    {   initval = m->cf[m->n] / m->n;
        end = m->cf[m->n] % m->n;
        for (i=0; i<end; i++)
            m->newf[i] = initval+1;
        for (; i<m->n; i++)
            m->newf[i] = initval;
    } else
        for(i=0; i<m->n; i++)
            m->newf[i] = init[i];
    dorescale(m);
}


/* deletion of qsmodel m                               */
void deleteqsmodel( qsmodel *m )
{   free(m->cf);
    free(m->newf);
    if (m->search != NULL)
        free(m->search);
}


/* retrieval of estimated frequencies for a symbol     */
/* m   qsmodel to be questioned                        */
/* sym  symbol for which data is desired; must be <n   */
/* sy_f frequency of that symbol                       */
/* lt_f frequency of all smaller symbols together      */
/* the total frequency is 1<<lg_totf                   */
void qsgetfreq( qsmodel *m, int sym, int *sy_f, int *lt_f )
{   *sy_f = m->cf[sym+1] - (*lt_f = m->cf[sym]);
}


/*same as previously but with maximal contextual value */
void qsgetfreqcont( qsmodel *m, int sym, int max_context, int *sy_f, int *lt_f )
{   
	if(sym >= max_context) printf(" ERROR in qsgetfreqcont!!!!!!!\n");
	double factor = (double)m->cf[m->n]/(double)m->cf[max_context];
	*lt_f = (int)(((double)m->cf[sym]*factor));
	*sy_f = (int)(((double)m->cf[sym+1]*factor)- *lt_f);
	if(sym+1 == max_context)
		*sy_f = m->cf[m->n]- *lt_f;
}

/*same as previously but with minimal and maximal contextual values */
void qsgetfreqcont( qsmodel *m, int sym, int min_context, int max_context, int *sy_f, int *lt_f )
{   
	if(sym >= max_context) printf(" ERROR in qsgetfreqcont!!!!!!!\n");
	if(sym < min_context) printf(" ERROR in qsgetfreqcont!!!!!!!\n");
	double factor = (double)(m->cf[m->n])/(double)(m->cf[max_context]-m->cf[min_context]);
	*lt_f = (int)(((double)(m->cf[sym]-m->cf[min_context])*factor));
	*sy_f = (int)(((double)(m->cf[sym+1]-m->cf[min_context])*factor)- *lt_f);
	if(sym+1 == max_context)
		*sy_f = m->cf[m->n]- *lt_f;
}

//get frequency for a combination encoding: only valid for 0 and 1
void qsgetfreqcomb( qsmodel *m, int sym, double proba0, int *sy_f, int *lt_f )
{   
    if(sym > 1)     printf(" ERROR in qsgetfreqcomb!!!!!!!\n");

	// separator divides the intervall according to proba0
    int separator = (int)(proba0*m->cf[2]); 

	if(sym == 0) 
	{
	  *lt_f = 0;
      *sy_f = separator; 
	}
	
    if(sym == 1)
	{
	  *lt_f = separator; 
	  *sy_f = m->cf[2]-(*lt_f);
	}
}


//void qsgetfreqcomb( qsmodel *m, int sym, double proba0, int *sy_f, int *lt_f )
void qsgetfreqcombmulti( qsmodel *m, double cum_proba0, double cum_proba1, int *sy_f, int *lt_f )
{   
	// separator divides the intervall according to proba0
   int inf_separator = (int)(cum_proba0*m->cf[m->n]); 
   int sup_separator = (int)(cum_proba1*m->cf[m->n]);
   
   if(sup_separator > m->cf[m->n])
   {
     fprintf(stderr,"PROBLEME DANS qsmodel.cpp: sup_separator > m->cf[m->n]\n");    
   }

 
   if(inf_separator == sup_separator)
   {
     fprintf(stderr,"PROBLEME DANS qsmodel.cpp: difference entre les proba trop petites, inf et sup_sep egaux\n");
     printf("inf_separator : %d, sup_separator : %d \n \n",inf_separator,sup_separator);
     printf("cum_proba0 : %g, cum_proba1 : %g \n \n",cum_proba0,cum_proba1);
     
   }

   *lt_f = inf_separator;
   *sy_f = sup_separator-(*lt_f); 
   
}



// find out symbol for a given cumulative frequency  
// m   qsmodel to be questioned                        
// lt_f  cumulative frequency     
int qsgetsym( qsmodel *m, int lt_f )
{   int lo, hi;
    uint2 *tmp;
    tmp = m->search+(lt_f>>m->searchshift);
    lo = *tmp;
    hi = *(tmp+1) + 1;
    while (lo+1 < hi )
    {   int mid = (lo+hi)>>1;
        if (lt_f < m->cf[mid])
            hi = mid;
        else
            lo = mid;
    }
    return lo;
}


//same as previous with context 
void qsgetsymcont( qsmodel *m, int lt_f, int max_context,int* Value )
{  
	int lo, hi;
    uint2 *tmp;
	double factor = (double)(m->cf[m->n])/(double)(m->cf[max_context]);
	int local_lt_f = (int)((double)(lt_f)/factor);
    tmp = m->search+(local_lt_f>>m->searchshift);
    lo = int(*tmp);
    hi = int(*(tmp+1) + 1);
    while (lo+1 < hi )
    {   int mid = (lo+hi)>>1;
        if (local_lt_f < m->cf[mid])
            hi = mid;
        else
            lo = mid;
    }
    //return lo;
	*Value = int(lo);
}


void qsgetsymcont( qsmodel *m, int lt_f, int min_context,int max_context,int* Value )
{  
	int lo, hi;
    uint2 *tmp;
	double factor = (double)(m->cf[m->n])/(double)(m->cf[max_context]-m->cf[min_context]);
	int local_lt_f = (int)((double)(lt_f)/factor) + m->cf[min_context];
    tmp = m->search+(local_lt_f>>m->searchshift) ;
    lo = int(*tmp);
    hi = int(*(tmp+1) + 1);
    while (lo+1 < hi )
    {   int mid = (lo+hi)>>1;
        if (local_lt_f < m->cf[mid])
            hi = mid;
        else
            lo = mid;
    }
    //return lo;
	*Value = int(lo);
}


void qsgetsymcomb( qsmodel *m, int lt_f, double proba0, int* Value )
{
  int separator = (int)(proba0*m->cf[2]); 

  if(lt_f < separator)
	*Value = 0;
  else
	*Value = 1;
}


void qsgetcombmult( qsmodel *m, int lt_f_dec, double* proba, int* Value,  int *st_f, int *lt_f)
{
  int inf_separator;
  int sup_separator;
  inf_separator = (int)(proba[0]*m->cf[m->n]);
    
  if(lt_f_dec < inf_separator)
  {
    *Value = 0;
    *lt_f = 0;
    *st_f = inf_separator;
  }
  else
  {
    for(int i=1;i<m->n;i++)
    {
      inf_separator = (int)(proba[i-1]*m->cf[m->n]);
      sup_separator = (int)(proba[i]*m->cf[m->n]);
      if(lt_f_dec < sup_separator && lt_f_dec >= inf_separator)
      {
        *Value = i;
        *lt_f = inf_separator;
        *st_f = sup_separator-*lt_f;
        break;
      }
    }
  }
  
}


// update model                      
// m   qsmodel to be updated                           
// sym  symbol that occurred (must be <n from init)    
void qsupdate( qsmodel *m, int sym )
{   if (m->left <= 0)
        dorescale(m);
    m->left--;
    m->newf[sym] += m->incr;
}
