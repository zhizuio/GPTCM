/* header file for arms function */

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
/* #endif */

#ifdef __cplusplus
extern "C" {
#endif

namespace ARMS{

int arms_simple (int ninit, double *xl, double *xr,
	         double (*myfunc)(double x, void *mydata), void *mydata,
                 int dometrop, double *xprev, double *xsamp);

int arms (double *xinit, int ninit, double *xl, double *xr,
	 double (*myfunc)(double x, void *mydata), void *mydata,
         double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
         int nsamp, double *qcent, double *xcent, int ncent,
         int *neval);

//double arms( );

double expshift(double y, double y0);

#define YCEIL 50.                /* maximum y avoiding overflow in exp(y) */


}    /*** end of the namespace ARS ***/

#ifdef __cplusplus
}
#endif

#endif