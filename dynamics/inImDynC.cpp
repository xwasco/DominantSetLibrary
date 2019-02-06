/*
 Infection-Immunization dynamics. 

 This mex file has been written by Dr. Samuel Rota Bulo
 
 If you use the Infection Immunization Dynamics please cite this work:
 [2] S. Rota Bulo, and I. M. Bomze.  Infection and immunization:  a new class of evolutionarygame dynamics.
 Games and Economic Behaviour, vol.  71, pp.  193–211, 2011.Special issue in honor of John F. Nash, jr.
 
 Released under MIT License
 2019, Samuel Rota Bulo
 
*/

#include <string.h>
#include "mex.h"

void get_integer_scalar(int &scalar,const mxArray *mat_buf) {
    // Check input
    if (!mxIsInt32(mat_buf))
        mexPrintf("Integer scalar is not int32, but that's OK.\n");

    if (mxGetNumberOfElements(mat_buf) == 1) {
        scalar = mxGetScalar(mat_buf);
    } else {
        mexErrMsgTxt("Integer scalar is not of size == [1 1].\n");
    }
}


void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
double *A=mxGetPr(prhs[0]),*pA;
int n=mxGetN(prhs[0]);
double *Ax=new double[n],*pAx;
double *r=new double[n], *pr;
double *x0=mxGetPr(prhs[1]);
double *ptoll=mxGetPr(prhs[2]);
double toll=*ptoll;

int maxIters;
get_integer_scalar(maxIters,prhs[3]);

toll=toll*toll;

plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
double *x = mxGetPr(plhs[0]),*px;
memcpy(x, x0, n*sizeof(double));

/* Computing Ax */
pAx=Ax; pA=A;
for(int i=0;i<n;++i,++pAx){
  *pAx=0;
  px=x;
  for(int j=0;j<n;++j,++pA,++px)
    *pAx+=*pA**px;
}

double error=0;

int count=0;
//while(true){
while(count < maxIters){    

  /* Computing xAx */
  double xAx=0;
  pAx=Ax; px=x;
  for(int i=0;i<n;++i,++px,++pAx)
    xAx+=*pAx**px;


  /* Computing r=Ax-xAx */
  pAx=Ax; pr=r;
  for(int i=0;i<n;++i,++pr,++pAx) *pr=*pAx-xAx;

  /* Check on Nash error */
  //double error=0;
  error=0;
  px=x; pr=r;
  for(int i=0;i<n;++i,++px,++pr)
    error+=*px>-*pr?*pr**pr:*px**px;
  if(error<toll) break;

//  printf("xAx: %e %e\n",xAx,error);
  
  /* Selecting infective strategy */
  double maxVal=0,minVal=0;
  int maxIdx=-1,minIdx=-1;
  px=x; pr=r;
  for(int i=0;i<n;++i,++px,++pr)
    if(maxVal<*pr){
      maxVal=*pr;
      maxIdx=i;
    }else if(minVal>*pr&&*px>0){
      minVal=*pr;
      minIdx=i;
    }
  int infective=maxVal>=-minVal?maxIdx:minIdx;


  double den=A[infective*(n+1)]-Ax[infective]-r[infective];
  bool do_remove=false;
  double mu;
  double optDelta;
  if(r[infective]>=0){
    mu=1;
    if(den<0){
      optDelta=-r[infective]/den;
      if(optDelta<mu) mu=optDelta;
      //if(mu<0) mu=0;
    }
  }else{
    do_remove=true;
    mu=x[infective]/(x[infective]-1);
    if(den<0){
      optDelta=-r[infective]/den;
      if(optDelta>=mu){
        mu=optDelta;
        do_remove=false;
      }
    }
  }

  px=x;
  for(int i=0;i<n;++i,++px)
    *px=mu*((i==infective)-*px)+*px;
    //*px=(1-mu)**px  
  //x[infective]+=mu;
  if(do_remove)x[infective]=0;
  
  pAx=Ax; pA=A+infective*n;
  for(int i=0;i<n;++i,++pAx,++pA)
    *pAx=mu*(*pA-*pAx)+*pAx;

  ++count;
}

//printf("iterations: %d\n",count);
plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
*mxGetPr(plhs[1])=count;

plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
*mxGetPr(plhs[2])=error;

delete[] Ax;
delete[] r;
}    
