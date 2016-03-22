
#include "mex.h"
#include "matrix.h"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Get the inputs */
    double *GPT, *W1, *W2;
    long m, mom_cur, mom_all;
    
    GPT = mxGetPr(prhs[0]);
    W1 = mxGetPr(prhs[1]);
    W2 = mxGetPr(prhs[2]);
    m = mxGetScalar(prhs[3]);
    mom_all = mxGetScalar(prhs[4]);
    mom_cur = mxGetScalar(prhs[5]);
    
    /* Allocate the output */
    mxArray *GW_slice_out;
    GW_slice_out = mxCreateDoubleMatrix(2*m-1, 2*m-1, mxREAL); /* the filt1,filt2,mom slice of GW */
    double *GW_slice;
    GW_slice = mxGetPr(GW_slice_out);
    
    /* Setup temporary variables */
    long s1, s2, i1;
    const mwSize *dims;
    mxArray *ind_list_pr;
    double *ind_list;
    
    for( s1=0; s1<(2*m-1); s1++){
        for( s2=0; s2<(2*m-1); s2++){
          ind_list_pr = mxGetCell(prhs[0],(mom_cur-1)*(2*m-1)*(2*m-1) + s2*(2*m-1) + s1);
          ind_list = mxGetPr(ind_list_pr);
          dims = mxGetDimensions(ind_list_pr);
          GW_slice[s1*(2*m-1)+s2] = 0;
          for( i1 =0; i1<dims[0]; i1++){
            GW_slice[s1*(2*m-1)+s2] += W1[(long)ind_list[i1]-1] * W2[(long)ind_list[dims[0]+i1]-1];
          }

        }
    }
    
    plhs[0] = GW_slice_out;
}

