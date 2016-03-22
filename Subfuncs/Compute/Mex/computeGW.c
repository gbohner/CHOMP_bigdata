
#include "mex.h"
#include "matrix.h"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //printf("called");
    /* Get the inputs */
    double *GPT, *W1, *W2;
    long m, mom_cur, mom_all;
    
    GPT = mxGetPr(prhs[0]); /* GPT */
    W1 = mxGetPr(prhs[1]); /* GPT */
    W2 = mxGetPr(prhs[2]); /* GPT */
    m = mxGetScalar(prhs[3]); /* GPT */
    mom_all = mxGetScalar(prhs[4]);
    mom_cur = mxGetScalar(prhs[5]);/* GPT */
    //mxArray *W1 = prhs[1]; /* Wcur1c */
    //mxArray *W2 = prhs[2]; /* Wcur2c */
    //mwIndex m = prhs[3]; /* m */
    //mwIndex mom = prhs[4] /* mom */
    
    //printf("input alloc");
    /* Allocate the output */
    mxArray *GW_slice_out;
    GW_slice_out = mxCreateDoubleMatrix(2*m-1, 2*m-1, mxREAL); /* the filt1,filt2,mom slice of GW */
    double *GW_slice;
    GW_slice = mxGetPr(GW_slice_out);
    
    //printf("output alloc");
    /* Setup temporary variables */
    long s1, s2, i1;
    const mwSize *dims;
    mxArray *ind_list_pr;
    double *ind_list;
    
    //printf("temp alloc");
    for( s1=0; s1<(2*m-1); s1++){
        for( s2=0; s2<(2*m-1); s2++){
          //printf("loop");
          ind_list_pr = mxGetCell(prhs[0],(mom_cur-1)*(2*m-1)*(2*m-1) + s2*(2*m-1) + s1);
          ind_list = mxGetPr(ind_list_pr);
          //printf("cell"); 
          dims = mxGetDimensions(ind_list_pr);
          if (s1==0 && s2<2){
            for (i1 = 0; i1<(dims[0]*dims[1]); i1++)
            {
              printf("%d \n", (long) ind_list[i1]); 
            }
            for (i1 = 0; i1<2; i1++)
            {
              printf("%d \n", (long) dims[i1]); 
            }
          }
          GW_slice[s1*(2*m-1)+s2] = 0;
          //printf("slicealloc"); 
          for( i1 =0; i1<dims[0]; i1++){
            if (s1==0 && s2<2){
              printf("\n %d %d", (long)ind_list[i1], (long)ind_list[dims[0] + i1]);
            }
            GW_slice[s1*(2*m-1)+s2] += W1[(long)ind_list[i1]-1] * W2[(long)ind_list[dims[0]+i1]-1];
            //printf("slicealloc_i1"); 
          }
          //printf("slicealloc_done"); 
        }
    }
    
    plhs[0] = GW_slice_out;
    //printf("out_returned"); 
}

