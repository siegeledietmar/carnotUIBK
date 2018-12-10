/*
 * pdepesim_c_v0.c: 
 *
 * Copyright 
 * $Revision: 0.1 $
 */

#define S_FUNCTION_NAME  pdepesim_hygro_c_v4
#define S_FUNCTION_LEVEL 2

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include <math.h>
#include "simstruc.h"

/*
 *   Defines for easy access to the parameters
 */
#define S_XMESH(S)      ssGetSFcnParam(S,0)         /* xmesh */
#define NX     *mxGetPr(ssGetSFcnParam(S,1))        /* xmesh length */
#define S_D(S)          ssGetSFcnParam(S,2)         /* d */
#define ND     *mxGetPr(ssGetSFcnParam(S,3))        /* d length */
#define S_LAMBDA(S)      ssGetSFcnParam(S,4)        /* lambda */
#define S_RHO(S)      ssGetSFcnParam(S,5)           /* rho */
#define S_CP(S)      ssGetSFcnParam(S,6)            /* cp */
#define S_MU(S)      ssGetSFcnParam(S,7)            /* mu */
#define S_UFS(S)      ssGetSFcnParam(S,8)           /* ufs */
#define S_FSF_u(S)      ssGetSFcnParam(S,9)         /* FSF_u */
#define S_FSF_phi(S)      ssGetSFcnParam(S,10)      /* FSF_phi */
#define S_FSF_dudphi(S)      ssGetSFcnParam(S,11)   /* FSF_dudphi */
#define S_FSF_k(S)      ssGetSFcnParam(S,12)        /* FSF_k */
#define S_K_kleff(S)      ssGetSFcnParam(S,13)      /* K_kleff */
#define S_K_u(S)      ssGetSFcnParam(S,14)          /* K_u */
#define S_K_kl(S)      ssGetSFcnParam(S,15)         /* K_kl */
#define S_FSF_length(S)      ssGetSFcnParam(S,16)   /* FSF_length */
#define S_K_length(S)      ssGetSFcnParam(S,17)     /* K_length */
#define S_T_ini(S)      ssGetSFcnParam(S,18)        /* T_ini */
#define S_Phi_ini(S)      ssGetSFcnParam(S,19)      /* Phi_ini */
#define S_T_source(S)      ssGetSFcnParam(S,20)     /* T_source */
#define S_Phi_source(S)      ssGetSFcnParam(S,21)   /* Phi_source */
#define S_Tdactive(S) *mxGetPr(ssGetSFcnParam(S,22))       /* Tdactive */
#define S_Phidactive(S) *mxGetPr(ssGetSFcnParam(S,23))     /* Phidactive */

#define NPDE            2
#define NN              (NX*NPDE)

#define NMAX            1000
#define NNMAX           500


/* Error handling
 * --------------
 *
 * You should use the following technique to report errors encountered within
 * an S-function:
 *
 *       ssSetErrorStatus(S,"Error encountered due to ...");
 *       return;
 *
 * Note that the 2nd argument to ssSetErrorStatus must be persistent memory.
 * It cannot be a local variable. For example the following will cause
 * unpredictable errors:
 *
 *      mdlOutputs()
 *      {
 *         char msg[256];         {ILLEGAL: to fix use "static char msg[256];"}
 *         sprintf(msg,"Error due to %s", string);
 *         ssSetErrorStatus(S,msg);
 *         return;
 *      }
 *
 * See matlabroot/simulink/src/sfuntmpl_doc.c for more details.
 */

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    
    /*
    printf("NX : %d  %d\n",(int) NX,__LINE__);
    printf("NN : %d  %d\n",(int) NN,__LINE__);
    */
    
    
    /* See sfuntmpl_doc.c for more details on the macros below */
    
    ssSetNumSFcnParams(S, 24);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }
    
    /* number of cont/disc states */
    ssSetNumContStates(S, (int)NN);
    ssSetNumDiscStates(S, 0);
    
    /* number of input ports */
    if (!ssSetNumInputPorts(S, 10)) return;
	
    /* dimension of input ports */
	// qright
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 0, 0);
	// qleft
    ssSetInputPortWidth(S, 1, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 0);
	// Tright
    ssSetInputPortWidth(S, 2, 1);
    ssSetInputPortDirectFeedThrough(S, 2, 0);
	// phiright
    ssSetInputPortWidth(S, 3, 1);
    ssSetInputPortDirectFeedThrough(S, 3, 0);
	// Tleft
    ssSetInputPortWidth(S, 4, 1);
    ssSetInputPortDirectFeedThrough(S, 4, 0);
	// phileft
    ssSetInputPortWidth(S, 5, 1);
    ssSetInputPortDirectFeedThrough(S, 5, 0);
	// alphae
    ssSetInputPortWidth(S, 6, 1);
    ssSetInputPortDirectFeedThrough(S, 6, 0);
	// alphai
    ssSetInputPortWidth(S, 7, 1);
    ssSetInputPortDirectFeedThrough(S, 7, 0);
	// Tsource
    ssSetInputPortWidth(S, 8, 1);
    ssSetInputPortDirectFeedThrough(S, 8, 0);
	// Phisource
    ssSetInputPortWidth(S, 9, 1);
    ssSetInputPortDirectFeedThrough(S, 9, 0);
        /*
         * Set direct feedthrough flag (1=yes, 0=no).
         * A port has direct feedthrough if the input is used in either
         * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
         * See matlabroot/simulink/src/sfuntmpl_directfeed.txt.
         */
    
    /* number of output ports */
    if (!ssSetNumOutputPorts(S, 10)) return;
    
    /* dimension of output ports */
    ssSetOutputPortWidth(S, 0, NX);
    ssSetOutputPortWidth(S, 1, NX);
    ssSetOutputPortWidth(S, 2, 1);
    ssSetOutputPortWidth(S, 3, 1);
    ssSetOutputPortWidth(S, 4, 1);
    ssSetOutputPortWidth(S, 5, 1);
    ssSetOutputPortWidth(S, 6, NX);
    ssSetOutputPortWidth(S, 7, 1);
    ssSetOutputPortWidth(S, 8, 1);
    ssSetOutputPortWidth(S, 9, 1);
    
    /* sample times */
    ssSetNumSampleTimes(S, 1);
    
    /* DWork */
    ssSetNumDWork(S, 20);
    ssSetDWorkWidth(S, 0, NMAX);      // xi
    ssSetDWorkWidth(S, 1, NMAX);      // xim
    ssSetDWorkWidth(S, 2, NMAX);      // zxmp1
    ssSetDWorkWidth(S, 3, NMAX);      // xzmp1
    ssSetDWorkWidth(S, 4, NMAX);      // lambda
    ssSetDWorkWidth(S, 5, NMAX);      // rho
    ssSetDWorkWidth(S, 6, NMAX);      // cp
    ssSetDWorkWidth(S, 7, NMAX);      // mu
    ssSetDWorkWidth(S, 8, NMAX);      // ufs
    ssSetDWorkWidth(S, 9, NMAX*NNMAX);      // FSF_u
    ssSetDWorkWidth(S, 10, NMAX*NNMAX);     // FSF_phi
    ssSetDWorkWidth(S, 11, NMAX*NNMAX);     // FSF_dudphi
    ssSetDWorkWidth(S, 12, NMAX*2);     // FSF_k
    ssSetDWorkWidth(S, 13, NMAX);     // K_kleff
    ssSetDWorkWidth(S, 14, NMAX*NNMAX);     // K_u
    ssSetDWorkWidth(S, 15, NMAX*NNMAX);     // K_kl
    ssSetDWorkWidth(S, 16, NMAX);     // FSF_length
    ssSetDWorkWidth(S, 17, NMAX);     // K_length
    ssSetDWorkWidth(S, 18, NMAX);     // T_source
    ssSetDWorkWidth(S, 19, NMAX);     // Phi_source
    
    ssSetDWorkDataType(S, 0, SS_DOUBLE);
    ssSetDWorkDataType(S, 1, SS_DOUBLE);
    ssSetDWorkDataType(S, 2, SS_DOUBLE);
    ssSetDWorkDataType(S, 3, SS_DOUBLE);
    ssSetDWorkDataType(S, 4, SS_DOUBLE);
    ssSetDWorkDataType(S, 5, SS_DOUBLE);
    ssSetDWorkDataType(S, 6, SS_DOUBLE);
    ssSetDWorkDataType(S, 7, SS_DOUBLE);
    ssSetDWorkDataType(S, 8, SS_DOUBLE);
    ssSetDWorkDataType(S, 9, SS_DOUBLE);
    ssSetDWorkDataType(S, 10, SS_DOUBLE);
    ssSetDWorkDataType(S, 11, SS_DOUBLE);
    ssSetDWorkDataType(S, 12, SS_DOUBLE);
    ssSetDWorkDataType(S, 13, SS_DOUBLE);
    ssSetDWorkDataType(S, 14, SS_DOUBLE);
    ssSetDWorkDataType(S, 15, SS_DOUBLE);
    ssSetDWorkDataType(S, 16, SS_DOUBLE);
    ssSetDWorkDataType(S, 17, SS_DOUBLE);
    ssSetDWorkDataType(S, 18, SS_DOUBLE);
    ssSetDWorkDataType(S, 19, SS_DOUBLE);
    
//     ssSetNumRWork(S, (int)NN);       /* real vector    */
//     ssSetNumIWork(S, 1);        /* integer vector */
//     ssSetNumPWork(S, 0);        /* pointer vector */
//     ssSetNumModes(S, 0);        /* mode vector    */
//     ssSetNumNonsampledZCs(S, 0);
    
    /* Dimension Modes of the Ports */
    /*ssSetOutputPortDimensionsMode(S, 2, INHERIT_DIMS_MODE);
    ssSetInputPortRequiredContiguous(S, 0, true);
    ssSetInputPortRequiredContiguous(S, 1, true);*/

    /* Specify the sim state compliance to be same as a built-in block */
    /* ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);
    ssSetOptions(S, 0); */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
    /* ssSetOptions(S,
                 SS_OPTION_WORKS_WITH_CODE_REUSE |
                 SS_OPTION_EXCEPTION_FREE_CODE |
                 SS_OPTION_USE_TLC_WITH_ACCELERATOR); */
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
//     ssSetSampleTime(S, 0, VARIABLE_SAMPLE_TIME);
//     ssSetOffsetTime(S, 0, 0.0);
//     ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
//     ssSetOffsetTime(S, 0, 0.0);
//     ssSetModelReferenceSampleTimeDefaultInheritance(S);
}



#define MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */
#if defined(MDL_INITIALIZE_CONDITIONS)
  /* Function: mdlInitializeConditions ========================================
   * Abstract:
   *    In this function, you should initialize the continuous and discrete
   *    states for your S-function block.  The initial states are placed
   *    in the state vector, ssGetContStates(S) or ssGetRealDiscStates(S).
   *    You can also perform any other initialization activities that your
   *    S-function may require. Note, this routine will be called at the
   *    start of simulation and if it is present in an enabled subsystem
   *    configured to reset states, it will be call when the enabled subsystem
   *    restarts execution to reset the states.
   */
  static void mdlInitializeConditions(SimStruct *S)
  {
    real_T  *x0     = ssGetContStates(S);
    
    real_T  *XMESH  = mxGetPr(S_XMESH(S));
    real_T  *D  = mxGetPr(S_D(S));
    real_T  *x_LAMBDA  = mxGetPr(S_LAMBDA(S));
    real_T  *x_RHO  = mxGetPr(S_RHO(S));
    real_T  *x_CP  = mxGetPr(S_CP(S));
    real_T  *x_MU  = mxGetPr(S_MU(S));
    real_T  *x_UFS  = mxGetPr(S_UFS(S));
    real_T  *x_FSF_u  = mxGetPr(S_FSF_u(S));
    real_T  *x_FSF_phi  = mxGetPr(S_FSF_phi(S));
    real_T  *x_FSF_dudphi  = mxGetPr(S_FSF_dudphi(S));
    real_T  *x_FSF_k  = mxGetPr(S_FSF_k(S));
    real_T  *x_K_kleff  = mxGetPr(S_K_kleff(S));
    real_T  *x_K_u  = mxGetPr(S_K_u(S));
    real_T  *x_K_kl  = mxGetPr(S_K_kl(S));
    real_T  *x_FSF_length  = mxGetPr(S_FSF_length(S));
    real_T  *x_K_length  = mxGetPr(S_K_length(S));
    
    real_T  *T_ini  = mxGetPr(S_T_ini(S));
    real_T  *Phi_ini  = mxGetPr(S_Phi_ini(S));
    
    real_T  *x_T_source  = mxGetPr(S_T_source(S));
    real_T  *x_Phi_source  = mxGetPr(S_Phi_source(S));
    
    real_T  *xi_    = ssGetDWork(S,0);
    real_T  *xim_   = ssGetDWork(S,1);
    real_T  *zxmp1_ = ssGetDWork(S,2);
    real_T  *xzmp1_ = ssGetDWork(S,3);
    
    real_T  *lambda = ssGetDWork(S,4);
    real_T  *rho = ssGetDWork(S,5);
    real_T  *cp = ssGetDWork(S,6);
    real_T  *mu = ssGetDWork(S,7);
    real_T  *ufs = ssGetDWork(S,8);
    real_T  *FSF_u = ssGetDWork(S,9);
    real_T  *FSF_phi = ssGetDWork(S,10);
    real_T  *FSF_dudphi = ssGetDWork(S,11);
    real_T  *FSF_k = ssGetDWork(S,12);
    real_T  *K_kleff = ssGetDWork(S,13);
    real_T  *K_u = ssGetDWork(S,14);
    real_T  *K_kl = ssGetDWork(S,15);
    real_T  *FSF_length = ssGetDWork(S,16);
    real_T  *K_length = ssGetDWork(S,17);
    
    real_T  *T_source = ssGetDWork(S,18);
    real_T  *Phi_source = ssGetDWork(S,19);
    
    
    int     L;
    int     j;
    int     jj;
    int     nnm  = ((int)NN-1);
    int     nnx  = ((int)NX-1);
    int     nnpde  = (int)NPDE;
    
    double  h[NMAX],
            xi[NMAX],
            xim[NMAX],
            zeta[NMAX],
            midpoint[NMAX],
            zxmp1[NMAX],
            xzmp1[NMAX+1];
    
    int     m = 0;
    int     n;
    
    
    // initialize phi, T
    for (n = 0; n < NX; n++)
    {
        x0[n*nnpde] = Phi_ini[n];      // phi
        x0[n*nnpde+1] = T_ini[n];      // teta
    }
    
    // save constant sources to DWork
    for (n = 0; n < NX; n++)
    {
        T_source[n] = x_T_source[n];
        Phi_source[n] = x_Phi_source[n];
    }
    
    // we calculate 
    //      h ... length of element
    //      midpoint ... center of element
    //      xi ... in this case center of element
    //      zeta ... in this case xi
    //      xim ... xi^m, in this case m=0 always =1
    for (L = 0; L<(nnx); L++)
    {
        h[L] = XMESH[L+1]-XMESH[L];
        midpoint[L] = XMESH[L] + 0.5 * h[L];
        /*printf("h : %f  %d\n",h[L],__LINE__);
        printf("midpoint : %f  %d\n",midpoint[L],__LINE__);*/
        
        xi[L] = midpoint[L];
        zeta[L] = xi[L];

//         xim[L] = pow(xi[L],m);
        xim[L] = 1;
    }
    
    // we calculate zxmp1 and xzmp1 for first element
//     TODO: m > 0, singularity
//     zxmp1[0] = pow(zeta[1],(m+1)) - pow(XMESH[1],(m+1));
    zxmp1[0] = zeta[1] - XMESH[1];
    xzmp1[0] = 0;
    // and now for all mid elements
    for (L = 1; L<(nnx); L++)
    {
//         zxmp1[L] = pow(zeta[L],(m+1)) - pow(XMESH[L],(m+1));
//         xzmp1[L] = pow(XMESH[L],(m+1)) - pow(zeta[L-1],(m+1));
        zxmp1[L] = zeta[L] - XMESH[L];
        xzmp1[L] = XMESH[L] - zeta[L-1];
    }
    // and for the last element
//     xzmp1[nnx] = pow(XMESH[nnx-1],(m+1)) - pow(zeta[nnx-2],(m+1));
    xzmp1[nnx] = XMESH[nnx-1] - zeta[nnx-2];
    zxmp1[nnx] = 0;       // fill the rang
    
    /*
    printf("XMESH : %f  %d\n",XMESH[nnx-1],__LINE__);
    printf("zeta : %f  %d\n",zeta[nnx-2],__LINE__);
    printf("nnx : %d  %d\n",nnx,__LINE__);
    printf("zxmp1 : %f  %d\n",zxmp1[nnx-1],__LINE__);
    printf("xzmp1 : %f  %d\n",xzmp1[nnx-1],__LINE__);
    */
    
//     // for case m>0: divide zxmp1 and xzmp1
//     for (L = 0; L<=(nnx); L++)
//     {
//         zxmp1[L] = zxmp1[L] / (m+1);
//         xzmp1[L] = xzmp1[L] / (m+1);
// //         printf("zxmp1 : %f  %d\n",xzmp1[L],L);
// //         printf("xzmp1 : %f  %d\n",zxmp1[L],L);
//     }
    
    // now we write the values to DWork
    for (L = 0; L<=(nnx); L++)
    {
        xi_[L] = xi[L];
        xim_[L] = xim[L];
        zxmp1_[L] = zxmp1[L];
        xzmp1_[L] = xzmp1[L];
    }
    
    // define materials
    // we create a loop for all nodes expect the last one
    for (L = 0; L<nnx; L++)
    {
        // first layer
        if (XMESH[L] < D[0])
        {
            lambda[L] = x_LAMBDA[0];
            rho[L] = x_RHO[0];
            cp[L] = x_CP[0];
            mu[L] = x_MU[0];
            ufs[L] = x_UFS[0];
            
            K_kleff[L] = x_K_kleff[0];
            
            FSF_length[L] = x_FSF_length[0];
            K_length[L] = x_K_length[0];
            
            // loop for all values FSF and K
            for (jj=0; jj<NNMAX; jj++)
            {
                FSF_u[L*NNMAX+jj] = x_FSF_u[jj];
                FSF_phi[L*NNMAX+jj] = x_FSF_phi[jj];
                FSF_dudphi[L*NNMAX+jj] = x_FSF_dudphi[jj];
                
                K_u[L*NNMAX+jj] = x_K_u[jj];
                K_kl[L*NNMAX+jj] = x_K_kl[jj];
                
//                 printf("x_K_u[jj] : %f  %d\n",x_K_u[jj],jj);
//                 printf("x_K_kl[jj] : %f  %d\n",x_K_kl[jj],jj);
            }
        
        // all other layers
        } else {
//             printf("XMESH[L] : %f  %d\n",XMESH[L],L);
            for (j=0; j<ND; j++)
            {
//                 printf("j : %d  %d\n",j,L);
//                 printf("D[j] : %f  %d\n",D[j],L);
                // in some cases this numeric cheat could be a problem
                if (XMESH[L] >= (D[j]-.0001) && XMESH[L] < (D[j+1]-.0001))
                {
                    lambda[L] = x_LAMBDA[j+1];
                    rho[L] = x_RHO[j+1];
                    cp[L] = x_CP[j+1];
                    mu[L] = x_MU[j+1];
                    ufs[L] = x_UFS[j+1];
                    
                    K_kleff[L] = x_K_kleff[j+1];
                    
                    FSF_length[L] = x_FSF_length[j+1];
                    K_length[L] = x_K_length[j+1];
                    
                    // loop for all values FSF and K
                    for (jj=0; jj<NNMAX; jj++)
                    {
                        FSF_u[L*NNMAX+jj] = x_FSF_u[(j+1)*NNMAX+jj];
                        FSF_phi[L*NNMAX+jj] = x_FSF_phi[(j+1)*NNMAX+jj];
                        FSF_dudphi[L*NNMAX+jj] = x_FSF_dudphi[(j+1)*NNMAX+jj];
                        
                        K_u[L*NNMAX+jj] = x_K_u[(j+1)*NNMAX+jj];
                        K_kl[L*NNMAX+jj] = x_K_kl[(j+1)*NNMAX+jj];
                        
//                         printf("x_K_u[jj] : %f  %d\n",x_K_u[(j+1)*NNMAX+jj],jj);
//                         printf("x_K_kl[jj] : %f  %d\n",x_K_kl[(j+1)*NNMAX+jj],jj);
                        
//                         printf("FSF_u : %f  %d\n",FSF_u[L*NNMAX+jj],L);
                    }
                    
                }
            }
        }
        
    }
    
    // values for the last node
    lambda[nnx] = x_LAMBDA[(int)ND-1];
    rho[nnx] = x_RHO[(int)ND-1];
    cp[nnx] = x_CP[(int)ND-1];
    mu[nnx] = x_MU[(int)ND-1];
    ufs[nnx] = x_UFS[(int)ND-1];
    
    K_kleff[nnx] = x_K_kleff[(int)ND-1];
    
    FSF_length[nnx] = x_FSF_length[(int)ND-1];
    K_length[nnx] = x_K_length[(int)ND-1];
    
    // loop for all values FSF and K
    for (jj=0; jj<NNMAX; jj++)
    {
        FSF_u[L*NNMAX+jj] = x_FSF_u[((int)ND-1)*NNMAX+jj];
        FSF_phi[L*NNMAX+jj] = x_FSF_phi[((int)ND-1)*NNMAX+jj];
        FSF_dudphi[L*NNMAX+jj] = x_FSF_dudphi[((int)ND-1)*NNMAX+jj];
        
        K_u[L*NNMAX+jj] = x_K_u[((int)ND-1)*NNMAX+jj];
        K_kl[L*NNMAX+jj] = x_K_kl[((int)ND-1)*NNMAX+jj];
    }
    
    // here we can check all values
//     for (L = 0; L<=nnx; L++)
//     {
//         printf("XMESH[L] : %f  %d\n",XMESH[L],L);
//         printf("lambda : %f  %d\n",lambda[L],L);
//         printf("rho : %f  %d\n",rho[L],L);
//         printf("cp : %f  %d\n",cp[L],L);
//         printf("mu : %f  %d\n",mu[L],L);
//         printf("ufs : %f  %d\n",ufs[L],L);
//         printf("FSF_length : %f  %d\n",FSF_length[L],L);
//     }
    
  }
#endif /* MDL_INITIALIZE_CONDITIONS */



#undef MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
  /* Function: mdlStart =======================================================
   * Abstract:
   *    This function is called once at start of model execution. If you
   *    have states that should be initialized once, this is the place
   *    to do it.
   */
  static void mdlStart(SimStruct *S)
  {
  }
#endif /*  MDL_START */



/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block.
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T   *u    = ssGetContStates(S);    
    real_T   *y0   = ssGetOutputPortRealSignal(S,0);
    real_T   *y1   = ssGetOutputPortRealSignal(S,1);
    real_T   *y2   = ssGetOutputPortRealSignal(S,2);
    real_T   *y3   = ssGetOutputPortRealSignal(S,3);
    real_T   *y4   = ssGetOutputPortRealSignal(S,4);
    real_T   *y5   = ssGetOutputPortRealSignal(S,5);
    real_T   *y6   = ssGetOutputPortRealSignal(S,6);
    real_T   *y7   = ssGetOutputPortRealSignal(S,7);
    real_T   *y8   = ssGetOutputPortRealSignal(S,8);
    real_T   *y9   = ssGetOutputPortRealSignal(S,9);
    
    int     nnm  = ((int)NN-1);
    int     nnx  = ((int)NX-1);
    int     nnpde  = (int)NPDE;
    int     n;
    int     j;
    int     jj;
    
    int    L;
    
    real_T  *XMESH    = mxGetPr(S_XMESH(S));
    
    real_T  *D  = mxGetPr(S_D(S));
    
    real_T  *ufs = ssGetDWork(S,8);
    real_T  *FSF_u = ssGetDWork(S,9);
    real_T  *FSF_phi = ssGetDWork(S,10);
    real_T  *FSF_length = ssGetDWork(S,16);
    
    real_T  *Tdactive  = (int)S_Tdactive(S);
    real_T  *Phidactive  = (int)S_Phidactive(S);
    
    double u_phi;
    
    double u_meas[NNMAX];
    double phi_meas[NNMAX];
    
    int ip_idx;
    double ip_steps;
    
    double DD;
    
    
    // write all nodes for phi
    j = 0;
    for (n = 0; n < NX; n++)
    {
        y0[j] = u[n*nnpde];
        j++;
    }
    
    // write all nodes for T
    j = 0;
    for (n = 0; n < NX; n++)
    {
        y1[j] = u[n*nnpde+1];
        j++;
    }
    
    // write right an left theta
    y2[0] = u[nnm]-273.15;
    y3[0] = u[1]-273.15;
    
    // write right an left phi
    y4[0] = u[nnm-1];
    y5[0] = u[0];
    
    
    // loop for all nodes
    DD = 0.0;
    for (L = 0; L < nnx; L++)
    {
        /* read phi, u */
        for (jj = 0; jj<NNMAX; jj++)
        {
            phi_meas[jj] = FSF_phi[L*NNMAX+jj];
            u_meas[jj] = ufs[L]*FSF_u[L*NNMAX+jj];
        }
        
        // interpolate sorption isotherm
        u_phi = 0;
        if (u[L*nnpde]>0.0 && u[L*nnpde]<1.0) 
        {
            ip_idx = 0;
            while(ip_idx < FSF_length[L] && u[L*nnpde] <= phi_meas[ip_idx])
            {
                ip_idx++;
            }
            
            ip_steps = (u_meas[ip_idx] - u_meas[ip_idx-1]) / (phi_meas[ip_idx] - phi_meas[ip_idx-1]);
            u_phi = u_meas[ip_idx] - (ip_steps * (phi_meas[ip_idx] - u[L*nnpde]) );
        }
        
        // write value to output
        y6[L] = u_phi;
        
        
        // calculate moisture content
//         printf("XMESH %d: %f \n", L, (XMESH[L+1]-XMESH[L]));
//         printf("%f \n", u_phi);
        DD = DD + (XMESH[L+1]-XMESH[L]) * u_phi;
    }
    
    // moisture content in whole structure
    y7[0] = DD;
    
    
    // phi of element with source
    n = Phidactive;
    if (n > 0)
    {
        y8[0] = (u[n*nnpde-2] + u[n*nnpde]) / 2;
    }
    
    // temperature of element with source
    n = Tdactive;
    if (n > 0)
    {
        y9[0] = (u[n*nnpde+1-2] + u[n*nnpde+1]) / 2 - 273.15;
    }
    
}



#undef MDL_UPDATE  /* Change to #undef to remove function */
#if defined(MDL_UPDATE)
  /* Function: mdlUpdate ======================================================
   * Abstract:
   *    This function is called once for every major integration time step.
   *    Discrete states are typically updated here, but this function is useful
   *    for performing any tasks that should only take place once per
   *    integration step.
   */
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
  }
#endif /* MDL_UPDATE */



#define MDL_DERIVATIVES  /* Change to #undef to remove function */
#if defined(MDL_DERIVATIVES)
  /* Function: mdlDerivatives =================================================
   * Abstract:
   *    In this function, you compute the S-function block's derivatives.
   *    The derivatives are placed in the derivative vector, ssGetdX(S).
   */
  static void mdlDerivatives(SimStruct *S)
  {
    real_T            *du = ssGetdX(S);
    real_T            *u  = ssGetContStates(S);
    InputRealPtrsType u0  = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType u1  = ssGetInputPortRealSignalPtrs(S,1);
    InputRealPtrsType u2  = ssGetInputPortRealSignalPtrs(S,2);
    InputRealPtrsType u3  = ssGetInputPortRealSignalPtrs(S,3);
    InputRealPtrsType u4  = ssGetInputPortRealSignalPtrs(S,4);
    InputRealPtrsType u5  = ssGetInputPortRealSignalPtrs(S,5);
    InputRealPtrsType u6  = ssGetInputPortRealSignalPtrs(S,6);
    InputRealPtrsType u7  = ssGetInputPortRealSignalPtrs(S,7);
    InputRealPtrsType u8  = ssGetInputPortRealSignalPtrs(S,8);
    InputRealPtrsType u9  = ssGetInputPortRealSignalPtrs(S,9);
    
    real_T            qright = (*u0[0]);
    real_T            qleft = (*u1[0]);
    real_T            Tright = (*u2[0]);
    real_T            phiright = (*u3[0]);
    real_T            Tleft = (*u4[0]);
    real_T            phileft = (*u5[0]);
    real_T            alphae = (*u6[0]);
    real_T            alphai = (*u7[0]);
    real_T            Tsource = (*u8[0]);
    real_T            Phisource = (*u9[0]);
    
    real_T  *XMESH    = mxGetPr(S_XMESH(S));
    
    real_T  *Tdactive  = (int)S_Tdactive(S);
    real_T  *Phidactive  = (int)S_Phidactive(S);
    
    real_T  *xi    = ssGetDWork(S,0);
    real_T  *xim   = ssGetDWork(S,1);
    real_T  *zxmp1 = ssGetDWork(S,2);
    real_T  *xzmp1 = ssGetDWork(S,3);
    
    real_T  *lambda = ssGetDWork(S,4);
    real_T  *rho = ssGetDWork(S,5);
    real_T  *cp = ssGetDWork(S,6);
    real_T  *mu = ssGetDWork(S,7);
    real_T  *ufs = ssGetDWork(S,8);
    real_T  *FSF_u = ssGetDWork(S,9);
    real_T  *FSF_phi = ssGetDWork(S,10);
    real_T  *FSF_dudphi = ssGetDWork(S,11);
    real_T  *FSF_k = ssGetDWork(S,12);
    real_T  *K_kleff = ssGetDWork(S,13);
    real_T  *K_u = ssGetDWork(S,14);
    real_T  *K_kl = ssGetDWork(S,15);
    real_T  *FSF_length = ssGetDWork(S,16);
    real_T  *K_length = ssGetDWork(S,17);
    
    real_T  *T_source = ssGetDWork(S,18);
    real_T  *Phi_source = ssGetDWork(S,19);
    
    
    int    L;
    int    nnm  = ((int)NN-1);
    int    nnx  = ((int)NX-1);
    int    nnpde  = (int)NPDE;
    
    int     m = 0;
    
//     double  up[NMAX*NPDE];
//     double  U[NMAX*NPDE];
//     double  Ux[NMAX*NPDE];
    
//     double tempUL;
    double tempUxL[NPDE];
//     double tempUR;
    double tempUxR[NPDE];
    
    double cL[NPDE];
    double fL[NPDE];
    double sL[NPDE];
    double cR[NPDE];
    double fR[NPDE];
    double sR[NPDE];
    
    double pL[NPDE];
    double qL[NPDE];
    double pR[NPDE];
    double qR[NPDE];
    
    int n;
    int j;
    int jj;
    
    double denom;
    
    double ps_l;
    double ps_s_l;
    double ps_r;
    double ps_s_r;
    
    double D_v;
    double R_v;
    double rho_w;
    double h_v;
    double cp_v;
    double cp_w;
    
    double k_v;
    double ps;
    double dps_dtheta;
    double rho_v;
    double dH_dtheta;
    double lambda_diff;
    double a_m;
    double a_f;
    
    double kl;
    int ende;
    
    double delta;
    double alpha;
    double beta;
    
    double u_phi;
    double du_dphi;
    
    double u_meas[NNMAX];
    double phi_meas[NNMAX];
    double du_dphi_meas[NNMAX];
    
    double u_kl_meas[NNMAX];
    double kl_meas[NNMAX];
    
    int ip_idx;
    double ip_steps;
    
    
    D_v = 0.00002662;
    R_v = 462.0;
    rho_w = 1000.0;
    h_v = 2445000.0;
    cp_v = 2050.0;
    cp_w = 4180.0;
    
//     printf("nnm : %d  %d\n",nnm,__LINE__);
//     printf("nnx : %d  %d\n",nnx,__LINE__);
//     printf("nnpde : %d  %d\n",nnpde,__LINE__);
    
    // initialize du[n]=0
//     for (n = 0; n < NN; n++)
//         du[n] = 0;
    
    
    /* PDE -> ODE */
    
    /* calculate boundary conditions */
    // phi
        ps_l = 610.0*exp(17.08*(Tleft)/(234.2+(Tleft)));
        ps_s_l = 610.0*exp(17.08*(u[1]-273.15)/(234.2+(u[1]-273.15)));
        ps_r = 610.0*exp(17.08*(Tright)/(234.2+(Tright)));
        ps_s_r = 610.0*exp(17.08*(u[nnm]-273.15)/(234.2+(u[nnm]-273.15)));
        
//         printf("Tleft : %f\n",Tleft);
//         printf("u[1] : %f\n",u[1]);
//         printf("u[0] : %f\n",u[0]);
//         printf("ps_l : %f\n",ps_l);
//         printf("ps_s_l : %f\n",ps_s_l);
//         
//         printf("Tright : %f\n",Tright);
//         printf("u[nnm] : %f\n",u[nnm]);
//         printf("u[nnm-1] : %f\n",u[nnm-1]);
//         printf("ps_r : %f\n",ps_r);
//         printf("ps_s_r : %f\n",ps_s_r);
    
    // phi
    pL[0] = -(0.000000007*alphai)*(u[0]*ps_s_l-phileft*ps_l);
    qL[0] = 1;
    pR[0] = +(0.000000007*alphae)*(u[nnm-1]*ps_s_r-phiright*ps_r);
    qR[0] = 1;
    // teta
    pL[1] = qleft;
    qL[1] = 1;
    pR[1] = -qright;
    qR[1] = 1;
//     printf("pR[0] : %f\n",pR[0]);
//     printf("pR[1] : %f\n",pR[1]);
    
    
    /* for the left point */
    /* calculate UL and dUL/dx */ 
//     tempUL = 1 + (u[1]-u[0]) * ( (xi[0] - XMESH[0]) / (XMESH[1] - XMESH[0]));
    for (j = 0;j<(nnpde);j++)
    {
        tempUxL[j] = (u[j+nnpde]-u[j]) * (1 / (XMESH[1] - XMESH[0]));
//         printf("tempUxL : %f\n",tempUxL[j]);
    }
        
    /* calculate cL, fL, sL */
    for (jj = 0; jj<NNMAX; jj++)
    {
        phi_meas[jj] = FSF_phi[jj];
        u_meas[jj] = ufs[0]*FSF_u[jj];
        du_dphi_meas[jj] = ufs[0]*FSF_dudphi[jj];
    }
    
    if (u[0]>0.0 && u[0]<1.0) 
    {
//         interp1(phi_meas,FSF_length[0],u_meas,u[0],1,u_phi);
//         interp1(phi_meas,FSF_length[0],du_dphi_meas,u[0],1,du_dphi);
        ip_idx = 0;
//         ip_steps = (u_meas[1] - u_meas[0]) / (phi_meas[1] - phi_meas[0]);
//         while(ip_idx < FSF_length[0] && u[0] >= phi_meas[ip_idx])
        while(ip_idx < FSF_length[0] && u[0] <= phi_meas[ip_idx])
        {
            ip_idx++;
        }
//         if(ip_idx != FSF_length[0])
//         {
            ip_steps = (u_meas[ip_idx] - u_meas[ip_idx-1]) / (phi_meas[ip_idx] - phi_meas[ip_idx-1]);
            u_phi = u_meas[ip_idx] - (ip_steps * (phi_meas[ip_idx] - u[0]) );
//         }
        
        ip_idx = 0;
//         ip_steps = (du_dphi_meas[1] - du_dphi_meas[0]) / (phi_meas[1] - phi_meas[0]);
//         while(ip_idx < FSF_length[0] && u[0] >= phi_meas[ip_idx])
        while(ip_idx < FSF_length[0] && u[0] <= phi_meas[ip_idx])
        {
            ip_idx++;
        }
//         if(ip_idx != FSF_length[0])
//         {
            ip_steps = (du_dphi_meas[ip_idx] - du_dphi_meas[ip_idx-1]) / (phi_meas[ip_idx] - phi_meas[ip_idx-1]);
            du_dphi = du_dphi_meas[ip_idx] - (ip_steps * (phi_meas[ip_idx] - u[0]) );
//         }
    } else if (u[0]>=0.99 && u[0] <=1.01) {
        u_phi = ufs[0] + (u[0]-1)/0.01*(ufs[0]*0.10*1000-ufs[0]);
        du_dphi = (ufs[0] + (u[0]-1)/0.01*(ufs[0]*0.10*1000-ufs[0]))*1000000;
    } else if (u[0]>=0.99) {
        u_phi = ufs[0]*100;
        du_dphi = ufs[0]*100000000;
    } else {
        u_phi = 0;
        du_dphi = 1;
    }
    
    k_v = D_v/(mu[0]*R_v*u[1]);
    ps = 610.0*exp(17.08*(u[1]-273.15)/(234.2+(u[1]-273.15)));
    dps_dtheta = 610.0*exp(17.08*(u[1]-273.15)/(234.2+(u[1]-273.15)))*(17.08/(234.2+(u[1]-273.15))-17.08*(u[1]-273.15)/((234.2+(u[1]-273.15))*(234.2+(u[1]-273.15))));
    rho_v = u[0]*ps/(R_v*(u[1]));
    dH_dtheta = rho[0]*cp[0]+cp_w*u_phi+rho_v*cp_v*(ufs[0]-u_phi)/rho_w+u[0]*(ufs[0]-u_phi)*(cp_v*u[1]+h_v)*(dps_dtheta-ps/u[1])/(rho_w*R_v*u[1]);
    
    lambda_diff = -h_v*k_v*u[0]*dps_dtheta;
    
    if (K_kleff[0] > 0)
    {
        for (jj = 0; jj<NNMAX; jj++)
        {
            u_kl_meas[jj] = ufs[0]*K_u[jj];
            kl_meas[jj] = K_kleff[0]*K_kl[jj];
        }
        
        kl = 0.0;
        
        if (u[0]>0.0 && u[0]<1.0) 
        {
            ip_idx = 0;
//             printf("u_kl_meas[ip_idx] : %f\n",u_kl_meas[ip_idx]);
            
            while(ip_idx < K_length[0] && u_phi >= u_kl_meas[ip_idx])
            {
                ip_idx++;
//                 printf("ip_idx : %d\n",ip_idx);
            }
            ip_steps = (kl_meas[ip_idx] - kl_meas[ip_idx-1]) / (u_kl_meas[ip_idx] - u_kl_meas[ip_idx-1]);
            kl = kl_meas[ip_idx] - (ip_steps * (u_kl_meas[ip_idx] - u_phi) );
//             printf("u_kl_meas[ip_idx] : %f\n",u_kl_meas[ip_idx]);
//             printf("kl : %f\n",kl);
        } else if (u[0]>=0.99) {
            ende = (int)K_length[0];
            kl = kl_meas[ende];
        } else {
            kl = 0;
        }
//         printf("kl : %f\n",kl);
//         printf("u_phi : %f\n",u_phi);
        
        a_f = -kl*(-rho_w*u[1]*R_v/u[0]);
    } else {
        a_f = 0.0;       // no fluid transport!
    }
    
    a_m = k_v*ps;
    
    delta = u[0]*dps_dtheta/ps;
    alpha = (lambda[0]-lambda_diff);
    beta = h_v;
    
    // phi
    cL[0] = du_dphi;
    fL[0] = (a_m+a_f)*tempUxL[0] + (a_m*delta)*tempUxL[1];
    sL[0] = Phi_source[0]/(2*zxmp1[0]);
//     sL[0] = 0.0/(2*zxmp1[0]);
    if (Phidactive == 1)
    {
        sL[0] = sL[0] + Phisource/(2*zxmp1[0]);
    }
    // teta
    cL[1] = dH_dtheta;
    fL[1] = alpha*tempUxL[1] + (a_m*beta)*tempUxL[0];
    sL[1] = T_source[0]/(2*zxmp1[0]);
    if (Tdactive == 1)
    {
        sL[1] = sL[1] + Tsource/(2*zxmp1[0]);
    }
    if (Phidactive == 1)
    {
        sL[1] = sL[1] - Phisource*2450*1000/(2*zxmp1[0]);
    }
    
//     printf("cL : %f\n",cL);
//     printf("fL : %f\n",fL);
//     printf("sL : %f\n",sL);
    
    for (j = 0;j<(nnpde);j++)
    {
        if (qL[j] == 0)
        {
            du[j] = pL[j];
        } else {
//             du[j] = pL[j] + (qL[j] / pow(XMESH[0],m)) * (xim[0] * fL[j] + zxmp1[0] * sL[j]);
            du[j] = pL[j] + (qL[j]) * (xim[0] * fL[j] + zxmp1[0] * sL[j]);
//             printf("du[1] : %f\n",du[1]);
//             printf("xzmp1[0] : %f\n",xzmp1[0]);
//             printf("zxmp1[0] : %f\n",zxmp1[0]);
//             printf("xim[0] : %f\n",xim[0]);
//             denom = (qL[j] / pow(XMESH[0],m)) * zxmp1[0] * cL[j];
            denom = (qL[j]) * zxmp1[0] * cL[j];
//             printf("denom : %f\n",denom);
            if (denom != 0) {
                du[j] = du[j] / denom;
//                 printf("du[1] : %f\n",du[1]);
            }
        }
    }
    
    
    /* for the interior points */
    for (L = 1; L<(nnx); L++)   // numbers of points
    {
        /* calculate U and dU/dx */ 
//         tempUR = 1 + (u[L+1]-u[L]) * ( (xi[L] - XMESH[L]) / (XMESH[L+1] - XMESH[L]));
        for (j = 0;j<(nnpde);j++)
        {
            tempUxR[j] = (u[nnpde*(L+1)+j]-u[nnpde*L+j]) * (1 / (XMESH[L+1] - XMESH[L]));
        }
        
        /* calculate cR, fR, sR */
        for (jj = 0; jj<NNMAX; jj++)
        {
            phi_meas[jj] = FSF_phi[L*NNMAX+jj];
            u_meas[jj] = ufs[L]*FSF_u[L*NNMAX+jj];
            du_dphi_meas[jj] = ufs[L]*FSF_dudphi[L*NNMAX+jj];
        }
        
        u_phi = 0;
        du_dphi = 1;
        
//         printf("L : %d\n",L);
                
        if (u[L*nnpde]>0.0 && u[L*nnpde]<1.0) 
        {
    //         interp1(phi_meas,FSF_length[0],u_meas,u[0],1,u_phi);
    //         interp1(phi_meas,FSF_length[0],du_dphi_meas,u[0],1,du_dphi);
            ip_idx = 0;
    //         ip_steps = (u_meas[1] - u_meas[0]) / (phi_meas[1] - phi_meas[0]);
//             while(ip_idx < FSF_length[L] && u[L*nnpde] >= phi_meas[ip_idx])
            while(ip_idx < FSF_length[L] && u[L*nnpde] <= phi_meas[ip_idx])
            {
                ip_idx++;
            }
            
//             printf("u : %f\n",u[L*nnpde]);
//             printf("ip_idx : %d\n",ip_idx);
            
//             if(ip_idx != FSF_length[L])
//             {
                ip_steps = (u_meas[ip_idx] - u_meas[ip_idx-1]) / (phi_meas[ip_idx] - phi_meas[ip_idx-1]);
                u_phi = u_meas[ip_idx] - (ip_steps * (phi_meas[ip_idx] - u[L*nnpde]) );
//             }

//             ip_idx = 0;
    //         ip_steps = (du_dphi_meas[1] - du_dphi_meas[0]) / (phi_meas[1] - phi_meas[0]);
//             while(ip_idx < FSF_length[L] && u[L*nnpde] >= phi_meas[ip_idx])
//             while(ip_idx < FSF_length[L] && u[L*nnpde] <= phi_meas[ip_idx])
//             {
//                 ip_idx++;
//             }
//             printf("ip_idx : %d\n",ip_idx);
//             if(ip_idx != FSF_length[L])
//             {
                ip_steps = (du_dphi_meas[ip_idx] - du_dphi_meas[ip_idx-1]) / (phi_meas[ip_idx] - phi_meas[ip_idx-1]);
                du_dphi = du_dphi_meas[ip_idx] - (ip_steps * (phi_meas[ip_idx] - u[L*nnpde]) );
                
//                 printf("ip_idx : %d\n",ip_idx);
//                 printf("du_dphi_meas[ip_idx] : %f\n",du_dphi_meas[ip_idx]);
//                 printf("phi_meas[ip_idx] : %f\n",phi_meas[ip_idx]);
//                 printf("u[L*nnpde]) : %f\n",u[L*nnpde]));
//             }

        } else if (u[L*nnpde]>=0.99 && u[L*nnpde] <=1.01) {
            u_phi = ufs[L] + (u[L*nnpde]-1)/0.01*(ufs[L]*0.10*1000-ufs[L]);
            du_dphi = (ufs[L] + (u[L*nnpde]-1)/0.01*(ufs[L]*0.10*1000-ufs[L]))*10000000;
//             printf("du_dphi_meas[L] 1.01 : %f %d\n",du_dphi_meas[L],L);
        } else if (u[L*nnpde]>=0.99) {
            u_phi = ufs[L]*100;
            du_dphi = ufs[L]*100000000;
//             printf("du_dphi_meas[L] : %f %d\n",du_dphi_meas[L],L);
        } else {
            u_phi = 0;
            du_dphi = 1;
        }
        
//         printf("L : %d\n",L);
//         printf("phi_meas[ip_idx] : %f\n",phi_meas[ip_idx]);
//         printf("phi_meas[ip_idx-1] : %f\n",phi_meas[ip_idx-1]);
//         printf("du_dphi_meas[ip_idx] : %f\n",du_dphi_meas[ip_idx]);
//         printf("du_dphi_meas[ip_idx-1] : %f\n",du_dphi_meas[ip_idx-1]);
//         printf("ip_steps : %f\n",ip_steps);
//         printf("u[L*nnpde] : %f\n",u[L*nnpde]);
//         printf("u_phi : %f\n",u_phi);
//         printf("du_dphi : %f\n",du_dphi);
        
        k_v = D_v/(mu[L]*R_v*u[L*nnpde+1]);
        ps = 610.0*exp(17.08*(u[L*nnpde+1]-273.15)/(234.2+(u[L*nnpde+1]-273.15)));
        dps_dtheta = 610.0*exp(17.08*(u[L*nnpde+1]-273.15)/(234.2+(u[L*nnpde+1]-273.15)))*(17.08/(234.2+(u[L*nnpde+1]-273.15))-17.08*(u[L*nnpde+1]-273.15)/((234.2+(u[L*nnpde+1]-273.15))*(234.2+(u[L*nnpde+1]-273.15))));
        rho_v = u[L*nnpde]*ps/(R_v*(u[L*nnpde+1]));
        dH_dtheta = rho[L]*cp[L]+cp_w*u_phi+rho_v*cp_v*(ufs[L]-u_phi)/rho_w+u[L*nnpde]*(ufs[L]-u_phi)*(cp_v*u[L*nnpde+1]+h_v)*(dps_dtheta-ps/u[L*nnpde+1])/(rho_w*R_v*u[L*nnpde+1]);

        lambda_diff = -h_v*k_v*u[L*nnpde]*dps_dtheta;
        
        if (K_kleff[L] > 0)
        {
            for (jj = 0; jj<NNMAX; jj++)
            {
                u_kl_meas[jj] = ufs[L]*K_u[jj];
                kl_meas[jj] = K_kleff[L]*K_kl[jj];
            }

            kl = 0.0;

            if (u[L*nnpde]>0.0 && u[L*nnpde]<1.0) 
            {
                ip_idx = 0;
    //             printf("u_kl_meas[ip_idx] : %f\n",u_kl_meas[ip_idx]);

                while(ip_idx < K_length[L] && u_phi >= u_kl_meas[ip_idx])
                {
                    ip_idx++;
    //                 printf("ip_idx : %d\n",ip_idx);
                }
                ip_steps = (kl_meas[ip_idx] - kl_meas[ip_idx-1]) / (u_kl_meas[ip_idx] - u_kl_meas[ip_idx-1]);
                kl = kl_meas[ip_idx] - (ip_steps * (u_kl_meas[ip_idx] - u_phi) );
    //             printf("u_kl_meas[ip_idx] : %f\n",u_kl_meas[ip_idx]);
            } else if (u[L*nnpde]>=0.99) {
                ende = (int)K_length[L];
                kl = kl_meas[ende];
            } else {
                kl = 0;
            }
            
//             printf("kl : %f\n",kl);
//             printf("u_phi : %f\n",u_phi);

            a_f = -kl*(-rho_w*u[L*nnpde+1]*R_v/u[L*nnpde]);
//             printf("a_m : %.20f\n",a_m);
//             printf("k_v : %.20f\n",k_v);
//             printf("ps : %.20f\n",ps);
//             printf("-kl*(-rho_w*u[L*nnpde+1]*R_v/u[L*nnpde]) : %.20f\n",-kl*(-rho_w*u[L*nnpde+1]*R_v/u[L*nnpde]));
//             printf("k_v*ps : %.20f\n",k_v*ps);
        } else {
            a_f = 0.0;       // no fluid transport!
        }
        
        a_m = k_v*ps; 

        delta = u[L*nnpde]*dps_dtheta/ps;
        alpha = (lambda[L]-lambda_diff);
//         printf("alpha : %f\n",alpha);
        beta = h_v;

        // phi
        cR[0] = du_dphi;
        fR[0] = (a_m+a_f)*tempUxR[0] + (a_m*delta)*tempUxR[1];
        sR[0] = Phi_source[L]/(2*zxmp1[L]);
//         printf("Phidactive : %d\n",Phidactive);
        if (Phidactive == (L+1))
        {
            sR[0] = sR[0] + Phisource/(2*zxmp1[L]);
//             printf("Phisource : %f\n",Phisource);
        }
        // teta
        cR[1] = dH_dtheta;
        fR[1] = alpha*tempUxR[1] + (a_m*beta)*tempUxR[0];
        sR[1] = T_source[L]/(2*zxmp1[L]);
        if (Tdactive == (L+1))
        {
            sR[1] = sR[1] + Tsource/(2*zxmp1[L]);
//             printf("Tsource %d : %f\n",L, Tsource);
        }
        if (Phidactive == (L+1))
        {
            sR[1] = sR[1] - Phisource*2450*1000/(2*zxmp1[L]);
        }
//         sR[1] = 0.0/(2*zxmp1[L]);
        
//         printf("L : %d\n",L);
//         printf("dH_dtheta : %f\n",dH_dtheta);
//         printf("rho[L] : %f\n",rho[L]);
//         printf("cp[L] : %f\n",cp[L]);
//         printf("ufs[L] : %f\n",ufs[L]);
//         printf("u_phi : %f\n",u_phi);
//         printf("u[L*nnpde] : %f\n",u[L*nnpde]);
//         printf("u[L*nnpde+1] : %f\n",u[L*nnpde+1]);
//         printf("dps_dtheta : %f\n",dps_dtheta);
//         printf("ps : %f\n",ps);
//         printf("rho_v : %f\n",rho_v);
//         printf("cp_v : %f\n",cp_v);
//         printf("rho_w : %f\n",rho_w);
//         printf("R_v : %f\n",R_v);
//         printf("h_v : %f\n",h_v);
        
//         printf("tempUxR[0] : %f\n",tempUxR[0]);
//         printf("tempUxR[1] : %f\n",tempUxR[1]);
        
//         printf("L : %d\n",L);
        for (j = 0;j<(nnpde);j++)
        {
            du[nnpde*L+j] = (xim[L] * fR[j] - xim[L-1] * fL[j]) + (zxmp1[L] * sR[j] + xzmp1[L] * sL[j]);
            denom = zxmp1[L] * cR[j] + xzmp1[L] * cL[j];
            if (denom != 0) {
                du[nnpde*L+j] = du[nnpde*L+j] / denom;
            }
//             printf("xzmp1[L] : %f\n",xzmp1[L]);
//             printf("zxmp1[L] : %f\n",zxmp1[L]);
//             printf("xim[L] : %f\n",xim[L]);
//             printf("j : %d\n",j);
//             printf("xim[L] : %f\n",xim[L]);
//             printf("xim[L-1] : %f\n",xim[L-1]);
//             printf("fR[j] : %f\n",fR[j]);
//             printf("fL[j] : %f\n",fL[j]);
//             printf("sR[j] : %f\n",sR[j]);
//             printf("sL[j] : %f\n",sL[j]);
//             printf("du[nnpde*L+j] : %f\n",du[nnpde*L+j]);
        }
        
        for (j = 0;j<(nnpde);j++)
        {
            cL[j] = cR[j];
            fL[j] = fR[j];
            sL[j] = sR[j];
        }
        
    }
    
    
    /* for the right point */
    for (j = 0;j<(nnpde);j++) 
    {
        if (qR[j] == 0)
        {
            du[nnpde*nnx+j] = pR[j];
        } else {
//             du[nnpde*nnx+j] = pR[j] + (qR[j] / pow(XMESH[nnx],m)) * (xim[nnx-1] * fL[j] - xzmp1[nnx] * sL[j]);
            du[nnpde*nnx+j] = pR[j] + (qR[j]) * (xim[nnx-1] * fL[j] - xzmp1[nnx] * sL[j]);
//             denom = -(qR[j] / pow(XMESH[nnx],m)) * xzmp1[nnx] * cL[j];
            denom = -(qR[j]) * xzmp1[nnx] * cL[j];
            if (denom != 0.0) {
                du[nnpde*nnx+j] = du[nnpde*nnx+j] / denom;
            }
        }
    }
    
  }
#endif /* MDL_DERIVATIVES */



/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
}


/*======================================================*
 * See sfuntmpl_doc.c for the optional S-function methods *
 *======================================================*/

/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
