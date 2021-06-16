/***********************************************************************
 *  M O D E L    O R    F U N C T I O N
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * multinode termal wall model according to the Beuken Modell
 *
 * Author list
 *  Thomas Wenzel -> tw
 *  Bernd Hafner -> hf
 *  Christian Winteler -> wic
 *  Arnold Wohlfeil -> aw
 *
 * Version  Author  Changes                                     Date
 * 0.9.0    hf      created                                     19jul1998
 *          tw      changed to level 2 s-function               11jun1999 
 *                  dyn. sized vector: power per node
 *                  conductivity + capacity in RWORK   
 * 1.0.0    tw      created new, following matwall.m            06jul1999
 * 1.0.1    hf      corrected power in active layers            16jul99
 * 6.1.0    hf      Power to active layer distributed between   03apr2014
 *                  nodes if active node is not exactly on a node.
 *                  Added define for iwork.
 * 6.1.1    wic     If active node is not exactly on a existing 07apr14
 *                  node, an additional node is created.
 *                  corrected if-clause in derivatives to 
 *                  prevent strange behaviour if active layer
 *                  is in first or last layer.
 * 6.1.2    hf      Corrected parameter check, layerwall can    11apr2014
 *                  also be used for walls, not only for floors
 *                  so active layers must be arranged in in-
 *                  creasing depth (not from top to bottom)
 * 6.1.3    aw      line 355, NUMACTIVE changed to numactive    11nov2014
 *
 ***********************************************************************
 * This file is part of the CARNOT Blockset.
 * Copyright (c) 1998-2015, Solar-Institute Juelich of the FH Aachen.
 * Additional Copyright for this file see list auf authors.
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are 
 * met:
 * 1. Redistributions of source code must retain the above copyright notice, 
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright 
 *    notice, this list of conditions and the following disclaimer in the 
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its 
 *    contributors may be used to endorse or promote products derived from 
 *    this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
 * THE POSSIBILITY OF SUCH DAMAGE.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  D E S C R I P T I O N
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The wall is devided into "NODES" nodes according to the Beuken-Model
 *
 *        ------------ W A L L -----------
 *        |    Layer 1    |   Layer 2    | 
 *        |               |              | 
 *        |   --------    |   --------   |
 *  Q1 -> T1 -|  R1  |-   T2 -|  R2  |-  T3  <- Q3
 *        |   --------    |   --------   |
 *        |               |              | 
 *      ______         ______          ______
 *      __C1__         __C2__          __C3__
 *        |               |              | 
 *        |               |              | 
 *
 * Energy-balance for every node with the differential equation:
 *
 * rho*cp*d_node * dT/dt =
 *        q_outside                             % only first node
 *      + q_inside                              % only last node
 *      + cond/d_node^2 * (Tnextnode - Tnode)   % not last node
 *      + cond/d_node^2 * (Tlastnode - Tnode)   % not first node
 *      + qdot_heating
 *
 *  symbol      used for                                        unit
 *	cond        effective axial heat conduction                 W/(m*K)
 *  cp          heat capacity of wall                           J/(kg*K)
 *  d_node      distance between two nodes                      m
 *  rho         density of layer material                      kg/m³
 *  T           temperature                                     K
 *  t           time                                            s
 *  q_          power per surface (positive for energy gain)    W/(m^2)
 *
 *
 *  Matrix representation
 *      R = d / cond
 *      C = d*c
 *
 *  |dT1/dt|   |-1/(R1*C1)     1/(R1*C1)            0    |   |T1|   |1/C1  0  |   |Q1|
 *  |dT2/dt| = | 1/(R1*C2)  (R1+R2)/(R1*R2*C2)    1/C2   | * |T2| + | 0    0  | * |Q1|
 *  |dT3/dt|   |   0           1/(R2*C3)        1/(R2*C3)|   |T3|   | 0   1/C3|   |Q1|
 *         
 * structure of u (input vector)
 *  see defines below
 *
 * Literature: 
 * Feist, W.: Thermische Gebaeudesimulation, Dissertation Uni Kassel, 
 *              Müller 2004
 * Wimmer, A.: Thermoaktive Bauteilsysteme, ein neuer simulationstechnischer
 *              Berechnungsansatz, Dissertation Uni Kassel, 2004
 */

#define S_FUNCTION_NAME     layerwall_60
#define S_FUNCTION_LEVEL    2

#include <math.h>
#include "simstruc.h"

/*
 *   Defines for easy access to the parameters
 */

#define NODES       *mxGetPr(ssGetSFcnParam(S,0)) /* number of nodes */
#define TAU         *mxGetPr(ssGetSFcnParam(S,1)) /* time-constant */
#define TINI        *mxGetPr(ssGetSFcnParam(S,2)) /* initial temperature [°C]  */
#define S_DNODE              ssGetSFcnParam(S,3)  /* thickness of node in m */
#define S_COND               ssGetSFcnParam(S,4)  /* conductivity [W/(m*K)] */
#define S_CWALL              ssGetSFcnParam(S,5)  /* capacity [J/(kg*K)] */
#define S_RHO                ssGetSFcnParam(S,6)  /* density [kg/m^3] */
#define NUMACTIVE   *mxGetPr(ssGetSFcnParam(S,7)) /* number of active layers */
#define S_DEPTH              ssGetSFcnParam(S,8)  /* depth of active layers [m] */
#define NPARAMS                               9

#define NDNODE              mxGetM(S_DNODE)      
#define NCOND               mxGetM(S_COND)
#define NCWALL              mxGetM(S_CWALL)
#define NRHO                mxGetM(S_RHO)
#define NDEPTH              mxGetM(S_DEPTH)

#define Q_OUTSIDE           (*u0[0])    /* power per surface outside node */
#define POWER_PER_NODE(n)   (*u1[n])    /* power per node */
#define Q_INSIDE            (*u2[0])    /* power per surface inside node */

#define DIMMAT              iwork[0]

#define MAXNODES                    20      /* maximum layer */
#define MAX_LAYERS                  10      /* Layers per Layer */
#define MAX_L                       200     /* has to be atleast MAXNODES * MAX_LAYERS !!! */
#define DMAX_ACTIVE_LAYER_TO_NODE   0.00999 /* maximum distance from active layer to node in m *
                                             * if bigger, an additional node is added          */



#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
  /* Function: mdlCheckParameters =============================================
   * Abstract:
   *    Validate our parameters to verify they are okay.
   */
  static void mdlCheckParameters(SimStruct *S)
  {
      /* Check 1st parameter: number of layers */
      {
          if (NODES <= 0.0) {
              ssSetErrorStatus(S,"Error in wall: number of layers must be > 0");
              return;
          }
      }
      
      {
        if (NODES > MAXNODES){
           printf("WARNING: Wall cannot have more than %i nodes. \n", MAXNODES);
           printf("         Disregarding all other nodes. \n");
           NODES = MAXNODES;
        }
      }
      
      /* Check 2nd parameter: time-constant */
      {
          if (TAU <= 0.0) {
              ssSetErrorStatus(S,"Error in wall: time-constant must be > 0 s");
              return;
          }
      }
      /* Check 4th parameter: thickness of layers in m */
      {
          int i;
          for (i = 0; i < NODES; i++)
            if (mxGetPr(S_DNODE)[i] <= 0.0) {
              ssSetErrorStatus(S,"Error in wall: thickness of layers must be > 0");
              return;
          }
      }
      /* Check 5th parameter: conductivity [W/(m*K)] */
      {
          int i;
          for (i = 0; i < NODES; i++)
            if (mxGetPr(S_COND)[i] < 0.0) {
              ssSetErrorStatus(S,"Error in wall: heat conductivity must be >= 0");
              return;
          }
      }
      /* Check 6th parameter: capacity [J/(kg*K)] */
      {
          int i;
          for (i = 0; i < NODES; i++)
            if (mxGetPr(S_CWALL)[i] <= 0.0) {
              ssSetErrorStatus(S,"Error in wall: heat capacity must be > 0");
              return;
          }
      }
      /* Check 7th parameter: density [kg/m^3] */
      {
          int i;
          for (i = 0; i < NODES; i++)
            if (mxGetPr(S_RHO)[i] <= 0.0) {
              ssSetErrorStatus(S,"Error in wall: density must be > 0");
              return;
          }
      }
      /* Check 8th parameter: number of active layers */
      {
        if (NUMACTIVE > MAXNODES){
              ssSetErrorStatus(S,"Error in wall: number of acitve nodes exceeds maximum number of nodes."
                " Recompile layerwall.c with higher number for MAXNODES.");
              return;
        }
      }
      /* Check 9th parameter: depth of active layers [m] */
      {
          int i;
          double sum = 0.0;

          for (i = 0; i < NODES; i++) {
            sum += mxGetPr(S_DNODE)[i];
          }
          i = (int_T)NUMACTIVE-1;
          if (sum < mxGetPr(S_DEPTH)[i]) {
              ssSetErrorStatus(S,"Error in wall: position of last active layer is outside of wall. "
                "Reduce depth of active layer or increase dimension of layers.");
              return;
          }
          for (i=1;i<NUMACTIVE;i++)
              if (mxGetPr(S_DEPTH)[i] < mxGetPr(S_DEPTH)[i-1]){
              ssSetErrorStatus(S,"Error in wall: depth of active layers must be monotonicaly increasing.");
              return;
          }
      }
  }
#endif /* MDL_CHECK_PARAMETERS */
 



/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *   Setup sizes of the various vectors.
 */
static void mdlInitializeSizes(SimStruct *S)
{
    int numactive = (int)(NUMACTIVE+0.5);

    ssSetNumSFcnParams(S, NPARAMS);
#if defined(MATLAB_MEX_FILE)
    if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S)) {
        mdlCheckParameters(S);
        if (ssGetErrorStatus(S) != NULL) {
            return;
        }
    } else {
        return; /* Parameter mismatch will be reported by Simulink */
    }
#endif

    ssSetNumContStates(    S, (int)NODES*MAX_LAYERS);  /* number of continuous states */
    /*    
    * printf("NODES : %d                 l. %d\n",(int) NODES,__LINE__);
    * printf("MAX_LAYERS %d              l. %d\n",MAX_LAYERS,__LINE__);
    * printf("MAX_L : %d                 l. %d\n",MAX_L,__LINE__);
    */
    ssSetNumDiscStates(    S, 0);      /* number of discrete states */

    if (!ssSetNumInputPorts(S, 3)) return;
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 0, 0);
    ssSetInputPortWidth(S, 1, numactive);
    ssSetInputPortDirectFeedThrough(S, 1, 0);
    ssSetInputPortWidth(S, 2, 1);
    ssSetInputPortDirectFeedThrough(S, 2, 0);

    if (!ssSetNumOutputPorts(S,3)) return;
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, numactive);
    ssSetOutputPortWidth(S, 2, 1);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, MAX_L);

    ssSetNumIWork(S, numactive+1);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* Take care when specifying exception free code - see sfuntmpl.doc */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}


/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy that we inherit our sample time from the driving block.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}


#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    Initialize both states to one
 */
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *x0   = ssGetContStates(S);
    real_T t0    = TINI;
    real_T *rwork = ssGetRWork(S);
    int    *iwork = ssGetIWork(S);

    real_T   *dnode    = mxGetPr(S_DNODE);
    real_T   *cond     = mxGetPr(S_COND);
    real_T   *cwall    = mxGetPr(S_CWALL);
    real_T   *rho      = mxGetPr(S_RHO);
    real_T   *depth    = mxGetPr(S_DEPTH);

    int_T    nodes = (int)NODES;
    int_T    n;
    int_T    numactive = (int)(NUMACTIVE+0.5);

    int_T numlayer = nodes,
        DimLambda,
        DimCap,
        StartCell,
        EndCell,
        j,
        i;
          
    real_T delx,
           celldepth[MAXNODES],
           m[MAX_L],
           k[MAX_L],
           c[MAX_L],
           MatLambda1[MAX_L],
           MatLambdaO[MAX_L],
           MatLambdaU[MAX_L],
           MatLambdaM[MAX_L],
           MatCap1[MAX_L],
           MatCap[MAX_L];            /* Matrix, nur auf Hauptdiagonaler besetzt, daher Vektor */

	/* no active layer if depth of active layer < 0*/
    if (depth[0]<0)
        numactive=0;
    
    /* layer capacities, layer resistances und layer conductivities:*/
    DIMMAT = 0;
    for (j = 0;j<numlayer;j++)
    {
        /* m is the number of nodes per layer, equation from Feist */
        m[j] = ceil(sqrt(rho[j] * cwall[j] / (2.0 * cond[j] * TAU)) * dnode[j]);
        if (m[j] < 1)
           m[j] = 1;
        else if (m[j] > MAX_LAYERS)
           m[j] = MAX_LAYERS;
        
        c[j] = rho[j]*cwall[j]*dnode[j]/(2.0*m[j]); // thermal capacity of node (half at each surface) C = d*rho*cp
        k[j] = cond[j]*m[j]/dnode[j];               // heat transfer by conduction
        DIMMAT = DIMMAT + (int)(m[j]+0.01);
    }
    DIMMAT = DIMMAT+1;
     
    /*%%%%% depth of cells from upper surface*/
    celldepth[0] = 0.0;
    n=1;
    for (i=0;i<numlayer;i++)
    {
         for (j=1;j<=m[i];j++)
         {
             celldepth[n]=celldepth[n-1]+dnode[i]/m[i];
             n++;
         }
    }
     
    /*%%%%% create help matrices "MatKap" and "MatLambda":*/
    DimCap    = DIMMAT+1;
    DimLambda = DimCap;

    /*%%%%% occupy elements of "MatCap"*/
    for (j=0; j<DimCap; j++)
    {
        MatCap1[j] = 0.0;
        MatCap[j] = 0.0;
        MatLambda1[j] = 0.0;
    }  
 
    /*%%%%% set initial MatCap for all layers */
    StartCell = 0;
    for (j=0; j<numlayer; j++)
    {
        if (j == 0)
            EndCell = StartCell + (int_T)m[j];
        else
            EndCell = StartCell + (int_T)m[j];
     
        for (i = StartCell; i<EndCell; i++)
        {
            MatCap1[i] = c[j];
        }

        StartCell = i;
    }

     /*%%%%% occupy elements of "MatLambda":*/
     StartCell = 1;
     MatLambda1[0] = 0.0;
     for (j=0; j<numlayer; j++)
     {
        EndCell = StartCell + (int_T)m[j];
        for (i = StartCell; i<EndCell; i++)
           MatLambda1[i] = k[j];

        StartCell = i;
     }
     MatLambda1[i+1] = 0.0;

     /*%%%%% which of the layers are active layers? store in iwork[1, ...] 
			 insert new nodes if necessary */
     j = 0;
     for (i=0;i<numactive;i++)
     {
        while (depth[i]>celldepth[j])
          j++;

        iwork[i+1] = j;
        
         if (depth[i] == celldepth[j])
             continue;
        
        if (depth[i]-celldepth[j-1] > DMAX_ACTIVE_LAYER_TO_NODE)
        {
            if (celldepth[j]-depth[i] > DMAX_ACTIVE_LAYER_TO_NODE)
            {
                /*%%%%% Insert a new node for the active layer at curent position*/
                
                /* shift elements after current position to the right */
                for (n=DIMMAT+1;n>=j;n--)
                {
                    celldepth[n+1] = celldepth[n];
                    MatCap1[n+1] = MatCap1[n];
                    MatLambda1[n+1] = MatLambda1[n];
                }
                DIMMAT++;
                
                celldepth[j] = depth[i]; /* insert element at current position*/
                delx = (depth[i]-celldepth[j-1])/(celldepth[j+1]-celldepth[j-1]);
                MatCap1[j] = MatCap1[j-1]*(1.0-delx); /*distribute capacity on two nodes*/
                MatCap1[j-1] = MatCap1[j-1]-MatCap1[j];
                MatLambda1[j] = 1.0/(delx/MatLambda1[j+1]); /* distribute lambda on two nodes */
                MatLambda1[j+1] = 1.0/(1.0/MatLambda1[j+1] - 1.0/MatLambda1[j]);
            }
        }
        else
            iwork[i+1]=j-1;
     }

     /*%%%%% reset dimensions of MatCap and MatLambda:*/
     DimCap    = DIMMAT+1;
     DimLambda = DimCap;
     
     /*%%%%% create final MatCap by summing up capacities of inner nodes*/
     /*%%%%% i.e. except first and last element*/
     MatCap[0] = MatCap1[0];
     for (j=1;j<DimCap-1;j++)
        MatCap[j] = MatCap1[j-1]+MatCap1[j];
     
    /* coefficient matrix "MatA" (consists of 4 single matrices),  
     * (differential equations, see W.Feist, S.136, [4-37]):        
     *  Matrix representation
     *      R = d / cond
     *      C = d*c
     *
     *  |dT1/dt|   |-1/(R1*C1)     1/(R1*C1)            0    |   |T1|   |1/C1  0  |   |Q1|
     *  |dT2/dt| = | 1/(R1*C2)  (R1+R2)/(R1*R2*C2)    1/C2   | * |T2| + | 0    0  | * |Q1|
     *  |dT3/dt|   |   0           1/(R2*C3)        1/(R2*C3)|   |T3|   | 0   1/C3|   |Q1|
     */
    for (j=0; j<DimLambda-1; j++)
    {
        MatLambdaM[j] = -(MatLambda1[j]+MatLambda1[j+1])/(MatCap[j]);   // main diagonal
        MatLambdaO[j] = MatLambda1[j+1] / MatCap[j];                    // upper diagonal
        MatLambdaU[j] = MatLambda1[j] / MatCap[j];                      // lower digagonal
    }

    /*
     Im rwork werden die Matrizen A und B zeilenweise gespeichert. Da A und B 
     beide schwach besetzt sind, brauchen nur die Werte auf den besetzten 
     Diagonalen gespeichert werden. Hierbei werden A und B zeilenweise 
     hintereinander in rwork geschrieben.
         [ A11 A12                 ]         [ B11                     ]
         [ A21 A22 A23             ]         [     B22                 ]
         [     A32 A33 A34         ]         [         B33             ]
     A = [         A43 A44 A45     ]     B = [             B44         ]
         [             A54  .   .  ]         [                  .      ]
         [                  .   .  ]         [                      .  ]
         [                      .  ]         [                         ]


     rwork = [ A11 A12 B11   A21 A22 A23 B22   A32 A33 A34 B33   A43 A44 A45 B44 ... ]   
    */

    i = 0;
    for (j=0;j<DIMMAT;j++)
    {
        if (j>0)
           rwork[i++] = MatLambdaU[j];
        rwork[i++] = MatLambdaM[j];
        if (j<DIMMAT-1)
           rwork[i++] = MatLambdaO[j];
        rwork[i++] = 1./MatCap[j];      /* MatB */
    }

    for (n = 0; n <NODES*MAX_LAYERS; n++)
        x0[n] = t0;             /* state-vector is initialized with TINI */             
}



/* Function: mdlOutputs =======================================================
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T   *y0   = ssGetOutputPortRealSignal(S,0);
    real_T   *y1   = ssGetOutputPortRealSignal(S,1);
    real_T   *y2   = ssGetOutputPortRealSignal(S,2);
    real_T   *Tn    = ssGetContStates(S);

    real_T  *rwork = ssGetRWork(S);    
    int     *iwork = ssGetIWork(S);
   
    int_T nodes = (int)NODES;
    int_T numactive = (int)(NUMACTIVE+0.5);
    int_T n;

    y0[0] = Tn[0];                   	/* temperature first node */
    for (n = 0; n < numactive; n++) 	/* all node temperatures for active layers */
        y1[n] = Tn[iwork[n+1]];  		/* = x[ pos[n] ];   */        
    y2[0] = Tn[DIMMAT-1];            	/* temperature last node */
   
}



#define MDL_DERIVATIVES 
/* Function: mdlDerivatives =================================================
 * Abstract:
 *      xdot = Ax + Bu
 */
static void mdlDerivatives(SimStruct *S)
{
    real_T          *dTdt = ssGetdX(S);
    real_T           *Tn  = ssGetContStates(S);
    InputRealPtrsType u0  = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType u1  = ssGetInputPortRealSignalPtrs(S,1);
    InputRealPtrsType u2  = ssGetInputPortRealSignalPtrs(S,2);

    real_T  *rwork = ssGetRWork(S);
    int_T   *iwork = ssGetIWork(S);
    
    real_T qinside   = Q_INSIDE;
    real_T qoutside  = Q_OUTSIDE;
    int_T  nodes     = (int_T)NODES;
    int_T  numactive = (int_T)(NUMACTIVE+0.5);
    int_T  n,i,j;

    /* loop over all layers and nodes */
    i = 0;
    j = 1;

    for (n=0; n<DIMMAT; n++)   
    {
        dTdt[n] = 0;
       
        if (n>0)
            dTdt[n] += rwork[i++] * Tn[n-1];    /* untere Diagonale von A */
          
        dTdt[n] += rwork[i++] * Tn[n];          /* Hauptdiagonale von A   */
       
        if (n < DIMMAT-1)
            dTdt[n] += rwork[i++] * Tn[n+1];    /* obere Diagonale von A  */

        if (n == 0)                             /* Hauptdiagonale von B   */
            dTdt[n] += rwork[i] * qoutside;     
        else if (n == DIMMAT-1)
            dTdt[n] += rwork[i] * qinside;
        
        if (j <= numactive) 
        {
            if (iwork[j]==n) 
            {
                dTdt[n] += (rwork[i]*POWER_PER_NODE(j-1));
                j++;
            } 
        } 
        i++;                                    /* next Matrix Index */
    } /* end fpr n */

    for (n = DIMMAT; n < NODES*MAX_LAYERS; n++)   
        dTdt[n] = 0;
}




/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
}


#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif