  #include "mex.h"
  #include <stdio.h>
  #include <math.h>
  #include <stdlib.h>


  /***********************************/
  /*********global variables*********/
  /**********************************/
  /*********time discretization*******/

   double dt=0.1;

  /*********space discretization********/
   double  h=1.0;

   /********parameter epsilon in the regularized versions of H and delta********/
   double epsilon=1.0;

   /*******coefficients of the fidelity terms****************/
   double a1=1.0;
   double a2=1.0;

   /**********nombre d'it�rations maximal*********************/
   int Itermax=400;

   /**********coefficient of the length term******************/
   double eta=0.1*255.0*255.0;


   void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){


   int i,j,k,m,n;
   int mrows,ncols;
   double eps=0.000001;
   double c1,c2;
   double aintum,aintun;
   double c01,c02,c03,c04;
   double phix,phiy,GradPhi;
   double c0,t;
   double fit1,fit2;

   double ** nu;
   double ** Phi;
   double * kcpy;

   /*mxArray * affichageS[4]; */

   /**************************************************************************/
   /*************declaration of the functions used in the main program********/
   /**************************************************************************/

   double dirac(double,double);

   if(nrhs!=2){
   mexErrMsgTxt("Two inputs required.");
   }
   if(nlhs>1){
   mexErrMsgTxt("One output argument.");
   }

   mrows=mxGetM(prhs[0]);
   ncols=mxGetN(prhs[0]);

   nu=(double**)malloc(mrows*sizeof(double *));
   Phi=(double**)malloc(mrows*sizeof(double *));

   for(i=0;i<mrows;i++){
   nu[i]=(double*)calloc(ncols,sizeof(double));
   Phi[i]=(double*)calloc(ncols,sizeof(double));

   }


   for(i=0;i<ncols; i++){
            for(j=0; j<mrows; j++){
   Phi[j][i]=(mxGetPr(prhs[1]))[i*mrows+j];
   nu[j][i]=(mxGetPr(prhs[0]))[i*mrows+j];

   }
   }




    /*affichageS[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    affichageS[1]=mxCreateDoubleMatrix(1,1,mxREAL);
    affichageS[2]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
    affichageS[3]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);*/


    printf("initialization ok");



  for(k=0;k<Itermax;k++){

      m=0;
      n=0;
      aintum=0.0;
      aintun=0.0;

      for(i=1;i<mrows-1; i++){
              for(j=1;j<ncols-1; j++){

              if(Phi[i][j]>=0.0){
              m=m+1;
              aintum=aintum+nu[i][j];
              }
              else{
              n=n+1;
              aintun=aintun+nu[i][j];
              }

              }
      }


      if(m>0){
      c1=aintum/m;
      }


      if(n>0){
      c2=aintun/n;
      }



      for(i=1;i<mrows-1; i++){
              for(j=1;j<ncols-1; j++){

              phix=(Phi[i+1][j]-Phi[i][j]);
              phiy=(Phi[i][j+1]-Phi[i][j-1])/2.0;
              GradPhi=sqrt(eps+phix*phix+phiy*phiy);
              c01=1.0/GradPhi;


              phix=(Phi[i][j]-Phi[i-1][j]);
              phiy=(Phi[i-1][j+1]-Phi[i-1][j-1])/2.0;
              GradPhi=sqrt(eps+phix*phix+phiy*phiy);
              c02=1.0/GradPhi;



              phix=(Phi[i+1][j]-Phi[i-1][j])/2.0;
              phiy=(Phi[i][j+1]-Phi[i][j]);
              GradPhi=sqrt(eps+phix*phix+phiy*phiy);
              c03=1.0/GradPhi;



              phix=(Phi[i+1][j-1]-Phi[i-1][j-1])/2.0;
              phiy=(Phi[i][j]-Phi[i][j-1]);
              GradPhi=sqrt(eps+phix*phix+phiy*phiy);
              c04=1.0/GradPhi;



              c0=1.0+dt*eta*dirac(epsilon,Phi[i][j])*(c01+c02+c03+c04);

              t=c01*Phi[i+1][j]+c02*Phi[i-1][j]+c03*Phi[i][j+1]+c04*Phi[i][j-1];



              fit1=a1*(nu[i][j]-c1)*(nu[i][j]-c1);
              fit2=a2*(nu[i][j]-c2)*(nu[i][j]-c2);
              Phi[i][j]=(1.0/c0)*(Phi[i][j]+dt*dirac(epsilon,Phi[i][j])*(eta*t-fit1+fit2));

            }



         }

      for(j=1; j<ncols-1; j++){
      Phi[0][j]=Phi[1][j];
      Phi[mrows-1][j]=Phi[mrows-2][j];
      }

      for(i=1;i<mrows-1; i++){
      Phi[i][0]=Phi[i][1];
      Phi[i][ncols-1]=Phi[i][ncols-2];

      }



      Phi[0][0]=Phi[1][1];
      Phi[0][ncols-1]=Phi[1][ncols-2];
      Phi[mrows-1][0]=Phi[mrows-2][1];
      Phi[mrows-1][ncols-1]=Phi[mrows-2][ncols-2];

     /* if(k%50==0&k!=0){



           for(i=0;i<1;i++){
                for(j=0;j<1;j++){

                    (mxGetPr(affichageS[0]))[i+j]=1.0*mxGetM(prhs[1]);
                    (mxGetPr(affichageS[1]))[i+j]=1.0*mxGetN(prhs[1]);
                }
            }


           affichageS[2]=prhs[1];

           for(i=0;i<ncols;i++){
                for(j=0;j<mrows;j++){

            (mxGetPr(affichageS[3]))[i*mrows+j]=Phi[j][i];

                }
           }





            mexCallMATLAB(0,NULL,4,affichageS,"affichage");
        }


       */




}


  plhs[0]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
  kcpy=mxGetPr(plhs[0]);

  for(i=0;i<ncols;i++){
            for(j=0;j<mrows;j++){

              kcpy[i*mrows+j]=Phi[j][i];

  }
  }

  for(i=0;i<mrows;i++){
  free(nu[i]);
  free(Phi[i]);


  }


  free(nu),free(Phi);
}





    /**********************************************************************/
    /**********C infinity regularization of the Dirac function*************/
    /**********************************************************************/

    double dirac(double epsi,double Phi){
    #define PI 3.1415926
    return((epsi/PI)/(epsi*epsi+Phi*Phi));
    }
