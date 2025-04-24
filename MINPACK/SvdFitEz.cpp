//#include <iostream>
//#include <iomanip>
#include <cmath>
//#include "nr.h"
#include "\NumRep\cpp210\cpp\other\nr.h"
using namespace std;

// Driver for routine svdfit

void PolyFit(const DP x, Vec_O_DP &p)
{
	int j;

	int np=p.size();
	p[0]=1.0;
	for (j=1;j<np;j++) p[j]=p[j-1]*x;
}

void PolyFit0(const DP x, Vec_O_DP &p)
{
	int j;

	int np=p.size();
	p[0]=x;
	for (j=1;j<np;j++) p[j]=p[j-1]*x;
}


int SvdFitEz(double *xData, double *yData, int nData, int nParms, double *parm, double *stdErr, double *chisqPas)
{
        //const int NPT=100,NPOL=5;
        const DP SPREAD=0.002;
        int i,idum=(-911);
		int errCode=0;
        DP chisq;
        Vec_DP x(nData),y(nData),sig(nData),a(nParms),w(nParms);
        Mat_DP cvm(nParms,nParms),u(nData,nParms),v(nParms,nParms);

        //cout << fixed << setprecision(6);
        for (i=0;i<nData;i++) {
          x[i]=xData[i];
          y[i]=yData[i];
          //y[i] *= (1.0+SPREAD*NR::gasdev(idum));
          //sig[i]=y[i]*SPREAD;
          sig[i]=1;
        }
        NR::svdfit(x,y,sig,a,u,v,w,chisq,PolyFit);
        NR::svdvar(v,w,cvm);
        for (i=0;i<nParms;i++) {
          //cout << setw(12) << a[i] << "  +-";
          //cout << setw(11) << sqrt(cvm[i][i]) << endl;
			parm[i]=a[i];
			stdErr[i]=cvm[i][i];
        }
        //cout << endl << "Chi-squared " << setw(12) << chisq << endl;
		*chisqPas=chisq;
        return errCode;
}

int SvdFitEz0(double *xData, double *yData, int nData, int nParms, double *parm, double *stdErr, double *chisqPas)
{
        //const int NPT=100,NPOL=5;
        const DP SPREAD=0.002;
        int i,idum=(-911);
		int errCode=0;
        DP chisq;
        Vec_DP x(nData),y(nData),sig(nData),a(nParms),w(nParms);
        Mat_DP cvm(nParms,nParms),u(nData,nParms),v(nParms,nParms);

        //cout << fixed << setprecision(6);
        for (i=0;i<nData;i++) {
          x[i]=xData[i];
          y[i]=yData[i];
          //y[i] *= (1.0+SPREAD*NR::gasdev(idum));
          //sig[i]=y[i]*SPREAD;
          sig[i]=1;
        }
        NR::svdfit(x,y,sig,a,u,v,w,chisq,PolyFit0);
        NR::svdvar(v,w,cvm);
        for (i=0;i<nParms;i++) {
          //cout << setw(12) << a[i] << "  +-";
          //cout << setw(11) << sqrt(cvm[i][i]) << endl;
			parm[i]=a[i];
			stdErr[i]=sqrt( cvm[i][i]*chisq/(nData-nParms) );
        }
        //cout << endl << "Chi-squared " << setw(12) << chisq << endl;
		*chisqPas=chisq;
        return errCode;
}
