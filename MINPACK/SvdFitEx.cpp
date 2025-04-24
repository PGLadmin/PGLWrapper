//#include <iomanip.h>
#include <cmath>
#include "c:\NumRep\cpp210\cpp\other\nr.h"
#include <iomanip>
using namespace std;
#include "c:\NumRep\cpp210\cpp\recipes\fleg.cpp"
#include "c:\NumRep\cpp210\cpp\recipes\fpoly.cpp"

// Driver for routine svdfit
//"svdcmp.cpp" "svbksb.cpp" "svdfit.cpp" "svdvar.cpp" "fpoly.cpp" "gasdev.cpp" "fleg.cpp"

	//SvdFitEz0(etaData,zFitData,nDensities,nPolyOrder,zRefCoeff,stdErr,&chisq);
int SvdFitEz0(double* xPas,double* yPas,int NPT,int NPOL,double* aPas,double* sigPas,double* chisqPas)
{
        //const int NPT=100,NPOL=5;
        //const DP SPREAD=0.002;
        int i,idum=(-911);
        DP chisq;
        Vec_DP x(NPT),y(NPT),sig(NPT),a(NPOL),w(NPOL);
        Mat_DP cvm(NPOL,NPOL),u(NPT,NPOL),v(NPOL,NPOL);
		for(i=0;i<NPT;i++) {
			x[i]=xPas[i];
			y[i]=yPas[i];
		}

        cout << fixed << setprecision(6);
        NR::svdfit(x,y,sig,a,u,v,w,chisq,NR::fpoly);
        NR::svdvar(v,w,cvm);
        cout << endl << "polynomial fit:" << endl << endl;
        for (i=0;i<NPOL;i++) {
          cout << setw(12) << a[i] << "  +-";
          cout << setw(11) << sqrt(cvm[i][i]) << endl;
        }
        cout << endl << "Chi-squared " << setw(12) << chisq << endl;
        NR::svdfit(x,y,sig,a,u,v,w,chisq,NR::fleg);
        NR::svdvar(v,w,cvm);
        cout << endl << "Legendre polynomial fit:" << endl << endl;
        for (i=0;i<NPOL;i++) {
          cout << setw(12) << a[i] << "  +-";
          cout << setw(11) << sqrt(cvm[i][i]) << endl;
        }
        cout << endl << "Chi-squared " << setw(12) << chisq << endl;
		*chisqPas=chisq;
		for(i=0;i<NPT;i++) {
			aPas[i]=a[i];
			sigPas[i]=sig[i];
		}
        return 0;
}
