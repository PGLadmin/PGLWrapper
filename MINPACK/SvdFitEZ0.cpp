#include <iostream>
#include <iomanip>
#include <cmath>
#include "nr.h"
using namespace std;


int SVDFITEZ0(Vec_DP x,Vec_DP y,int ndata, int nParms,Vec_DP &parm,DP &chisq)
{


    int NPT=ndata;
	const int NPOL=nParms;
    int i;
    //Vec_DP x(NPT),y(NPT);
	Vec_DP sig(NPT),a(NPOL),w(NPOL);
    Mat_DP cvm(NPOL,NPOL),u(NPT,NPOL),v(NPOL,NPOL);
	
    for (i=0;i<NPT;i++) 
	{
		sig[i]=1;
		//sig[i]=1.0/y[i] // gives excel result for unweighted least squares.
		//sig[i]=y[i]*0.02
		//sig[i]=y[i]*0.05 // gives %err MINimum	and err estimates based on 5% err
    }
	
	NR::svdfit(x,y,sig,a,u,v,w,chisq,NR::fpoly);
    NR::svdvar(v,w,cvm);
    for (i=0;i<NPOL;i++) 
	{
		parm[i]=a[i];

	}
    //cout << endl << "polynomial fit:" << endl << endl;
    //for (i=0;i<NPOL;i++) {
    //  cout << setw(12) << a[i] << "  +-";
    //  cout << setw(11) << sqrt(cvm[i][i]) << endl;
    //}
    // cout << endl << "Chi-squared " << setw(12) << chisq << endl;

	return 0;

}


void NR::svdfit(Vec_I_DP &x, Vec_I_DP &y, Vec_I_DP &sig, Vec_O_DP &a,Mat_O_DP &u, Mat_O_DP &v, Vec_O_DP &w, DP &chisq,void funcs(const DP, Vec_O_DP &))
{
	int i,j;
	const DP TOL=1.0e-13;
	DP wmax,tmp,thresh,sum;

	//int ndata=x.size();
	int ndata=sig.size();
	int ma=a.size();
	Vec_DP b(ndata),afunc(ma);
	for (i=0;i<ndata;i++) 
	{
		funcs(x[i],afunc);
		tmp=1.0/sig[i];
		for (j=0;j<ma;j++) 
		{
			u[i][j]=afunc[j]*tmp;
		}
		b[i]=y[i]*tmp;
	}
	svdcmp(u,w,v);
	wmax=0.0;
	for (j=0;j<ma;j++)
		if (w[j] > wmax) wmax=w[j];
	thresh=TOL*wmax;
	for (j=0;j<ma;j++)
		if (w[j] < thresh) w[j]=0.0;
	svbksb(u,w,v,b,a);
	chisq=0.0;
	for (i=0;i<ndata;i++) {
		funcs(x[i],afunc);
		sum=0.0;
		for (j=0;j<ma;j++) sum += a[j]*afunc[j];
		chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
	}
}


void NR::svdcmp(Mat_IO_DP &a, Vec_O_DP &w, Mat_O_DP &v)
{
	bool flag;
	int i,its,j,jj,k,l,nm;
	DP anorm,c,f,g,h,s,scale,x,y,z;

	int m=a.nrows();
	int n=a.ncols();
	Vec_DP rv1(n);
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) 
			{
				scale += fabs(a[k][i]);
			}
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i != n) {
			for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l-1]=f-g;
				for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=false;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l-1;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if (fabs(f)+anorm == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) nrerror("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
}


void NR::svdvar(Mat_I_DP &v, Vec_I_DP &w, Mat_O_DP &cvm)
{
	int i,j,k;
	DP sum;

	int ma=w.size();
	Vec_DP wti(ma);
	for (i=0;i<ma;i++) {
		wti[i]=0.0;
		if (w[i] != 0.0) wti[i]=1.0/(w[i]*w[i]);
	}
	for (i=0;i<ma;i++) {
		for (j=0;j<i+1;j++) {
			sum=0.0;
			for (k=0;k<ma;k++)
				sum += v[i][k]*v[j][k]*wti[k];
			cvm[j][i]=cvm[i][j]=sum;
		}
	}
}


void NR::svbksb(Mat_I_DP &u, Vec_I_DP &w, Mat_I_DP &v, Vec_I_DP &b, Vec_O_DP &x)
{
	int jj,j,i;
	DP s;

	int m=u.nrows();
	int n=u.ncols();
	Vec_DP tmp(n);
	for (j=0;j<n;j++) {
		s=0.0;
		if (w[j] != 0.0) {
			for (i=0;i<m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=0;j<n;j++) {
		s=0.0;
		for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
}

DP NR::pythag(const DP a, const DP b)
{
	DP absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


void NR::fpoly(const DP x, Vec_O_DP &p)
{
	int j;
	int np=p.size();
	p[0]=x;//The lowest term is x because the y-intercept is zero
	for (j=1;j<np;j++) 
		p[j]=p[j-1]*x;
	return;
}

