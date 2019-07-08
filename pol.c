/*--------------------------------------------------

Polarizability calculation for a Rydberg state ns1/2 given by
the sum over all n->n2 dipole allowed matrix elements 

|<n,l,j,m|z|n2,k,j2,m2>|²/(En-En2)

------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf.h>
#define max 50000
//#define pc 6.626070040E-34 //planck constant J*s
#define pc 4.135667662E-15 //planck constant eV*s 
#define kb 1.38064852E-23 //Boltzmann constant
#define PI 3.14159265359
#define c 137.036 //speed of light in cm/s au
#define fs 7.2973525698E-3 // fine structure constant
#define aB 5.2917721092E-11 //Bohr radius m
#define Rc 3.289841960355E15 //Rydberg constant in Hz c*R
#define Rr 1.09736605E7 //Rb Rydberg constant in m-1
#define csi 2.99792458E8 //speed of light in m/s
#define e0 8.854187817620E-12 //vacuum permitivitty SI units
#define ec 1.60217662E-19//electric charge in C


double potential(double En, double mu, double xvec[max], int l, int i, double j, int a);
void funguarda(double vec[max], double vec2[max], char *nombre, int N);
double norm(double X[max], double xvec[max], int N);
double rydritz(int nn, double matd[][2], int l);
double dme(double X[max], double Xk[max], double xvec[max], double nc, double nck, int N, int N2, int k);
double lev(double j2, int n2);

int main(void)
{
	/*---Variables definition---*/
	int i, l, k, s, p, n, n2, m, m2, N, N2; //Integers for condition, loop, energy and momentum numbers
	double d0s=3.1311804, d2s=0.1784, d0p1=2.6548849, d2p1=0.2900, d0p3=2.6416737, d2p3= 0.2950, d0d3=1.34809171, d2d3=-0.60286, 		d0d5=1.34646572, d2d5= -0.59600, d0f5= 0.0165192, d2f5=-0.085, d0f7= 0.0165437, d2f7= -0.086; //Quantum defects for ns to nf states
	double matd[7][2]; //Matrix for quantum defects; 
	double j, j2, mj, mj2, qd, qd2, En,En2, nqd, nqd2, dip, den; //Tot angular mom, Energy of state with principal quantum number n
	char name[30];//Name for text file where eigenenergies will be saved
	double x, xs, xs2, xmin, mon, h, Xk[max], X[max], xvec[max], xveck[max], f1, f2, f3, f1k, f2k, f3k; //Vectors and variables for radial wavefunction calculation
	double rbm=87*1836+36;//Rb core mass in au, 87 nucleons +36 core e
	double mu=rbm/(1+rbm);// Reduced mass of Rb in au
	double R=1/(4*PI*c); //Rydberg constant in au
	double nc, nck; //Normalization constant
	double Ei= 3369080.48; //87Rb Ionization limit m-1 from Steck 0.00002 uncertainty
	double ED2=csi*1281654.938993, ED1=csi*1257895.0985, ED5, deg; //d2 line, d1 line degeneracy

	
	//---Parameters defined by user---\\

	printf("\nPrincipal quantum number:\n");
	scanf("%d", &n);

	printf("\nFor l=0, write 0; l=1, write 1; l=2, write 2; l=3, write 3:\n");
	scanf("%d", &l);
	
	printf("\nFor j=1/2, write 0.5; j=3/2, write 1.5; j=5/2, write 2.5; j=7/2, write 3.5:\n");
	scanf("%lf", &j);

	printf("\nFor mj=1/2, write 0.5; j=3/2, write 1.5; j=5/2, write 2.5; j=7/2, write 3.5:\n");
	scanf("%lf", &mj);
	
	
	//---Conditions for choosing the quantum defect matrix matd[m][i=1,2] element to compute quantum defect dn---\\
				
	if(l==0&&j==0.5)
		m=0;
	else if(l==1&&j==0.5)
		m=1;
	else if(l==1&&j==1.5)
		m=2;
	else if(l==2&&j==1.5)
		m=3;
	else if(l==2&&j==2.5)
		m=4;
	else if(l==3&&j==2.5)
		m=5;
	else if(l==3&&j==3.5)
		m=6;

	//---Quantum defect matrix matd[m][i]---\\	
	
	matd[0][0]=d0s;  matd[0][1]=d2s;
	matd[1][0]=d0p1; matd[1][1]=d2p1;
	matd[2][0]=d0p3; matd[2][1]=d2p3;
	matd[3][0]=d0d3; matd[3][1]=d2d3;
	matd[4][0]=d0d5; matd[4][1]=d2d5;
	matd[5][0]=d0f5; matd[5][1]=d2f5;
	matd[6][0]=d0f7; matd[6][1]=d2f7;

	double am1, am2,am3, am4, ang, avec[max], asum,ae; //Angular momentum constants and polarizability array
	int q, v=0, t=1;

	k=l+1; //Angular momentum for the second state, in particualr for l=s, l2=p
	
	int a;
	printf("\nFor 1/x2, write 0; For mc, write 1; For qc, write 2:\n"); //Calculation method chosen by user: 0 Coulomb, 1 Marinescu, 2 Quasiclassical
	scanf("%d", &a);
	
	qd=rydritz(n, matd, m);// Quantum defect function call for n
	nqd=n-qd;
	En=-1/(2*pow(nqd,2)); //Ionization energy w qd from energy units (hcRy) to au
	
	if(a==0) //Innermost point for radial wavefunction calculation
		xmin=2.08;
	else 
		xmin=0.05;
	
	//---Radial wavefunction calculation for n---\\
	
	h=0.01; //step size
	xs=sqrt(2*n*(n+15)); //outermost value of r
	N=(xs-xmin)/h; //Number of points
	xvec[0]=xs; //specify first 2 values of xvec
	xvec[1]=xs-h;
	X[0]=1.0E-10; //specify first 2 values of transformed radial function X
	X[1]=-1.1E-10;
	x=xvec[1];
	xvec[2]=xvec[1]-h;
	i=0;
	while(x>=xmin)
	{
		xvec[i+2]=x=xvec[i+1]-h;
		//mon=potential(En, mu, xvec, l,i+2, j, a);	//monitor
		f1=potential(En, mu, xvec, l, i+1, j, a)-(12/pow(h,2));
		f2=10*potential(En, mu, xvec, l, i+1, j, a)+(24/pow(h,2));
		f3=(12/pow(h,2))-potential(En, mu, xvec, l, i+2, j, a);
		X[i+2]=(f1*X[i]+f2*X[i+1])/f3; //Transformed radial wavefunction evaluated at xvec[i]
		i++;	
	}
	
	nc=norm(X, xvec, N); //Normalization constant for wavefunction 

	ae=-(aB*aB)/(Rr*pc*pc*csi)*10000;// Constant to convert atomic untis to MHz (cm/V)²

	
	//--Dipole Matrix Elements calculation for polarizability---\\

	j2=0; //Cicle for angular momentum j=1/2, 3/2 initialization

	for(p=0; p<2; p++)
	{
		j2=p+0.5; //Angular momentum j2 assignment
	
		//Conditions for choosing element from the quantum defect matrix matd[m][i=1,2] and compute quantum defect dn2
		if(k==0&&j2==0.5)
			m2=0;
		else if(k==1&&j2==0.5)
			m2=1;
		else if(k==1&&j2==1.5)
			m2=2;
		else if(k==2&&j2==1.5)
			m2=3;
		else if(k==2&&j2==2.5)
			m2=4;
		else if(k==3&&j2==2.5)
			m2=5;
		else if(k==3&&j2==3.5)
			m2=6;	
		
		mj2=mj; // For <z> mj2=mj due to angular momentum conservation
		q=0;		

		//Angular momentum constants
		am1=sqrt(k*(2*j+1)*(2*j2+1));
		am2=gsl_sf_coupling_6j(2*j2, 2, 2*j, 2*l, 1, 2*k);
		am3=gsl_sf_coupling_3j(2*j2, 2, 2*j, -2*mj, 2*q, 2*mj2);
		ang=am1*am2*am3;

		//printf("j2=%lf, ang=%lf, 3j=%lf \n",j2, am1, am3);
		
		printf("\nContribution for np%1.1lf", j2);

		//---Loop for computing |<n|z|n2>|²/(En-En2)---\\

		nck=0;
		qd2=0;
		nqd2=0;
		En2=0;
		
		for(n2=5; n2<n; n2++) //Polarizability calculation for n2<n
		{
			if(n2<=10)  //Condition for Energy, for n<10 quantum defect theory is not very accurate, energy levels are taken from NIST database
				En2=lev(j2, n2);
			else
			{
				qd2=rydritz(n2, matd, m2); //Quantum defect for n2
				nqd2=n2-qd2;
				En2=-1/(2*pow(nqd2,2)); //Ionization energy w qd from energy units (hcRy) in au
			}
		
			xs2=sqrt(2*n2*(n2+15)); //outermost value of r for n2
			N2=(xs2-xmin)/h; //step size
			xveck[0]=xs2; //specify first 2 values of xveck
			xveck[1]=xs2-h;
			Xk[0]=1.0E-10; //specify first 2 values of transformed radial function Xk
			Xk[1]=-1.1E-10;
			x=xveck[1];
			xveck[2]=xveck[1]-h;
			i=0;
			while(x>=xmin)
			{
				xveck[i+2]=x=xveck[i+1]-h;
				f1k=potential(En2, mu, xveck, k,i+2, j2, a)-(12/pow(h,2));
				f2k=10*potential(En2, mu, xveck, k, i+1, j2, a)+(24/pow(h,2));
				f3k=(12/pow(h,2))-potential(En2, mu, xveck, k, i+2, j2, a);
				Xk[i+2]=(f1k*Xk[i]+f2k*Xk[i+1])/f3k; //Transformed radial wavefunction evaluated at xvec[i]
				i++;	
			}
			nck=norm(Xk, xveck, N2); //Normalization constant for wavefunction 
			den=1/(En-En2); //Polarizability denominator
			am4=pow(-1, j2-mj2+j-0.5+k-l); //Sign of <n|z|n2> not important at the end
			avec[n2-5]=ae*pow(ang,2)*pow((dme(X, Xk, xveck, nc, nck, N, N2, t)),2)*den; //Polarizability contribution from <n|z|n2>
			//printf("<%d|r|%d>=%lf \t %3.6le\n", n, n2, dme(X, Xk, xveck, nc, nck, N, N2, t), avec[n2-5]);
		}
				

		for(n2=n; n2<n*2+6; n2++)  //Polarizability calculation for n2<2n+6
		{
			qd2=rydritz(n2, matd, m2);
			nqd2=n2-qd2;
			En2=-1/(2*pow(nqd2,2)); //Ionization energy w qd from energy units (hcRy) in au
			xs2=sqrt(2*n2*(n2+15)); //outermost value of r for n2
			h=0.01; //step size
			N2=(xs2-xmin)/h; //number of vector entries
			xveck[0]=xs2; //specify first 2 values of xveck
			xveck[1]=xs2-h;
			Xk[0]=1.0E-10; //specify first 2 values of transformed radial function Xk
			Xk[1]=-1.1E-10;
			x=xveck[1];
			xveck[2]=xveck[1]-h;
			i=0;
			while(x>=xmin)
			{
				xveck[i+2]=x=xveck[i+1]-h;
				f1k=potential(En2, mu, xveck, k,i+2, j2, a)-(12/pow(h,2));
				f2k=10*potential(En2, mu, xveck, k, i+1, j2, a)+(24/pow(h,2));
				f3k=(12/pow(h,2))-potential(En2, mu, xveck, k, i+2, j2, a);
				Xk[i+2]=(f1k*Xk[i]+f2k*Xk[i+1])/f3k; //Transformed radial wavefunction evaluated at xvec[i]
				i++;	
			}
			nck=norm(Xk, xveck, N2); //Normalization constant for wavefunction 
			den=1/(En-En2); //Polarizability denominator
			am4=pow(-1, j2-mj2+j-0.5+k-l);  //Sign of <n|z|n2> not important at the end
			avec[n2-5]=ae*pow(am1*am2*am3,2)*pow((dme(Xk, X, xvec, nc, nck, N2, N, t)),2)*den; //Polarizability contribution from <n|z|n2>
			//printf("<%d|r|%d>=%lf \t %3.6le\n", n, n2, dme(Xk, X, xvec, nc, nck, N2, N, t), avec[n2-5]);
		}
	
	
	//---Final polarizability calculation---\\
	
	asum=0; //Sum of every dipole <n|z|n2> contribution from avec[n2-5] vector element 
	for(i=0; i<n*2+1; i++)
	{
		asum=asum+avec[i];
		//printf("%d\t%lf\n", i+5, asum);
	}
	
	printf("\nTotalMDE=%d\talfa=%3.5le\tConst=%3.5le\n", i+5, asum, ae);

	//printf("decay=%1.20lf, decay1=%1.20lf, decay4=%1.20lf\n", decay, decay1, decay4);
	
	sprintf(name, "%ds%1.1lfPolarizability_mj%lf", n, j2, mj2); //Name assignment 
	
	//funguarda(avec, avec, name, n*2+1); //Saves avec in a file

	}	
	
}


/*------------------------------------------------------------------
This function calculates the potential energy to compute the radial wavefunction given the input arguments:
	En: Energy of n2
	mu: Rb reduced mass
	l: orbital angular momentum
	i: position of radial wavefunction
	j: angular momentum J
	a: integer to choose which potential to use 
Returns gx as the potential energy depending whether a=0 Coulomb, 1 Marinescu, 2 Quasiclassical
------------------------------------------------------------------*/
double potential(double En, double mu, double xvec[max], int l, int i, double j, int a)
{
	double zl, zl1, zl2, cp, pot, gx, tot, soc, cp2;
	double zld, zl1d, zl2d, cpd, potd, gjl, soct, soc2;
	double x=xvec[i];
	double r=x*x;
	double matp[5][5]; //Matrix for model parameters
	double ac=9.0760; //Rb core polarizability
	double a1l0=3.69628474, a2l0=1.64915255, a3l0=-9.86069196, a4l0=0.19579987, rcl0=1.66242117; //l-dependent parameters for potential model from Marinescu l=0
	double a1l1=4.44088978, a2l1=1.92828831, a3l1=-16.79597770, a4l1=-0.81633314, rcl1=1.50195124; //l-dependent parameters for potential model from Marinescu l=1
	double a1l2=3.78717363, a2l2=1.57027864, a3l2=-11.65588970, a4l2=0.52942835, rcl2=4.86851938; //l-dependent parameters for potential model from Marinescu l=2
	double a1l3=2.39848933, a2l3=1.76810544, a3l3=-12.07106780, a4l3=0.77256589, rcl3=4.79831327; //l-dependent parameters for potential model from Marinescu l=3
	
	matp[0][0]=rcl0; matp[0][1]=a1l0; matp[0][2]=a2l0; matp[0][3]=a3l0; matp[0][4]=a4l0;
	matp[1][0]=rcl1; matp[1][1]=a1l1; matp[1][2]=a2l1; matp[1][3]=a3l1; matp[1][4]=a4l1;
	matp[2][0]=rcl2; matp[2][1]=a1l2; matp[2][2]=a2l2; matp[2][3]=a3l2; matp[2][4]=a4l2;
	matp[3][0]=rcl3; matp[3][1]=a1l3; matp[3][2]=a2l3; matp[3][3]=a3l3; matp[3][4]=a4l3;
	
	//Marinescu potential calculation: compensates for valence electron penetration and core polarization
	zl1=1+(36/exp(matp[l][1]*r)); //Marinescu potential sub part 1a
	zl2=r*(matp[l][3]+(matp[l][4]*r))/exp(matp[l][2]*r);//Marinescu potential sub part 2a
	cp2=pow((r/matp[l][0]), 6);//Marinescu potential sub sub part b
	cp=0.5*ac*(1-(1/(exp(cp2))));//Marinescu potential sub part 2b
	zl=zl1-zl2;//Marinescu potential part a
	pot=-(zl/r)-(cp/pow(r,4));//Marinescu total potential part a+b

	//Quasiclassical potential calculation Marinescu potential + spin orbit coupling (the latter is negligible) 	doi: 10.1103/PhysRevA.91.032509
	zl1d=(1+((36/exp(matp[l][1]*r))*(1+matp[l][1]*r)))/(pow(r,2));//zl1 derivative
	zl2d=(matp[l][4]-(matp[l][2]*(matp[l][3]+(matp[l][4]*r))))/exp(matp[l][2]*r);//zl2 derivative
	cpd=2*fs*(1-((1-1.5*cp2)/exp(cp2)))/pow(r,5);//derivative of cp
	potd=zl1+zl2+cpd;
	if(l==0)
		gjl=0;  //spin orbit coupling potential for l=0
	else
		gjl=0.5*pow(fs,2)*(j*(j+1)-l*(l+1)-0.75)/pow(r, 3); //spin orbit coupling potential
	soc=(pow(fs,2)/r)*potd*gjl;
	soc2=1-pow(fs,2)*pot;
	soct=soc/pow(soc2,2); //total spin orbit coupling potential	

	if(a==0) // given int a, return the total potential of:
		tot=-1/r; //Coulomb potential
	else if(a==1)
		tot=pot+gjl;//Marinestu potential
	else if(a==2)
		tot=pot+soct;//Quasiclassical potential

	gx=-((8*mu*pow(x,2))*(En-tot))+((2*l+1.5)*(2*l+0.5)/(pow(x,2)));

	return gx;
}

/*------------------------------------------------------------------
This function calculates the quantum defect using the Rydberg-Ritz formula given
	nn: Principal quantum number
	matd: Quantum defect matrix
	l: Matrix element matd[l][i]
------------------------------------------------------------------*/
double rydritz(int nn, double matd[][2], int l)
{
	double dn;
	dn=matd[l][0]+matd[l][1]/(pow(nn-matd[l][0], 2));
	return dn;
}


/*------------------------------------------------------------------
This function calculates the normalization constant of the radial wavefunction using transformation x=sqrt(r), see SA Bhatti et al 94 
	X: transformed radial wavefunction vector
	xvec: position vector
	N: number of vector elements
------------------------------------------------------------------*/
double norm(double X[max], double xvec[max], int N)
{
	double nc, nc2, aux[max], R, r;
	int i;
	nc=0;
	r=0;
	for(i=0;i<=N;i++)
	{
		r=xvec[i]*xvec[i];
		aux[i]=pow(X[i],2)*r;
		nc=nc+aux[i];
	}
	nc2=sqrt(nc);
	return nc2;
}

/*------------------------------------------------------------------
This function calculates the dipole matrix element <n1|r^k|n2> given
	X: transformed radial wavefunction vector for the largest n
	Xk: transformed radial wavefunction vector for the lowest n
	xvec: position vector of the lowest n
	nc: n1 normalization constant
	nck: n2 normalization constant
	N: number of points of the largest vector which corresponds to largest n
	n: number of points of the shortest vector which corresponds to lowest n
	k: Power of r
------------------------------------------------------------------*/
double dme(double X[max], double Xk[max], double xvec[max], double nc, double nck, int N, int n, int k)
{
	long double aux[max],hch;
	int i;
	int p=2*k+2;
	long double dipt, dip=0;
	for(i=0;i<=n;i++)
	{
		aux[i]=X[i+(N-n)]*Xk[i]*pow(xvec[i],p);
		dip=dip+aux[i];
	}
	dipt=dip/(nc*nck);
	return dipt;
}

/*------------------------------------------------------------------
This function returns the energy for n=5-10 with l=p, to get more accurate results 
	j2: angular momentum J
	n2: principal quantum number
------------------------------------------------------------------*/
double lev(double j2, int n2)
{
	double energy, Et;
	double Ei=3369080.48; //87Rb Ionization limit m-1 from Steck 0.00002 uncertainty
	if(j2==0.5&&n2==5)
		energy=1257895.0;
	else if(j2==1.5&&n2==5)
		energy=1281654.5;
	if(j2==0.5&&n2==6)
		energy=2371508.1;
	else if(j2==1.5&&n2==6)
		energy=2379259.1;
	if(j2==0.5&&n2==7)
		energy=2783502;
	else if(j2==1.5&&n2==7)
		energy=2787011;
	if(j2==0.5&&n2==8)
		energy=2983494;
	else if(j2==1.5&&n2==8)
		energy=2985379;
	if(j2==0.5&&n2==9)
		energy=3095891;
	else if(j2==1.5&&n2==9)
		energy=3097019;
	if(j2==0.5&&n2==10)
		energy=3165385;
	else if(j2==1.5&&n2==10)
		energy=3166116;
	Et=(energy-Ei)/(2*Rr);
	return Et;
}

/*------------------------------------------------------------------
This function saves two vectors in a file of the form 

	i	vec[i] 		vec2[i]

	vec: vector 1 to save
	vec2: vector 1 to save
	nombre: file name
	N: number of array elements
------------------------------------------------------------------*/
void funguarda(double vec[max], double vec2[max], char *nombre, int N)
{
	FILE *p_vector;
	int i;
	double r, r2;
	double R[max];
	p_vector=fopen(nombre, "wt");
	if(p_vector!=NULL)
	{
		for(i=0;i<=N-5;i++)
		{
			fprintf(p_vector, "%d\t%lf\t%lf\n", i+5, vec[i], vec[i]*vec2[i]);
		fprintf(p_vector, "\n");
		}
		fclose(p_vector);	
	}
}

