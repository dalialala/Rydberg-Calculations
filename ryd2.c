/*Program to compute Rydberg eigenenergies using 
quantum defect values of Gallagher PRA 67, 052502 (2003)*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#define max 200
#define PI 3.14159265359
#define c 2.99792458E8 //speed of light in m/s
#define R 10973731.5685 //Rydberg constant in m-1
#define mm 10E-1
#define um 10E-4
#define nm 10E-7

long double rydritz(int nn, long double matd[][2], int l);
void save(long double nmat[][max], char *nombre, int nmax, int nmin, int dm);
long double error(int nn, long double ns, long double matd[][2], long double matde[][2], long double Rr, int l);
void saveexp(long double nmat[][max], char *nombre, int nmax, int nmin, int dm);

int main(void)
{
	/*---Variables definition---*/
	int i, l, k, nmin, nmax; //Integers for condition and loop
	long double d0s=3.1311804, d2s=0.1784, d0p1=2.6548849, d2p1=0.2900, d0p3=2.6416737, d2p3= 0.2950, d0d3=1.34809171, d2d3=-0.60283, 		d0d5=1.34646572, d2d5= -0.59600, d0f5= 0.0165192, d2f5=-0.085, d0f7= 0.0165437, d2f7= -0.086; //Quantum defects for ns to nf states
	long double de0s=0.0000010, de2s=0.0006, de0p1=0.0000010, de2p1=0.0006, de0p3=0.000001, de2p3= 0.0007, de0d3=0.00000040, de2d3=0.00026, 	de0d5=0.00000030, de2d5= 0.00018, de0f5= 0.0000009, de2f5=-0.009, de0f7=  0.0000007, de2f7= -0.007; //Quantum defects errors for ns to nf 		states
	long double matd[7][2], matde[7][2], vns[max][max], En[max][max], Enf[max][max], Errw[max][max], Errf[max][max]; //Matrix for quantum defects, quantum defects errors, 	eigen energies and their errors; 
	long double me=9.10938188E-31, mR=1.44316060E-25; //mass of electron and mass of 87Rb in kg
	long double Rr= 1.09736605E7, Rr2=R/(1+me/mR); //87Rb Rydberg constant taken from Gallagher/Calculated from effective mass
	long double Ei= 3369080.48; //87Rb Ionization limit m-1 from Steck 0.00002 uncertainty
	long double D2= 1281654.938993; //87Rb D2 line m-1 from Steck 0.0000000021 uncertainty
	long double EF3= 193.7407E6/c; //Hyperfine energy splitting for 5p3/2 F=3 from Steck 0.000046 uncertainty in MHz/c
	long double ED2F3= D2+EF3; //D2 line energy plus hyperfine splitting from F=3 state in cm-1
	long double D1= 1257895.0985; //87Rb D1 line m-1 from Steck 0.0000000013 uncertainty
	long double EF2= 305.4324E6/c; //Hyperfine energy splitting for 5p1/2 F=2 from Steck 0.00011 uncertainty in MHz/c
	long double ED1F2= D1+EF2; //D1 line energy plus hyperfine splitting from F=2 state in m-1
	long double Ef, ns, mon; //Energy of either D1 or D2 transition, quantum defect, error of energies
	char name[50], name2[50], name3[50], name4[50], name5[50]; //Name for text file where eigenenergies and errors will be saved
		
	//printf("%7.10Lf\t%7.10Lf\n", EF3, ED2F3);
	//printf("%7.10Lf\t%7.10Lf\n", Rr, Rr2);

	//---Parameters defined by user---\\

	printf("\nFor D1 transition write 1, D2 line transition write 2:\n");
	scanf("%d", &k);

	if(k==1)
		Ef=ED1F2;
	else if(k==2)
		Ef=ED2F3;
	
	printf("\nMinimum principal quantum number:\n");
	scanf("%d", &nmin);

	printf("\nMaximum principal quantum number:\n");
	scanf("%d", &nmax);
	
	sprintf(name, "D%d_n*_%d-%d", k, nmin, nmax); //Name assignment 

	sprintf(name2, "D%d_Freq_%d-%d", k, nmin, nmax); //Name assignment 
	
	sprintf(name3, "D%d_Wave_%d-%d", k, nmin, nmax); //Name assignment 

	sprintf(name4, "D%d_Freq_err_%d-%d", k, nmin, nmax); //Name assignment 
	
	sprintf(name5, "D%d_Wave_err_%d-%d", k, nmin, nmax); //Name assignment 

	
	//---Quantum defect matrix matd[m][i]---\\	

	matd[0][0]=d0s;  matd[0][1]=d2s;
	matd[1][0]=d0p1; matd[1][1]=d2p1;
	matd[2][0]=d0p3; matd[2][1]=d2p3;
	matd[3][0]=d0d3; matd[3][1]=d2d3;
	matd[4][0]=d0d5; matd[4][1]=d2d5;
	matd[5][0]=d0f5; matd[5][1]=d2f5;
	matd[6][0]=d0f7; matd[6][1]=d2f7;

	//---Quantum defect matrix uncertanties matd[m][i]---\\		
	
	matde[0][0]=de0s;  matde[0][1]=de2s;
	matde[1][0]=de0p1; matde[1][1]=de2p1;
	matde[2][0]=de0p3; matde[2][1]=de2p3;
	matde[3][0]=de0d3; matde[3][1]=de2d3;
	matde[4][0]=de0d5; matde[4][1]=de2d5;
	matde[5][0]=de0f5; matde[5][1]=de2f5;
	matde[6][0]=de0f7; matde[6][1]=de2f7;

	/*for(i=0; i<7; i++)
	{
		printf("%Lf\t\t%Lf\t\n", matd[i][0], matd[i][1]);
	}*/
	

	//---Loop for computing n*, energies and their uncertanties vectors in units of s-1 and m---\\

	for(i=nmin; i<nmax+1; i++)
	{
		//printf("%d\t", i);
		for(l=0; l<7; l++)
			{
				ns=rydritz(i, matd, l);
				vns[i-nmin][l]=i-ns;
				En[i-nmin][l]=1/(Ei-(Rr/(pow(i-ns,2)))-Ef);
				Enf[i-nmin][l]=c/En[i-nmin][l];
				Errw[i-nmin][l]=error(i, ns, matd, matde, Rr, l)/pow(Ei-Rr/(pow(i-ns,2))-Ef,2);
				Errf[i-nmin][l]=c*Errw[i-nmin][l]/pow(En[i-nmin][l],2);
				//printf("%1.5Le   ",  Ei-(Rr/(pow(i-ns,2))));
					
			}
			
			//printf("\n");
	}

	//Saves all vectors calculated in the loop above in different files
	save(vns, name, nmax, nmin, 7);
	saveexp(Enf, name2, nmax, nmin, 7);
	saveexp(En, name3, nmax, nmin, 7);
	saveexp(Errf, name4, nmax, nmin, 7);
	saveexp(Errw, name5, nmax, nmin, 7);
}


/*------------------------------------------------------------------
This function calculates the quantum defect using the Rydberg-Ritz formula given
	nn: Principal quantum number
	matd: Quantum defect matrix
	l: Matrix element matd[l][i]
------------------------------------------------------------------*/
long double rydritz(int nn, long double matd[][2], int l)
{
	long double dn;
	dn=matd[l][0]+matd[l][1]/(pow(nn-matd[l][0], 2));
	return dn;
}

//Function to calculate error in s-1 using propagation of error from d0 and d2
long double error(int nn, long double ns, long double matd[][2], long double matde[][2], long double Rr, int l)
{
	long double se;
	se=(2*Rr/pow(nn-ns,3))*sqrt(pow((1-(2*matd[l][1]/pow(nn-matd[l][0],3))),2)*pow(matde[l][0],2)+(pow(matde[l][1],2)/pow(nn-matd[l][0],4)));
	return se;
}


//Function to save the Eigenenergies in a table in a text file
void save(long double nmat[][max], char *nombre, int nmax, int nmin, int dm)
{
	FILE *p_matriz;
	int i,j;
	p_matriz=fopen(nombre, "wt");
	if(p_matriz!=NULL)
	{
		fprintf(p_matriz, "n\t\ts1/2\t\tp1/2\t\tp3/2\t\td3/2\t\td5/2\t\tf5/2\t\tf7/2\n");
		for(i=0;i<nmax-nmin+1;i++)
		{
			fprintf(p_matriz, "%d\t", i+nmin);
			for(j=0; j<dm; j++)
				fprintf(p_matriz, "    %3.6Lf\t", nmat[i][j]);
			fprintf(p_matriz, "\n");
		}
		fclose(p_matriz);	
	}
}

//Function to save the error in exponential form in a table in a text file
void saveexp(long double nmat[][max], char *nombre, int nmax, int nmin, int dm)
{
	FILE *p_matriz;
	int i,j;
	p_matriz=fopen(nombre, "wt");
	if(p_matriz!=NULL)
	{
		fprintf(p_matriz, "n\t\ts1/2\t\tp1/2\t\tp3/2\t\td3/2\t\td5/2\t\tf5/2\t\tf7/2\n");
		for(i=0;i<nmax-nmin+1;i++)
		{
			fprintf(p_matriz, "%d\t", i+nmin);
			for(j=0; j<dm; j++)
				fprintf(p_matriz, "  %3.6Le  ", nmat[i][j]);
			fprintf(p_matriz, "\n");
		}
		fclose(p_matriz);	
	}
}
