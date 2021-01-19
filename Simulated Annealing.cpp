//#include <stdafx.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>

//////////////////////////////////////////////////////////////////////////////////
///////////////////Simulated Annealing////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

# define	m_N		22        //the number of variables in the wild type model
# define    M_N     19        //the number of variables in the mutant model
# define	k_N		39        //the number of parameter
# define	t_N		8         //time points of experimental data 
 	
double 	m[m_N];  //arry to save concentration values in the wild type model
double  M[M_N];  //arry to save concentration values in the mutant model
double 	k[k_N];  //arry to save parameters
//arry to save concentration values of various proteins in the model
double 	pro_com[t_N];  
double 	pro_comd[t_N];
double  pro_com1[t_N];
double  pro_comd1[t_N];
double  pro_mRNA[t_N];
double  pro_mRNA1[t_N];
double  pro_c1[t_N];
double  pro_cd1[t_N];
double  pro_c[t_N];
double  pro_cd[t_N];
double 	pro_rrp[t_N]; 
double 	pro_rrp1[t_N]; 
//specific initial parameters
double 	k_ini[k_N]={0.057,   0.0093,	    67.44,	    17.52,	    1.28,	    0.015,	    1.7,	    21.62,       0.000105,	 51.39,	     0.047,	     0.0044,	2.08,	    0.037, 	     0.00155,	   0.000403,     0.0082,    0.018,	     0.0024,	    0.00764,     0.4,	  0.7,	   0.3,  	1.2,     1.2,0.36,0.26,0.00058,0.033,0.00044,      0.0126,     0.00089,      0.162,	 0.072,	0.00006,	  0.001,	0.021,	  0.20,  0.20 };
		
//or randomized initial parameters    
/*double* random_kini (double k_ini[k_N])
{
	int i;
	double *y;
	for(i=0; i<k_N; i++)
	{
		k_ini[i]=(rand()%(100000000+1))/100000.0;
		printf("%lf \t",k_ini[i]);
	}

	printf("\n");
		y=&k_ini[0];
	return y;
}*/


double myRand()
{  //Generate random number between 0 and 1
    return rand()/(RAND_MAX+1);
}

double*	random_k (double kk[k_N])
{
	int i; //,flag=1;
	double	r1=0.6; // rate to detemine the neighbourhood range of the parameter
    double	r2=0.6; // rate to detemine the neighbourhood range of the parameter
         double r3=0.6;
//	double 	kkk[k_N],k_up[k_N-2],k_down[k_N-2];
	double *x;

	/////////////////////set k_up and k_down for the k known
 
	/////////////////////////////

    /////////////////////generate two last ones of k  
//	for(i=0; i<4; i++)
//	{
//	    kk[i]=kk[i]*(1-r2/2+r2*myRand());
//			printf("%lf \t",kk[i]);
//	}
//	      kk[4]=kk[4]*(1-r3/2+r3*myRand());
//		  printf("%lf \t",kk[4]);

		  for(i=0; i<k_N; i++)
	{
	    kk[i]=kk[i]*(1-r2/2+r2*myRand());
			printf("%lf \t",kk[i]);
	}
  //            kk[12]=kk[12]*(1-r1/2+r1*myRand());
  //               printf("%lf \t",kk[12]);
		printf(" \n");
  //  /////////////////////////

	/////////////////////generate the first 13 values of k
  

	x=&kk[0];
	return x;

}





double	sum_difference (double kk[k_N])
{
	int i,ii,tt,flag=1;
    double lw,t,SD,SD1; // time, sum of difference
	double h=0.0001; // time step
	//wild type model
    double x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21;
    double a[m_N][4];
	double dm[m_N];
	double m[m_N]={0.492, 0,  2.6, 0, 975,  0, 0, 0, 125, 0, 0, 0, 0, 25, 0.6, 0, 0.4, 0, 0.86, 0, 0.14, 0};
	//mutant model
    double X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18;
	double A[M_N][4];
	double dM[M_N];
	double M[M_N]={0.85, 0,  14, 0, 975,  0, 0, 0, 125, 0, 0, 0, 0, 25,  0.55, 0.54, 0, 0.46, 0};

    double time[t_N+1]={0,10,15,20,30,45,60,75,90}; // record time (s)


	//  nomalized experimental results of various proteins in the complex in different time points:

   	double lab_rrp[5]={100,  50.11554,  19.30678,  14.29006,  12.92186};
	double lab_com[4]={0.8826 ,   0.94734 ,   0.88606 ,   0.75434 };
	double lab_com1[4]={0.84929 ,  0.84545 ,  0.85593 ,   0.75575 };
	double lab_mRNA[5]={0.65063,   0.846,  1,  0.3831,   0.43239};
	double lab_mRNA1[5]={0.35406,  0.34614,  0.42798,  0.35934,  0.33206};

	double *kkk;
   // FILE	*fp;


	while (flag==1) // flag loop: do while any of the concentration value is negative
	{
		SD=0;
        flag=0;
		

		for(lw=0; lw<=1000; lw+=h)  //the model reach to steady state at 10 mM Mg2+
			{
				x0=m[0];x1=m[1];x2=m[2];x3=m[3];x4=m[4];x5=m[5];x6=m[6];x7=m[7];x8=m[8];x9=m[9];x10=m[10];x11=m[11];x12=m[12];x13=m[13];x14=m[14];x15=m[15];x16=m[16];x17=m[17];x18=m[18];x19=m[19];x20=m[20];x21=m[21];
				X0=M[0];X1=M[1];X2=M[2];X3=M[3];X4=M[4];X5=M[5];X6=M[6];X7=M[7];X8=M[8];X9=M[9];X10=M[10];X11=M[11];X12=M[12];X13=M[13];X14=M[14];X15=M[15];X16=M[16];X17=M[17];X18=M[18];
  
for(i = 0;i < 4;i++)  
        { 
		 
        //wild-type model
        dm[0] = -0*x0*x4+kk[1]*x5+kk[9]*x10-kk[13]*x0*x3+kk[14]*x12+kk[15]*x12-0*x0*x8+kk[17]*x7+kk[27]*x17+kk[29]-kk[31]*x0; //PhoQ   
        dm[1] = kk[3]*x6-kk[7]*x1*x2+kk[8]*x10-kk[31]*x1;  //P-PhoQ
        dm[2] = -kk[4]*x6*x2+kk[5]*x9-kk[7]*x1*x2+kk[8]*x10+kk[12]*x11+kk[15]*x12+kk[26]*x17+kk[28]-kk[30]*x2;  //PhoP
        dm[3] = kk[6]*x9+kk[9]*x10-kk[10]*x7*x3+kk[11]*x11-kk[13]*x0*x3+kk[14]*x12-kk[30]*x3;  //P-PhoP
        dm[4] = -0*x0*x4+kk[1]*x5+kk[18]*x8*x13-kk[19]*x4;  //ATP
        dm[5] = 0*x0*x4-kk[1]*x5-kk[2]*x5-kk[31]*x5;  //PhoQ-ATP
        dm[6] = kk[2]*x5-kk[3]*x6-kk[4]*x6*x2+kk[5]*x9-kk[31]*x6;  //P-PhoQ-ADP
        dm[7] = kk[6]*x9-kk[10]*x7*x3+kk[11]*x11+kk[12]*x11+0*x0*x8-kk[17]*x7-kk[31]*x7;  //PhoQ-ADP
        dm[8] = kk[3]*x6-0*x0*x8+kk[17]*x7-kk[18]*x8*x13+kk[19]*x4;  //ADP
        dm[9] = kk[4]*x6*x2-kk[5]*x9-kk[6]*x9-kk[31]*x9;  //P-PhoQ-ADP-PhoP
        dm[10] = kk[7]*x1*x2-kk[8]*x10-kk[9]*x10-kk[31]*x10;  //P-PhoQ-PhoP
        dm[11] = kk[10]*x7*x3-kk[11]*x11-kk[12]*x11-kk[31]*x11;  //PhoQ-ADP-P-PhoP
        dm[12] = kk[13]*x0*x3-kk[14]*x12-kk[15]*x12-kk[31]*x12;  //PhoQ-P-PhoP
		dm[13] = kk[12]*x11+kk[15]*x12-kk[18]*x8*x13+kk[19]*x4;  //Pi
        dm[15] = kk[20]*x3*x14-kk[21]*x15;  //p_1
        dm[16] = kk[22]*x2*x14-kk[23]*x16;  //p_2
		dm[17] = kk[24]*x15-kk[25]*x17;     //mRNA
		dm[19] = kk[32]*x3*x18-kk[33]*x19;  //p_1'
		dm[20] = kk[34]*x2*x18-kk[35]*x20;  //p_2'
		dm[21] = kk[36]*x19-kk[37]*x21;     //mRNA_((pmrD))
     
           
        //mutant model
        dM[0] = -kk[0]*X0*X4+kk[1]*X5+kk[9]*X10-kk[13]*X0*X3+kk[14]*X12+kk[15]*X12-kk[16]*X0*X8+kk[17]*X7+kk[27]*X14+kk[29]-kk[31]*X0;  //PhoQ 
        dM[1] = kk[3]*X6-kk[7]*X1*X2+kk[8]*X10-kk[31]*X1;  //P-PhoQ
        dM[2] = -kk[4]*X6*X2+kk[5]*X9-kk[7]*X1*X2+kk[8]*X10+kk[12]*X11+kk[15]*X12+kk[26]*X14+kk[28]-kk[30]*X2;  //PhoP
        dM[3] = kk[6]*X9+kk[9]*X10-kk[10]*X7*X3+kk[11]*X11-kk[13]*X0*X3+kk[14]*X12-kk[30]*X3;  //P-PhoP
        dM[4] = -kk[0]*X0*X4+kk[1]*X5+kk[18]*X8*X13-kk[19]*X4;  //ATP
        dM[5] = kk[0]*X0*X4-kk[1]*X5-kk[2]*X5-kk[31]*X5;   //PhoQ-ATP
        dM[6] = kk[2]*X5-kk[3]*X6-kk[4]*X6*X2+kk[5]*X9-kk[31]*X6;  //P-PhoQ-ADP
        dM[7] = kk[6]*X9-kk[10]*X7*X3+kk[11]*X11+kk[12]*X11+kk[16]*X0*X8-kk[17]*X7-kk[31]*X7;  //PhoQ-ADP
        dM[8] = kk[3]*X6-kk[16]*X0*X8+kk[17]*X7-kk[18]*X8*X13+kk[19]*X4;  //ADP
        dM[9] = kk[4]*X6*X2-kk[5]*X9-kk[6]*X9-kk[31]*X9;  //P-PhoQ-ADP-PhoP
        dM[10] = kk[7]*X1*X2-kk[8]*X10-kk[9]*X10-kk[31]*X10;  //P-PhoQ-PhoP
        dM[11] = kk[10]*X7*X3-kk[11]*X11-kk[12]*X11-kk[31]*X11;  //PhoQ-ADP-P-PhoP
        dM[12] = kk[13]*X0*X3-kk[14]*X12-kk[15]*X12-kk[31]*X12;  //PhoQ-P-PhoP
		dM[13] = kk[12]*X11+kk[15]*X12-kk[18]*X8*X13+kk[19]*X4;  //Pi
		dM[14] = kk[38]-kk[25]*X14;   //mRNA
		dM[16] = kk[32]*X3*X15-kk[33]*X16;  //p_1'
		dM[17] = kk[34]*X2*X15-kk[35]*X17;  //p_2'
		dM[18] = kk[36]*X16-kk[37]*X18;  //mRNA_((pmrD))
		
            
            a[0][i] = dm[0];  
            a[1][i] = dm[1];  
            a[2][i] = dm[2];  
            a[3][i] = dm[3];  
            a[4][i] = dm[4];  
            a[5][i] = dm[5];
			a[6][i] = dm[6];  
            a[7][i] = dm[7];  
            a[8][i] = dm[8];  
            a[9][i] = dm[9];  
            a[10][i] = dm[10];  
            a[11][i] = dm[11]; 
			a[12][i] = dm[12];  
            a[13][i] = dm[13]; 
            a[15][i] = dm[15]; 
			a[16][i] = dm[16];  
            a[17][i] = dm[17]; 
			a[19][i] = dm[19];  
			a[20][i] = dm[20];  
            a[21][i] = dm[21]; 
            
            A[0][i] = dM[0];  
            A[1][i] = dM[1];  
            A[2][i] = dM[2];  
            A[3][i] = dM[3];  
            A[4][i] = dM[4];  
            A[5][i] = dM[5];
			A[6][i] = dM[6];  
            A[7][i] = dM[7];  
            A[8][i] = dM[8];  
            A[9][i] = dM[9];  
            A[10][i] = dM[10];  
            A[11][i] = dM[11]; 
			A[12][i] = dM[12];  
            A[13][i] = dM[13]; 
            A[14][i] = dM[14]; 
			A[16][i] = dM[16];  
            A[17][i] = dM[17]; 
			A[18][i] = dM[18];  
		              
            if(i==0||i==1)  
            {  
    
				
                x0 = m[0]+h*a[0][i]/2;  
                x1 = m[1]+h*a[1][i]/2;  
                x2 = m[2]+h*a[2][i]/2;   
				x3 = m[3]+h*a[3][i]/2;  
                x4 = m[4]+h*a[4][i]/2;  
                x5 = m[5]+h*a[5][i]/2; 
                x6 = m[6]+h*a[6][i]/2;  
                x7 = m[7]+h*a[7][i]/2;  
                x8 = m[8]+h*a[8][i]/2;   
				x9 = m[9]+h*a[9][i]/2;  
                x10 = m[10]+h*a[10][i]/2;  
                x11 = m[11]+h*a[11][i]/2; 
				x12 = m[12]+h*a[12][i]/2; 
				x13 = m[13]+h*a[13][i]/2; 
                x15 = m[15]+h*a[15][i]/2; 
				x16 = m[16]+h*a[16][i]/2;  
                x17 = m[17]+h*a[17][i]/2; 
				x19 = m[19]+h*a[19][i]/2; 
            	x20 = m[20]+h*a[20][i]/2;  
                x21 = m[21]+h*a[21][i]/2; 
				
                X0 = M[0]+h*A[0][i]/2;  
                X1 = M[1]+h*A[1][i]/2;  
                X2 = M[2]+h*A[2][i]/2;   
				X3 = M[3]+h*A[3][i]/2;  
                X4 = M[4]+h*A[4][i]/2;  
                X5 = M[5]+h*A[5][i]/2; 
                X6 = M[6]+h*A[6][i]/2;  
                X7 = M[7]+h*A[7][i]/2;  
                X8 = M[8]+h*A[8][i]/2;   
				X9 = M[9]+h*A[9][i]/2;  
                X10 = M[10]+h*A[10][i]/2;  
                X11 = M[11]+h*A[11][i]/2; 
				X12 = M[12]+h*A[12][i]/2; 
				X13 = M[13]+h*A[13][i]/2; 
                X14 = M[14]+h*A[14][i]/2; 
				X16 = M[16]+h*A[16][i]/2;  
                X17 = M[17]+h*A[17][i]/2; 
				X18 = M[18]+h*A[18][i]/2; 
			
            }  
              
            if(i==2)  
            {  
     
								
				x0 = m[0]+h*a[0][i];  
                x1 = m[1]+h*a[1][i];  
                x2 = m[2]+h*a[2][i];   
				x3 = m[3]+h*a[3][i];  
                x4 = m[4]+h*a[4][i];  
                x5 = m[5]+h*a[5][i]; 
                x6 = m[6]+h*a[6][i];  
                x7 = m[7]+h*a[7][i];  
                x8 = m[8]+h*a[8][i];   
				x9 = m[9]+h*a[9][i];  
                x10 = m[10]+h*a[10][i];  
                x11 = m[11]+h*a[11][i]; 
				x12 = m[12]+h*a[12][i];
				x13 = m[13]+h*a[13][i];
                x15 = m[15]+h*a[15][i];
				x16 = m[16]+h*a[16][i];  
                x17 = m[17]+h*a[17][i]; 
				x19 = m[19]+h*a[19][i]; 
				x20 = m[20]+h*a[20][i];  
                x21 = m[21]+h*a[21][i]; 
			
                X0 = M[0]+h*A[0][i];  
                X1 = M[1]+h*A[1][i];  
                X2 = M[2]+h*A[2][i];   
				X3 = M[3]+h*A[3][i];  
                X4 = M[4]+h*A[4][i];  
                X5 = M[5]+h*A[5][i]; 
                X6 = M[6]+h*A[6][i];  
                X7 = M[7]+h*A[7][i];  
                X8 = M[8]+h*A[8][i];   
				X9 = M[9]+h*A[9][i];  
                X10 = M[10]+h*A[10][i];  
                X11 = M[11]+h*A[11][i]; 
				X12 = M[12]+h*A[12][i];
				X13 = M[13]+h*A[13][i];
                X14 = M[14]+h*A[14][i];
				X16 = M[16]+h*A[16][i];  
                X17 = M[17]+h*A[17][i]; 
				X18 = M[18]+h*A[18][i]; 
              
            }  
      
        }  

m[0] = m[0]+h*(a[0][0]+2*a[0][1]+2*a[0][2]+a[0][3])/6;  
m[1] = m[1]+h*(a[1][0]+2*a[1][1]+2*a[1][2]+a[1][3])/6;  
m[2] = m[2]+h*(a[2][0]+2*a[2][1]+2*a[2][2]+a[2][3])/6;
m[3] = m[3]+h*(a[3][0]+2*a[3][1]+2*a[3][2]+a[3][3])/6;  
m[4] = m[4]+h*(a[4][0]+2*a[4][1]+2*a[4][2]+a[4][3])/6;  
m[5] = m[5]+h*(a[5][0]+2*a[5][1]+2*a[5][2]+a[5][3])/6;
m[6] = m[6]+h*(a[6][0]+2*a[6][1]+2*a[6][2]+a[6][3])/6;  
m[7] = m[7]+h*(a[7][0]+2*a[7][1]+2*a[7][2]+a[7][3])/6;  
m[8] = m[8]+h*(a[8][0]+2*a[8][1]+2*a[8][2]+a[8][3])/6;
m[9] = m[9]+h*(a[9][0]+2*a[9][1]+2*a[9][2]+a[9][3])/6;  
m[10] = m[10]+h*(a[10][0]+2*a[10][1]+2*a[10][2]+a[10][3])/6;  
m[11] = m[11]+h*(a[11][0]+2*a[11][1]+2*a[11][2]+a[11][3])/6;
m[12] = m[12]+h*(a[12][0]+2*a[12][1]+2*a[12][2]+a[12][3])/6; 
m[13] = m[13]+h*(a[13][0]+2*a[13][1]+2*a[13][2]+a[13][3])/6;
m[15] = m[15]+h*(a[15][0]+2*a[15][1]+2*a[15][2]+a[15][3])/6;  
m[16] = m[16]+h*(a[16][0]+2*a[16][1]+2*a[16][2]+a[16][3])/6;
m[14] = 1-m[15]-m[16];  //p_0
m[17] = m[17]+h*(a[17][0]+2*a[17][1]+2*a[17][2]+a[17][3])/6;
m[19] = m[19]+h*(a[19][0]+2*a[19][1]+2*a[19][2]+a[19][3])/6; 
m[20] = m[20]+h*(a[20][0]+2*a[20][1]+2*a[20][2]+a[20][3])/6;
m[18] = 1-m[19]-m[20];  //p_0'
m[21] = m[21]+h*(a[21][0]+2*a[21][1]+2*a[21][2]+a[21][3])/6;

M[0] = M[0]+h*(A[0][0]+2*A[0][1]+2*A[0][2]+A[0][3])/6;  
M[1] = M[1]+h*(A[1][0]+2*A[1][1]+2*A[1][2]+A[1][3])/6;  
M[2] = M[2]+h*(A[2][0]+2*A[2][1]+2*A[2][2]+A[2][3])/6;
M[3] = M[3]+h*(A[3][0]+2*A[3][1]+2*A[3][2]+A[3][3])/6;  
M[4] = M[4]+h*(A[4][0]+2*A[4][1]+2*A[4][2]+A[4][3])/6;  
M[5] = M[5]+h*(A[5][0]+2*A[5][1]+2*A[5][2]+A[5][3])/6;
M[6] = M[6]+h*(A[6][0]+2*A[6][1]+2*A[6][2]+A[6][3])/6;  
M[7] = M[7]+h*(A[7][0]+2*A[7][1]+2*A[7][2]+A[7][3])/6;  
M[8] = M[8]+h*(A[8][0]+2*A[8][1]+2*A[8][2]+A[8][3])/6;
M[9] = M[9]+h*(A[9][0]+2*A[9][1]+2*A[9][2]+A[9][3])/6;  
M[10] = M[10]+h*(A[10][0]+2*A[10][1]+2*A[10][2]+A[10][3])/6;  
M[11] = M[11]+h*(A[11][0]+2*A[11][1]+2*A[11][2]+A[11][3])/6;
M[12] = M[12]+h*(A[12][0]+2*A[12][1]+2*A[12][2]+A[12][3])/6; 
M[13] = M[13]+h*(A[13][0]+2*A[13][1]+2*A[13][2]+A[13][3])/6;
M[14] = M[14]+h*(A[14][0]+2*A[14][1]+2*A[14][2]+A[14][3])/6;
M[16] = M[16]+h*(A[16][0]+2*A[16][1]+2*A[16][2]+A[16][3])/6; 
M[17] = M[17]+h*(A[17][0]+2*A[17][1]+2*A[17][2]+A[17][3])/6;
M[15] = 1-M[16]-M[17];  //p_0'
M[18] = M[18]+h*(A[18][0]+2*A[18][1]+2*A[18][2]+A[18][3])/6;

			}
		
		
		
		
		for(tt=0;tt<t_N;tt++)   //whole time loop
		{
		   		
			
		    for(t=time[tt]; t<time[tt+1]; t+=h) // part time loop
			{
	
  
      x0=m[0];x1=m[1];x2=m[2];x3=m[3];x4=m[4];x5=m[5];x6=m[6];x7=m[7];x8=m[8];x9=m[9];x10=m[10];x11=m[11];x12=m[12];x13=m[13];x14=m[14];x15=m[15];x16=m[16];x17=m[17];x18=m[18];x19=m[19];x20=m[20];x21=m[21];
				X0=M[0];X1=M[1];X2=M[2];X3=M[3];X4=M[4];X5=M[5];X6=M[6];X7=M[7];X8=M[8];X9=M[9];X10=M[10];X11=M[11];X12=M[12];X13=M[13];X14=M[14];X15=M[15];X16=M[16];X17=M[17];X18=M[18];
  
for(i = 0;i < 4;i++)  
        { 
		 
        //wild-type model
        dm[0] = -kk[0]*x0*x4+kk[1]*x5+kk[9]*x10-kk[13]*x0*x3+kk[14]*x12+kk[15]*x12-kk[16]*x0*x8+kk[17]*x7+kk[27]*x17+kk[29]-kk[31]*x0; //PhoQ   
        dm[1] = kk[3]*x6-kk[7]*x1*x2+kk[8]*x10-kk[31]*x1;  //P-PhoQ
        dm[2] = -kk[4]*x6*x2+kk[5]*x9-kk[7]*x1*x2+kk[8]*x10+kk[12]*x11+kk[15]*x12+kk[26]*x17+kk[28]-kk[30]*x2;  //PhoP
        dm[3] = kk[6]*x9+kk[9]*x10-kk[10]*x7*x3+kk[11]*x11-kk[13]*x0*x3+kk[14]*x12-kk[30]*x3;  //P-PhoP
        dm[4] = -kk[0]*x0*x4+kk[1]*x5+kk[18]*x8*x13-kk[19]*x4;  //ATP
        dm[5] = kk[0]*x0*x4-kk[1]*x5-kk[2]*x5-kk[31]*x5;  //PhoQ-ATP
        dm[6] = kk[2]*x5-kk[3]*x6-kk[4]*x6*x2+kk[5]*x9-kk[31]*x6;  //P-PhoQ-ADP
        dm[7] = kk[6]*x9-kk[10]*x7*x3+kk[11]*x11+kk[12]*x11+kk[16]*x0*x8-kk[17]*x7-kk[31]*x7;  //PhoQ-ADP
        dm[8] = kk[3]*x6-kk[16]*x0*x8+kk[17]*x7-kk[18]*x8*x13+kk[19]*x4;  //ADP
        dm[9] = kk[4]*x6*x2-kk[5]*x9-kk[6]*x9-kk[31]*x9;  //P-PhoQ-ADP-PhoP
        dm[10] = kk[7]*x1*x2-kk[8]*x10-kk[9]*x10-kk[31]*x10;  //P-PhoQ-PhoP
        dm[11] = kk[10]*x7*x3-kk[11]*x11-kk[12]*x11-kk[31]*x11;  //PhoQ-ADP-P-PhoP
        dm[12] = kk[13]*x0*x3-kk[14]*x12-kk[15]*x12-kk[31]*x12;  //PhoQ-P-PhoP
		dm[13] = kk[12]*x11+kk[15]*x12-kk[18]*x8*x13+kk[19]*x4;  //Pi
        dm[15] = kk[20]*x3*x14-kk[21]*x15;  //p_1
        dm[16] = kk[22]*x2*x14-kk[23]*x16;  //p_2
		dm[17] = kk[24]*x15-kk[25]*x17;     //mRNA
		dm[19] = kk[32]*x3*x18-kk[33]*x19;  //p_1'
		dm[20] = kk[34]*x2*x18-kk[35]*x20;  //p_2'
		dm[21] = kk[36]*x19-kk[37]*x21;     //mRNA_((pmrD))
     
           
        //mutant model
        dM[0] = -kk[0]*X0*X4+kk[1]*X5+kk[9]*X10-kk[13]*X0*X3+kk[14]*X12+kk[15]*X12-kk[16]*X0*X8+kk[17]*X7+kk[27]*X14+kk[29]-kk[31]*X0;  //PhoQ 
        dM[1] = kk[3]*X6-kk[7]*X1*X2+kk[8]*X10-kk[31]*X1;  //P-PhoQ
        dM[2] = -kk[4]*X6*X2+kk[5]*X9-kk[7]*X1*X2+kk[8]*X10+kk[12]*X11+kk[15]*X12+kk[26]*X14+kk[28]-kk[30]*X2;  //PhoP
        dM[3] = kk[6]*X9+kk[9]*X10-kk[10]*X7*X3+kk[11]*X11-kk[13]*X0*X3+kk[14]*X12-kk[30]*X3;  //P-PhoP
        dM[4] = -kk[0]*X0*X4+kk[1]*X5+kk[18]*X8*X13-kk[19]*X4;  //ATP
        dM[5] = kk[0]*X0*X4-kk[1]*X5-kk[2]*X5-kk[31]*X5;   //PhoQ-ATP
        dM[6] = kk[2]*X5-kk[3]*X6-kk[4]*X6*X2+kk[5]*X9-kk[31]*X6;  //P-PhoQ-ADP
        dM[7] = kk[6]*X9-kk[10]*X7*X3+kk[11]*X11+kk[12]*X11+kk[16]*X0*X8-kk[17]*X7-kk[31]*X7;  //PhoQ-ADP
        dM[8] = kk[3]*X6-kk[16]*X0*X8+kk[17]*X7-kk[18]*X8*X13+kk[19]*X4;  //ADP
        dM[9] = kk[4]*X6*X2-kk[5]*X9-kk[6]*X9-kk[31]*X9;  //P-PhoQ-ADP-PhoP
        dM[10] = kk[7]*X1*X2-kk[8]*X10-kk[9]*X10-kk[31]*X10;  //P-PhoQ-PhoP
        dM[11] = kk[10]*X7*X3-kk[11]*X11-kk[12]*X11-kk[31]*X11;  //PhoQ-ADP-P-PhoP
        dM[12] = kk[13]*X0*X3-kk[14]*X12-kk[15]*X12-kk[31]*X12;  //PhoQ-P-PhoP
		dM[13] = kk[12]*X11+kk[15]*X12-kk[18]*X8*X13+kk[19]*X4;  //Pi
		dM[14] = kk[38]-kk[25]*X14;   //mRNA
		dM[16] = kk[32]*X3*X15-kk[33]*X16;  //p_1'
		dM[17] = kk[34]*X2*X15-kk[35]*X17;  //p_2'
		dM[18] = kk[36]*X16-kk[37]*X18;  //mRNA_((pmrD))
		
            
            a[0][i] = dm[0];  
            a[1][i] = dm[1];  
            a[2][i] = dm[2];  
            a[3][i] = dm[3];  
            a[4][i] = dm[4];  
            a[5][i] = dm[5];
			a[6][i] = dm[6];  
            a[7][i] = dm[7];  
            a[8][i] = dm[8];  
            a[9][i] = dm[9];  
            a[10][i] = dm[10];  
            a[11][i] = dm[11]; 
			a[12][i] = dm[12];  
            a[13][i] = dm[13]; 
            a[15][i] = dm[15]; 
			a[16][i] = dm[16];  
            a[17][i] = dm[17]; 
			a[19][i] = dm[19];  
			a[20][i] = dm[20];  
            a[21][i] = dm[21]; 
            
            A[0][i] = dM[0];  
            A[1][i] = dM[1];  
            A[2][i] = dM[2];  
            A[3][i] = dM[3];  
            A[4][i] = dM[4];  
            A[5][i] = dM[5];
			A[6][i] = dM[6];  
            A[7][i] = dM[7];  
            A[8][i] = dM[8];  
            A[9][i] = dM[9];  
            A[10][i] = dM[10];  
            A[11][i] = dM[11]; 
			A[12][i] = dM[12];  
            A[13][i] = dM[13]; 
            A[14][i] = dM[14]; 
			A[16][i] = dM[16];  
            A[17][i] = dM[17]; 
			A[18][i] = dM[18];  
		              
            if(i==0||i==1)  
            {  
    
				
                x0 = m[0]+h*a[0][i]/2;  
                x1 = m[1]+h*a[1][i]/2;  
                x2 = m[2]+h*a[2][i]/2;   
				x3 = m[3]+h*a[3][i]/2;  
                x4 = m[4]+h*a[4][i]/2;  
                x5 = m[5]+h*a[5][i]/2; 
                x6 = m[6]+h*a[6][i]/2;  
                x7 = m[7]+h*a[7][i]/2;  
                x8 = m[8]+h*a[8][i]/2;   
				x9 = m[9]+h*a[9][i]/2;  
                x10 = m[10]+h*a[10][i]/2;  
                x11 = m[11]+h*a[11][i]/2; 
				x12 = m[12]+h*a[12][i]/2; 
				x13 = m[13]+h*a[13][i]/2; 
                x15 = m[15]+h*a[15][i]/2; 
				x16 = m[16]+h*a[16][i]/2;  
                x17 = m[17]+h*a[17][i]/2; 
				x19 = m[19]+h*a[19][i]/2; 
            	x20 = m[20]+h*a[20][i]/2;  
                x21 = m[21]+h*a[21][i]/2; 
				
                X0 = M[0]+h*A[0][i]/2;  
                X1 = M[1]+h*A[1][i]/2;  
                X2 = M[2]+h*A[2][i]/2;   
				X3 = M[3]+h*A[3][i]/2;  
                X4 = M[4]+h*A[4][i]/2;  
                X5 = M[5]+h*A[5][i]/2; 
                X6 = M[6]+h*A[6][i]/2;  
                X7 = M[7]+h*A[7][i]/2;  
                X8 = M[8]+h*A[8][i]/2;   
				X9 = M[9]+h*A[9][i]/2;  
                X10 = M[10]+h*A[10][i]/2;  
                X11 = M[11]+h*A[11][i]/2; 
				X12 = M[12]+h*A[12][i]/2; 
				X13 = M[13]+h*A[13][i]/2; 
                X14 = M[14]+h*A[14][i]/2; 
				X16 = M[16]+h*A[16][i]/2;  
                X17 = M[17]+h*A[17][i]/2; 
				X18 = M[18]+h*A[18][i]/2; 
			
            }  
              
            if(i==2)  
            {  
     
								
				x0 = m[0]+h*a[0][i];  
                x1 = m[1]+h*a[1][i];  
                x2 = m[2]+h*a[2][i];   
				x3 = m[3]+h*a[3][i];  
                x4 = m[4]+h*a[4][i];  
                x5 = m[5]+h*a[5][i]; 
                x6 = m[6]+h*a[6][i];  
                x7 = m[7]+h*a[7][i];  
                x8 = m[8]+h*a[8][i];   
				x9 = m[9]+h*a[9][i];  
                x10 = m[10]+h*a[10][i];  
                x11 = m[11]+h*a[11][i]; 
				x12 = m[12]+h*a[12][i];
				x13 = m[13]+h*a[13][i];
                x15 = m[15]+h*a[15][i];
				x16 = m[16]+h*a[16][i];  
                x17 = m[17]+h*a[17][i]; 
				x19 = m[19]+h*a[19][i]; 
				x20 = m[20]+h*a[20][i];  
                x21 = m[21]+h*a[21][i]; 
			
                X0 = M[0]+h*A[0][i];  
                X1 = M[1]+h*A[1][i];  
                X2 = M[2]+h*A[2][i];   
				X3 = M[3]+h*A[3][i];  
                X4 = M[4]+h*A[4][i];  
                X5 = M[5]+h*A[5][i]; 
                X6 = M[6]+h*A[6][i];  
                X7 = M[7]+h*A[7][i];  
                X8 = M[8]+h*A[8][i];   
				X9 = M[9]+h*A[9][i];  
                X10 = M[10]+h*A[10][i];  
                X11 = M[11]+h*A[11][i]; 
				X12 = M[12]+h*A[12][i];
				X13 = M[13]+h*A[13][i];
                X14 = M[14]+h*A[14][i];
				X16 = M[16]+h*A[16][i];  
                X17 = M[17]+h*A[17][i]; 
				X18 = M[18]+h*A[18][i]; 
              
            }  
      
        }  

m[0] = m[0]+h*(a[0][0]+2*a[0][1]+2*a[0][2]+a[0][3])/6;  
m[1] = m[1]+h*(a[1][0]+2*a[1][1]+2*a[1][2]+a[1][3])/6;  
m[2] = m[2]+h*(a[2][0]+2*a[2][1]+2*a[2][2]+a[2][3])/6;
m[3] = m[3]+h*(a[3][0]+2*a[3][1]+2*a[3][2]+a[3][3])/6;  
m[4] = m[4]+h*(a[4][0]+2*a[4][1]+2*a[4][2]+a[4][3])/6;  
m[5] = m[5]+h*(a[5][0]+2*a[5][1]+2*a[5][2]+a[5][3])/6;
m[6] = m[6]+h*(a[6][0]+2*a[6][1]+2*a[6][2]+a[6][3])/6;  
m[7] = m[7]+h*(a[7][0]+2*a[7][1]+2*a[7][2]+a[7][3])/6;  
m[8] = m[8]+h*(a[8][0]+2*a[8][1]+2*a[8][2]+a[8][3])/6;
m[9] = m[9]+h*(a[9][0]+2*a[9][1]+2*a[9][2]+a[9][3])/6;  
m[10] = m[10]+h*(a[10][0]+2*a[10][1]+2*a[10][2]+a[10][3])/6;  
m[11] = m[11]+h*(a[11][0]+2*a[11][1]+2*a[11][2]+a[11][3])/6;
m[12] = m[12]+h*(a[12][0]+2*a[12][1]+2*a[12][2]+a[12][3])/6; 
m[13] = m[13]+h*(a[13][0]+2*a[13][1]+2*a[13][2]+a[13][3])/6;
m[15] = m[15]+h*(a[15][0]+2*a[15][1]+2*a[15][2]+a[15][3])/6;  
m[16] = m[16]+h*(a[16][0]+2*a[16][1]+2*a[16][2]+a[16][3])/6;
m[14] = 1-m[15]-m[16];  //p_0
m[17] = m[17]+h*(a[17][0]+2*a[17][1]+2*a[17][2]+a[17][3])/6;
m[19] = m[19]+h*(a[19][0]+2*a[19][1]+2*a[19][2]+a[19][3])/6; 
m[20] = m[20]+h*(a[20][0]+2*a[20][1]+2*a[20][2]+a[20][3])/6;
m[18] = 1-m[19]-m[20];  //p_0'
m[21] = m[21]+h*(a[21][0]+2*a[21][1]+2*a[21][2]+a[21][3])/6;

M[0] = M[0]+h*(A[0][0]+2*A[0][1]+2*A[0][2]+A[0][3])/6;  
M[1] = M[1]+h*(A[1][0]+2*A[1][1]+2*A[1][2]+A[1][3])/6;  
M[2] = M[2]+h*(A[2][0]+2*A[2][1]+2*A[2][2]+A[2][3])/6;
M[3] = M[3]+h*(A[3][0]+2*A[3][1]+2*A[3][2]+A[3][3])/6;  
M[4] = M[4]+h*(A[4][0]+2*A[4][1]+2*A[4][2]+A[4][3])/6;  
M[5] = M[5]+h*(A[5][0]+2*A[5][1]+2*A[5][2]+A[5][3])/6;
M[6] = M[6]+h*(A[6][0]+2*A[6][1]+2*A[6][2]+A[6][3])/6;  
M[7] = M[7]+h*(A[7][0]+2*A[7][1]+2*A[7][2]+A[7][3])/6;  
M[8] = M[8]+h*(A[8][0]+2*A[8][1]+2*A[8][2]+A[8][3])/6;
M[9] = M[9]+h*(A[9][0]+2*A[9][1]+2*A[9][2]+A[9][3])/6;  
M[10] = M[10]+h*(A[10][0]+2*A[10][1]+2*A[10][2]+A[10][3])/6;  
M[11] = M[11]+h*(A[11][0]+2*A[11][1]+2*A[11][2]+A[11][3])/6;
M[12] = M[12]+h*(A[12][0]+2*A[12][1]+2*A[12][2]+A[12][3])/6; 
M[13] = M[13]+h*(A[13][0]+2*A[13][1]+2*A[13][2]+A[13][3])/6;
M[14] = M[14]+h*(A[14][0]+2*A[14][1]+2*A[14][2]+A[14][3])/6;
M[16] = M[16]+h*(A[16][0]+2*A[16][1]+2*A[16][2]+A[16][3])/6; 
M[17] = M[17]+h*(A[17][0]+2*A[17][1]+2*A[17][2]+A[17][3])/6;
M[15] = 1-M[16]-M[17];  //p_0'
M[18] = M[18]+h*(A[18][0]+2*A[18][1]+2*A[18][2]+A[18][3])/6;
       
             
        for(i=0; i<m_N; i++)
				{
		            if (m[i]<0)
					{
                        flag=1;
						break;
					}
				}
				
				
		for(i=0; i<M_N; i++)
				{
		            if (M[i]<0)
					{
                        flag=1;
						break;
					}
				}



		} // exit the part time loop


            //concentration values of various proteins in the model
           
	         pro_rrp[tt]=m[3];    //P-PhoP
             pro_com[tt]=(m[19]+m[20])/(m[18]);   //occupancy of the pmrD promoter
			 pro_com1[tt]=(M[16]+M[17])/(M[15]);
			 pro_mRNA[tt]=m[21];   //mRNA expression of pmrD
		     pro_mRNA1[tt]=M[18];
		     

			if (flag==1)
			{
               // generate random reaction rates in the neighbourhood of the previous ones
			   kkk=random_k(kk);
	 	       for(i=0; i<k_N; i++)
			   {
		           kk[i]=kkk[i];
			   }
               break;
			} 

		} // exit the whole time loop	
           
	  

	for(tt=0;tt<t_N;tt++)

		{
       	pro_rrp1[tt]=(pro_rrp[tt]/pro_rrp[1])*100;     
        pro_comd[tt]=pro_mRNA[tt]/pro_mRNA[3];          
		pro_comd1[tt]=pro_mRNA1[tt]/pro_mRNA[3];       
		
		} 
      
        //this can use pow( ,2) or fabs()
		SD=fabs(pro_rrp[1]-lab_rrp[0])+fabs(pro_rrp[3]-lab_rrp[1])+fabs(pro_rrp[4]-lab_rrp[2])+fabs(pro_rrp[5]-lab_rrp[3])+fabs(pro_rrp[6]-lab_rrp[4]);   //sum of P-PhoP difference
        //sum of occupancy of the pmrD promoter difference
        SD+=fabs(pro_com[0]-lab_com[0])+fabs(pro_com[2]-lab_com[1])+fabs(pro_com[3]-lab_com[2])+fabs(pro_com[5]-lab_com[3]);                
		SD+=fabs(pro_com1[0]-lab_com1[0])+fabs(pro_com1[2]-lab_com1[1])+fabs(pro_com1[3]-lab_com1[2])+fabs(pro_com1[5]-lab_com1[3]);        
        //sum of mRNA expression of pmrD promoter difference
		SD+=fabs(pro_comd[0]-lab_mRNA[0])+fabs(pro_comd[2]-lab_mRNA[1])+fabs(pro_comd[3]-lab_mRNA[2])+fabs(pro_comd[5]-lab_mRNA[3])+fabs(pro_comd[7]-lab_mRNA[4]); 
        SD+=fabs(pro_comd1[0]-lab_mRNA1[0])+fabs(pro_comd1[2]-lab_mRNA1[1])+fabs(pro_comd1[3]-lab_mRNA1[2])+fabs(pro_comd1[5]-lab_mRNA1[3])+fabs(pro_comd1[7]-lab_mRNA1[4]); 

	} // exit the flag loop

	for (i=0; i<k_N; i++) // Change the global variable k[k_N] to be kk[k_N]
	{
        k[i]=kk[i];
    }
	return SD;
}




int main()
{
	int i,j, ii,n,l,  NN=100;
	int	iterations,   threshold=100;
	double	temperature, SD_best,  cooling_rate=1;
	double  k_previous[k_N], k_current[k_N], k_best[k_N],m_best[m_N],pro_com_best[t_N];
	double previous_SD,current_SD,diff;
	double *kk;
//	double *k_ini;
    FILE	*fp;
	
	srand( (unsigned)time( NULL ) ); // generate the random seed

   
for(l=1;l<=10;l++)
{


//    	k_ini=random_kini(k);

        temperature=50.00001;
	    iterations = 1;
     //   SD_best=10000000000.0;
		previous_SD=sum_difference(k_ini);
		SD_best=previous_SD;

			///////////////////// initial value of k
        
		for(i=0; i<k_N; i++)
		{
			k_previous[i]=k_ini[i];
		    k[i]=k_ini[i];
			k_best[i]=k_ini[i];
		}

       ////////////////////
    
 
	



       for(n=1;n<=50;n++) // trial times loop
	{
	
        
	

        //while (iterations < threshold) // threshold loop
 	   	for (j=1; j<=NN; j++)
		{
	        	for(i=0; i<k_N; i++)
				{
		 //	      k_previous[i]=k_best[i];
				  printf("%lf\t",k_previous[i]);
				}
		
	            printf("\n");


			kk=random_k(k);
	 	    for(i=0; i<k_N; i++)
			{
		        k[i]=kk[i];
			}
			current_SD = sum_difference(k);

			for(i=0; i<k_N; i++)
			{
		 	    k_current[i]=k[i];
			}

			printf("%d % lf\t  %lf\n",j,previous_SD,current_SD);


			  if (current_SD < SD_best)
			{
				SD_best=current_SD;
			
                for(i=0; i<k_N; i++)
				{
		 	      
				   k_best[i]=k[i];
				}
               
					
			}////////////////////////////////////////////
             
			diff = fabs(current_SD - previous_SD);
			//////////////////////////////////////////// record the best solution
 	        if (current_SD < previous_SD)
			{
		//		SD_best=current_SD;
			    previous_SD=current_SD;
                for(i=0; i<k_N; i++)
				{
		 	      k[i]= k_current[i];
				  k_previous[i]=k_current[i];
		//		   k_best[i]=k[i];
				}
               
					
			}////////////////////////////////////////////

		    // absolute difference between current_SD and previous_SD
            

            //////////////////////////////////// simulated annealing process

		    if ((current_SD > previous_SD) && (myRand() < exp(-diff/(temperature))))
			{
			//	SD_best=current_SD;
			   previous_SD=current_SD;
                for(i=0; i<k_N; i++)
				{
		 	       k[i]= k_current[i];
				   k_previous[i]=k_current[i];
		//		   k_best[i]=k[i];
				}
			}

	
			else 
			{
			//	SD_best=previous_SD;
		//		previous_SD=previous_SD;
				 for(i=0; i<k_N; i++)
				{
		 	       k[i]= k_previous[i];
	//			   k_best[i]=k[i];
				}
			}
			//


			
	         for(i=0; i<k_N; i++)
				{
		      
				  printf("%lf\t",k_best[i]);
				}
		
	            printf("\n");
		  
		
		   printf("%d % lf\t %lf\t %lf\n",j,previous_SD,current_SD,SD_best);


			 //////////////////// generate random reaction rates in the neighbourhood of the previous ones
			
            ////////////////////////
			

		} 
		
		temperature=temperature-cooling_rate;  ///cooling

		
		printf("%d %lf\n",n,temperature);
		previous_SD=SD_best;
		for(i=0; i<k_N; i++)
				{
		 	      k_previous[i]=k_best[i];
				   k[i]=k_best[i];
		//		  printf("%lf\t",k_best[i]);
				}
		
	     printf("\n");
		
		
	} // exit the threshold loop

        ///////////////////////////////////////// output
		
       

        // output of k_best
        fp=fopen("KK.txt","a");
	    fprintf(fp,"%d %.6f\t",l,SD_best);
	    for(i=0; i<k_N; i++)
		{    
		    fprintf(fp,"%.10f\t",k_best[i]);
		}
	    fprintf(fp,"\n");

        fclose(fp);
}
      
   
	return 0 ;
}
