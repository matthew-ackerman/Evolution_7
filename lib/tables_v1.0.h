#ifndef _TABLES_HPP
#define _TABLES_HPP

#include "MersenneTwister.h"
#include <cstdlib>
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

#define MAXINT 		(unsigned int)(0-1)

using namespace std;

long double A, B;

//sets the first parameter of a distribution.

void set_A (long double a){
	A=a;
}

//sets the second parameter of a distribution. 

void set_B (long double b){
	B=b;
}

//double tables 

class DOUBLE_TABLE
{	
	private:

	unsigned int Allocate, sAllocate;

	int allocate(){

		element = (double *) calloc(Allocate, sizeof(double));

		if(element == NULL) {
			fprintf(stderr, "error allocating memory\n");
			return 1;
		}
		return 0;
	}

	public:

	double *element;

	DOUBLE_TABLE(){
		element=NULL;
	}

	void openXorMakeTable(char * str, double ( *return_function  ) (unsigned int, unsigned int max), int size){
		fstream file;
		file.open(str,ios::in);
		if( file.is_open() ) {
			file.close();
			makeTableFromFile(str);
		}
		else {
			file.close();
			makeTable(return_function, size);
			saveTableToFile(str);
		}
	}

	int makeTable(double ( *return_function  ) (unsigned int, unsigned int), int size){
		Allocate=size;
		sAllocate=Allocate-1;
		if(allocate()!=0){
			fprintf(stderr, "error allocating memory\n");
			return 1;
		}

		for(unsigned int i=0; i < Allocate; i++){
			element[i]=return_function(i+1,Allocate);
  		}
	}

	int saveTableToFile(char *str){
		FILE* TableFile;
		TableFile=fopen(str, "wb");
		if(TableFile==NULL){
			fprintf(stderr, "cannot open file %s\n", str);
			return 1;
		}

		fwrite (&Allocate, sizeof(unsigned int), 1, TableFile);

		for(int i=0; i<Allocate; i++) fwrite (&element[i], sizeof(double), 1, TableFile);

		fclose(TableFile);
	}
	int makeTableFromFile(char *str){
		FILE* TableFile;
		TableFile=fopen(str, "rb");
		if(TableFile==NULL){
			fprintf(stderr, "cannot open file %s\n", str);
			return 1;
		}

		fread(&Allocate, sizeof(Allocate),1,TableFile);
		sAllocate=Allocate-1;

		if(allocate()!=0){
			fprintf(stderr, "error allocating memory\n");
			return 1;
		}

		for(int i=0; i<Allocate; i++) fread (&element[i], sizeof(double), 1, TableFile);

		fclose(TableFile);
		return 0;
	}
	
	~DOUBLE_TABLE(){
		if(element!=NULL) free(element);
	}
	void clear(void){
		if(element!=NULL) free(element);
		element=NULL;
	}

	void print(void){
		for(int x=0;x<Allocate;x++) cout << element[x] << " ";
		cout << endl;
	}

	inline double getElement(unsigned int x){
		return *(element+x);
	}
	inline double getRand(MTRand *mtrand){
		return *(element+mtrand->randInt(sAllocate));
	}

	double mean(){
		double mean=0;
		for(int X=0;X<Allocate;X++) mean+=element[X];
		return mean/Allocate;
	}

	unsigned int min(){
		double min=0;
		min--;
		for(int X=0;X<Allocate;X++) if(element[X]<min) min=element[X];
		return min;
	}

	unsigned int maximum(){
		double max=0;
		for(int X=0;X<Allocate;X++) if(element[X]>max) max=element[X];
		return max;
	}
};


/* Function_tables are made to store lookup tables of mathmatical functions. */

class FUNCTION_TABLE
{	
	private:

	unsigned int Allocate, sAllocate;

	int allocate(){

		element = (unsigned int *) calloc(Allocate, sizeof(unsigned int));

		if(element == NULL) {
			fprintf(stderr, "error allocating memory\n");
			return 1;
		}
		return 0;
	}

	public:

	unsigned int *element;

	FUNCTION_TABLE(){
		element=NULL;
	}

	void openXorMakeTable(char * str, unsigned int ( *return_function  ) (unsigned int, unsigned int max), int size){
		fstream file;
		file.open(str,ios::in);
		if( file.is_open() ) {
			file.close();
			makeTableFromFile(str);
		}
		else {
			file.close();
			makeTable(return_function, size);
			saveTableToFile(str);
		}
	}

	int makeTable(unsigned int ( *return_function  ) (unsigned int, unsigned int), int size){
		Allocate=size;
		sAllocate=Allocate-1;
		if(allocate()!=0){
			fprintf(stderr, "error allocating memory\n");
			return 1;
		}

		for(unsigned int i=0; i < Allocate; i++){
			element[i]=return_function(i+1,Allocate);
  		}
	}

	int saveTableToFile(char *str){
		FILE* TableFile;
		TableFile=fopen(str, "wb");
		if(TableFile==NULL){
			fprintf(stderr, "cannot open file %s\n", str);
			return 1;
		}

		fwrite (&Allocate, sizeof(unsigned int), 1, TableFile);

		for(int i=0; i<Allocate; i++) fwrite (&element[i], sizeof(unsigned int), 1, TableFile);

		fclose(TableFile);
	}
	int makeTableFromFile(char *str){
		FILE* TableFile;
		TableFile=fopen(str, "rb");
		if(TableFile==NULL){
			fprintf(stderr, "cannot open file %s\n", str);
			return 1;
		}

		fread(&Allocate, sizeof(Allocate),1,TableFile);
		sAllocate=Allocate-1;

		if(allocate()!=0){
			fprintf(stderr, "error allocating memory\n");
			return 1;
		}

		for(int i=0; i<Allocate; i++) fread (&element[i], sizeof(unsigned int), 1, TableFile);

		fclose(TableFile);
		return 0;
	}
	
	~FUNCTION_TABLE(){
		if(element!=NULL) free(element);
	}
	void clear(void){
		if(element!=NULL) free(element);
		element=NULL;
	}

	void print(void){
		for(int x=0;x<Allocate;x++) cout << element[x] << " ";
		cout << endl;
	}

	inline unsigned int getElement(unsigned int x){
		return *(element+x);
	}
	inline unsigned int getRand(MTRand *mtrand){
		return *(element+mtrand->randInt(sAllocate));
	}

	double mean(){
		double mean=0;
		for(int X=0;X<Allocate;X++) mean+=element[X];
		return mean/Allocate;
	}

	unsigned int min(){
		unsigned int min=0;
		min--;
		for(int X=0;X<Allocate;X++) if(element[X]<min) min=element[X];
		return min;
	}

	unsigned int maximum(){
		unsigned int max=0;
		for(int X=0;X<Allocate;X++) if(element[X]>max) max=element[X];
		return max;
	}
};

/* Probability tabels store lookup tables of mathmatical functions, and are similar to function tables with a few additional constraints.
 * the function passed should be the pdf of the desired distribution. It should set the parameters acourding to the set(A) and set(B) functions 
 * of this  library....
 *  
 */

class PROBABILITY_TABLE{
	private:
	unsigned int xr, randnum;

	int allocate(){
	
		element = (unsigned int **) calloc(max, sizeof(unsigned int *));

		if(element == NULL) {
			fprintf(stderr, "error allocating memory\n");
			return 1;
		}

		for(int i = 0; i < max; i++) {
			element[i] = (unsigned int *) calloc(width, sizeof(unsigned int));
			if(element[i] == NULL) {
				fprintf(stderr, "out of memory\n");
				return 1;
			}
		}

		for(unsigned int i=0; i < max; i++){
			for(unsigned int p=0; p<width;p++){
				element[i][p]=i*width+p;
			}
  		}

		for(unsigned int i=0; i < max; i++){
			for(unsigned int p=0; p<width;p++){
				if(element[i][p]!=i*width+p) {printf("Memory managment error"); return 1;}
			}
  		}
		return 0;
	}

	public:

	unsigned int width;
	unsigned int max;

	unsigned int **element;

	PROBABILITY_TABLE(){
		element=NULL;
	}

	int saveTableToFile(char *str){
		FILE* TableFile;
		TableFile=fopen(str, "wb");
		if(TableFile==NULL){
			fprintf(stderr, "cannot open file %s\n", str);
			return 1;
		}
		fwrite (&max, sizeof(unsigned int), 1, TableFile);
		fwrite (&width, sizeof(unsigned int), 1, TableFile);

		for(int i=0; i<max; i++) for(int j=0;j<width;j++) fwrite (&element[i][j], sizeof(unsigned int), 1, TableFile);
		fclose(TableFile);
	}

	void openXorMakeTable(char *str, unsigned int ( *function  ) (unsigned int, unsigned int), void (*clear_function) (void), int Q, int R){
		max=Q;
		width=R;

		fstream file;
		file.open(str,ios::in);
		if( file.is_open() )
		{
			file.close();
			if(makeTableFromFile(str)!=0){
				makeTable(function, clear_function, Q, R);
				saveTableToFile(str);
			}
		}
		else
		{
			file.close();
			makeTable(function, clear_function,Q, R);
			saveTableToFile(str);
		}
	}

	int makeTableFromFile(char *str){
		FILE* TableFile;
		TableFile=fopen(str, "rb");
		if(TableFile==NULL){
			fprintf(stderr, "cannot open file %s\n", str);
			return 1;
		}

		unsigned int Q;
		unsigned int R;

		fread(&Q, sizeof(unsigned int),1,TableFile);
		fread(&R, sizeof(unsigned int),1,TableFile);

		if(max!=Q||width!=R){
			printf("%s :file does not match parameters\n", str);
			return 1;
		}

		if(allocate()!=0){
			fprintf(stderr, "error allocating memory\n");
			return 1;
		}
 
		for(int i=0; i<max; i++) for(int j=0;j<width;j++) fread (&element[i][j], sizeof(unsigned int), 1, TableFile);

		fclose(TableFile);
		return 0;
	}
	
	int makeTable(unsigned int ( *function  ) (unsigned int, unsigned int), void ( *clear_function) (void) , int Q, int R){
		//cout << "make Table ()\n";
		max=Q;
		width=R;
		double local_A=A;

		if(allocate()!=0){
			fprintf(stderr, "error allocating memory\n");
			return 1;
		}

		for(unsigned int i=0; i < max; i++){
			set_A(local_A*((double)(max-i)/double(max)));
			set_B(sqrt(local_A*((double)(max-i)/double(max))));

			clear_function();
			for(unsigned int p=0; p<width;p++){
				element[i][p]=function(p+1,width);
			}
  		}

	}

	void clear(){
		if(element!=NULL){
			for(int i = 0; i < (max); i++) free(element[i]);
			free(element);
			element=NULL;
		}
	}

	~PROBABILITY_TABLE(){
		if(element!=NULL){
			for(int i = 0; i < (max); i++) free(element[i]);
			free(element);
			element=NULL;
		}
	}

	unsigned int getRand(unsigned int x, MTRand *mtrand){
		return *(element[x]+mtrand->randInt(width-1));
	}

	unsigned int getElement(unsigned int x, unsigned int y){
		return element[x][y];
	}

	double mean(unsigned int x){
		double mean=0;
		for(int X=0;X<width;X++) mean+=element[x][X];
		return mean/width;
	}

	unsigned int min(unsigned int x){
		unsigned int min=0;
		min--;
		for(int X=0;X<width;X++) if(element[x][X]<min) min=element[x][X];
		return min;
	}

	unsigned int maximum(unsigned int x){
		unsigned int max=0;
		for(int X=0;X<width;X++) if(element[x][X]>max) max=element[x][X];
		return max;
	}
};

unsigned int getMutation(unsigned int a, int x, unsigned int max){

		long long unsigned int value;

		long double lamda, x_factorial, e_neg_lamda, lambda_x, C;

		lamda=(long double)(A*B);
		
		value=0;
		x_factorial=1;
		for(long double y=1; y<=x; y++) x_factorial=(x_factorial)*y;
		e_neg_lamda=exp((lamda*-1));

		lambda_x=pow(lamda,x);
		C=(lambda_x*e_neg_lamda)/x_factorial;
		return (unsigned int)((long double)(C)*(long double)(MAXINT));
};

/*
unsigned int getRecombination(unsigned int a, unsigned int x, unsigned int max){

		long long unsigned int value;

		long double lamda, x_factorial, e_neg_lamda, lambda_x, C;

		lamda=2.0;
		
		value=0;
		x_factorial=1;
		for(long double y=1; y<=x; y++) x_factorial=(x_factorial)*y;
		e_neg_lamda=exp((lamda*-1));

		lambda_x=pow(lamda,x);
		C=(lambda_x*e_neg_lamda)/x_factorial;
		return (unsigned int)((long double)(C)*(long double)(MAXINT));
}

unsigned int getLn (unsigned int a, unsigned int b, unsigned int scaleA, unsigned int scaleB){		
		
  		double base=(double)(a)/(double)(scaleA);
		double X=(double)(b)/(double)(scaleB);
		return (unsigned int)((log(X)/log(base))*double(scaleA));
};

unsigned int getPow (unsigned int a, unsigned int b, unsigned int scaleA, unsigned int scaleB){

  		double base=(double)(a)/(double)(scaleA);
		double X=(double)(b)/(double)(scaleB);
		return (unsigned int)((pow(base,X))*(double)(scaleA));
};*/

unsigned int getCrossOver(unsigned int a, unsigned int b, unsigned int max)
{
//	cout << "HI!\n";//c << " : " << P << endl;
	double P=0, c=0;
	//m=((double)(max)-
	c=(double)(log(1.0-.01/.5)/(-2.0*A));
	P=(double)(.5*(1-(exp(-2*(double)(a)*c))));
//	cout << c << " : " << P << endl;
	if(b==0)
	  	return (unsigned int)(P*(double)(MAXINT));
	else return (unsigned int)((double)(MAXINT)*(1.0-P))-1;
};

unsigned int getExponential2(unsigned int a, unsigned int b, unsigned int max){
		double X=(double)(a)/(double)(b+1);
		return (unsigned int)(-log(1.0-X)*A);
};

unsigned int getExponential (unsigned int a, unsigned int b){
		double X=(double)(a)/(double)(b+1);
		return (unsigned int)(-log(1.0-X)*A+1);
};

long double last_p, last_x, last_X;

void clear_Poisson (void){
	last_p=0;
	last_x=0;
	last_X=0;
}

double getNormal (unsigned int a, unsigned int b){
	long double p=(double)(a)/(double)(b+1);
	long double z=(2*p-1);
	//cout << "p:" << p << "z:" << z << "=" <<A+B*sqrt(2.0)*0.5*sqrt(M_PI)*(z+M_PI/12.0*pow(z,3)+7*pow(M_PI,2)/480*pow(z,5)+127*pow(M_PI,3)/40320*pow(z,7)+4369*pow(M_PI,4)/5806080*pow(z,9)+34807*pow(M_PI,5)/182476800*pow(z,11)) << endl;
	return A+B*sqrt(2.0)*0.5*sqrt(M_PI)*(z+M_PI/12.0*pow(z,3)+7*pow(M_PI,2)/480*pow(z,5)+127*pow(M_PI,3)/40320*pow(z,7)+4369*pow(M_PI,4)/5806080*pow(z,9)+34807*pow(M_PI,5)/182476800*pow(z,11));
	
	/*const double b0= 0.2316419, b1 = 0.319381530, b2 = -0.356563782, b3 = 1.781477937, b4 = -1.821255978, b5 = 1.330274429;
	double X=0,S;

	if(x>last_x){
		p=last_p;
		X=last_X;
	}
	S=(X-A)/B;
	double t=1.0/(1+S*b0), k=1.0/sqrt(2*M_PI)*exp(-0.5*pow(S,2));
	while (p<x){
		t=1.0/(1+S*b0); k=1.0/sqrt(2*M_PI)*exp(-0.5*pow(S,2));
		if (S>-1.904765 && S<6) p=(1-k*(b1*t+b2*pow(t,2)+b3*pow(t,3)+b4*pow(t,4)+b5*pow(t,5)));
		X++;
		S=(X-A)/B;
	}

	last_p=p;
	last_X=X;
	last_x=x;

	return X-1; */
}

unsigned int getPoisson(unsigned int a, unsigned int b){


		long double x=(double)(a)/(double)(b+1);

		long double p=0;

		long double lambda, k_term=1, e=exp(1.0), pe;

		unsigned int X=0;		

		lambda=A;

		if(x>last_x){
			p=last_p;
			X=last_X;
		}

/*		if((lambda-800)>X){
			X=lambda-800;
			p=0;
		//	last_p=0;
		//	last_X=lambda-101;
		}*/
		if(lambda<500){
			while(p<x){
			k_term=1;
			if(X<lambda/100){
				pe=pow(e,100);
				for (long double y=1;(y<=X);y++) k_term*=lambda/(y*pe);
				p+=(k_term*exp((lambda-X*100)*-1));
			}
			else if(X<lambda/10){
				pe=pow(e,10);
				for (long double y=1;(y<=X);y++) k_term*=lambda/(y*pe);
				p+=(k_term*exp((lambda-X*10)*-1));
			}
			else if(X<lambda/5){
				pe=pow(e,5);
				for (long double y=1;(y<=X);y++) k_term*=lambda/(y*pe);
				p+=(k_term*exp((lambda-X*5)*-1));
			}
			else if(X<lambda/4){
				pe=pow(e,4);
				for (long double y=1;(y<=X);y++) k_term*=lambda/(y*pe);
				p+=(k_term*exp((lambda-X*4)*-1));
			}
			else if(X<lambda/3){
				pe=pow(e,3);
				for (long double y=1;(y<=X);y++) k_term*=lambda/(y*pe);
				p+=(k_term*exp((lambda-X*3)*-1));
			}
			else if(X<lambda/2){
				pe=pow(e,2);
				for (long double y=1;(y<=X);y++) k_term*=lambda/(y*pe);
				p+=(k_term*exp((lambda-X*2)*-1));
			}
			else {
				for (long double y=1;(y<=X);y++) k_term*=lambda/(y*e);
				p+=(k_term*exp((lambda-X)*-1));
			}
			X++;
		}

		//cout << "y:" << 1 << "X:" << X << "dp:" << p-last_p << endl;

		last_p=p;
		last_X=X;
		last_x=x;

		return X-1;
		}
		else {
			B=sqrt(A);
			return getNormal(a,b);
		}
};

/*unsigned int getPoisson2(unsigned int a, unsigned int b, unsigned int max){

		unsigned int X=b;
		long double lambda, k_term=1, p, e=exp(1.0);

		lambda=(1.0-(long double)(a)/(long double)(max))*(long double)(A);
		if(lambda==0){
			if (X==0) return (MAXINT);
			else return 0;
		}

		if(X<lambda/100){
			for (long double y=1;(y<=X&&y<20);y++) k_term*=lambda/(y*pow(e,100));
			p=(k_term*exp((lambda-X*100)*-1));
		}
		else if(X<lambda/10){
			for (long double y=1;(y<=X&&y<20);y++) k_term*=lambda/(y*pow(e,10));
			p=(k_term*exp((lambda-X*10)*-1));
		}
		else if(X<lambda/5){
			for (long double y=1;(y<=X&&y<20);y++) k_term*=lambda/(y*pow(e,5));
			p=(k_term*exp((lambda-X*5)*-1));
		}
		else if(X<lambda/4){
			for (long double y=1;(y<=X&&y<20);y++) k_term*=lambda/(y*pow(e,4));
			p=(k_term*exp((lambda-X*4)*-1));
		}
		else if(X<lambda/3){
			for (long double y=1;(y<=X&&y<20);y++) k_term*=lambda/(y*pow(e,3));
			p=(k_term*exp((lambda-X*3)*-1));
		}
		else if(X<lambda/2){
			for (long double y=1;(y<=X&&y<);y++) k_term*=lambda/(y*pow(e,2));
			p=(k_term*exp((lambda-X*2)*-1));
		}
		else {
			for (long double y=1;(y<=X&&y<T);y++) k_term*=lambda/(y*e);
			p=(k_term*exp((lambda-X)*-1));
		}
		return (unsigned int)(p*(long double)(MAXINT));

};*/
long double choose(unsigned int N, unsigned int i){
	long double term=1;
	for (int x=0;x<i;x++) term*=(double)(N-x)/(double)(i-x);
	return term;
}

unsigned int getBinomial(unsigned int a, unsigned int b){
	long double P=(double)(a)/(double)(b+1);

	int X=0;
	int k=B;
	double p=A;
	double cp=0;
	while (cp<P) if (X<=k) p+=choose(X,k)*pow(p,X)*pow(1-p,X-k);
	return 0;
}

#endif