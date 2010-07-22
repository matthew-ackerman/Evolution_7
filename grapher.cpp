//  GLUT CODE BASED ON Song Ho Ahn's demo (song.ahn@gmail.com);

#define ESCAPE 27

#include <stdlib.h> //Needed for "exit" function
#include <sstream>
#include <iomanip>
#include <unistd.h>//Needed for _exit

#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include "lib/tensor.h"
#include "lib/graphics.h"
#include <gsl/gsl_randist.h>
#include <sys/wait.h>

#include <math.h>
#include <time.h>
#include <iostream>
#include "lib/tensor.h"
#include <vector>
#include "lib/tables_v1.0.h"
#include "lib/MersenneTwister.h"
#include <algorithm>
#include <string.h>
#include "time.h"
#include <gsl/gsl_randist.h>

using std::stringstream;
using std::cout;
using std::endl;
using std::ends;

using namespace std;


void init(void);
void display();
void keyboard(unsigned char, int, int);

double Umin=0.5;
int omega=20;
int n=20;

class PAIRS {
	public:
	double * number, weight;
	PAIRS (){};

	PAIRS (double *n, double w){
		number=n;
		weight=w;
	};
	void operator() (double *n, double w){
		number=n;
		weight=w;
	}
	~PAIRS(){};
};

void reduce (tensor <double> *p, vector <vector <PAIRS> > function){
		double *p_it, *p_end;

		vector < vector <PAIRS> >::iterator outer_it, start_outer, end_outer;
	 	vector <PAIRS>::iterator inner_it, start_inner, end_inner;

		start_outer=function.begin();
		end_outer=function.end();
		outer_it=start_outer;

		p_it=p->begin;
		p_end=p->end;

		while(p_it!=p_end){
			inner_it=outer_it->begin();
			end_inner=outer_it->end();
			(*p_it)=0;
			while (inner_it!=end_inner){
				(*p_it)+=(*(inner_it->number))*(inner_it->weight);
				inner_it++;
			}
			p_it++;
			outer_it++;
		}
};

tensor <double> P_even, P_odd;
tensor <unsigned int> N_tensor;
gsl_rng *gslrand;


vector <vector <PAIRS> > even_to_odd;
vector <vector <PAIRS> > odd_to_even;

void drift(tensor <double> *t, tensor <unsigned int> *K, unsigned int N){
	double *it=t->begin, *end=t->end;
	unsigned int *Ns=K->begin;
	//gsl_ran_multinomial(gslrand, t->get_size(), N, it, Ns);
	while(it!=end){
		*it=1/t->get_size();//(double)(*Ns)/(double)(N);
		it++;
		Ns++;
	}
}

void time_step(){
	reduce(&P_odd,even_to_odd);
	P_odd.normalize();
	drift(&P_odd, &N_tensor, 1000);
	reduce(&P_even,odd_to_even);
	P_even.normalize(); 	
	drift(&P_even, &N_tensor, 1000);
};

void update(int value){
	reduce(&P_odd,even_to_odd);
	P_odd.normalize();
	reduce(&P_even,odd_to_even);
	P_even.normalize(); 	
	glutPostRedisplay();
	glutTimerFunc(2, tick, 0);
};

double S(int x){
	return pow(1-0.1,x);
};

double U(int x){
	return (double)(x)/double(omega)*1.0+Umin;
}

double UofW(int x, int y, int z){
	if(y<x) return 0;
	else return gsl_ran_poisson_pdf(double(y-x), U(z));
};

double UofU(int x, int y){
	double p=(1-pow(0.5,x+Umin))/2;
	if(abs(y-x)>1) return 0;
	if(x!=0 && x!=omega){
		if(x!=y) {
			return p;
		}
		else return (1-2*p);
	}
	else{
		if(x!=y) {
			return p;
		}
		else return (1-p);
	}
};

double get_Wbar(tensor <double> *t){
	double SUM=0;
	for(int x=0;x<n;x++) for(int y=0;y<omega;y++) SUM+=(*t)(x,y)*S(x);
	return SUM;
}

double get_Ubar(tensor <double> *t){
	double SUM=0;
	for(int x=0;x<n;x++) for(int y=0;y<omega;y++) SUM+=(*t)(x,y)*U(y);
	return SUM;
}

int main (int argc, char **argv){
	#define X_SIZEc	101
	#define Y_SIZEc	1

	gsl_rng_env_setup();
	const gsl_rng_type * T;     
       	T = gsl_rng_default;
       	gslrand = gsl_rng_alloc (T);

	INDEX index_2(n,omega);
	INDEX index_4(n,omega,n,omega);

	P_even.resize(&index_2);
	P_odd.resize(&index_2);
	N_tensor.resize(&index_2);


	unsigned int size=P_even.get_size();

	even_to_odd.resize(size);
	odd_to_even.resize(size);

	double R;
	PAIRS Q;
	int q;
	vector <PAIRS>::iterator it, end;
	double *test;

	setDrawTensor(&P_even);

	cout << "starting " << size << "\n";

	for(int w=0;w<n;w++) for(int x=0;x<omega;x++) for(int y=0;y<n;y++) for(int z=0;z<omega;z++){
			R=S(w)*UofW(w,y,x)*UofU(x,z);
			if(R>0) {
				q=y+z*(omega);

				it=even_to_odd[q].begin();
				end=even_to_odd[q].end();

				test=&P_even(w,x);

				Q(test,R);

				while(it!=end && it->number != test) it++;
				if(it==end) even_to_odd[q].push_back(Q);
				else (it->weight)+=R;

				it=odd_to_even[q].begin();
				end=odd_to_even[q].end();

				test=&P_odd(w,x);

				Q(test,R);

				while(it!=end && it->number != test) it++;
				if(it==end) odd_to_even[q].push_back(Q);
				else (it->weight)+=R;

			}
	}

	cout << "done\n";
	P_even(0,0)=1;

	glutInit(&argc, argv);
	setUpdateTensor(time_step);
	glutTimerFunc(2,update, 0);
	Start();
}
