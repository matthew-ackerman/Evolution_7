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

double Umin=0.001;
int omega=10;
int n=150;

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

double S(int x){
	return pow(1-0.1,x);
};

double U(int x){
	return 0.01;//(double)(x)/double(omega)*10.0+Umin;
}

double UofW(int x, int y, int z){
	if(y<x) return 0;
	else return gsl_ran_poisson_pdf(double(y-x), U(z));
};

double p=0.1;

double UofU(int x, int y){
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

int main (void){
	INDEX index_2(n,omega);
	INDEX index_4(n,omega,n,omega);

	tensor <double> P_even(&index_2);
	tensor <double> P_odd(&index_2);

	unsigned int size=P_even.get_size();

	vector <vector <PAIRS> > even_to_odd(size);
	vector <vector <PAIRS> > odd_to_even(size);

	double R;
	PAIRS Q;
	int q;
	vector <PAIRS>::iterator it, end;
	double *test;

	cout << "starting " << size << "\n";

	for(int w=0;w<n;w++){
		for(int x=0;x<omega;x++) for(int y=w;y<n;y++) for(int z=0;z<omega;z++){
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
	}

	cout << "done\n";
	P_even(0,0)=1;

	for(int t=0;t<1000;t+=2){
		cout << get_Wbar(&P_even) << " : " << get_Ubar(&P_even) << endl;
		reduce(&P_odd,even_to_odd);
		P_odd.normalize();
		cout << get_Wbar(&P_odd) << " : " << get_Ubar(&P_odd) << endl;
		reduce(&P_even,odd_to_even);
		P_even.normalize(); 
	}
	cout << P_odd.print_quite() << endl;
}

