/* A very general population genetics simulation, by Matt Ackerman.
 */

#include "stdafx.h"
#include <math.h>
#include <time.h>
#include <iostream>
#include "lib/tensor.h"
#include <vector>
#include "lib/tables_v1.0.h"
#include "lib/MersenneTwister.h"
#include <algorithm>
#include <string.h>
#include "lib/bin.h"
#include "lib/evo7.h"

double const_function (double H, double h, vector <double> arg){
	return arg[0];
};

double generic_fitness_function (double H, double h, vector <double> arg) {
	return 1-(pow(arg[0], H)*pow(arg[1], h));
};

double generic_one_way_mutation_function (double A, double Adot, vector <double> arg) {
	if( (A+1==Adot || A==Adot ) && A!=arg[1] ) return arg[0];
	else return 0;
};

double generic_two_way_mutation_function (double A, double Adot, vector <double> arg) {
	if( (A+1==Adot || A==Adot+1) && (A!=arg[1] || A!=0 ) ) return arg[0];
	if (A==Adot && (A!=arg[1] || A!=0 ) ) return arg[0]*2;
	if (A==Adot && (A==arg[1] || A==0 ) ) return arg[0];
	else return 0;
};


double multiplicative (double A, double B){
	return 1.0-(1.0-A)*(1.0-B);
};

double wiered_multiplicative (double A, double B){
	if (A!=0){
		if (B!=0) return A+B;
		else return A;
	}
	else return B;
};

using namespace std;
using std::cout;

double rate (locus *A, locus *B){
	if(A->chromosome!=B->chromosome){
		return 0.5;
	}
	else {
		return 0.5*(1-exp(-2*(double)(fabs((double)(A->cM_from_p0)-(double)(B->cM_from_p0))/100.0)));
	}
}

class population{
	private:
	/*recombination is a vector of vectors. recombination[X] is a vector of pointers to the numbers which will...*/

	vector <vector <PAIRS> > transition;
	double *it, *end_minus_one, *begin, *end;

	public:

	INDEX *gamete_index, *zygote_index;

	tensor *haplotype_density;
	tensor *genotype__density;
	tensor *selection_tensor;
	tensor *mutation_tensor;

	double mean_W;
	double N;
	unsigned int t;

	unsigned int number_of_genotypes;
	
	MTRand *mtrand; DOUBLE_TABLE *normal;

	population(vector <locus*> loci, double _N, tensor *t_selection_tensor, tensor *t_mutation_tensor, MTRand *_mtrand, DOUBLE_TABLE *_normal ){

		selection_tensor=t_selection_tensor;
		mutation_tensor=t_mutation_tensor;

		sort(loci.begin(), loci.end(), compare_loci);

		int sum=0;

		N=_N;

		mtrand=_mtrand;
		normal=_normal;

		t=0;

		gamete_index=new INDEX (loci.size());
		zygote_index=new INDEX(loci.size()*2);

		for(int x=0; x<loci.size(); x++){
			(*gamete_index)[x]=(loci[x]->number_of_alleles);
			(*zygote_index)[x]=(loci[x]->number_of_alleles);
			(*zygote_index)[x+loci.size()]=(loci[x]->number_of_alleles);
		}

		gamete_index->make_step_array();
		zygote_index->make_step_array();

		haplotype_density=new tensor(gamete_index);
		end=haplotype_density->end;
		end_minus_one=end-1;
		begin=haplotype_density->begin;		

		genotype__density=new tensor(zygote_index);

		INDEX *paternal_index=new INDEX(loci.size()), *maternal_index= new INDEX(loci.size()), *pre_mutation_index=new INDEX(loci.size()), *post_mutation_index=new INDEX(loci.size());

		(*paternal_index)=(*gamete_index);
		(*maternal_index)=(*gamete_index);
		(*pre_mutation_index)=(*gamete_index);
		(*post_mutation_index)=(*gamete_index);


		number_of_genotypes=haplotype_density->get_size();

		for(int x=0;x<number_of_genotypes; x++) transition.push_back(*new vector <PAIRS>);

		unsigned int number_of_loci=loci.size();

		for(int X=0;X<number_of_genotypes; X++) {
		for(int y=0;y<number_of_genotypes; y++) {

		for(int x=0;x<number_of_genotypes; x++) {
		for(int Y=0;Y<number_of_genotypes; Y++) {

			paternal_index->setFromOffset(x);
			maternal_index->setFromOffset(y);
			pre_mutation_index->setFromOffset(X);
			post_mutation_index->setFromOffset(Y);
			zygote_index->setFromOffset(x+y*number_of_genotypes);

			double R=1; /*the rate at which genotype y,z gives rise to haplotype x */
				    /*R is mendelian rate of production of n, M is rate which n becomes x considering mutation */

			/*first "we need to calculate the rate at which gamete X is made from parent x/y before mutaiton*/;
			for(int w=0;w<number_of_loci;w++){
				if(R!=0){
					bool from_mom;
					unsigned int H;					
					if((*maternal_index)[w]==(*pre_mutation_index)[w] || (*paternal_index)[w]==(*pre_mutation_index)[w]){
						if((*maternal_index)[w] != (*paternal_index)[w]){
							if(R==1){
								R*=0.5;
								if((*maternal_index)[w]==(*pre_mutation_index)[w]) from_mom=true;
								else from_mom=false;
								H=w;
							}
							else {
								if((*maternal_index)[w]==(*pre_mutation_index)[w]){
									if(from_mom) {R*=(1-rate(loci[w], loci[H])); H=w;}
									else { R*=rate(loci[w],loci[H]); H=w; from_mom=false;}
								}
								else {
									if(!from_mom) {R*=(1-rate(loci[w], loci[H])); H=w;}
									else { R*=rate(loci[w],loci[H]); H=w; from_mom=true;}
								}
							}
						}
					}
					else{								
						R=0;
					}
				}		
			}
			double M=1, K=1+(*(*mutation_tensor)(zygote_index).it);

			/*now, we need to calculate the rate at which the pre_mutaion is converted to a post_mutation*/

			for(int w=0;w<number_of_loci;w++){
				if( (*pre_mutation_index)[w]==(*post_mutation_index)[w] ) M*=(1-K*loci[w]->mutation_rate[(*pre_mutation_index)[w]][(*post_mutation_index)[w]]);
				else M*=loci[w]->mutation_rate[(*pre_mutation_index)[w]][(*post_mutation_index)[w]]*K;
			}

			double S=1+(*(*selection_tensor)(zygote_index).it);

		 	vector <PAIRS>::iterator it, end;
			double * g_it;

			if(R*M*S!=0) {
				it=transition[Y].begin();
				end=transition[Y].end();
				g_it=(*genotype__density)(zygote_index).it;
				while(it!=end && it->number != g_it) it++;
				if(it==end) transition[Y].push_back(*new PAIRS(g_it,R*M*S));
				else (it->weight)+=(R*M*S);
			}
		}
		}
		}
		}
		skip=false;
	}

	/* needs to be called if the haplotype tensor chages, but currently doesn't reset */

	void reset_ptr(void){
		end=haplotype_density->end;
		end_minus_one=end-1;
		begin=haplotype_density->begin;		
	}

	long double sum;
	double *skip_to;
	bool skip;

	void finite_iterate(void){
			genotype__density->product(haplotype_density,haplotype_density);
			reduce(haplotype_density,transition);
			haplotype_density->normalize();

			it=begin;
			sum=0;
			long double K;
			long double Q;
			/* If only one haplotype exist, the time to the next mutation is exponentially distributed. Since this is a markov
			 * process, the transition probability depends only on the current state of the population, so that we can use the 
			 * waiting time (which is always exponentially distributed) to tell us the time of the next transition, and then see
			 * which state it transitions to.   
			 */
			if(skip) {
				t-=(log((*mtrand)())*pow((long double)(10),8)/N);
				skip=false;
				sum+=*skip_to;
				K=1;
				while(it!=end_minus_one){
					if(it!=skip_to){
						if (K>0) {
							Q=(*it)/(long double)((long double)(1.0)-sum);
							sum+=(*it);
							//this is kinda stupid. Why do I devide by (2*N) then multiply by (2*N)? ...
							(*it)=ignbin(K, Q ,mtrand)/(2*N);
							K-=(*it)*(2*N);
						}
						else *it=0;
					}
					it++;
				}
				(*it)=K/(N*2);
			}
			else if(!skip){
				K=(N*2);
				while(it!=end_minus_one){
					/* The commented code that follows is attempts to use the normal and poisson aproximations to speed up the 						 * simulation. Unfortunately I've been having problems with this, so it is currently commented out.
					*/

					/*if(K > 500){
						if (K-Q< 500 ){
							sum+=(*it);
							cout << "(L poiss) E(left-x)=" << K-(*it)*N << " K: " << K << " Q: " <<  Q;
							long double L=exp(-(K-Q)), k=0, p=1;
							while(p>L){p*=(*mtrand)(); k++;}
							if(k==0) k++;
							//cout << k << 
							(*it)=(double)(K-(k-1))/(double)(N);
							K-=(*it)*N;
							cout << " O(x)=" << (*it)*N << " , left " << K << "\n";
						}
						else if(Q<500){
							sum+=(*it);
							cout << "(R poiss) E(x)=" << (*it)*N;
							long double L=exp(-Q), k=0, p=1;
							while(p>L){p*=(*mtrand)(); k++;}
							if(k==0) k++;
							(*it)=(k-1)/(long double)(N);
							K-=(*it)*N;
							cout << " ~= " << K-Q << " O(x)=" << (*it)*N << " , left " << K << "\n";
						}
						else {
							if(Q>N/2){
								sum+=(*it);
								cout << "(L normal) E(x)=" << (*it)*N;
								(*it)+=(long double)(normal->getRand(mtrand) *sqrt( N-Q ) / N);
								K-=(*it)*N;
								cout << " ~= " << Q << " O(x)=" << (*it)*N << " , left " << K << "\n";
							}
							else{
								sum+=(*it);
								cout << "(S normal) E(x)=" << (*it)*N;
								(*it)+=(long double)(normal->getRand(mtrand) *sqrt( Q ) / N);
								K-=(*it)*N;
								cout << " ~= " << Q << " O(x)=" << (*it)*N << " , left " << K << "\n";
							}
						}
					}*/

					if (K>=0) {
						/* K here is simply the number of samples we will be taking from a binomial distribution (i.e. N minus the 							 * draws we have alread made. If K is large, we can use an normal or poisson distriubiton (depending on 
						 * whether or not K*p is large). Q is simply the probability of X, given that we do not draw Y	 							 * (where Y is every haplotype we have already iterated through).   
						 */
	
						Q=(*it)/(long double)((long double)(1.0)-sum);
						sum+=(*it);
						(*it)=(double)(ignbin(K, Q, mtrand))/(2*N);
						K-=(*it)*(N*2);
						if(*it==1) {skip=true;skip_to=it;}
					}
					else *it=0;
					it++;
				}
				/* K may be sligtly less (or more) than a whole number (because of floating point magic). So, we add a small amount to it 
				 * and round down. */
				(*it)=(int)(K+0.01)/(2*N);
			}
			t++;
	}
	double get_w(void){
		double *d_it=genotype__density->begin;
		double *d_end=genotype__density->end;
		double *s_it=selection_tensor->begin;
		sum=0;
		while(d_it!=d_end){
			sum+=(*d_it)*(*s_it);
			d_it++;
			s_it++;
		}
		return sum;
	}
	/*locus, allele: This function is quite slow.*/
	double get_allele_freq(int x, int y){
		it=begin;
		sum=0;
		for(int z=0;z<number_of_genotypes;z++){
			gamete_index->setFromOffset(z);
			if((*gamete_index)[x]==y) sum+=(*it);
			it++;
		}
		return sum;
	}
	/* iterates deterministically */
	void infinite_iterate(void){
		genotype__density->product(haplotype_density,haplotype_density);
		reduce(haplotype_density,transition);
		haplotype_density->normalize();
		t++;
	}
	~population(){	
//		gamete_counts->kill(); 
	}
};

int main(int argc, char *argv[]){

	/* the syntax (class)(varaible) is a c++ style 'cast' operation, which tells the compiler to convert the variable to the class 
	 * (class).
	 */

	long double N=pow((long double)(10),3);

	DOUBLE_TABLE normal;

	/* sets the mean of a statstical distribution */
	set_A(0);

	/* set the variance, if not determined by the mean */
	set_B(1);

	/* makes a file in memory for generation normally distributed RV with mean set_A and variance set_B; 	*
	 * here I am making a standard normal distribution 							*/
	normal.makeTable(getNormal, 10000);

	/* a variable for quickly generating uniformly distributed RV */ 
	MTRand mtrand;

	/* a vector of all loci */ 
	vector <locus*> loci;

	double s1=0.0, s2=0.0, s3=0.00;
	double h1=0.5, h2=0.5, h3=0.5;
	double mh1=0.5, mh2=0.5, mh3=0.0;
	double m1=0.0, m2=0*pow((long double)(10),2), m3=0.0;

	/* 
	 * Declaring tensors: tensors are simply an N dimensional extension of the concept of a matrix. They can be declared in a number of ways. Either 	  * with an index, or with the style used here, where the rank (i.e. dimensions) of the tensor is declared, the dimension (i.e. width) of each rank, 		 * and the value of each cell. The seletion and mutation tensors list should always have a rank (i.e. dimensions) of number_of_loci * 2,
	 * and dimensions of the number of alleles at the X th or X-rank th locus (iff x>=rank).     
	 */


	//Empty tensor declared. These guys will segfault if you do almost anything to them except assign a tensor; 
	tensor selection_tensor;
	tensor mutation_tensor;

	/* 
	 * Adding a new locus: When a locus is declared you choose how many alleles exist at that locus, which chromosome the locus resides on, where on 
	 * the chromosome it resides, the base mutation rates to other alleles at the loci, the fitness of those alleles, and the effect those alleles have
	 * on the mutation rate.
	 *
	 * locus (number_of_loci, chromosome, cM, allelic_mutation_rates) -or- locus (number_of_loci, chromosome, cM, allelic_mutation_rates,
	 * selection_matirx, mutation_matrix) -or- locus (number_of_loci, chromosome, cM, allelic_mutation_rates, &selection_function, &mutation_function); 
	 * 
	 * The matrix for allelic_mutation_rates is not a stocastic matrix. Instead, the element [A][B] should be the rate per generation that x->y,
	 * for x=A y=B and A!=B. iff A==B, then element [A][A] or [B][B] or whatever should be the sum of the rates of x->y for all x=A and all y!=x.
	 * 
	 * The values of each element can be set manually, as in the fist example below, or using a function pointer, as bellow. 
	 * To pass a function pointer take a look at the function pointers declared at the begining of this file.
	 */

	/* push_back adds an element to the end of a vector. */

	/* functors are just functions which store some local variables, in this case the local variables are actually saved in an argument vector 
	 * in the functor class, and we just need to pass the functor class a function pointer */

	functor allelic_mutation_functor, selection_functor, mutator_functor;

	vector <double> args;

	/* each functor takes 1,2,3 or 4 doubles, and a vector of arguments. The vector arguments is a private static member of the functor class 
	 * which is passed to the function pointed to be the function pointer each time it is called. This may not be the fastet way to implement 
	 * the functors, but they are only called durring set up, so I'm not that worried. Technically, I think it should be possible to call the 
	 * constructors w/o an explicit args vector, but I haven't figured out the syntax yet. */

	args.push_back(pow((long double)(10),-8));
	args.push_back(2.0);


	allelic_mutation_functor.set_2(&generic_one_way_mutation_function, args);
	args.clear();

	args.push_back(1-s1);
	args.push_back(1-s1*h1);
	selection_functor.set_2(&generic_fitness_function,args);
	args.clear();

	args.push_back(0);
	mutator_functor.set_2(&const_function,args);
 	args.clear();

	loci.push_back (new locus ("locus_1", 2, 0, 0.1, 
			/*allelic_mutation_rate*/
				allelic_mutation_functor,
			/*selection matrix*/ 
				selection_functor,
			/*mutator matrix*/
				mutator_functor
								) );
	loci.push_back (new locus ("locus_2", 2, 0, 0.1, 
			/*allelic_mutation_rate*/
				allelic_mutation_functor,
			/*selection matrix*/ 
				selection_functor,
			/*mutator matrix*/
				mutator_functor
								) );

	/* code like this can be used to save all of your loci to a file*/

	/*
	ofstream myfile_out;
	myfile_out.open ("loci.txt");
	myfile_out << *(loci[0]) << *(loci[1]) << *(loci[2]) << endl;
	myfile_out.close();
	*/
	
	/*and code like this can be used to read your loci in from a file*/

	/*
	ifstream myfile_in;
	myfile_in.open ("loci.txt");

	locus new_locus_1("locus_4", 2, 0, 0.1, (double []) {0.1, 0.1, 0.1, 0.1}, (double []) {0.1, 0.1, 0.1, 0.1},(double []) {0.1, 0.1, 0.1, 0.1} );
	locus new_locus_2("locus_5", 2, 0, 0.1, (double []) {0.1, 0.1, 0.1, 0.1}, (double []) {0.1, 0.1, 0.1, 0.1},(double []) {0.1, 0.1, 0.1, 0.1} );
	locus new_locus_3("locus_6", 2, 0, 0.1, (double []) {0.1, 0.1, 0.1, 0.1}, (double []) {0.1, 0.1, 0.1, 0.1},(double []) {0.1, 0.1, 0.1, 0.1} );

	myfile_in >> new_locus_1 >> new_locus_2 >> new_locus_3;
	*/

	/*This function takes your list of loci and an pointer to an epistatic function as arguments
	 *to auto-magicly make the neccisary tensors*/

	/*I had to switch to the somewhat more awkward declearation of function(tensor *destination, args) rather than 
	 * destination=function(args) because the vissual c++ compiler seems to be automagically deconstructing the return 
	 * tensor before I have had a change to assign the values from it to our destination tensor. I have no idea why it 
	 * does this, nor how to make it stop. */
	make_mutation_tensor(&mutation_tensor, loci, wiered_multiplicative);
	make_selection_tensor(&selection_tensor, loci, multiplicative);

	/* populations are called with a vector of all loci and the population size */
	population pop(loci, N, &selection_tensor, &mutation_tensor, &mtrand, &normal);
	population pop_inf(loci, N, &selection_tensor, &mutation_tensor, &mtrand, &normal);


	/* ?? */
	//printf("T,abc,abC,aBc,aBC,Abc,AbC,ABc,ABC\n");
	double P=0, F=0, WTF=0;

	(*pop_inf.haplotype_density).clear();
	*(*pop_inf.haplotype_density)(0,0).it=1-1/(2*N);
	*(*pop_inf.haplotype_density)(0,1).it=1/(2*N);
	
	for(int j=0;j<5000;j++){
		pop_inf.infinite_iterate();
	}

	for(int x=0;x<10000;x++){

		//*(pop.haplotype_density)=*(pop_inf.haplotype_density);

		(*pop.haplotype_density).clear();
		(*pop.haplotype_density)(0,0)=1-1/(2*N);
		(*pop.haplotype_density)(0,1)=1/(2*N);

		pop.t=0;
		double Q=1/(2*N);
		//(*(*pop.haplotype_density)(6,0,0).it)+(*(*pop.haplotype_density)(6,0,1).it)+(*(*pop.haplotype_density)(6,1,0).it)+(*(*pop.haplotype_density)(6,1,1).it) < 0.99;

		while(Q>0 && Q<1){
		//	pop_inf.infinite_iterate();
			pop.finite_iterate();
		//	if(pop_inf.t==10000) 
			Q=pop.get_allele_freq(1,1);
		//	cout << pop.t << "," << Q << endl;
		//	cout << "---===" << pop.t << "===---" << endl <<(pop.haplotype_density)->print(N*2) << endl;
		//	cout << endl;
		}
		if(Q==0) P++;
		else if(Q==1) F++;
		else {
			WTF++;
			cout << "WTF " << WTF << "/" << F+P+WTF << " ?!?! " << Q << endl;
		}
	//	cout << pop.t << "\t " << sum/(double)(x) << "\t " << (sum/x-(sum-pop.t)/(x-1)) << endl;
	}
	cout << 1/(2*N) << "->" << F/(F+P) << endl;
	//myfile_out.open ("data.txt");
	return 0;
}
