#ifndef EVO7_H
#define EVO7_H

#ifndef DATA
#include "data_record.h"
#endif

#include <list>

#ifndef _TENSOR_HPP
#include "tensor.h"
#endif


#ifdef USE_EVO7_DEFS

#define _number_of_arguments	6
#define _number_of_alleles	0
#define	_dominance		1
#define	_selection_coefficient	2
#define	_mutation_rate		3
#define	_mutator_effect		4
#define	_user_argument		5

#endif

using namespace std;

class functor{
	private:
		double (* function_4)(double, double, double, double, vector <double>);
		double (* function_3)(double, double, double, vector <double>);
		double (* function_2)(double, double, vector <double>);
		double (* function_1)(double, vector <double>);
		vector <double> args;
	public:
		functor(){};
		~functor(){};
		double operator() (double a){
			return function_1(a, args);
		};
		double operator() (double a, double b){
			return function_2(a, b, args);
		};		
		double operator() (double a, double b, double c){
			return function_3(a, b, c, args);
		};		
		double operator() (double a, double b, double c, double d){
			return function_4(a, b, c, d, args);
		};
		void set_1(double (* function) (double, vector <double>), vector <double> _args){
			args=_args;
			function_1=function;
		};
		void set_2(double (* function) (double,double, vector <double>), vector <double> _args){
			args=_args;
			function_2=function;
		};
		void set_3(double (* function) (double,double,double, vector <double>), vector <double> _args){
			args=_args;
			function_3=function;
		};
		void set_4(double (* function) (double, double, double, double, vector <double>), vector <double> _args){
			args=_args;
			function_4=function;
		};
};

/*arg:
	0 	locally used variable,
	1	locally used variable,
	2	# of alleles,
	3	S
	4	s
	5	h
	6 
*/

class locus{
	private :

	public :
	char *name;
	int number_of_alleles;
	/*mutation rate is a conditional rate*/
	double **mutation_rate;
	double **selection;
	double **mutation_kicker;

	int chromosome;
	double cM_from_p0;

	locus(char _name [], int N, int C, double M ,double mutation_array []){
		name=_name;
		number_of_alleles=N;
		chromosome=C;
		cM_from_p0=M;

		mutation_rate=new double *[N];
		selection=new double *[N]; 
		mutation_kicker=new double *[N]; 

		for(int x=0;x<N;x++){
			mutation_rate[x]=new double[N];
		}
		for(int x=0;x<N;x++){
			for(int y=0;y<N;y++){
				mutation_rate[x][y]=mutation_array[x*N+y];
			}
		}
	}

	locus(char _name [], int N, int C, double M ,double mutation_array [], double trait_array [], double kicker_array []){
		name=_name;
		number_of_alleles=N;
		chromosome=C;
		cM_from_p0=M;

		mutation_rate=new double *[N];
		selection=new double *[N]; 
		mutation_kicker=new double *[N]; 

		for(int x=0;x<N;x++){
			mutation_rate[x]=new double[N];
			selection[x]=new double[N];
			mutation_kicker[x]=new double[N];
		}
		for(int x=0;x<N;x++){
			for(int y=0;y<N;y++){
				mutation_rate[x][y]=mutation_array[x*N+y];
				selection[x][y]=trait_array[x*N+y];
				mutation_kicker[x][y]=kicker_array[x*N+y];
			}
		}
	}

	locus(char _name [], int N, int C, double M, functor mutation_function, functor trait_function, functor kicker_function){
		name=_name;
		number_of_alleles=N;
		chromosome=C;
		cM_from_p0=M;

		mutation_rate=new double *[N];
		selection=new double *[N]; 
		mutation_kicker=new double *[N]; 

		for(int x=0;x<N;x++){
			mutation_rate[x]=new double[N];
			selection[x]=new double[N];
			mutation_kicker[x]=new double[N];
		}
		for(int x=0;x<N;x++){
			for(int y=0;y<N;y++){
				mutation_rate[x][y]=mutation_function(x,y);
				selection[x][y]=trait_function(x,y);
				mutation_kicker[x][y]=kicker_function(x,y);
			}
		}
	}
	locus(locus &l){
		name=l.name;
		number_of_alleles=l.number_of_alleles;
		chromosome=l.chromosome;
		cM_from_p0=l.cM_from_p0;

		mutation_rate=new double *[number_of_alleles];
		selection=new double *[number_of_alleles]; 
		mutation_kicker=new double *[number_of_alleles]; 

		for(int x=0;x<number_of_alleles;x++){
			mutation_rate[x]=new double[number_of_alleles];
			selection[x]=new double[number_of_alleles];
			mutation_kicker[x]=new double[number_of_alleles];
		}
		for(int x=0;x<number_of_alleles;x++){
			for(int y=0;y<number_of_alleles;y++){
				mutation_rate[x][y]=l.mutation_rate[x][y];
				selection[x][y]=l.selection[x][y];
				mutation_kicker[x][y]=l.mutation_kicker[x][y];
			}
		}
	}
	~locus(){
		if (mutation_rate != NULL) {
			for(int x=0;x<number_of_alleles;x++) delete [] mutation_rate[x];
			delete [] mutation_rate;
			mutation_rate=NULL;
		}	
		if (selection != NULL) {
			for(int x=0;x<number_of_alleles;x++) delete [] selection[x];
			delete [] selection;
			selection=NULL;
		}	
		if (mutation_kicker != NULL) {
			for(int x=0;x<number_of_alleles;x++) delete [] mutation_kicker[x];
			delete [] mutation_kicker;
			mutation_kicker=NULL;
		}	
	}
	friend ostream &operator<<(ostream &stream, const locus &l);
	friend istream &operator>>(istream &stream, locus &l);
};

ostream &operator<<(ostream &stream, const locus &l){
	stream << "(locus:" << endl << endl;
	stream << "(name: "<< l.name << ")" << endl << endl;
	stream << "(number: " << l.number_of_alleles << ")" << endl << endl;
	if (l.mutation_rate != NULL){
		stream << "(mutation_rate:" << endl;
		for(int x=0;x<l.number_of_alleles;x++){
			stream << l.mutation_rate[x][0]; 
			for(int y=1;y<l.number_of_alleles;y++) stream << ", " << l.mutation_rate[x][y];
			if (x!=l.number_of_alleles-1) stream << endl;
		}
		stream << ")" << endl << endl;
	}
	if (l.selection != NULL){
		stream << "(selection_matrix:" << endl;
		for(int x=0;x<l.number_of_alleles;x++){
			stream << l.selection[x][0]; 
			for(int y=1;y<l.number_of_alleles;y++) stream << ", " << l.selection[x][y];
			if (x!=l.number_of_alleles-1) stream << endl;
		}
		stream << ")" << endl << endl;
	}
	if (l.mutation_kicker != NULL){
		stream << "(mutation_kicker:" << endl;
		for(int x=0;x<l.number_of_alleles;x++){
			stream << l.mutation_kicker[x][0]; 
			for(int y=1;y<l.number_of_alleles;y++) stream << ", " << l.mutation_kicker[x][y];
			if (x!=l.number_of_alleles-1) stream << endl;
		}
		stream << ")" << endl << endl;
	}
	stream << "end_locus)" << endl << endl;
	return stream; // must return stream
}

ifstream &operator>>(ifstream &stream, locus &l){
	char buffer[256], buffer_2[256];
	char *start, *stop;
	std::streamoff g;
	int size, x, y;

	int state=0;
	while (state!=7){
		switch ( state ){
			case 0: 
				stream.getline(buffer,256);
				if(strstr(buffer,"locus")!=NULL) state=1;
			break;
			case 1: 
				stream.getline(buffer,256);
				if(strstr(buffer,"name")!=NULL){
					start=strstr(buffer,":");
					start++;
					while(*start==' ' && *start!= '\0' && *start!='\n') start++;
					stop=start;
					while(*stop!=' ' && *stop!='\0' && *stop!='\n' && *stop!=')') stop++;
			//		delete [] l.name;
					size=stop-start;
					l.name=new char [size+1];
					memcpy (l.name,start,size+1);
					l.name[size]='\0';

					state=2;
				}
			break;
			case 2:
				stream.getline(buffer,256);
				if(strstr(buffer,"number")!=NULL){
					start=strstr(buffer,":");
					start++;
					while(*start==' ' && *start!= '\0' && *start!='\n') start++;
					stop=start;
					while(*stop!=' ' && *stop!='\0' && *stop!='\n' && *stop!=')') stop++;
					size=stop-start;
					if(size<256) memcpy (buffer_2,start,size+1);
					buffer_2[size]='\0';
					l.number_of_alleles=atoi(buffer_2);

					state=3;
					l.mutation_rate=new double *[l.number_of_alleles];
					l.selection=new double *[l.number_of_alleles];
					l.mutation_kicker=new double *[l.number_of_alleles];

					for(int q=0;q<l.number_of_alleles;q++){
						l.mutation_rate[q]=new double [l.number_of_alleles];
						l.selection[q]=new double [l.number_of_alleles];
						l.mutation_kicker[q]=new double [l.number_of_alleles];
					}
				}
			break;
			case 3:
				g=stream.tellg();
				stream.getline(buffer,256);
				if(strstr(buffer,"mutation_rate")!=NULL){
					x=stream.gcount();
					y=0;
					state=4;
					//stream.seekg(g);
				}
				if(strstr(buffer,"selection_matrix")!=NULL){
					x=stream.gcount();
					y=0;
					state=5;
//					stream.seekg(g);
				}
				if(strstr(buffer,"mutation_kicker")!=NULL){
					x=stream.gcount();
					y=0;
					state=6;
//					stream.seekg(g);
				}
				if(strstr(buffer,"end_locus")!=NULL) state=7;
			break;
			case 4:
			case 5:
			case 6:
				if (x==stream.gcount()) {
					stream.getline(buffer,256);
					x=0;
				}
				start=buffer+x;
				while(*start==' ' || *start==',') start++;
				if(*start==',') start++;
				stop=start;
				while(*stop!=' ' && *stop!='\0' && *stop!='\n' && *stop!=')' && *stop!=',') stop++;
				size=stop-start;
				if(size<256){
					memcpy (buffer_2,start,size+1);
					buffer_2[size]='\0';
				}
				x+=(stop-buffer-x);
				if(size==0) x++;
				int q=y/l.number_of_alleles;
				int z=y-l.number_of_alleles*q;

				if(state==4 && size!=0) l.mutation_rate[q][z]=atof(buffer_2);
				else if (state==5 && size!=0) l.selection[q][z]=atof(buffer_2);
				else if (state==6 && size!=0) l.mutation_kicker[q][z]=atof(buffer_2);
				if(size!=0) y++;
				if(*stop==')') state=3;

			break;

		}
		cout << "read {" << buffer << "} state=" << state << endl;
	}
//		size=stream.gcount();
//		l.name=new char [size];
/*	
		sprintf(l.name, "%s", buffer);	

	stream.getline(buffer,256);
	stream.getline(buffer,256);

	l.number_of_alleles=atoi(buffer);

	l.mutation_rate=new double *[l.number_of_alleles];
	l.selection=new double *[l.number_of_alleles];
	l.mutation_kicker=new double *[l.number_of_alleles];

	for(int x=0;x<l.number_of_alleles;x++){
		l.mutation_rate[x]=new double[l.number_of_alleles];
		l.selection[x]=new double[l.number_of_alleles];
		l.mutation_kicker[x]=new double[l.number_of_alleles];
	}

	stream.getline(buffer,256);
	stream.getline(buffer,256);

	for(int x=0;x<l.number_of_alleles;x++){
		for(int y=0;y<l.number_of_alleles-1;y++) {
			stream.getline(buffer,256, ',');
			l.mutation_rate[x][y]=atof(buffer);
		}
		stream.getline(buffer,256);
		l.mutation_rate[x][l.number_of_alleles-1]=atof(buffer);
	}

	stream.getline(buffer,256);
	stream.getline(buffer,256);

	for(int x=0;x<l.number_of_alleles;x++){
		for(int y=0;y<l.number_of_alleles-1;y++) {
			stream.getline(buffer,256, ',');
			l.selection[x][y]=atof(buffer);
		}
		stream.getline(buffer,256);
		l.selection[x][l.number_of_alleles-1]=atof(buffer);
	}

	stream.getline(buffer,256);
	stream.getline(buffer,256);

	for(int x=0;x<l.number_of_alleles;x++){
		for(int y=0;y<l.number_of_alleles-1;y++) {
			stream.getline(buffer,256, ',');
			l.mutation_kicker[x][y]=atof(buffer);
		}
		stream.getline(buffer,256);
		l.mutation_kicker[x][l.number_of_alleles-1]=atof(buffer);
	}*/

	return stream; // must return stream
};

void make_selection_tensor (tensor <double> *p, vector <locus*> loci, double (*epistatic_function) (double, double) ) {	
	unsigned int number_of_loci=loci.size();
	INDEX hap_index (number_of_loci), dip_index (number_of_loci*2);
	for(int x=0; x<number_of_loci;x++) {
		hap_index[x]=loci[x]->number_of_alleles;
		dip_index[x]=loci[x]->number_of_alleles;
		dip_index[x+number_of_loci]=loci[x]->number_of_alleles;
	}
	hap_index.make_step_array();
	
	tensor <double> t(&dip_index);

	unsigned int number_of_genotypes=hap_index.step[hap_index.D];

	double fitness;
	INDEX hap1(0), hap2(0);
	hap2=hap_index;
	hap1=hap_index; 

	for(int x=0; x<number_of_genotypes;x++) {

		hap1.setFromOffset(x);
		for(int y=0; y<number_of_genotypes;y++) {

			hap2.setFromOffset(y);

			fitness=loci[0]->selection[hap1[0]][hap2[0]];
			for(int z=1;z<number_of_loci;z++){
				fitness=epistatic_function(fitness, loci[z]->selection[hap1[z]][hap2[z]]);
			}
			t(x+y*number_of_genotypes)=fitness;
		}
	}
	*p=t;
};

void make_mutation_tensor (tensor <double> *p, vector <locus*> loci, double (*epistatic_function) (double, double) ) {
	unsigned int number_of_loci=loci.size();
	INDEX hap_index (number_of_loci), dip_index (number_of_loci*2);
	for(int x=0; x<number_of_loci;x++) {
		hap_index[x]=loci[x]->number_of_alleles;
		dip_index[x]=loci[x]->number_of_alleles;
		dip_index[x+number_of_loci]=loci[x]->number_of_alleles;
	}
	dip_index.make_step_array();
	hap_index.make_step_array();
	
	tensor <double> t(&dip_index);

	unsigned int number_of_genotypes=hap_index.step[hap_index.D];

	double kicker, temp;
	INDEX hap1(0), hap2(0);
	hap2=hap_index;
	hap1=hap_index; 

	for(int x=0; x<number_of_genotypes;x++) {

		hap1.setFromOffset(x);
		for(int y=0; y<number_of_genotypes;y++) {

			hap2.setFromOffset(y);

			kicker=loci[0]->mutation_kicker[hap1[0]][hap2[0]];
			for(int z=1;z<number_of_loci;z++){
				temp=loci[z]->mutation_kicker[hap1[z]][hap2[z]];
				kicker=epistatic_function(kicker, loci[z]->mutation_kicker[hap1[z]][hap2[z]]);
			}
			t(x+y*number_of_genotypes)=kicker;
		}
	}
	*p=t;
};

bool compare_loci (locus *first, locus *second)
{
	if (first->chromosome < second->chromosome) return true; 		
	else if (first->chromosome < second->chromosome) return false;
	else if (first->cM_from_p0 < second->cM_from_p0) return true;
	else if (first->cM_from_p0 < second->cM_from_p0) return false;
	else return false;
}

class PAIRS {
	public:
	double * number, weight;
	PAIRS (double *n, double w){
		number=n;
		weight=w;
	};
	~PAIRS(){};
};

tensor <double> & reduce (tensor <double> *p, vector <vector <PAIRS> > function){
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
		return *p;
};

class ENV_VAR{
	private:
	public:
	//used in calculation, there must be a seperate record for each thread!

	double *mean_W, *mean_U;
	double extinction_time;

	MTRand mtrand;

	list <DATA_RECORD> data;

	list <DATA_RECORD>::iterator it, begin, end;

	unsigned int records;

	//static durring thread execution.

	int number_of_mutator_loci;
	int number_of_nuetral_loci;

	int N, F, ID;
	double ms,s,max_ms;
	double W, U, U_base;

	unsigned int T;  
	unsigned int replicates;  

	ENV_VAR(){
		mean_W=NULL;
		mean_U=NULL;
		records=0;
		begin=data.begin();
		end=data.end();
	}

	ENV_VAR& operator= (const ENV_VAR& E){
		list <DATA_RECORD>::iterator temp;
		number_of_mutator_loci=E.number_of_mutator_loci;
		number_of_nuetral_loci=E.number_of_nuetral_loci;
		N=E.N;
		F=E.F;
		ID=E.ID;
		ms=E.ms;
		s=E.s;
		max_ms=E.max_ms;
		W=E.W;
		U=E.U;
		U_base=E.U_base;
		T=E.T;
		replicates=E.replicates;

		records=0;
		data.clear();
		it=data.begin();
		temp=E.begin;

		while(temp!=E.end){
			data.push_back(*temp);
			records++;
			temp++;
		}

		mean_W=new double[T];
		mean_U=new double[T];
	}

	~ENV_VAR(){
		if(mean_W!=NULL) {delete [] mean_W; mean_W=NULL;}
		if(mean_U!=NULL) {delete [] mean_U; mean_U=NULL;}
		if(!data.empty()) data.clear();
	}

	list <DATA_RECORD> ::iterator set_name (string name){
		it=data.begin();
		list <DATA_RECORD>::iterator end=data.end();
		while(it!=end){
			if(it->name.compare(name)==0) return it;
			it++;
		}
		return it;
	}

	void push_back(double _double, string name){
		set_name(name);
		if(it!=data.end()){
			it->push_back(_double);
		}
		else{
			DATA_RECORD temp(name, _double);
			data.push_back(temp);
			records++;
		}
		//end=data.end();
	}

	void confuse(double *_double, unsigned int size, string name){
		set_name(name);
		if(it!=data.end()){
			it->confuse(_double, size);
		}
		else{
			DATA_RECORD temp(name, _double, size);
			data.push_back(temp);
			records++;
		}
		//end=data.end();

	}

	void push_back(vector <double> _double, string name){
		set_name(name);
		if(it!=data.end()){
			it->push_back(_double);
		}
		else{
			DATA_RECORD temp(name, _double);
			data.push_back(temp);
			records++;
		}
		//end=data.end();
	}

	void confuse(vector <double> _double, string name){
		set_name(name);
		if(it!=data.end()){
			it->confuse(_double);
		}
		else{
			DATA_RECORD temp(name, _double);
			data.push_back(temp);
			records++;
		}
		//end=data.end();
	}

	void confuse(list <DATA_RECORD> DATA){
		list <DATA_RECORD>::iterator d_it, d_end;
		d_end=DATA.end();
		d_it=DATA.begin();
		while(d_it!=d_end){
			confuse(*d_it);
			d_it++;
		}
		//end=data.end();
	}

	void confuse(DATA_RECORD DATA){
		set_name(DATA.name);
		if(it!=data.end()){
			it->confuse(DATA);
		}
		else{
			data.push_back(DATA);
			records++;
		}		
	}

	void pop_back(){
		data.pop_back();
		records--;
		//end=data.end();
	}

	vector <double> get_double (unsigned int x){
		it=data.begin();
		for(int y=0;y<x;y++) it++;
		return (it->datum);		
	}

	DATA_RECORD get_data_record (unsigned int x){
		it=data.begin();
		for(int y=0;y<x;y++) it++;
		return *it;		
	}

	string get_name (unsigned int x){
		it=data.begin();
		for(int y=0;y<x;y++) it++;
		return (it->name);			
	}

	void set_begin (void){
		it=data.begin();
	}

	void clear (void){
		if(mean_W!=NULL) {delete [] mean_W; mean_W=NULL;}
		if(mean_U!=NULL) {delete [] mean_U; mean_U=NULL;}
		if(!data.empty()) data.clear();
	}
	void print (void){
		cout << "number_of_mutator_loci " << number_of_mutator_loci << endl;
		cout << "number_of_nuetral_loci " << number_of_nuetral_loci << endl;
		cout << "N " << N << endl;
		cout << "F " << F << endl;
		cout << "ID " << ID << endl;
		cout << "ms=" << ms << endl;
		cout << "s=" << s << endl;
		cout << "max_ms=" << max_ms << endl;
		cout << "W=" << W << endl << "U=" << U << endl << "U_base=" << U_base << endl << "T " << T << endl << "replicates=" << replicates << endl;		
		cout << "records=" << records << endl;
	}
};
#endif

