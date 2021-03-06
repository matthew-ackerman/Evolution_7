#ifndef _TENSOR_HPP
#define _TENSOR_HPP

#include <iostream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
using std::cout;

class INDEX
{
	public:
	int D;

	unsigned int *cell;
	unsigned int *step;

	int size;

	INDEX(){
		cell=NULL;
		step=NULL;
	}

	INDEX(int d){
		cell=NULL;
		step=NULL;
		D=d;
		cell=new unsigned int [d];
		step = new unsigned int [D+1];
	}

	INDEX(int a, int b){
		cell=NULL;
		step=NULL;

		cell=new unsigned int [2];
		cell[0]=a; cell[1]=b;

		D=2;

		step = new unsigned int [D+1];
		step[0]=1;
		for(int x=1;x<D+1;x++){
			step[x]=step[x-1]*cell[x-1];
		}
	}

	INDEX(int a, int b, int c){
		cell=NULL;
		step=NULL;

		cell=new unsigned int [3];
		cell[0]=a; cell[1]=b; cell[2]=c;

		D=3;

		step = new unsigned int [D+1];
		step[0]=1;
		for(int x=1;x<D+1;x++){
			step[x]=step[x-1]*cell[x-1];
		}
	}

	INDEX(int a, int b, int c, int d){
		cell=NULL;
		step=NULL;

		cell=new unsigned int [4];
		cell[0]=a; cell[1]=b; cell[2]=c; cell[3]=d;

		D=4;

		step = new unsigned int [D+1];
		step[0]=1;
		for(int x=1;x<D+1;x++){
			step[x]=step[x-1]*cell[x-1];
		}
	}

	~INDEX(){
		if(cell!=NULL) {delete [] cell; cell=NULL;}
		if(step!=NULL) {delete [] step; step=NULL;}
	}	

	void kill(void){
		if(cell!=NULL) {delete [] cell; cell=NULL;}
		if(step!=NULL) {delete [] step; step=NULL;}
	}	

	bool operator==(INDEX &a) const {
		if(D==a.D){
			for(int x=0;x<D;x++){
				if(cell[x]!=a[x]){
					return false;
				}
			}
			return true;
		}
		else{
			 return false;
		}
	}

	unsigned int operator[] (unsigned int _cell) const{
  		return cell[_cell];
	}
	unsigned int & operator[] (unsigned int _cell){
  		return cell[_cell];
	}
	
	unsigned int  * make_step_array(void){
		if(step!=NULL) delete [] step;
		step = new unsigned int [D+1];
		step[0]=1;
		for(int x=1;x<D+1;x++){
			step[x]=step[x-1]*cell[x-1];
		}
		return step;
	}

	void setFromOffset(unsigned int off){
		int total=1;
		int R=0;
		for (int x=D-1;x>=0;x--){
			if(step[x]!=0){
				cell[x]=off/step[x];
				R=off%(step[x]);
				off=R;
			}
		}
	}

	std::string print (void) const{
		std::stringstream out(" ");
		if(cell!=NULL){
			out << "(";
			for(int x=0;x<D;x++){
				out << cell[x] << ",";
			}
			out << '\b' <<")";
		}
		else out << "(cell is NULL)";
		
		if(step!=NULL){
			out << ":(";
			for(int x=0;x<D;x++){
				out << step[x] << ",";
			}
			out << '\b' <<")";
		}
		else out << ":(step is NULL)";
		
		return out.str();
	}

	INDEX & operator = (const INDEX & ind){
		if (this!=&ind) {
			if (step!=NULL) delete [] step;
			if (cell!=NULL) delete [] cell;

			D=ind.D;
			size=ind.size;

			step = new unsigned int [D+1];
			cell=new unsigned int [D];
	
			for(int x=0;x<D;x++) {
				if(ind.cell!=NULL) cell[x]=ind.cell[x];
				if(ind.step!=NULL) step[x]=ind.step[x];
			}
			if(ind.step!=NULL) step[D]=ind.step[D];
		}
		return *this;
	}
	void clear (void){
		for(int x=0;x<D;x++) cell[x]=0;
	}
};

template  <typename type>
class tensor 
{
	private:

	type *cell; 	//Data member.
	INDEX *width;	//number of data items in a row. This should be private.

	unsigned int rank;	//numbr of dimensions.
	unsigned int size;	//total number of data items.

	public: 

	unsigned int *step;	//step cell.

	type *it;		//iterator.
	type *begin, *end;	//start and stop

	tensor (){
		width=NULL;
		cell=NULL;
		step=NULL;
	};

	/* constructors */
	tensor (unsigned int X){

		width=NULL;
		cell=NULL;
		step=NULL;

		INDEX *w;

		w=new INDEX(X);
		for (unsigned int x=0;x<X;x++){
			(*w)[x]=1;
		}

		rank=(*w).D;

		width=new INDEX (rank);

		for(unsigned int x=0;x<rank;x++) (*width).cell[x]=(*w).cell[x];

		step = (unsigned int *) calloc ((rank+1), sizeof(unsigned int));
		step[0]=1;
		for(unsigned int x=1;x<rank+1;x++){
			step[x]=step[x-1]*(*width).cell[x-1];
		}

		cell=new type [step[rank]];

		begin=cell;
		end=&cell[step[rank]];
		size=step[rank];

		clear();

		it=begin;
		width->make_step_array();
	}

	tensor (unsigned int X, unsigned int * int_cell, type * _cell){
		width=NULL;
		cell=NULL;
		step=NULL;

		INDEX *w;

		w=new INDEX(X);
		for (unsigned int x=0;x<X;x++){
			(*w)[x]=*(int_cell+x);
		}

		rank=(*w).D;

		width=new INDEX (rank);

		for(unsigned int x=0;x<rank;x++) (*width).cell[x]=(*w).cell[x];

		step = (unsigned int *) calloc ((rank+1), sizeof(unsigned int));
		step[0]=1;
		for(unsigned int x=1;x<rank+1;x++){
			step[x]=step[x-1]*(*width).cell[x-1];
		}

		cell=new type [step[rank]];

		begin=cell;
		end=&cell[step[rank]];
		size=step[rank];

		clear();

		for(unsigned int x=0;x<size;x++){
			cell[x]=_cell[x];
		}

		it=begin;
		width->make_step_array();
	}

	tensor (unsigned int X, unsigned int * int_cell){
		width=NULL;
		cell=NULL;
		step=NULL;

		INDEX *w;

		w=new INDEX(X);
		for (unsigned int x=0;x<X;x++){
			(*w)[x]=*(int_cell+x);
		}

		rank=(*w).D;

		width=new INDEX (rank);

		for(unsigned int x=0;x<rank;x++) (*width).cell[x]=(*w).cell[x];

		step = (unsigned int *) calloc ((rank+1), sizeof(unsigned int));
		step[0]=1;
		for(unsigned int x=1;x<rank+1;x++){
			step[x]=step[x-1]*(*width).cell[x-1];
		}

		cell=new type [step[rank]];
		begin=cell;
		end=&cell[step[rank]];
		size=step[rank];

		clear();
		it=begin;
		width->make_step_array();
	}

	tensor (INDEX *w)
	{  
		width=NULL;
		cell=NULL;
		step=NULL;

		rank=(*w).D;

		width=new INDEX (rank);

		for(unsigned int x=0;x<rank;x++) (*width).cell[x]=(*w).cell[x];		

		step = (unsigned int *) calloc ((rank+1), sizeof(unsigned int));
		step[0]=1;

		for(unsigned int x=1;x<rank+1;x++){
			step[x]=step[x-1]*(*width).cell[x-1];
		}

		cell=new type [step[rank]];
		begin=cell;
		end=&cell[step[rank]];
		size=step[rank];

		clear();
		it=begin;
		width->make_step_array();
	}

	tensor (char *str)
	{  
		unsigned int a;
		width=NULL;
		cell=NULL;
		step=NULL;

		FILE* TensorFile;
		TensorFile=fopen(str, "rb");
		if(TensorFile==NULL){
			fprintf(stderr, "cannot open file %s\n", str);
		}
		fread (&rank, sizeof(int), 1, TensorFile);
		width=new INDEX(rank);

		for(int x=0;x<rank;x++){
			fread (a, sizeof(unsigned int), 1, TensorFile);
			(*width)[x]=a;
		}

		step = (unsigned int *) calloc ((rank+1), sizeof(unsigned int));
		step[0]=1;
		for(int x=1;x<rank+1;x++){
			step[x]=step[x-1]*(*width).cell[x-1];
		}

		cell=(type *) calloc(step[rank], sizeof(type));

		begin=&cell[0];
		end=&cell[step[rank]];
		size=step[rank];
		it=begin;
		while(it<end){
			fread (it, sizeof(type),1,TensorFile);
			it++;
		}
		fclose(TensorFile);
	}

	void resize (INDEX *w){
		
		kill();

		width=NULL;
		cell=NULL;
		step=NULL;

		rank=(*w).D;

		width=new INDEX (rank);

		for(unsigned int x=0;x<rank;x++) (*width).cell[x]=(*w).cell[x];		

		step = (unsigned int *) calloc ((rank+1), sizeof(unsigned int));
		step[0]=1;

		for(unsigned int x=1;x<rank+1;x++){
			step[x]=step[x-1]*(*width).cell[x-1];
		}

		cell=new type [step[rank]];
		begin=cell;
		end=&cell[step[rank]];
		size=step[rank];

		clear();
		it=begin;
		width->make_step_array();		
	}
	
	void clear (void){
		it=begin;
		while ((it < end)){(*it)=0; it++;}
	}

	~tensor()
	{
		if (width != NULL){
			delete width;
			width=NULL;
		}
		if(cell!=NULL) {
			cell=NULL;
			delete [] cell;
		}
		if(step!=NULL) {
			free(step);
			step=NULL;
		}
	}

        void kill(void)
	{
		if (width != NULL){
			delete width;
			width=NULL;
		}
		if(cell!=NULL) {
			cell=NULL;
			delete [] cell;
		}
		if(step!=NULL) {
			free(step);
			step=NULL;
		}
	}

	std::string print (type N) const {
		std::stringstream out(" ");

		INDEX temp(0);
		temp=*width;

		int x=0;
		while (x<size){
			temp.setFromOffset(x);
			out << temp.print() << "x" << x << " " << (begin+x) << " " << (*(begin+x))*N << endl;
			x++;
		}
		return out.str();
	}

	std::string print_quite (void) const{
		std::stringstream out(" ");
		int x=0;
		while (x<size){
			out << ", " << (*(begin+x));
			x++;
		}
		return out.str();
	}

	std::string print (void) const{
		std::stringstream out(" ");
		INDEX index(0);
		index=*width;
		int x=0;
		while (x<size){
			index.setFromOffset(x);
			out << index.print() << " " << (begin+x) << " " << (*(begin+x)) << endl;
				x++;
		}
		return out.str();
	}

	void operator= (type z)
	{
		(*it)=z;
	}

	tensor & operator= (const tensor &z)
	{
		//cout << z.print() << endl;

		if(width!=NULL && *(width)==*(z.width)){
			//cout << this << " fast eq\n";

			it=begin;
			type *z_it=z.begin;

			while (it!=end){
			//	cout << *it << "/" << *z_it << endl;
				*it=*z_it;
			//	cout << *it << "\\" << *z_it << endl;
				it++;
				z_it++;
			}
		}
		else{
			//cout << this << " " <<" slow eq\n";
			kill();
			//cout << " set INDEX ";
			width=new INDEX(z.width->D);
			*(width)=*(z.width);
			//cout << "done." << endl;
			rank=z.rank;

			step = (unsigned int *) calloc ((rank+1), sizeof(unsigned int));
			step[0]=1;
			for(int x=1;x<rank+1;x++){
				step[x]=step[x-1]*(*width).cell[x-1];
			}
	
			cell=new type [step[rank]];
			begin=cell;
	
			end=&cell[step[rank]];
			size=step[rank];
	
			clear();
			it=begin;

			type *z_it=z.begin;

			while (it!=end){
				*it=*z_it;
				it++;
				z_it++;
			}
		}
		return *this;
	}

	void operator+= (type z) const
	{
		(*it)+=z;
	}

	/*absolute element access operators*/

	type & operator() (INDEX *index) const
	{
		double * t_it=begin;
		for(int x=0;x<rank;x++){
			t_it+=(step[x])*(*index)[x];
		}

		return *t_it;
	}

	type & operator() (unsigned int x) 
	{
		return *(begin+x);
	}

	type & operator() (unsigned int x, unsigned int y)
	{
		return *(begin+x+(step[1])*y);
	}


	type & operator() (unsigned int w, unsigned int x, unsigned int y)
	{
		return *(begin+w+step[1]*x+step[2]*y);
	}

	type & operator() (unsigned int w, unsigned int x, unsigned int y, unsigned int z)
	{
		return *(begin+w+step[1]*x+step[2]*y+step[3]*z);
	}

	/*relative element access operators*/

	type & relative (INDEX *index)
	{
		for(int x=0;x<rank;x++){
			it+=(step[x])*(*index)[x];
		}
		return *(it);
	}

	type & relative (unsigned int x) 
	{
		return *(it+x);
	}

	type & relative (unsigned int x, unsigned int y)
	{
		return *(it+x+(step[1])*y);
	}

	type & relative (unsigned int w, unsigned int x, unsigned int y)
	{
		return *(it+w+step[1]*x+step[2]*y);
	}

	type & relative (unsigned int w, unsigned int x, unsigned int y, unsigned int z)
	{
		return *(it+w+step[1]*x+step[2]*y+step[3]*z);
	}

	/*depricated functions*/

	tensor & go (INDEX *index)
	{
		it=begin;
		for(int x=0;x<rank;x++){
			it+=(step[x])*(*index)[x];
		}

		return *this;
	}

	type read(void) 
	{ 
		return (*it);
	}

	type readpp()
	{
		it++;
#if (SANE)
		if (it > end) it=begin;
#endif
		return *(it-1);
	}

	void setpp(type z)
	{
		(*it)=z;
		it++;
#if (SANE)
		if (it > end) it=begin;
#endif
	}

	/*private data access*/

	bool start(void){
		it=begin;
		return true;
	}

	unsigned int size_of(void) const {
		return step[rank];
	}

	unsigned int get_order(void) const {
		return rank;
	}

	unsigned int get_rank (void) const {
		return rank;
	}

	unsigned int get_size (void) const {
		return size;
	}

	type * get_cell(void) const {
		return cell;
	}

	unsigned int get_width(int x) const {
		return (*width)[x];
	}
	
	INDEX get_width() const{
		return *width;
	}

	unsigned int *get_step(void) const {
		return step;
	}
	unsigned int get_step(int x) const {
		return step[x];
	}

	/* mathmatical operators*/

	/* A[i,j]=Sum[k] Sum[l] Sum [M] B[k,l,m]*C[k,l,m,i,j];*/

	tensor & transition (tensor * t1, tensor * t2){
		type sum;
		type *t1_it, *t2_it, *t1_start, *t1_end;
	
		t1_it=t1->begin;
		t2_it=t2->begin;
		t1_start=t1_it;
		t1_end=t1->end;
		it=begin;

		while(it!=end){
			t1_it=t1_start;
			sum=0;
			while (t1_it!=t1_end){
				sum+=(*t1_it)*(*t2_it);
				t1_it++;
				t2_it++;
			}
			(*it)=sum;
			it++;
		}
		return *this;
	}

	/* A[i,j,k,l]=B[i,j]*C[k,l];*/

	tensor & product (tensor * t1, tensor * t2){
		type *t1_it, *t2_it, *t2_start, *t2_end;
	
		t1_it=t1->begin;
		t2_it=t2->begin;
		t2_start=t1_it;
		t2_end=t2->end;
		it=begin;

		while(it!=end){
			t2_it=t2_start;
			while (t2_it!=t2_end){
				(*it)=(*t1_it)*(*t2_it);
				t2_it++;
				it++;
			}
			t1_it++;
		}
		return *this;
	}

	/* A[i,j,k,l]=A[i,j,k,l]/density(A);*/

	void normalize(void){
		type SUM=0;
		it=begin;
		while(it!=end){
			SUM+=(*it);
			it++;
		}
		it=begin;
		while(it!=end){
			(*it)=(*it)/SUM;
			it++;
		}
	}

	/*density=sum [i] sum [j] sum [k] sum [l] A[i,j,k,l];*/ 
	
	type density(void){
		type SUM=0;
		it=begin;
		while(it<end){
			SUM+=(*it);
			it++;
		}
		return SUM;
	}

	/*edge_sum=sum [j] sum [k] sum [l] A[0,j,k,l];*/

	type edge_sum(void){
		type SUM=0;

		it=begin;
		it+=step[1];
		it-=step[0];
		while(it<end){
			SUM+=(*it);
			it+=step[1];
		}
		it=begin;
		it+=step[rank];
		it-=step[1];
		while(it<end){
			SUM+=(*it);
			it++;
		}		
		return SUM;
	}

	/* i/o */

	int saveToFile(char *str){
		FILE* TensorFile;
		TensorFile=fopen(str, "wb");
		if(TensorFile==NULL){
			fprintf(stderr, "cannot open file %s\n", str);
			return 1;
		}
		fwrite (&rank, sizeof(int), 1, TensorFile);
		for(int x=0;x<rank;x++){
			fwrite (&(*width)[x], sizeof(unsigned int), 1, TensorFile);
		}
		it=begin;
		while(it<end){
			fwrite (it, sizeof(type),1,TensorFile);
			it++;
		}
		fclose(TensorFile);
		return 0;
	}

	std::string printXBarGivenY (void){
		std::stringstream out(" ");
		for(unsigned y=0;y<(*width)[1];y++){
			it=(begin+(step[1])*y);
			type XBar=0;
			type Y=0;
			for(type x=0;x<(*width)[0];x++){
				XBar+=(*it)*x;
				Y+=(*it);
				it++;
			}
			out << XBar/Y << ","; //<< (*it) << ",";
		}
		
		return out.str();
	};

	std::string printYBarGivenX (void){
		std::stringstream out(" ");
		if (width!=NULL){
			for(unsigned x=0;x<(*width)[0];x++){
				it=(begin+x);
				type YBar=0;
				type X=0;
				for(type y=0;y<(*width)[1];y++){
					YBar+=(*it)*y;
					X+=(*it);
					it+=step[1];
				}
				out << (YBar/X) << "," ; //<< (*it) << ",";
			}
		}
		return out.str();
	};

	/* boolain operators */

	bool operator< (type a) const
	{
		if(*it<a) return true;
		return false;
	}

	bool operator<= (type a) const
	{
		if(*it<=a) return true;
		return false;
	}

	bool operator> (type a) const
	{
		if(*it>a) return true;
		return false;
	}

	bool operator>= (type a) const
	{
		if(*it>=a) return true;
		return false;
	}

	bool operator== (type a) const
	{
		if(*it==a) return true;
		return false;
	}

	bool operator!= (type a) const
	{
		if(*it!=a) return true;
		return false;
	}

	/* casts */ 
	/*type operator (tensor){
		return tensor.it;
	}*/
};

struct PAIR;

class ELEMENT{
	public:
	double x,y,z,R,B,G,A;
	ELEMENT *next;
	PAIR *it, *start, *end;
	bool wire, quad;
	ELEMENT(double X, double Y, double Z){
		x=X;
		y=Y;
		z=Z;
		start=NULL;
		it=NULL;
		end=NULL;
		wire=false;
		quad=false;
		R=0.0;
		B=0.0;
		G=0.0;
		A=1.0;
	}
};

struct PAIR{
	ELEMENT *element;
	PAIR *next;
	double R,B,G,A, w;
	bool wire, quad;
};

class LIST
{
	public:
	ELEMENT *it, *end, *start;
	int size;
	double max_x, max_y, max_z, min_x, min_y, min_z;
	double rel_x, rel_y, rel_z, x_unit, y_unit, z_unit;
	LIST(double X, double Y, double Z){
		start=new ELEMENT(X,Y,Z);
		end=start;
		size=1;
		max_x=X;
		max_y=Y;
		max_z=Z;
		min_x=X;
		min_y=Y;
		min_z=Z;
	}
	ELEMENT * add_element(double X, double Y, double Z){
		it=end;
		it->next=new ELEMENT(X,Y,Z);
		size++;
		end=it->next;

		if(X<min_x) min_x=X;
		else if(X>max_x) max_x=X;

		if(Y<min_y) min_y=Y;
		else if(Y>max_y) max_y=Y;

		if(Z<min_z) min_z=Z;
		else if(Z>max_z) max_z=Z;

		return end;
	}
	void join_element(ELEMENT *from, ELEMENT *to, bool wire, bool quad, double red, double green, double blue, double width){
		if (from->start==NULL){
			from->start=new PAIR;
			from->start->element=to;
			from->start->next=NULL;
			from->it=from->start;
			from->end=from->start;
			from->end->wire=wire;
			from->end->quad=quad;
			from->start->R=red;
			from->start->G=green;
			from->start->B=blue;
			from->start->w=width;
		}
		else{
			from->it=from->end;
			from->it->next=new PAIR;
			from->end=from->it->next;
			from->it->next->next=NULL;
			from->it->next->element=to;
			from->end->wire=wire;
			from->end->quad=quad;
			from->it->R=red;
			from->it->G=green;
			from->it->B=blue;
			from->it->w=width;
		}
	}
	void join_element(ELEMENT *from, ELEMENT *to, bool wire, bool quad, double red, double green, double blue){
		if (from->start==NULL){
			from->start=new PAIR;
			from->start->element=to;
			from->start->next=NULL;
			from->it=from->start;
			from->end=from->start;
			from->end->wire=wire;
			from->end->quad=quad;
			from->start->R=red;
			from->start->G=green;
			from->start->B=blue;
			from->start->w=1.0;
		}
		else{
			from->it=from->end;
			from->it->next=new PAIR;
			from->end=from->it->next;
			from->it->next->next=NULL;
			from->it->next->element=to;
			from->end->wire=wire;
			from->end->quad=quad;
			from->it->R=red;
			from->it->G=green;
			from->it->B=blue;
			from->start->w=1.0;
		}
	}
	void scale(double FOV_size){
		it=start;
		x_unit=(1.0/(max_x-min_x))*FOV_size;
		y_unit=(1.0/(max_y-min_y))*FOV_size;
		z_unit=(1.0/(max_z-min_z))*FOV_size;
		if(max_x==min_x) x_unit=0;
		if(max_y==min_y) y_unit=0;
		if(max_z==min_z) z_unit=0;

		while (it!=NULL){
			it->x=it->x*x_unit-(FOV_size/2.0);
			it->y=it->y*y_unit-(FOV_size/2.0);
			it->z=it->z*z_unit-(FOV_size/2.0);
			it=it->next;
		}
	}
	void scale_keepRatio(double FOV_size){
		it=start;
		x_unit=(1.0/(max_x-min_x))*FOV_size;
		y_unit=(1.0/(max_y-min_y))*FOV_size;
		z_unit=(1.0/(max_z-min_z))*FOV_size;
		if(max_x==min_x) x_unit=0;
		if(max_y==min_y) y_unit=0;
		if(max_z==min_z) z_unit=0;

		double max_unit=x_unit, min_unit=x_unit;
		if (y_unit>max_unit) max_unit=y_unit;
		else if (y_unit<min_unit) min_unit=y_unit;

		if (z_unit>max_unit) max_unit=z_unit;
		else if (z_unit<min_unit) min_unit=z_unit;

		rel_x=x_unit/max_unit;
		rel_y=y_unit/max_unit;
		rel_z=z_unit/max_unit;

		while (it!=NULL){
			it->x=it->x*min_unit-(FOV_size/2.0);
			it->y=it->y*min_unit-(FOV_size/2.0);
			it->z=it->z*min_unit-(FOV_size/2.0);
			it=it->next;
		}
	}

	void reticulate_from_vector (ELEMENT *** vector_list, int xsize, int ysize){
		double A;
		double B;
		double C;
		A=(max_z-min_z)/15.0;
		B=A/10.0;
		C=B/10.0;
		for(int x=0;x<xsize;x++){
			for (int y=0;y<ysize;y++){
				if(x<xsize-1) join_element(vector_list[x][y],vector_list[x+1][y],true,true,0,0,0);
				if(x<xsize-1&&y<ysize-1){
					double meanValue=(vector_list[x][y]->y+vector_list[x+1][y]->y+vector_list[x][y+1]->y+vector_list[x+1][y+1]->y)/4.0;
					if(meanValue>0.0){
						if ((meanValue*C) < 1.0) {
							vector_list[x][y]->R=1.0-meanValue*C;
							vector_list[x][y]->B=1.0-meanValue*B;
							vector_list[x][y]->G=1.0-meanValue*A;
							vector_list[x][y]->A=0.9;
						}
						else if ((meanValue*B)<1.0) {
							vector_list[x][y]->R=0.0;
							vector_list[x][y]->G=1.0-meanValue*B;
							vector_list[x][y]->G=1.0-meanValue*A;
							vector_list[x][y]->A=0.9;
						}
						else if ((meanValue*A)<1.0) {
							vector_list[x][y]->R=0.0;
							vector_list[x][y]->B=0.0;
							vector_list[x][y]->G=1.0-meanValue*A;
							vector_list[x][y]->A=0.9;
						}
						else {
							vector_list[x][y]->R=0.0;
							vector_list[x][y]->B=0.0;
							vector_list[x][y]->G=0.0;
							vector_list[x][y]->A=0.9;
						}
					}
					else {
						if ((meanValue*A)>-1.0) {vector_list[x][y]->R=1.0+meanValue*A; vector_list[x][y]->B=1.0+meanValue*B; vector_list[x][y]->G=1.0+meanValue*C; vector_list[x][y]->A=0.7;}
						else if ((meanValue*B)>-1.0) {vector_list[x][y]->R=0.0; vector_list[x][y]->B=1.0+meanValue*B; vector_list[x][y]->G=1.0+meanValue*C; vector_list[x][y]->A=0.7;}
						else if ((meanValue*C)>-1.0) {vector_list[x][y]->R=0.0; vector_list[x][y]->G=0.0; vector_list[x][y]->G=1.0+meanValue*C; vector_list[x][y]->A=0.7;}
						else {vector_list[x][y]->R=0.0; vector_list[x][y]->B=0.0;vector_list[x][y]->G=0.0; vector_list[x][y]->A=0.7;}
					}
					join_element(vector_list[x][y],vector_list[x+1][y+1],false,true,0,0,0);
					vector_list[x][y]->quad=true;
				}
				if(y<ysize-1) join_element(vector_list[x][y],vector_list[x][y+1],true,true,0,0,0);
			}
		}	
	}
	void reticulate_without_vector (){
		double A;
		double B;
		double C;
		A=(max_z-min_z)/15.0;
		B=A/10.0;
		C=B/10.0;
	}
};
#endif


///This does not calculate a general tensor product, but instead calculates ....

#if(0==1)
void fast_product(tensor *v1, tensor *v2, tensor *p)
{
	PRECISION sum;

	type *v1_it, *v2_it, *v1_start, *v1_end, *p_it, *p_end;

	v1_it=v1->begin;
	v2_it=v2->begin;
	v1_start=v1_it;
	v1_end=v1->end;
	p_it=p->begin;
	p_end=p->end;

	while(p_it!=p_end){
		v1_it=v1_start;
		sum=0;
		while (v1_it!=v1_end){
			sum+=(*v1_it)*(*v2_it);
			v1_it++;
			v2_it++;
		}
		(*p_it)=sum;
		p_it++;
	}
}

tensor & fast_product(tensor *v1, tensor *v2)
{
	type sum;
	type *v1_it, *v2_it, *v1_start, *v1_end, *p_it, *p_end;

	v1_it=v1->begin;
	v2_it=v2->begin;
	v1_start=v1_it;
	v1_end=v1->end;

//to do: finish
	INDEX p_index(2,2);

	tensor P(&p_index), *p;

	p=&P;

	p_it=p->begin;
	p_end=p->end;

	while(p_it!=p_end){
		v1_it=v1_start;
		sum=0.0;
		while (v1_it!=v1_end){
			sum+=(*v1_it)*(*v2_it);
			v1_it++;
			v2_it++;
		}
		(*p_it)=sum;
		p_it++;
	}
	return P;
}
#endif

/*
void partial_product(tensor *v1, tensor *v2, tensor *p, )	
{

	int size=v1->size_of();
	int fork_return, fork_count=0, size_over_forks=size/forks;

	type *v1_begin=v1->begin, *v1_end=v1->end, *p_begin=p->begin, *v2_begin=v2->begin;

	unsigned int rank1=v1->get_order();
	unsigned int rank2=v2->get_order();

	int pid=-1;
	int status=1;
	for(int x=0;x<(forks-1);x++){
		if(pid!=0){
			pid=fork();
			if(pid!=0)fork_count++;
		}
		else {
			x=forks;
		}
	}

	type sum;


	type * columnEnd, *v1_it, *v2_it, *p_it, *p_end;

	p_it=p_begin+fork_count*size_over_forks;
	v2_it=v2_begin+(fork_count*size_over_forks)*size;
	p_end=p_it+size_over_forks;

	v1_it=v1_begin;

	while(p_it!=p_end){
		v1_it=v1_begin;
		sum=0;
		while (v1_it!=v1_end){
			sum+=(*v1_it)*(*v2_it);
			v1_it++;
			v2_it++;
		}
		(*p_it)=sum;
		p_it++;
	}

	_exit(0);
}*/

#if (OLD)
/*void product2(newtensor *v1, tensor *v2, newtensor *p, double error)
{
	double sum;
	unsigned int pos;
	unsigned int width=(*v1).get_width();
	unsigned int size=(*v1).size_of();

	unsigned int rank1=(*v1).get_rank();
	unsigned int rank2=(*v2).get_rank();

	INDEX index ((*v2).get_rank());

	for(int x=0; x<index.D;x++) index[x]=0;
	
	int indexFlip; 
			
	(*v1).start();
	(*v2).start();
	(*p).start();
	(*p).purgeSum();

	sum=0;
	
	for (int z=0; z<size;z++){	
		if((*v1).readSum()<=error){
			(*v1).it+=width;
			(*v2).it+=width;
			z+=(width-1);
		}
		else{
			 for(int zz=0; zz<width;zz++){
				sum+=(*v1).read()*(*v2).read();
				(*v1).it++;
				(*v2).it++;
			}
			z+=(width-1);
		}
	}

	(*p)=sum;
	(*p).it++;

	index[rank2-rank1]++;
	indexFlip=1;

	while(index[rank2-rank1]>=width){
		if(index[rank2-rank1+indexFlip]>width) indexFlip++;
		else {
			index[rank2-rank1+indexFlip]++;
			indexFlip--;
			while (indexFlip>=0){
				index[rank2-rank1+indexFlip]=0;
				indexFlip--;
			}
		}
	}	


	while ((*p).it<(*p).end){

		(*v1).start();
		(*v2).go(&index);

		sum=0;
		for (int z=0; z<size;z++){	
			if((*v1).readSum()<=error){
				(*v1).it+=width;
				(*v2).it+=width;
				z+=(width-1);
			}
			else{
				 for(int zz=0; zz<width;zz++){
					sum+=(*v1).read()*(*v2).read();
					(*v1).it++;
					(*v2).it++;
				}
				z+=(width-1);
			}
		}
		(*p)=sum;
		(*p).it++;
		 
		index[rank2-rank1]++;
		indexFlip=1;
		while(index[rank2-rank1]>=width){
			if(index[rank2-rank1+indexFlip]>=width) indexFlip++;
			else {
				index[rank2-rank1+indexFlip]++;
				indexFlip--;
				while (indexFlip>=0){
					index[rank2-rank1+indexFlip]=0;
					indexFlip--;
				}
			}
		}
	}
}*/

/*void product (tensor * tg, sparce_tensor * v2, tensor *gp)
{

	tensor_wrapper tp(tg), pp(gp), *p, *v1;

	v1=&tp;
	p=&pp;

	double sum;
	unsigned int pos, skip, end, size=(*v1).size_of(), width=v1->get_width();
	
	unsigned int * index2;
	unsigned int rank=(*v1).get_rank();
	index2=(unsigned int *) calloc(rank, sizeof(unsigned int));

	INDEX index ((*v2).get_rank());
	bool first=false, empty=true,last=false;

	for(int x=0; x<index.D;x++) index[x]=0;

	(*v1).start();
	(*v2).start();
	(*p).start();

	for(;index2[rank-1]<width;index2[0]++){
		for(int x=0; x<rank-1; x++){
			 if(index2[x]==width){
				index2[x]=0;
				index2[x+1]++;
			}
		}
		if (index2[rank-1]<width){
		
			for(int y=0; y<rank; y++) index[(index.D-rank)+y]=index2[y];	

			(*v1).start();
			(*v2).go(&index);

			sum=0;
			empty=true;
		
			for (int z=0; z<=size;z++){	
				sum+=(*v1).read()*(*v2).read();
				skip=(*v2).get_next();
				(*v1).it+=skip;
				(*v2).vinc(0,skip);
				z+=(skip);
				(*v1).it++;
				(*v2).inc(0);
			}

			if(sum!=0) empty=false;
			if(!first&&!empty)first=true;
			else if(first&&empty&&!last) index2[rank-1]=width;
			(*p)=sum;
			(*p).it++;
		}		
	}
	free(index2);
}*/
/*
class buffer {
    int head;
    int tail;
    int store[10];
    _Mutex_simple mutex;
    _Condition not_full;
    _Condition not_empty;
public:
    buffer() : head( 0 ) , tail( 0 ) { }
    void insert( int arg ) {
        _Wait( not_full; (head+1)%10 != tail ) {
            store[head] = arg;
            head = (head+1)%10;
            _Notify( not_empty );
        }
    }
    int remove() {
        int result;
        _Wait( not_empty; head != tail ) {
            result = store[tail];
            tail = (tail+1)%10;
            _Notify( not_full );
        }
    }
} */

#endif
