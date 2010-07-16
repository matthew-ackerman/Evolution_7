#ifndef _TENSOR_HPP
#define _TENSOR_HPP

#include <iostream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef PRECISION
#define PRECISION double
#endif
#endif

#define OLD	0

using namespace std;
using std::cout;

class INDEX
{
	public:
	int D;
	unsigned int *cell;
	unsigned int *step;

	int size;
	int off_set;

	INDEX(int d){
		cell=NULL;
		step=NULL;

		cell=new unsigned int [d];
		D=d;
	}
	INDEX(int a, int b, int c, int d){
		cell=NULL;
		step=NULL;

		cell=new unsigned int [4];
		cell[0]=a; cell[1]=b; cell[2]=c; cell[3]=d;
		D=4;

		step = (unsigned int *) calloc ((D+1), sizeof(unsigned int));
		step[0]=1;
		for(int x=1;x<D+1;x++){
			step[x]=step[x-1]*cell[x-1];
		}
	}
	INDEX(int a, int b){
		cell=NULL;
		step=NULL;

		cell=new unsigned int [2];
		cell[0]=a; cell[1]=b;
		D=2;

		step = (unsigned int *) calloc ((D+1), sizeof(unsigned int));
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
	INDEX & operator=(INDEX &a){
		kill();
		D=a.D;
		cell=new unsigned int [D];
		for(int x=0;x<D;x++) cell[x]=a.cell[x];
		make_step_array();
		return *this;
	} 

	unsigned int operator[] (unsigned int _cell) const{
  		return cell[_cell];
	}
	unsigned int & operator[] (unsigned int _cell){
  		return cell[_cell];
	}
	
	unsigned int  * make_step_array(void){
		step = (unsigned int *) calloc ((D+1), sizeof(unsigned int));
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
		out << "(";
		for(int x=0;x<D;x++){
			out << cell[x] << ",";
		}
		out << '\b' <<")";
		
		out << ":(";
		for(int x=0;x<D;x++){
			out << step[x] << ",";
		}
		out << '\b' <<")";
		
		return out.str();
	}

	INDEX & operator = (const INDEX & ind){
		if (this!=&ind) {
			if (step!=NULL) delete step;
			if (cell!=NULL) delete cell;

			D=ind.D;
			size=ind.size;
			off_set=ind.off_set;

			step = (unsigned int *) calloc ((D+1), sizeof(unsigned int));
			cell=new unsigned int [D];
	
			for(int x=0;x<D;x++) cell[x]=ind.cell[x];
			for(int x=0;x<D+1;x++) step[x]=ind.step[x];
		}

		return *this;
	}
};

class tensor 
{
	private:

	PRECISION *cell; 	//Data member.

	unsigned int rank;	//numbr of dimensions.
	unsigned int size;	//total number of data items.
	INDEX *width;		//number of data items in a row. This should be private.

	public: 

	unsigned int *step;	//step cell.
        int a;			//shmid is the shared memeory ID, and should be called by functions whanting to know where the begin pointer is.

	PRECISION *it;		//iterator.
	PRECISION *begin, *end;	//start and stop

	tensor (unsigned int X, unsigned int * int_cell, PRECISION * _cell){
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

		cell=new PRECISION [step[rank]];

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

		cell=new PRECISION [step[rank]];
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

		cell=new PRECISION [step[rank]];
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

	tensor (char *str)
	{  
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
			fread (&a, sizeof(unsigned int), 1, TensorFile);
			(*width)[x]=a;
		}

		step = (unsigned int *) calloc ((rank+1), sizeof(unsigned int));
		step[0]=1;
		for(int x=1;x<rank+1;x++){
			step[x]=step[x-1]*(*width).cell[x-1];
		}

		cell=(PRECISION *) calloc(step[rank], sizeof(PRECISION));

		begin=&cell[0];
		end=&cell[step[rank]];
		size=step[rank];
		it=begin;
		while(it<end){
			fread (it, sizeof(PRECISION),1,TensorFile);
			it++;
		}
		fclose(TensorFile);
	}
	tensor (){
		width=NULL;
		cell=NULL;
		step=NULL;
	};

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

	std::string print (PRECISION N) const {
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

/*	std::string printRowSum (void){
		std::stringstream out(" ");
//		char buffer[20];
		it=begin;
		int x=0;
		PRECISION * sum;
		sum=new PRECISION [(*width)[0]];
		for (int x=0;x<(*width)[0];x++) sum[x]=0;
		while ((it < end)){
			sum[x]+=(*it);
			it++;
			x++;
			if(x==(*width)[0]) {x=0;}
		}
		for (int x=0;x<(*width)[0];x++){
			out << sum[x] << ", ";
//			sprintf(buffer,"%f\",sum[x]);
//			out+=buffer;
//			out+=", ";
		}
		return out.str();
	}*/

	PRECISION read(void) 
	{ 
		return (*it);
	}

	void operator= (PRECISION z)
	{
		(*it)=z;
	}

	tensor & operator= (const tensor &z)
	{
	//	cout << z.print() << endl;

		if(width!=NULL && *(width)==*(z.width)){
		//	cout << this << " fast eq\n";
			it=begin;
			PRECISION *z_it=z.begin;

			while (it!=end){
				*it=*z_it;
				it++;
				z_it++;
			}
		}
		else{
		//	cout << this << " " <<" slow eq\n";
			kill();
			width=new INDEX(z.width->D);
			*(width)=*(z.width);
			rank=z.rank;

			step = (unsigned int *) calloc ((rank+1), sizeof(unsigned int));
			step[0]=1;
			for(int x=1;x<rank+1;x++){
				step[x]=step[x-1]*(*width).cell[x-1];
			}
	
			cell=new PRECISION [step[rank]];
			begin=cell;
	
			end=&cell[step[rank]];
			size=step[rank];
	
			clear();
			it=begin;

			PRECISION *z_it=z.begin;

			while (it!=end){
				*it=*z_it;
				it++;
				z_it++;
			}
		}
		return *this;
	}

	bool operator!= (PRECISION z)
	{
		if((*it)==z) return false;
		else return true;
	}

	bool operator== (PRECISION z)
	{
		if((*it)==z) return true;
		else return false;
	}

	void operator+= (PRECISION z)
	{
		(*it)+=z;
	}

	tensor & go (INDEX *index)
	{
		it=begin;
		for(int x=0;x<rank;x++){
			it+=(step[x])*(*index)[x];
		}
		if(it>end) {cout << "out of bounds ";
		index->print();
		exit(1);}
#if (SANE)
		while(it>end) it-=step[rank];
#endif
		return *this;
	}
	
	tensor & operator() (INDEX *index)
	{
		it=begin;
		for(int x=0;x<rank;x++){
			it+=(step[x])*(*index)[x];
		}
#if (SANE)
		while(it>end) it-=step[rank];
#endif
		return *this;
	}

	tensor operator() (INDEX *index) const
	{
		PRECISION *temp;
		temp=it;
		for(int x=0;x<rank;x++){
			temp+=(step[x])*(*index)[x];
		}
#if (SANE)
		while(it>end) it-=step[rank];
#endif
		return *this;
	}

	tensor & operator() (unsigned x)
	{
		it=(begin+x);
#if (SANE)
		while(it>end) it-=step[rank];
#endif
		return *this;
	}

	PRECISION operator() (unsigned int x) const
	{
		return (*(it+x));
	}

	tensor & operator() (unsigned int x, unsigned int y)
	{
		it=(begin+x+(step[1])*y);
#if (SANE)
		while(it>end) it-=step[rank];
#endif
		return *this;
	}

	PRECISION operator() (unsigned int x, unsigned int y) const
	{
		return *(it+x+(step[1])*y);
	}

	std::string printXBarGivenY (void){
		std::stringstream out(" ");
		for(unsigned y=0;y<(*width)[1];y++){
			it=(begin+(step[1])*y);
			PRECISION XBar=0;
			PRECISION Y=0;
			for(PRECISION x=0;x<(*width)[0];x++){
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
		for(unsigned x=0;x<(*width)[0];x++){
			it=(begin+x);
			PRECISION YBar=0;
			PRECISION X=0;
			for(PRECISION y=0;y<(*width)[1];y++){
				YBar+=(*it)*y;
				X+=(*it);
				it+=step[1];
			}
			out << (YBar/X) << "," ; //<< (*it) << ",";
		}
		
		return out.str();
	};


	tensor & operator() (unsigned int w, unsigned int x, unsigned int y)
	{
		it=begin+w+(step[1])*x+(step[2])*y;
#if (SANE)
		while(it>end){
			it-=step[rank];
		}
#endif
		return (*this);
	}

	PRECISION operator() (unsigned int w, unsigned int x, unsigned int y) const
	{
		return *(it+w+step[1]*x+step[2]*y);
	}

	tensor & operator() (unsigned int w, unsigned int x, unsigned int y, unsigned int z)
	{
		it=begin+w+(step[1])*x+(step[2])*y+(step[3])*z;
#if (SANE)
		while(it>end){
			it-=step[rank];
		}
#endif
		return (*this);
	}

	PRECISION operator() (unsigned int w, unsigned int x, unsigned int y, unsigned int z) const
	{
		return *(it+w+step[1]*x+step[2]*y+step[3]*z);
	}

	PRECISION readpp()
	{
		it++;
#if (SANE)
		if (it > end) it=begin;
#endif
		return *(it-1);
	}

	void setpp(PRECISION z)
	{
		(*it)=z;
		it++;
#if (SANE)
		if (it > end) it=begin;
#endif
	}

	bool start(void){
		it=begin;
		return true;
	}

	unsigned int size_of(void){
		return step[rank];
	}

	unsigned int get_order(void){
		return rank;
	}

	unsigned int get_rank(void){
		return rank;
	}

	unsigned int get_size(void){
		return size;
	}

	PRECISION * get_cell(void){
		return cell;
	}

	unsigned int get_width(int x){
		return (*width)[x];
	}
	
	INDEX * get_width(){
		return width;
	}

	unsigned int *get_step(void){
		return step;
	}
	unsigned int get_step(int x){
		return step[x];
	}

	void normalize(void){
		PRECISION SUM=0;
		it=begin;
		while(it<end){
			SUM+=(*it);
			it++;
		}
		it=begin;
		while(it<=end){
			(*it)=(*it)/SUM;
			it++;
		}
	}

	PRECISION density(void){
		PRECISION SUM=0;
		it=begin;
		while(it<end){
			SUM+=(*it);
			it++;
		}
		return SUM;
	}

	PRECISION edge_sum(void){
		PRECISION SUM=0;

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
			fwrite (it, sizeof(PRECISION),1,TensorFile);
			it++;
		}
		fclose(TensorFile);
		return 0;
	}

	/* A[i,j]=Sum[k] Sum[l] Sum [M] B[k,l,m]*C[k,l,m,i,j];*/

	tensor & transition (tensor * t1, tensor * t2){
		PRECISION sum;
		PRECISION *t1_it, *t2_it, *t1_start, *t1_end;
	
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
		PRECISION *t1_it, *t2_it, *t2_start, *t2_end;
	
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
};

///This does not calculate a general tensor product, but instead calculates ....

void fast_product(tensor *v1, tensor *v2, tensor *p)
{
	PRECISION sum;

	PRECISION *v1_it, *v2_it, *v1_start, *v1_end, *p_it, *p_end;

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
	PRECISION sum;
	PRECISION *v1_it, *v2_it, *v1_start, *v1_end, *p_it, *p_end;

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

/*
void partial_product(tensor *v1, tensor *v2, tensor *p, )	
{

	int size=v1->size_of();
	int fork_return, fork_count=0, size_over_forks=size/forks;

	PRECISION *v1_begin=v1->begin, *v1_end=v1->end, *p_begin=p->begin, *v2_begin=v2->begin;

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

	PRECISION sum;

	PRECISION * columnEnd, *v1_it, *v2_it, *p_it, *p_end;

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
