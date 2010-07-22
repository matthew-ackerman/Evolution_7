#ifndef DATA_H
#define DATA_H

/* DATA_RECORD working as of 5/21/10, do not alter with-out checking*/

class DATA_RECORD {
	private:
	public:
	vector <double> datum;
	vector <double> weight;
	string name;

	DATA_RECORD(string NAME){
		name=NAME;
	}

	DATA_RECORD(string NAME, vector <double> _double){
		name=NAME;
		weight.assign(_double.size(), 1);
		datum=_double;
	}

	DATA_RECORD(string NAME, double _double){
		name=NAME;
		weight.push_back(1);
		datum.push_back(_double);
	}

	DATA_RECORD(string NAME, double* _double, unsigned int size){
		name=NAME;
		weight.assign(size, 1);
		datum.assign(_double,_double+size);
	}

	DATA_RECORD(string NAME, vector <double> _double, unsigned int length){
		name=NAME;
		weight.assign(length, 1);
		datum.assign(_double.begin(),_double.begin()+length);
	}

	void push_back(vector <double> _double){
		datum.insert(datum.end(),_double.begin(),_double.end());
		weight.insert(weight.end(),_double.size(), 1);
	};

	void push_back(double _double){
		datum.push_back(_double);
		weight.push_back(1);
	};

	void confuse(vector <double> _double){
	//	cout << "data:" << _double.size() << " me " << datum.size() << " and " << weight.size() << endl;
		if(datum.size()>=_double.size()){
			//cout << "should be safe\n";
			for(int x=0; x<_double.size(); x++) {
	//			cout << (datum[x]*weight[x]+_double[x])/(weight[x]+1.0) << " -";
				datum[x]=(datum[x]*weight[x]+_double[x])/(weight[x]+1.0);
				weight[x]++;
			}
		}
		else{
			//cout << "should be safe\n";
			for(int x=0; x<datum.size(); x++) {
	//			cout << (datum[x]*weight[x]+_double[x])/(weight[x]+1.0) << " -";
				datum[x]=(datum[x]*weight[x]+_double[x])/(weight[x]+1.0);
				weight[x]++;
			//	cout << "end " << x << endl;
			}
			//cout << "about to segfault\n";
			datum.insert(datum.end(),_double.begin()+datum.size(), _double.end());
			//for(int x=weight.size(); x<_double.size();x++) cout << _double[x] << " *";
			weight.insert(weight.end(), (unsigned int)(_double.size()-weight.size()), 1.0);
			//cout << "totally didn't happen\n";
		}
	//	cout << endl;
	};

	void confuse(double *_double, unsigned int size){
		if(datum.size()>=size){
			for(int x=0; x<size; x++) {
				datum[x]=(datum[x]*weight[x]+_double[x])/(weight[x]+1.0); 
				weight[x]++;
			}
		}
		else{
			for(int x=0; x<datum.size(); x++) {
				datum[x]=(datum[x]*weight[x]+_double[x])/(weight[x]+1.0); 
				weight[x]++;
			}
			datum.insert(datum.end(),_double+datum.size(), _double+size);
			weight.insert(weight.end(),size-weight.size(), 1.0);
		}
	}
	
	DATA_RECORD confuse(DATA_RECORD data){
		if(datum.size()>=data.datum.size()){
			for(int x=0; x<data.datum.size(); x++) {
				datum[x]=(datum[x]*weight[x]+data.datum[x]*data.weight[x])/(weight[x]+data.weight[x]); 
				weight[x]+=data.weight[x]; 
			}
		}
		else{
			for(int x=0; x<datum.size(); x++) {
				datum[x]=(datum[x]*weight[x]+data.datum[x]*data.weight[x])/(weight[x]+data.weight[x]);
				weight[x]+=data.weight[x];
			}
			datum.insert(datum.end(), data.datum.begin()+datum.size(), data.datum.begin()+data.datum.size());
			weight.insert(weight.end(), data.weight.begin()+weight.size(), data.weight.begin()+data.weight.size());
		}
		return *this;
	}

	void write_to_binary(FILE * pFile){
		//if (!pFile) throw std::runtime_error( "Could not write file!" );  // don't halt your program, just complain 

		unsigned int size=datum.size();
		fwrite(&size, sizeof(unsigned int), 1, pFile);

		//TODO: write much better error check;
		size%=13;
		fwrite(&size, sizeof(unsigned int), 1, pFile);
		size=datum.size()%11;
		fwrite(&size, sizeof(unsigned int), 1, pFile);

		size=datum.size();
  		for (unsigned n = 0; n < size; n++) fwrite(&(datum[n]), sizeof(double), 1, pFile);
  		for (unsigned n = 0; n < size; n++) fwrite(&(weight[n]), sizeof(double), 1, pFile);		

	}

	void read_from_binary(FILE * pFile){
		//if (!pFile) throw std::runtime_error( "Could not read file!" );  // don't halt your program, just complain 
		//read hash_function;	
		unsigned int size, E13, E11;
		fread(&size, sizeof(unsigned int), 1, pFile);
		fread(&E13, sizeof(unsigned int), 1, pFile);
		fread(&E11, sizeof(unsigned int), 1, pFile);

		if(size%13==E13 && size%11==E11){
			double r;
			datum.clear();
			weight.clear();
  			for (unsigned n = 0; n < size; n++){
				fread(&r, sizeof(double),1, pFile);
				datum.push_back(r);
			}
  			for (unsigned n = 0; n < size; n++){
				fread(&r, sizeof(double),1, pFile);
				weight.push_back(r);
			}
		}
		//else throw std::runtime_error( "Data appears to be corrupt, or reading out of frame.\n");
	}

	void print_to_csv(FILE * pFile){
		//if (!pFile) throw std::runtime_error( "Could not write file!" );  // don't halt your program, just complain
 
		if (datum.size()!=weight.size()) {cout << "Data corrupt.\n"; exit(0);}
		if (datum.size()==0||weight.size()==0) {cout << "Data empty.\n"; exit(0);}


		fprintf (pFile, "\'%s value: \', %f", name.c_str(), datum[0]);
		for(int y=1;y<datum.size();y++){
			fprintf (pFile, ",%f", datum[y]);
		}
		fprintf (pFile, "\n\'%s weight:\', %f", name.c_str(), weight[0]);
		for(int y=1;y<weight.size();y++){
			fprintf (pFile, ",%f", weight[y]);
		}
		fprintf (pFile, "\n");
	}
	void print_header_to_csv(FILE * pFile){
		if (datum.size()!=weight.size()) {cout << "Data corrupt.\n"; exit(0);}
		if (datum.size()==0||weight.size()==0) {cout << "Data empty.\n"; exit(0);}

		fprintf (pFile, "\'%s cell #: \', %d", name.c_str(), 0);
		for(int y=1;y<datum.size();y++){
			fprintf (pFile, ",%d", y);
		}
		fprintf (pFile, "\n");
	}
	
	void read_from_csv(FILE *pFile){
	}
	
	void clear (void){
		datum.clear();
		weight.clear();
	}

	~DATA_RECORD (){
		datum.clear();
		weight.clear();
	}	
};

#endif
