
class RunStatistics{
	double *mean_W=NULL, *mean_U=NULL;
	unsigned int *count;
	list <unsigned int> extinction_time;
	RunStatistics(unsigned int length){
		mean_W=new double[length];
		mean_U=new double[length];
	};
	~RunStatistics(){
		if(mean_W!=NULL) delete [] mean_W;
		if(mean_U!=NULL) delete [] mean_U;
	}
}


