#include "lib/tensor.h"

int main(void){
	double a=0;
	unsigned int index []={2,2};	
	tensor <double> A(2, index);
	tensor <double> B(2, index);

	cout << a << " : " << A(0,0) << endl;
	A(0,0)=1;
	a=A(0,0);
	cout << a << " : "  << A(0,0) << endl;
	a=2;
	A(0,0)=a;
	cout << a << " : "  << A(0,0) << endl;

	cout << A(0,0) << " : " << A(0,1) << " : " << A(1,0) << " : " << A(1,1) << endl;
	cout << B(0,0) << " : " << B(0,1) << " : " << B(1,0) << " : " << B(1,1) << endl;
	B=A;
	cout << A(0,0) << " : " << A(0,1) << " : " << A(1,0) << " : " << A(1,1) << endl;
	cout << B(0,0) << " : " << B(0,1) << " : " << B(1,0) << " : " << B(1,1) << endl;

}
