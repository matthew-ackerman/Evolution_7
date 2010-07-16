#ifndef MERSENNETWISTER_H
#include "lib/MersenneTwister.h"
#endif

#ifndef UTILITIES_H
#define UTILITIES_H

inline void randomizeArray(ORGANISM **a, unsigned int size, MTRand *mtrand){
	ORGANISM *swap;
	int rand;
//	cout << "size " << size << endl;
	for (int x=size;(--x)>-1;){
		rand=mtrand->randInt(x);
		swap=a[rand];
		a[rand]=a[x];
		a[x]=swap;
	} 
}

#endif
