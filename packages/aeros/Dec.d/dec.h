#ifndef DEC_H
#define DEC_H

namespace Dec
{
int dec(int numProcessors,             // aka -p <number of processors>
	int numThreads,                // aka -n <number of threads>
	int numSubdomains,             // aka -s <number of subdomains> 
	//string decompositionFileName,  // aka -d <decomposition file name>
	int topFlag                // -1 --> do nothing
	                               // 0   ~   aka -t
	                               // 1   ~   aka -T
	                               // 2   ~   aka -m
	);

}

#endif
