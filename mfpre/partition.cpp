#include <metis.h>
#include <typeinfo>

#include "point.h"

using namespace std;

// the main metis partitioning
void Graph::partition(int n){

    /*  partitioning routine  */
    nParts    = n; // number of parts to partition the graph
    idx_t nWeights  = 1; // number of balancing constraints, should be at least 1
    part = new idx_t[nvtxs]; // stores the partition vector of the graph
    idx_t objval; 
    idx_t *xadj = &xadjVec[0]; // Indexes of starting points in adjacent array
    idx_t *adjncy = &adjncyVec[0]; // Adjacent vertices in consecutive index order
    if(n==1){
        for(int i=0; i<nvtxs; i++){
            part[i] = 0;
        }
    } else {
        int out = METIS_PartGraphKway(&nvtxs, &nWeights, xadj, adjncy,
                    NULL, NULL, NULL, &nParts, NULL,
                    NULL, NULL, &objval, part);
        std::cout << "output code: " << out << std::endl;
        std::cout << "total edge cuts: " << objval << std::endl;
    }

    /*  output processing routine  */
    // creates sets of ghost nodes for each partition
    ghosts = new set<int>[nParts];
    totalPoints = new int[nParts];
    fill_n(totalPoints,nParts,0);
    for(int i = 0; i < nvtxs; i++){
	    totalPoints[part[i]]++;
	    int partn = part[i];
	    for(idx_t j = xadjVec[i]; j < xadjVec[i+1]; j++){
		    int nbr = adjncyVec[(int)(j)];
		    if(part[nbr] != partn){
			    if(ghosts[part[i]].find(nbr) == ghosts[part[i]].end()){
				    ghosts[part[i]].insert(nbr);
			    }
		    }
	    }
    }

    
    // creates mapping of local and input numberings for each partition
    inputToLoc = new map<int, int>[nParts];
    int currLocalIndex[nParts];
    for(int i=0; i<nParts; i++)
	    currLocalIndex[i] = 1;
    for(int i = 0; i < nvtxs; i++){
	    inputToLoc[ part[i] ][ ptVec[i].id ] = currLocalIndex[ part[i] ];
	    currLocalIndex[ part[i] ]++;
    }
    
    for(int i=0; i < nParts; i++) {
	    set<int>::iterator itr;
	    for(itr = ghosts[i].begin(); itr!=ghosts[i].end(); itr++){
		    inputToLoc[i][ ptVec[*itr].id ] = currLocalIndex[i];
		    currLocalIndex[i]++;
	    }
    }

    // creates mapping of local and global numberings over all partitions
    currGlobalNum = new int[nParts];
    int globNumCounter = 1;
    for(int i=0; i<nParts; i++){
	    currGlobalNum[i] = globNumCounter;
	    globNumCounter += totalPoints[i];
    }
    for(int i = 0; i < nvtxs; i++){
	    inputToGlob[i]= currGlobalNum[ part[i] ];
	    currGlobalNum[ part[i] ]++;
    }

}
