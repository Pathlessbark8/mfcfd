#include <iostream>  // for input output

#include "point.h"

using namespace std;


int main(int argc, char *argv[]){
    int numPart;
    Graph g;
    if (argc == 1) {
	    numPart = 1;
    }
    else{
	    numPart = atoi(argv[1]);
	    if (argc == 3) {
		    g.format = atoi(argv[2]);	    
	    }
    }

    // choose the point format
    if (g.format == 0){ // Quad tree format
	    cout << " Quadtree input format " << endl;
	    g.read_point_create_graph_quad();
    }else if(g.format == 1){ // Legacy 1
	    cout << " Legacy 1 input format " << endl;
	    g.read_point_create_graph_legacy();
    }else if(g.format == 2){ // Legacy 2
	    cout << " Legacy 2 input format " << endl;
	    g.read_point_create_graph_legacy1();
    }else if(g.format == 3){ // Restart
	    cout << " Restart input format " << endl;
	    g.read_point_create_graph_restart();
    }

    g.cal_min_dist();
    g.partition(numPart);

    // Choose output format
    if (g.format == 1){ // Legacy output with normals
	    g.write_output_legacy();
    }else { // All other formats
	    g.write_output(); 
    }


    return 0;
}

