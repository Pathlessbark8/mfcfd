#include <iostream> 
#include <cstring>

#include "point.h"

using namespace std;


int main(int argc, char *argv[]){
    int numPart;
    Graph g;
    if (argc == 1)
       numPart = 1;
    else
    {
       for(int i=1; i<argc; ++i)
       {
          if(strcmp(argv[i],"-quadtree") == 0)
             g.format = 1;
          else if(strcmp(argv[i],"-legacy") == 0)
             g.format = 2;
          else if(strcmp(argv[i],"-restart") == 0)
             g.format = 3;
          else if(strcmp(argv[i],"-cuda") == 0)
             g.gpu = true;
          else
          {
            assert(atoi(argv[i]) > 0);
            numPart = atoi(argv[i]);
          }
       }
    }

    // choose the point format
    if (g.format == 1){ // Quad tree format
	    cout << " Quadtree input format " << endl;
	    g.read_point_create_graph_quad();
    }else if(g.format == 2){ // Legacy 2
	    cout << " Legacy input format " << endl;
	    g.read_point_create_graph_legacy();
    }else if(g.format == 3){ // Restart
	    cout << " Restart input format " << endl;
	    g.read_point_create_graph_restart();
    }else{
	    cout << "No input format chosen or invalid input format" << endl;
	    exit(0);
    }

    g.cal_min_dist();
    g.partition(numPart);

    // Choose output format
    if (g.gpu){ // gpu output
	    g.write_output_gpu();
    }else { // All other formats
	    g.write_output(); 
    }


    return 0;
}

