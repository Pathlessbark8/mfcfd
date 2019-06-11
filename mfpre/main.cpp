#include <iostream> 
#include <cstring>
#include <sys/stat.h>
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
          if(strcmp(argv[i],"-quadtree") == 0){
		  g.format = 1;
	  }
          else if(strcmp(argv[i],"-legacy") == 0){
		  g.format = 2;
	  }
          else if(strcmp(argv[i],"-restart") == 0){
		  g.format = 3;
	  }
          else if(strcmp(argv[i],"-cuda") == 0){
		  if(g.format == 1){
			  g.gpu = 1;
		  }
		  else if(g.format == 2){
			  g.gpu = 2;
		  }
	  }
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
    }else if(g.format == 2){ // Legacy
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

    // Creating a directory 
    if (mkdir("point", 0777) == -1) 
        cerr << "Error :  " << strerror(errno) << endl; 
  
    else
        cout << " Directory created " << endl; 

    // Choose output format
    if (g.gpu == 1){ // gpu output for quadtree
	    cout << " Writing quadtree cuda output " << endl;
	    g.write_output_gpu_quad();
    }
    else if (g.gpu == 2){ // gpu output for legacy
	    cout << " Writing legacy cuda output " << endl;
	    g.write_output_gpu_legacy();
    }else if (g.format == 1) { // quadtree mpi
	    cout << " Writing quadtree mpi output " << endl;
	    g.write_output_quad(); 
    }else if (g.format == 2) { // Legacy mpi
       	    cout << "Writing legacy mpi output" << endl;
	    g.write_output_legacy();
    }
    

    return 0;
}

