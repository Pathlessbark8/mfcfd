#include <metis.h>
#include <iostream>
#include <iomanip> 
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <set>
#include <map>
#include <climits>
#include <cassert>

using namespace std;

typedef struct point {
    idx_t id;
    double x, y;
    int left, right;
    int qtdepth;
    int nx, ny;
    double min_dist;
    idx_t flag1, flag2;
} point;

class Graph{
	public:
		idx_t nvtxs;
		idx_t nParts;
		int * totalPoints;
		int * currGlobalNum;
		idx_t * part;
		set<int> * ghosts;
		map<int, int> * inputToLoc;
		map<int, int> inputToGlob;
		vector<point> ptVec;
		vector<idx_t> xadjVec;
		vector<idx_t> adjncyVec;
		int format;
		int gpu;
		
		// read file and build graph
		Graph(){ format = 0; gpu = 0; };

		void read_point_create_graph_legacy();
		void read_point_create_graph_quad();
		void read_point_create_graph_restart();
		void cal_min_dist();
		void partition(int);
		string getfileName(int, int );
		void write_output_gpu_legacy();
		void write_output_gpu_quad();
		void write_output_legacy();
		void write_output_quad();
		 
};
