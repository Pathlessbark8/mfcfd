#include <iostream>
#include <fstream>

#include "point.h"

using namespace std;


// generates filename in string format from partition number
string Graph::getfileName(int partID, int nParts){
	string s;
        stringstream out;
        out << partID;
        s = out.str();
        if(nParts == 1)
		return "partGrid";
        else {
		if(partID > -1 && partID < 10)
			return "partGrid0" + s;
		else 
			return "partGrid0" + s;
	}
	
}


void Graph::write_output_legacy(){
	
	/*  writing output routine  */
        // create files and write basic info
	for(int i=0; i<nParts; i++){
		string filename = getfileName(i, nParts);
		ofstream outfile (filename.c_str());
		outfile << totalPoints[i]+ghosts[i].size() << " "<< totalPoints[i] << " " << ghosts[i].size() << endl;
	}

        // local nodes
        for(int i = 0; i < nvtxs; i++){
		int partnId = part[i];
		string filename = getfileName(partnId, nParts);
		fstream write;
		write.open(filename.c_str(), fstream::app);
		write << fixed;
		write << ptVec[i].id << " " << setprecision(20)  << " " << ptVec[i].x << " " << ptVec[i].y 
			<< " " << ptVec[i].nx << " " << ptVec[i].ny<< " " << ptVec[i].flag1 << " " << ptVec[i].flag2 << " " << ptVec[i].min_dist << " ";
		write << xadjVec[i+1]-xadjVec[i] << " ";
		for(idx_t j=xadjVec[i]; j<xadjVec[i+1]; j++){
			write << inputToLoc[partnId][adjncyVec[j]+1] << " ";
		}          
		currGlobalNum[partnId]++;   
		write << endl;
		write.close(); 
	}

        // ghost nodes
        for(int i=0; i < nParts; i++) {
		string filename = getfileName(i, nParts);
		fstream write;
		write.open(filename.c_str(), fstream::app);
		set<int>::iterator itr;
		write << fixed << setprecision(20);
		for(itr = ghosts[i].begin(); itr!=ghosts[i].end(); itr++){
			write  << inputToGlob[*itr] << " " << ptVec[*itr].x << " " << ptVec[*itr].y << " " << ptVec[*itr].min_dist << "  ";
			write << endl;
		}
		write.close(); 
	}
	
}

void Graph::write_output(){
	
	/*  writing output routine  */
        // create files and write basic info
	for(int i=0; i<nParts; i++){
		string filename = getfileName(i, nParts);
		ofstream outfile (filename.c_str());
		outfile << totalPoints[i]+ghosts[i].size() << " "<< totalPoints[i] << " " << ghosts[i].size() << endl;
	}

        // local nodes
        for(int i = 0; i < nvtxs; i++){
		int partnId = part[i];
		string filename = getfileName(partnId, nParts);
		fstream write;
		write.open(filename.c_str(), fstream::app);
		write << fixed;
		write << ptVec[i].id << " " << setprecision(20)  << " " << ptVec[i].x << " " << ptVec[i].y 
			<< " " << inputToLoc[partnId][ptVec[i].left] << " " << inputToLoc[partnId][ptVec[i].right]<< " " << ptVec[i].flag1 << " " << ptVec[i].flag2 << " " << ptVec[i].min_dist << " ";
		write << xadjVec[i+1]-xadjVec[i] << " ";
		for(idx_t j=xadjVec[i]; j<xadjVec[i+1]; j++){
			write << inputToLoc[partnId][adjncyVec[j]+1] << " ";
		}          
		currGlobalNum[partnId]++;   
		write << endl;
		write.close(); 
	}

        // ghost nodes
        for(int i=0; i < nParts; i++) {
		string filename = getfileName(i, nParts);
		fstream write;
		write.open(filename.c_str(), fstream::app);
		set<int>::iterator itr;
		write << fixed << setprecision(20);
		for(itr = ghosts[i].begin(); itr!=ghosts[i].end(); itr++){
			write << inputToGlob[*itr] << " " << ptVec[*itr].x << " " << ptVec[*itr].y << " " << ptVec[*itr].min_dist << "  ";
			write << endl;
		}
		write.close(); 
	}
	
}
