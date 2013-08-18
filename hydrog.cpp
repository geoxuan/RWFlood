/** RWFlood 1.0. Date of release: July 25th, 2013. 
 *
 * 
 * Copyright (C) 2013 Salles V. G. Magalhães, Marcus V. A. Andrade, W. Randolph Franklin, 
 * and Guilherme C. Pena.
 * 
 * This program is free software: you can redistribute it and/or modify it 
 * under the terms of the GNU General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * You should have received a copy of the GNU General Public License along 
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * RWFlood - A new and faster internal memory method to compute the drainage network, 
 * that is, the flow direction and accumulation on terrains represented by raster elevation matrix.
 * 
 * For further information, consider the following reference:
 * 
 * [1] Magalhães, S. V., Andrade, M. V., Franklin, W. R., & Pena, G. C. (2012). A new method for computing the drainage 
 * network based on raising the level of an ocean surrounding the terrain. In Bridging the Geographic Information Sciences 
 * (pp. 391-407). Springer Berlin Heidelberg
 * 
 */


#include <iostream>
#include <algorithm>
#include <queue>
#include <cmath>
#include <fstream>
using namespace std;

// http://en.wikipedia.org/wiki/ANSI_escape_code
const string redtty("\033[1;31m");   // tell tty to switch to red
const string greentty("\033[1;32m");   // tell tty to switch to bright green
const string bluetty("\033[1;34m");   // tell tty to switch to bright blue
const string magentatty("\033[1;35m");   // tell tty to switch to bright magenta
const string yellowbgtty("\033[1;43m");   // tell tty to switch to bright yellow background
const string underlinetty("\033[4m");   // tell tty to switch to underline
const string deftty("\033[0m");      // tell tty to switch back to default color

// Exec a string then print its time, e.g. TIME(init());
#define TIME(arg) { arg;  } ; ptime(#arg); 

time_t tStart;

void ptime(const char *const msg) {
  float t= ((float)clock())/CLOCKS_PER_SEC;
  float tr = difftime(time(NULL),tStart);
  cerr << magentatty << "Cumulative CPU time thru " << msg << "= " << t << "   |  Real time= " << tr << deftty << endl;
  cout  << "Cumulative CPU time thru " << msg << "= " << t << "   |  Real time= " << tr  << endl;
}

// Create a new array of size var of type.  Report error.
#define NEWA(var,type,size) { try  { if(0==(var=new type [(size)])) throw;} catch (...)  { cerr << "NEWA failed on " #var "=new " #type "[" #size "=" << (size) << "]" << endl; exit(1); }}


typedef  pair<int,int> pii;

class Point {
public:
  unsigned short int x,y;
 // short int z;
  
  Point(const short int y1, const short int x1) : 
    x(x1),y(y1) {}
};

/*bool operator<(const Point &a, const Point &b) {
  return a.z > b.z;
}*/




//Uses Young-He's algorithm to flood a terrain...
//Ref: http://www.google.com.br/url?sa=t&source=web&cd=3&ved=0CDIQFjAC&url=http%3A%2F%2Fwww.sciencedirect.com%2Fscience%2Farticle%2Fpii%2FS0098300405001937&ei=6FD2Tc-8JYSGhQfKxtjgBg&usg=AFQjCNF4EmCGPPtO6jCH1AMXyCzZbQcrNQ&sig2=oG8MIndi0fXNdixh_nXXvw

//Allocates matrix: dirs[nrows][nrows]

queue<Point> seaLevelQueues[65536];
void flood(short int **elevs,unsigned char **&dirs,const int nrows) {
	dirs = new unsigned char*[nrows]; 
  	for(int i=0;i<nrows;i++) {
		NEWA(dirs[i],unsigned char,nrows);
		for(int j=0;j<nrows;j++)
			dirs[i][j] = 0;
  	}

	//priority_queue<Point> seaLevelPQ;

	for (int i=0;i<nrows;i++) {		
		//seaLevelPQ.push( Point(i,0,elevs[i][0]) );
		//seaLevelPQ.push( Point(i,nrows-1,elevs[i][nrows-1]) );
		seaLevelQueues[elevs[i][0] + 32768].push(Point(i,0));
		seaLevelQueues[elevs[i][nrows-1]+ 32768].push(Point(i,nrows-1));

		dirs[i][0] = 128;
		dirs[i][nrows-1] = 8;
	}
	for (int i=1;i<nrows-1;i++) {
		//seaLevelPQ.push( Point(0,i,elevs[0][i]) );
		//seaLevelPQ.push( Point(nrows-1,i,elevs[nrows-1][i]) );
		seaLevelQueues[elevs[0][i]+ 32768].push(Point(0,i));
		seaLevelQueues[elevs[nrows-1][i]+ 32768].push(Point(nrows-1,i));

		dirs[0][i] =  2;
		dirs[nrows-1][i] = 32;
	}	

	
	int seaLevel = -32768;
	for(;seaLevel<=32767;seaLevel++) {
		while (!seaLevelQueues[seaLevel+32768].empty()) {	
			Point p = seaLevelQueues[seaLevel+32768].front(); seaLevelQueues[seaLevel+32768].pop();
		
			//seaLevel = p.z;

			if (p.y > 0 ) {
				if (p.x >0) {
					if (dirs[p.y-1][p.x-1]==0) {		
						dirs[p.y-1][p.x-1] = 16;		

						if ( elevs[ p.y-1][p.x-1 ] <= seaLevel )
							 elevs[ p.y-1][p.x-1 ] = seaLevel;
						seaLevelQueues[ 32768 +elevs[ p.y-1][p.x-1 ]].push( Point(p.y-1,p.x-1 ) );	
					}
				}

				if (dirs[p.y-1][p.x]==0)  {	
					dirs[p.y-1][p.x] = 32;		
	
					if ( elevs[ p.y-1][p.x ] <= seaLevel )
						 elevs[ p.y-1][p.x] = seaLevel;
					seaLevelQueues[ 32768 +elevs[ p.y-1][p.x ]].push( Point(p.y-1,p.x ) );	
				}

				if (p.x < nrows-1) {
					if (dirs[p.y-1][p.x+1]==0) {
						dirs[p.y-1][p.x+1] = 64;

						if ( elevs[ p.y-1][p.x+1 ] <= seaLevel )
							 elevs[ p.y-1][p.x+1 ] = seaLevel;
						seaLevelQueues[ 32768 +elevs[ p.y-1][p.x+1 ]].push( Point(p.y-1,p.x+1) );
					}	
				}
			} 
			if (p.x >0) {
				if (dirs[p.y][p.x-1]==0) {
					dirs[p.y][p.x-1] = 8;

					if ( elevs[ p.y][p.x-1 ] <= seaLevel )
						elevs[ p.y][p.x-1 ] = seaLevel;
					seaLevelQueues[ 32768 +elevs[p.y][p.x-1 ]  ].push( Point(p.y,p.x-1 ) );	
				}
			}
			if (p.x < nrows-1) {
				if (dirs[p.y][p.x+1]==0) { 
					dirs[p.y][p.x+1] = 128;

					if ( elevs[ p.y][p.x+1 ] <= seaLevel )
						elevs[ p.y][p.x+1 ] = seaLevel;
					seaLevelQueues[ 32768 +elevs[ p.y][p.x+1 ] ].push( Point(p.y,p.x+1 ) );	
				}
			}
	
			if (p.y < nrows-1 ) {
				if (p.x >0) {
					if (dirs[p.y+1][p.x-1]==0)  {
						dirs[p.y+1][p.x-1] = 4;

						if ( elevs[ p.y+1][p.x-1 ] <= seaLevel )
							 elevs[ p.y+1][p.x-1 ] = seaLevel;
						seaLevelQueues[ 32768 +elevs[ p.y+1][p.x-1 ] ].push( Point(p.y+1,p.x-1 ) );	
					}
				}

				if (dirs[p.y+1][p.x]==0)  {
					dirs[p.y+1][p.x] = 2;

					if ( elevs[ p.y+1][p.x ] <= seaLevel )
						 elevs[ p.y+1][p.x] = seaLevel;
					seaLevelQueues[ 32768 +elevs[ p.y+1][p.x ] ].push( Point(p.y+1,p.x ) );	
				}

				if (p.x < nrows-1) {
					if (dirs[p.y+1][p.x+1]==0) {
						dirs[p.y+1][p.x+1] = 1;

						if ( elevs[ p.y+1][p.x+1 ] <= seaLevel )
							 elevs[ p.y+1][p.x+1 ] = seaLevel;
						seaLevelQueues[ 32768 + elevs[ p.y+1][p.x+1 ] ].push( Point(p.y+1,p.x+1 ) );	
					}
				}
			} 
		}

	}
}










//Uses a "topological sorting" to compute flow accumulation
//Allocates matrix: flow[nrows][nrows]
void computeFlow(unsigned char **dirs,int **&flow,const int nrows) {
	flow = new int*[nrows];
  	for(int i=0;i<nrows;i++) {
		NEWA(flow[i],int,nrows);		
  	}	

	unsigned char **inputDegree = new unsigned char*[nrows]; //Matrix that represents the input degree of each point in the terrain
                      //the input degree of a point c represents the number of points in the 
                      //terrain that flows directly to c
  	for(int i=0;i<nrows;i++) {
		NEWA(inputDegree[i],unsigned char,nrows);
		for(int j=0;j<nrows;j++) inputDegree[i][j] = 0;
  	}


	int directionToDyDx[129][2]; //Represents the deltaY and deltaX of each direction code...
		         	     //For example, directionToDyDx[64][0] == 1 (deltaY) and directionToDyDx[64][1] == -1 (deltaX)
	directionToDyDx[1][0] = -1;
	directionToDyDx[1][1] = -1;
	directionToDyDx[2][0] = -1;
	directionToDyDx[2][1] = 0;
	directionToDyDx[4][0] = -1;
	directionToDyDx[4][1] = +1;

	directionToDyDx[128][0] = 0;
	directionToDyDx[128][1] = -1;
	directionToDyDx[8][0] = 0;
	directionToDyDx[8][1] = 1;


	directionToDyDx[64][0] = 1;
	directionToDyDx[64][1] = -1;
	directionToDyDx[32][0] = 1;
	directionToDyDx[32][1] = 0;
	directionToDyDx[16][0] = 1;
	directionToDyDx[16][1] = +1;

	#define PROCESSED_POINT_FLAG 255


	for(int i=0;i<nrows;i++)
		for(int j=0;j<nrows;j++) {	
			flow[i][j] = 1;		
			unsigned int ptDir = (unsigned int) dirs[i][j];
			if (ptDir!=0) {			
				if (	i+directionToDyDx[ptDir][0] <0  || j+directionToDyDx[ptDir][1] <0 || i+directionToDyDx[ptDir][0] >=nrows  || j+directionToDyDx[ptDir][1] >=nrows ) 
					continue;	
				inputDegree[ i+directionToDyDx[ptDir][0] ][ j+directionToDyDx[ptDir][1] ]++; //Increases the input degree of the point that receives the flow from point (i,j)
			}
		}



	for(int i=0;i<nrows;i++)
		for(int j=0;j<nrows;j++) {	
				/*if (inputDegree[i][j]==PROCESSED_POINT_FLAG) //If the point has already been processed.... 
					continue;*/
	
				pii point(i,j);
				//Distributes the flow from point (i,j) to the neighbor of (i,j) that will receive its flow				
				while(inputDegree[ point.first ][ point.second ]==0) {
					inputDegree[ point.first ][ point.second ] = PROCESSED_POINT_FLAG; //Mark this point as processed	
					int dirPt = dirs[point.first][point.second];	
					pii neighbor( point.first+directionToDyDx[dirPt][0], point.second+directionToDyDx[dirPt][1] );
					if (neighbor.first <0  || neighbor.second <0 || neighbor.first >=nrows  || neighbor.second >=nrows ) 
						break;

					
					flow[ neighbor.first ][ neighbor.second ]+= flow[point.first][point.second];
					inputDegree[ neighbor.first ][ neighbor.second ]--;
					point = neighbor; //Repeat the process with the neighbor (if it has input degree equal to 0).
				}			
		}


	for(int i=0;i< nrows;i++) delete []inputDegree[i];
	delete inputDegree;
}



//Allocates matrix: elevs[nrows][nrows]
void readElevs(short int **&elevs,const int nrows) {
  elevs = new short int*[nrows];
  for(int i=0;i<nrows;i++) {
	NEWA(elevs[i],short int,nrows);
  }
  for (int i=0;i<nrows;i++) {
    cin.read(reinterpret_cast<char*>(elevs[i]), 2*nrows);    
  }
}

//Allocates matrix: elevs[nrows][nrows] 
void init(const int argc, const char **argv,short int **&elevs,int &nrows) {
	if (argc!=2) {
		cerr << "Error, use: hydrog nrows <input" << endl;
		exit(1);
	}

	tStart = time(NULL);
	nrows = atoi(argv[1]);
	readElevs(elevs,nrows);
}

/*
void writeFlowPgm(int **flow,const int nrows) { 
  ofstream fcomp("flow.pgm");

  int maxflow = 0;
  for(int i=0;i<nrows;i++)
	for(int j=0;j<nrows;j++)
		maxflow = max(maxflow,flow[i][j]);
  if (fcomp) { 
    const float maxgrey(sqrt((float)maxflow));
    fcomp << "P2" << endl;
    fcomp << nrows << ' ' << nrows << endl;
    fcomp << (int(maxgrey)) << endl;
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<nrows; j++) fcomp << (int(sqrt(float(flow[i][j])))) << ' ';
      fcomp << '\n';
    }
    fcomp << flush;
  } else {
	cerr << "Error writing flow.pgm file" << endl;
	exit(1);
  }
}*/

void writeFlowPgm(int **flow,const int nrows) { 
  ofstream fcomp("flow.pgm");

  int maxflow = 0;
  for(int i=0;i<nrows;i++)
	for(int j=0;j<nrows;j++)
		maxflow = max(maxflow,flow[i][j]);
  if (fcomp) { 
    const float maxgrey(sqrt((float)maxflow));
    fcomp << "P5" << endl;
    fcomp << nrows << ' ' << nrows << endl;
    fcomp << (int(maxgrey)) << endl;
    
    if ((int(maxgrey)) <=255) {
	    for (int i=0; i<nrows; i++) {
	      for (int j=0; j<nrows; j++) {
			unsigned char tmp = (sqrt(float(flow[i][j])));
			fcomp.write(reinterpret_cast<char *>(&tmp),1);
	      }
	    }

    }
    else {
	    for (int i=0; i<nrows; i++) {
	      for (int j=0; j<nrows; j++) {
			unsigned short tmp = (sqrt(float(flow[i][j])));
			tmp = (tmp&255)<<8 | ((tmp>>8)&255);
			fcomp.write(reinterpret_cast<char *>(&tmp),2);

	      }
	    }
    }
    fcomp << flush;
  } else {
	cerr << "Error writing flow.pgm file" << endl;
	exit(1);
  }



}




int main(const int argc, const char **argv) {
	int nrows; //Number of rows in the terrain
	short int **elevs; //Elevation matrix, (we use signed shorts to represent elevation)

	
	TIME(init(argc,argv,elevs,nrows));


	unsigned char **dirs; //Flow direction matrix
		             // The direction of a cell c may be 1,2,4,8,16,32,64 or 128
		             // 1 means c points to the upper left cell, 2 -> upper cell, 
		             // 4 -> upper right , ...
		             // 1   2   4
		             // 8   c  16
		             // 32 64 128

  
	
	//Flood depressions and Compute directions
	TIME(flood(elevs,dirs,nrows));

	for(int i=0;i< nrows;i++) delete []elevs[i];
	delete elevs;	

	//Compute flow...
	int **flow;
	TIME(computeFlow(dirs,flow,nrows));

	for(int i=0;i< nrows;i++) delete []dirs[i];
	delete dirs;

	//Write results...
	TIME(writeFlowPgm(flow,nrows));
	for(int i=0;i< nrows;i++) delete []flow[i];
	delete flow;
	
}
