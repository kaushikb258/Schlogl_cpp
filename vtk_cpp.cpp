#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include "vtk_cpp.h"
#include "classes.h"
using namespace std;

void vtk_write(string filename,string title,int imax, int jmax, double* xyz, int nprim, double* prim);

void output_vtk(int i_max,int j_max,const double dx, GS* gs)
{
 string filename = "output_paraview.vtk";
 string title = "Schlogl2d";
 double coord[i_max][j_max][3];
 int i, j, l;
 int nprim = 2;
 double q[i_max][j_max][nprim];

// initialize coordinates
 for (i=0; i<i_max; i++)
 {
 for (j=0; j<j_max; j++)
 {
  coord[i][j][0] = ((double) (i) + 0.5)*dx; //x
  coord[i][j][1] = ((double) (j) + 0.5)*dx; //y
  coord[i][j][2] = 0.0; //z = 0 as 2D grid
 }
 }

 for (i=0; i<i_max; i++)
 {
 for (j=0; j<j_max; j++)
 {
  l = j + i*j_max;
  q[i][j][0] = gs[l].get_u();
  q[i][j][1] = gs[l].get_v();
 }
 }

 vtk_write(filename,title,i_max,j_max,&coord[0][0][0],nprim,&q[0][0][0]);
}


void vtk_write(string filename,string title,int imax, int jmax, double* xyz, int nprim, double* prim)
{
string cell_size_string, node_num_string;
int node;
int i, j, k, n;
stringstream s_node_num, s_cells, s_imax, s_jmax, s_kmax;
ofstream f;
int kmax = 1;
s_node_num << (imax*jmax);
s_cells << (imax-1)*(jmax-1);
s_imax << imax;
s_jmax << jmax;
s_kmax << kmax;
f.open (filename.c_str(),ios_base::out);
f<< "# vtk DataFile Version 2.0\n";
f<< title<<"\n";
f<< "ASCII\n";
f<< "DATASET STRUCTURED_GRID\n";
f<< "DIMENSIONS "<<"\t"<<s_imax.str()<<"\t\t"<<s_jmax.str()<<"\t\t"<<s_kmax.str()<<"\n";
f<< "POINTS "<<"\t"<<s_node_num.str()<<"\t"<<"double\n";
n = 0;
for (i=0; i<imax; i++)
{
for (j=0; j<jmax; j++)
{
for (k=0; k<3; k++)
{
f<<xyz[n]<<" ";
n++;
}
f<<"\n";
}
}
f<< "CELL_DATA "<<"\t"<<s_cells.str()<<"\n";
f<< "POINT_DATA "<<"\t"<<s_node_num.str()<<"\n";
f<< "SCALARS u double \n";
f<< "LOOKUP_TABLE default \n";
n = 0;
for (i=0; i<imax; i++)
{
for (j=0; j<jmax; j++)
{
f<< prim[n]<<" ";
n += nprim;
}
}
f<<"\n";
f<< "SCALARS v double \n";
f<< "LOOKUP_TABLE default \n";
n = 1;
for (i=0; i<imax; i++)
{
for (j=0; j<jmax; j++)
{
f<< prim[n]<<" ";
n += nprim;
}
}
f.close();
}

