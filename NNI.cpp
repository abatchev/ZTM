/// No non-STL dependencies. Language standard: C++14.
/// This library is released for public use as material supplementary to PhD thesis: Development of Specialized Non-Linear Inversion Algorithms, Basis Functions, Eikonal Solvers, and Their Integration for Use in Joint Seismic and Gravitational Tomographic Inversion
/// Author: Zagid Abatchev - abatchev@ucla.edu ; PI: Paul Davis. Department of Earth, Planetary, and Space Sciences; University of California, Los Angeles; Date of publication: March 2019.///
/// This library is published as is and may have unresolved bugs. Use at own discretion. Users are encouraged to adapt some or all of this code without restriction.  Refer to Github page for updated source code, located at: https://github.com/abatchev/ZTM .

/// NNI.cpp is a self contained C++14  Natural Neighbor Interpolation function library using C++ STL.///
/// NNI Implements Natural Neighbor/Discrete Sibson Interpolation (Park, 2006) to render scalar fields through a Cartesian grid using a set of input Voronoi nodes with defined coordinates and scalar field amplitudes. ///
/// Includes overloaded 2D and 3D user callable functions NNI(), and multi threaded 3D variant NNIm().
/// No non-STL dependencies. Language standard: C++14. Code may be compiled without multi threaded functionality by commenting out NNIm. Multi threaded compilation on Ubuntu can be done using -lpthread flag in GCC.
/// Format: 3D single threaded: vector<double> out = NNI(vector<double> b, int di, int dj, int dk);
/// Format: 3D multi threaded: vector<double> out = NNIm(vector<double> b, int di, int dj, int dk);
/// Format: 2D single threaded: vector<double> out = NNI(vector<double> b, int di, int dj);
/// Variables:
///     int di,dj,dk: i axis dimension, j axis dimension, k axis dimension (when applicable);
///     vector<double> b : flattened Voronoi node array with the following indexing:
///         3D case: b[n%4==0] = i coordinate; b[n%4==1] = j coordinate; b[n%4==2] = k coordinate; b[n%4==3] = voronoi node scalar field value;
///         2D case: b[n%3==0] = i coordinate; b[n%3==1] = j coordinate; b[n%3==2] = voronoi node scalar field value;
///     Axes {x,y,z} correspond to indices {i,j,k} with dimensions {di,dj,dk} and spacing 'h'. {x=i*h; y=j*h; z=k*h}; Mapping to 1D: n[i,j,k]=i+j*di+k*di*dj; to 3D: i=(n%di), j=((n/di)%dj)/h , k=(n/(di*dj)).
///     Rendered field indexing: out[n]=out[i+j*di+k*di*dj] where: x[n]=(n%di)/h, y[n]=((n/di)%dj)/h, z[n]=(n/(di*dj))/h.

#include <math.h>
#include <vector>
#include <algorithm>
using std::vector; using std::max; using std::min;

///3D multi threaded implementation of NNI
vector<double> NNIm(vector<double> b, int di, int dj, int dk){                 //Multi Threaded Natural Neighbor Interpolation function NNImt(): Identical in I/O to NNI(), runs on multiple CPU threads.
    int dn=di*dj*dk, db=b.size(),N,K; if(db<4){return vector<double>(dn);}     //Record array sizes dn(grid points) and db(b vector size). If b is empty, return zero array.
    vector<double> U(dn),C(dn); vector<vector<double>>u(dk,U),c(dk,C);         //Initialize vectors U and C to hold Cartesian grid of slownesses, vector u of vectors U and vector c of vectors C to hold cartesian grid of slownesses and hit counts from each thread.
    auto F=[&](int K){for(int J=0;J<dj;J++){for(int I=0;I<di;I++){             //Iterate through every grid/raster index in set 'u', computing its nearest neighbor field value 'vf' and distance, then distributing.
        double f=1,d=1e+18,t,s; int n,m,i,j,k,l;                               //Declare local variables, Initialize associated Voronoi node field value 'f' at 1, distance from point [n] to the Voronoi node 'd', and square of distance 't'.
        ///For each raster point n, determine the field value of, and distance to the nearest Voronoi node in b[].///
        for(m=0; m<db; m+=4){                                                  //Iterate through each Voronoi node to find grid point n's distance to and scalar value of nearest neighbor.
            t=pow(b[m]-I,2)+pow(b[m+1]-J,2)+pow(b[m+2]-K,2);                   //Compute grid index's square of distance to each Voronoi node in b as 't'.
            if(t<d){d=t; f=b[m+3];}} l=sqrt(d)+1;                              //If t<d, set Voronoi node distance square d=t, 'f' to cell field value. After loop, set distance s=sqrt(d).
        ///For each raster point n, distribute the field value of its associated (nearest) Voronoi node to all raster points within the distance between point 'n' and its associated Voronoi node 'd'.///
        for(i=max(0,I-l); i<=min(I+l,di-1); i++){                              //Iterate through axis i indices within l indices of contributing grid point.
            for(j=max(0,J-l); j<=min(J+l,dj-1); j++){                          //Iterate through axis j indices within l indices of contributing grid point.
                for(k=max(0,K-l); k<=min(K+l,dk-1); k++){                      //Iterate through axis k indices within l indices of contributing grid point.
                    if(d>=pow(I-i,2)+pow(J-j,2)+pow(K-k,2)){                   //Check if distance to nearest voroni node is less than "d".
                        n=i+j*di+k*di*dj; u[K][n]+=f; c[K][n]++;}}}}}}};       //Increment field value u[ii,jj,kk] by the field value at the originating point u[i,j,k], contributing points counter at c[ii,jj,kk] by 1.
    thread R[dk]; for(K=0;K<dk;K++){R[K]=thread(F,K);} for(auto&T:R){T.join();}//Execute k threads corresponding to k slices through the grid spanned by 'u', then join them.
    for(N=0;N<dn;N++){for(K=0;K<dk;K++){U[N]+=u[K][N]; C[N]+=c[K][N];}}        //Merge grids from every thread by summation.
    for(N=0;N<dn;N++){U[N]/=C[N];} return U;}                                  //Normalize the slowness field by dividing each vapue in 'u[n]' by the number of contributing Voronoi nodes 'c[n]'(if any). Return resultant field U.

///3D implementation of NNI
vector<double> NNI(vector<double> b, int di, int dj, int dk){
    int dn=di*dj*dk,db=b.size(); vector<double>u(dn),c(dn); if(db<4){return u;}//Record array sizes dn(grid points) and db(b vector size). If b is empty, return zero array.
    for(int K=0;K<dk;K++){for(int J=0;J<dj;J++){for(int I=0;I<di;I++){         //Iterate through every grid/raster index in set 'u', computing its nearest neighbor field value 'vf' and distance, then distributing.
        double t,f=1,d=1e+18; int l;                                           //Declare local variables, Initialize associated Voronoi node field value 'f' at 1, distance from point [n] to the Voronoi node 'd', and square of distance 't'.
        for(int m=0;m<db;m+=4){                                                //Iterate through each Voronoi node to find grid point n's distance to and scalar value of nearest neighbor.
            t=pow(b[m]-I,2)+pow(b[m+1]-J,2)+pow(b[m+2]-K,2);                   //Compute grid index's square of distance to each Voronoi node in b as 't'.
            if(t<d){d=t; f=b[m+3];}} l=sqrt(d)+1;                              //If t<d, set Voronoi node distance square d=t, 'f' to cell field value. After loop, set distance s=sqrt(d).
        ///For each raster point n, distribute the field value of its associated (nearest) Voronoi node to all raster points within the distance between point 'n' and its associated Voronoi node 'd'.///
        for(int i=max(0,I-l); i<=min(I+l,di-1); i++){                          //Iterate through axis i indices within l indices of contributing grid point.
            for(int j=max(0,J-l); j<=min(J+l,dj-1); j++){                      //Iterate through axis j indices within l indices of contributing grid point.
                for(int k=max(0,K-l); k<=min(K+l,dk-1); k++){                  //Iterate through axis k indices within l indices of contributing grid point.
                    if(d>=pow(I-i,2)+pow(J-j,2)+pow(K-k,2)){                   //Check if distance to nearest voroni node is less than "d".
                        int n=i+j*di+k*di*dj; u[n]+=f; c[n]++;}}}} }}}         //Increment field value u[ii,jj,kk] by the field value at the originating point u[i,j,k], contributing points counter at c[ii,jj,kk] by 1.
    for(int n=0;n<dn;n++){u[n]/=c[n];} return u;}                              //Divide each element u[n] by grid counter c[n] to normalize. return resultant field u.

///2D implementation of NNI
vector<double> NNI(vector<double> b, int di, int dj){
    int dn=di*dj,db=b.size(); vector<double>u(dn),c(dn); if(db<3){return u;}   //Record array sizes dn(grid points) and db(b vector size). If b is empty, return zero array.
    for(int J=0;J<dj;J++){for(int I=0;I<di;I++){                               //Iterate through every grid/raster index in set 'u', computing its nearest neighbor field value 'vf' and distance, then distributing.
        double t,f=1,d=1e+18; int l,n;                                         //Declare local variables, Initialize associated Voronoi node field value 'f' at 1, distance from point [n] to the Voronoi node 'd', and square of distance 't'.
        for(int m=0;m<db;m+=3){                                                //Iterate through each Voronoi node to find grid point n's distance to and scalar value of nearest neighbor.
            t=pow(b[m]-I,2)+pow(b[m+1]-J,2); if(t<d){d=t; f=b[m+2];}}          //Compute grid index's square of distance to each Voronoi node in b as 't'.
            l=sqrt(d)+1;                                                       //Set grid index offset range l to equal distance to nearest Voronoi node+1
        for(int i=max(0,I-l); i<=min(I+l,di-1); i++){                          //Iterate through axis i indices within l indices of contributing grid point.
            for(int j=max(0,J-l); j<=min(J+l,dj-1); j++){                      //Iterate through axis j indices within l indices of contributing grid point.
                if(d>=pow(I-i,2)+pow(J-j,2)){n=i+j*di; u[n]+=f; c[n]++;}}}}}   //Check if distance to nearest voroni node is less than "d".
    for(int n=0;n<dn;n++){u[n]/=c[n];} return u;}                              //Divide each element u[n] by grid counter c[n] to normalize. return resultant field u.
