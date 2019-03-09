/// No non-STL dependencies. Language standard: C++14. Code may be compiled without multi threaded functionality by commenting out MTM() on lines 268-270. Multi threaded compilation on Ubuntu can be done using -lpthread flag in GCC.
/// This library is released for public use as material supplementary to PhD thesis: Development of Specialized Non-Linear Inversion Algorithms, Basis Functions, Eikonal Solvers, and Their Integration for Use in Joint Seismic and Gravitational Tomographic Inversion
/// Author: Zagid Abatchev - abatchev@ucla.edu ; PI: Paul Davis. Department of Earth, Planetary, and Space Sciences; University of California, Los Angeles; Date of publication: March 2019.///
/// This library is published as is and may have unresolved bugs. Use at own discretion. Users are encouraged to adapt some or all of this code without restriction. Refer to Github link in readme file for updated source code.

/// zmath.cpp is a self contained C++14 Linear algebra and math library using C++ STL. ///

#include<vector>                                                               //Vector class.
#include<math.h>                                                               //Standard math library.
#include<thread>                                                               //Multithreading library. g++ requires flag: -lpthread to compile.
#include<time.h>                                                               //Timer: struct timespec aa,bb; clock_gettime(CLOCK_MONOTONIC,&aa); clock_gettime(CLOCK_MONOTONIC,&bb); print(bb.tv_sec-aa.tv_sec+(bb.tv_nsec-aa.tv_nsec)/1e+9,"timer");
#include<random>                                                               //Random number library.
#include<algorithm>                                                            //Algorithm library.
#include<iostream>                                                             //Console print library.
#include<fstream>                                                              //File i/o stream library.
#include<iomanip>                                                              //Library containing std::setw().
#include<cstdlib>                                                              //Standard C library.
#include<chrono>                                                               //STL time library.
typedef std::vector<double> vec ; typedef std::vector<int> ivec;               //Define template 'vec' for 1D array of doubles using std::vector. Assignment: vec v(array+start_ind, array+end_ind+1)
typedef std::vector<std::vector<double>> mat ;                                 //Define template 'mat' for 2D array of doubles using a vector of nested std::vector<double>.
using std::vector; using std::thread; using std::ref; using std::string;       //Use std::vector namespace.
using std::cout; using std::fstream; using std::ofstream; using std::to_string;//Use selected std namespace functions.
const double pi=3.14159265359, rad=57.29577951;                                //Define global variables pi, rad.

///arr(): Overloaded array generator function for types: vec, mat, indices i,j filled with value s. ///
mat arr(int r, int c, double s=0){ mat m(r,vec(c,s)); return m;}                                                                 //Generate/return i x j type mat filled with value 's'.
mat arr(double *a, int r, int c){mat m; for(int i=0;i<r;i++){m.push_back(vec(a+i*r,a+i*r+c));} return m;}                        //Convert double *a to 2D array mat using indexing: m[i][j]=*((a+i*n)+j)
mat arr(vec v, int r, int c){mat m(r,vec(c)); for(int i=0;i<r;i++){for(int j=0;j<c;j++){m[i][j]=v[i*c+j];}} return m;}           //Convert double *a to 2D array mat using indexing: m[i][j]=*((a+i*n)+j)
vec arr(vector<float> f){vec v; for(auto a:f){v.push_back(a);} return v;}                                                        //Convert vector<float> to new vector<double>
vec arr(vector<int  > i){vec v; for(auto a:i){v.push_back(a);} return v;}                                                        //Convert vector<int> to new vector<double>

///range()/arange(): range int/double vector generation functions.///
ivec range(int h){ivec v; for(int n=0;n<h;n++){v.push_back(n);} return v;}                                                       //Format: vector<int> range(int h): h=upper limit (non inclusive). First value starts at 0, increments by 1.
ivec range(int l, int h, int d){ivec v; for(int n=l;n<h;n+=d){v.push_back(n);} return v;}                                        //Format: vector<int> range(int l, int h): h=upper limit (non inclusive). First value starts at "l", increments by 1.
vec arange(double h){vec v; for(double n=0;n<h;n++){v.push_back(n);} return v;}                                                  //Format: vector<double> arange(double h): h=upper limit (non inclusive). First value starts at 0, increments by 1.
vec arange(double l, double h, double d){vec v; for(double n=l;n<h;n+=d){v.push_back(n);} return v;}                             //Format: vector<double> arange(double l, double h, double d): h=upper limit (non inclusive). First value starts at"l", increments by 'd'.

///Random double/vec functions. Set seed using DRE.seed(seed) in main. Format: vec u=rand(double min, double max, int rows); vec n=nrand(double mean, double variance, int rows); double u=rand(double min, double max); double n=nrand(double mean, double variance) ///
std::default_random_engine DRE; int seed=(std::chrono::system_clock::now().time_since_epoch().count());                          //Generate uniform random distribution DRE with seed set by system clock.
double randint(int h){return std::rand()%h;}                                                                                     //Format: double i=randint(max_int).
double rand (double l, double h){std::uniform_real_distribution<double>R(l,h); return R(DRE);}                                   //Format: double u=uniform(double min=0, double max=1)
double nrand(double a=0, double s=1){std::normal_distribution<double>R(a,s); return R(DRE);}                                     //Format: double n=normal(double mean=0, double variance=1)
vec rand (double l, double h, int r){std::uniform_real_distribution<double>R(l,h); vec v(r); for(auto& d:v){d=R(DRE);} return v;}//Format: vec u=uniform(double min, double max, int rows)
vec nrand(double a, double s, int r){std::normal_distribution<double>R(a,s); vec v(r); for(auto& d:v){d=R(DRE);} return v;}      //Format: vec n=normal(double mean, double variance, int rows)

///print(): Overloaded vec/mat/double/int console print function. String name is optional. Format: print(variable, string name); ///
void print(string s="checkpoint"){cout<<s<<"\n";}                                                                                //Checkpoint/string print function.
void print(int i=0 , string s="checkpoint"){cout<<s<<"= "<<i<<"\n";}                                                             //Checkpoint/int print function.
void print(double d, string s="d"){cout<<s<<"= "<<d<<"\n" ;}                                                                     //Double print function.
void print(vec v, string s="v"){cout<<s<<"=[ "; for(auto d:v){cout<<std::setw(8)<<d<<" ";} cout<<"]\n";}                         //Vec print function.
void print(mat m, string s="m"){cout<<s<<"=\n"; for(vec v:m){for(auto d:v){cout<<std::setw(8)<<d<<" ";} cout<<"\n";}}            //Mat print function.

///savetxt(): Overloaded .txt file save function. Saves vec or mat data using space delimited text files and file names supplied. Format: void savetxt(vec v, string s="v.txt"); void savetxt(mat m, string s="m.txt");///
void savetxt(vec v, string s="v.txt"){ofstream o; o.open(s); for(auto d:v){o<<d<<" ";} o.close();}                               //Save vec as space delimited 1D .txt array Format: void savetxt(vec v,"filename");
void savetxt(mat m, string s="m.txt"){ofstream o; o.open(s); for(vec v:m){for(auto d:v){o<<d<<" ";} o<<"\n";} o.close();}        //Save mat as space delimited 2D .txt array Format: void savetxt(vec v,"filename");

///loadtxt(): Overloaded .txt file load function. Loads vec or mat data from space delimited .txt files. Format: vec v=loadtxt(string fname); mat m=loadtxt(string fname, int columns); ///
vec loadtxt(string F="v.txt"){vec v; double d; fstream f(F); while(!f.eof()){f>>d; v.push_back(d);} v.pop_back(); f.close(); return v;}                                                //Load space delim 1D txt arr.
mat loadtxt(string F, int c){mat m; double d; fstream f(F); while(!f.eof()){vec v={}; for(int j=0;j<c;j++){f>>d; v.push_back(d);} m.push_back(v);} m.pop_back(); f.close(); return m;} //Load space delim 2D txt arr.

///mkdir(): Creates directory in root path "path", with name "name", returns string with full path inside created directory.
string mkdir(string path="", string name=""){
    string cmd, dir;                                                                                                                                //Declare strings cmd, dir.
    if(name==""){                                                                                                                                   //If argument string with name is not passed, or passed as empty.
        time_t rawtime; time(&rawtime); struct tm *ptm=gmtime(&rawtime);                                                                            //Look up raw system time.
        int y=ptm->tm_year+1900, m=ptm->tm_mon, d=ptm->tm_mday, hh=ptm->tm_hour, mm=ptm->tm_min, ss=ptm->tm_sec;                                    //Extract date and time elements as ints.
        name=std::to_string(y)+"_"+std::to_string(m+1)+"_"+std::to_string(d)+"_"+std::to_string(hh)+"_"+std::to_string(mm)+"_"+std::to_string(ss);} //Set directory name to be the system UTC date and time.
    dir=path+name+"/"; cmd="mkdir -p "+dir; int stm=system(cmd.c_str());                                                                            //Create directory in specified path (locally if none passed), and specified or created name.
    std::cout<<"\ncreated directory: "<<dir<<"\n\n"; return dir;}                                                                                   //Print new directory name and path, return created directory path as string.

///sort(): In-place matrix row quicksort by key at column index 'c'. Format: void sort(mat m, int c=0). ///
void sort(mat &m, int c=0, int l=0, int r=-1){
	if(r==-1){r=m.size()-1;}                                                   //If last row index is not specified, use the last index of X.
	int i=l, j=r; double p=m[(l+r)/2][c];                                      //Set sort rows:  i=L=upper, j=R=lower, set pivot value p= m[pivot][c] at midpoint row between rows L and R.
	while(i<=j) {                                                              //While valid swaps remain, iterate partitioning.
	    while(m[i][c]<p){i++;}  while(m[j][c]>p){j--;}                         //Increment i value while m[i][c] is less than the pivot, Decrement j value while m[j][c] is more than the pivot.
	    if(i<=j){m[i].swap(m[j]); i++; j--;}}                                  //If row i is less than row j, swap m[i] and m[j].
	if(l<j){sort(m,c,l,j);} if(i<r){sort(m,c,i,r);}}                           //If i<r: Recursively sort m[i] to m[r]; if l<j , recursively sort m from m[l] to m[j].

///argsort()/sort(): Vector index/value sorting functions, argsort() returns vector<int> with sorted inds of v. sort() returns sorted vec. Format: vec v_sorted=sorted(vec v) ; vector<int> sorted=argsort(vec v).////
vector<int> argsort(vec v){int l=v.size(); mat m(l,vec(2)); for(int i=0;i<l;i++){m[i]={v[i],(double)i};} sort(m); vector<int> a(l); for(int i=0;i<l;i++){a[i]=m[i][1]+.1;} return a;}
vec sort(vec v){std::sort(v.begin(),v.end()); return v;}                                           //Use <algorithm> std::sort(v.begin(),v.end());

///argmin()/argmax()/min()/max(): Return index/value of minimum/maximum value in a vector. Format: int ind=argmin(vec v), int ind=argmax(vec v). float m=min(vec v), float m=max(vec v)////
double min(vec v){auto f=v[0]; for(auto d:v){f=fmin(f,d);} return f;}                              //Return minimum element in vec v.
double max(vec v){auto f=v[0]; for(auto d:v){f=fmax(f,d);} return f;}                              //Return maximum element in vec v.
int argmin(vec v){int a=0; for(unsigned i=1;i<v.size();i++){if(v[i]<v[a]){a=i;}} return a;}        //Return index of minimum element in vec v.
int argmax(vec v){int a=0; for(unsigned i=1;i<v.size();i++){if(v[i]>v[a]){a=i;}} return a;}        //Return index of maximum element in vec v.

///mean() , sigma(): statistical sum, mean and standard deviation calculation functions, returns the mean/standard deviation of a vector. Format: double m=mean(vec v) , double s=sigma(vec A).////
double sum   (vec v){double s=0; for(auto f:v){s+=f;} return s;}                                    //Compute and return the sum of elements in vec v.
double mean  (vec v){double m=0; for(auto f:v){m+=f;} return m/v.size();}                           //Compute and return the mean of vec v.
double median(vec v){std::sort(v.begin(),v.end()); return v[(int)(v.size()/2)];}                    //Compute and return the median of vec v.
double sigma(vec v){double m=0, s=0; for(auto f:v){m+=f;} m/=v.size(); for(auto f:v){s+=pow(f-m,2);} return sqrt(s/(v.size()-1));} //Compute and return the standard deviation of elements in vec v.

///I(): Identity matrix generator, creates mat m=s*I[i,i], where 's' is a scalar (default s=1) Format: mat m=I(int i,double s=1)  ////
mat I(int l, double s=1){mat m(l,vec(l)); for(int i=0;i<l;i++){m[i][i]=s;} return m;}              //Generate m x m zeros matrix, Fill diagonals M[i][i] with 's', return.

///T() matrix transpose function. sets mat m[i,j]=a[j,i] and returns new mat. Format: mat m=T(mat a)////
mat T(mat m){int r=m.size(), c=m[0].size(); mat t(c,vec(r)); for(int i=0;i<r;i++){for(int j=0;j<c;j++){t[j][i]=m[i][j];}} return t;}//Get size m[r,c], create transpose matrix t[c,r], swap elements, return t.

///append(): mat/vec concatenation function. Format: vec o=append(vec a, vec b);
void append(vec &v, vec a){v.insert(v.end(),a.begin(),a.end());}                                   //Append vec a to vec v by reference (first vec is modified)
void append(vec &v, double s){v.push_back(s);}                                                     //Append double 's' to end of vec v.
void append(mat &m, vec v){m.push_back(v);}                                                        //Append vec v to bottom of mat m.
void append(vec v, mat &m){m.insert(m.begin(),v);}                                                 //Append vec v to top of mat m.

///diag(): Overloaded matrix diagonal to vector or vector to diagonal matrix conversion function: Format: vec/mat b=diag(mat/vec a). ///
vec diag(mat m){int l=m.size(); vec v(l);        for(int i=0;i<l;i++){v[i]=m[i][i];} return v;}    //Get l=dim(m), create vec v(l), assign v[i]=m[i][i], return v.
mat diag(vec v){int l=v.size(); mat m(l,vec(l)); for(int i=0;i<l;i++){m[i][i]=v[i];} return m;}    //Get l=dim(v), create output matrix m[l,l], assign m[i][i]=v[i], return m.

///trace(): Return trace of 2d array: sum(M[i][i]) ///
double trace(mat m){int l=m.size(); double d=0; for(int i=0;i<l;i++){d+=m[i][i];} return d;}       //Get sum of diagonal elements of matrix m, return.

///col(): overloaded mat column extraction/setting function. Return all elements in mat m[:][c] as a vec, or return mat m= m[:][c]=vec. Format: vec v=col(mat m, int c)  mat m=col(mat m, vec v); ///
vec  col(mat  m, int c){int r=m.size(); vec v(r); for(int i=0;i<r;i++){v[i]=m[i][c];} return v;}   //Find rows in m, define an output vector v, set v[i]=m[i][j], return v.
void col(mat &m, int c, vec v){for(unsigned i=0; i<m.size(); i++){m[i][c]=v[i];}}                  //Set all elements in m[i][c]=v[i].

///det(): Recursive square mat determinant function. Format: double d=det(mat a)///
double det(mat m){
	int l=m.size(); double d=0; mat mi=m, mb;                                  //Get size m[l,l], initialize determinant d, copy m to intermediate matrix mi.
	if(l>2){                                                                   //If dim(m)>2, decompose into cofactors and minors, call determinants of minors recursively.
		for(int i=0;i<l;i++){mi[i].erase(mi[i].begin(), mi[i].begin()+1);}     //Remove 0th column from matrix mi to give the columns of the minor mi.
		for(int i=0;i<l;i++){mb=mi; mb.erase(mb.begin()+i,mb.begin()+(i+1));   //Remove 'i'th row to give minor mi of m with row i and column 0 removed.
			d+=(1-2*(i%2))*m[i][0]*det(mb);}}                                  //Sum cofactors m[i][0]*(-1)^i*det(m) to compute determinant value. Note: (-1)^i=1-2*(i%2)
	else{d=m[0][0]*m[1][1]-m[1][0]*m[0][1];} return d;}                        //If dim(m)=2 return d=m[0][0]*m[1][1]-m[0][1]*m[1][0]; Return d.

///inv(): Matrix inversion function using Gauss-Jordan method. Stable against pivot instabilities. Format: mat mi=inv(mat m). ///
mat inv(mat a){
	int l=a.size(); double c; mat m(l,vec(l)); for(int i=0;i<l;i++){m[i][i]=1;}//Get size l=dim(a), initialize inverse matrix m[l][l] as identity. All operations on 'a' are duplicated on 'm'.
	for(int i=0; i<l; i++){                                                    ///Iterate through every row to transform 'a' to upper-right triangle form normalized along diagonal.
		int im=i; for(int n=i;n<l;n++){if(fabs(a[n][i])>fabs(a[im][i])){im=n;}}//Find the largest element in column i below row i: im=argmax(abs(a[ii>i][i])).
		a[i].swap(a[im]); m[i].swap(m[im]); if(a[i][i]==0){return arr(l,l);}   //Swap a[i] <-> a[im] to avoid pivot instabilities. If a[im][i]=0 (singular matrix), return zero.
		double f=a[i][i]; for(int j=0;j<l;j++){a[i][j]/=f; m[i][j]/=f;}        //Normalize diagonal element: divide a[i] by f=a[i][i] to make a[i][i]=1.
		for(int n=i+1;n<l;n++){ c=a[n][i];                                     //Set c=a[n][i], iterate through each row n>i.
            for(int j=0;j<l;j++){a[n][j]-=c*a[i][j]; m[n][j]-=c*m[i][j];}}}    //Subtract c*a[i] from a[n] so a[n][i]=0, bringing a to upper-right triangular form.
	for(int i=1;i<l;i++){                                                      ///Row reduce 'a' to identity matrix, duplicating operations on m.
		for(int n=0;n<i;n++){c=a[n][i];                                        //Subtract multiple c=a[n][i] of each row a[i] iteratively from each row above a[n<i] to bring 'a' to identity form.
            for(int j=0;j<l;j++){m[n][j]-=c*m[i][j]; a[n][j]-=c*a[i][j];}}}    //Set c, subtract c*a[i] from a[n] so a[n][i]=0: a[n]=a[n]-c*a[i].
	return m;}                                                                 //Return now inverted matrix m.

///Overloaded element wise basic math functions for vec/mat.
vec fabs(vec v){for(auto &d:v){d=fabs(d);} return v;}                          //Generalize fabs() to elements of a vec.
mat fabs(mat m){for(auto &v:m){for(auto &d:v){d=fabs(d);}} return m;}          //Generalize fabs() to elements of a mat.
vec sqrt(vec v){for(auto &d:v){d=sqrt(d);} return v;}                          //Generalize sqrt() to elements of a vec.
mat sqrt(mat m){for(auto &v:m){for(auto &d:v){d=sqrt(d);}} return m;}          //Generalize sqrt() to elements of a mat.
vec pow(vec v, int n){for(auto&d:v){d=pow(d,n);} return v;}                    //Generalize pow() to elements of a vec.
mat pow(mat m, int n){for(auto&v:m){for(auto&d:v){d=pow(d,n);}} return m;}     //Generalize pow() to elements of a mat.
vec exp(vec v){for(auto&d:v){d=exp(d);} return v;}                             //Generalize exp() to elements of a vec.
mat exp(mat m){for(auto&v:m){for(auto&d:v){d=exp(d);}} return m;}              //Generalize exp() to elements of a mat.
vec rint(vec v){for(auto &d:v){d=fabs(d);} return v;}                          //Generalize rint() to elements of a vec.
mat rint(mat m){for(auto &v:m){for(auto &d:v){d=rint(d);}} return m;}          //Generalize rint() to elements of a mat.
vec mult(vec a, vec b){for(unsigned n=0;n<a.size();n++){a[n]*=b[n];} return a;}//Element wise multiplication for equal sized vec a, vec b.
mat mult(mat a, mat b){for(unsigned n=0;n<a.size();n++){for(unsigned m=0; n<a[0].size(); n++){a[n][m]*=b[n][m];}} return a;}  //Element wise multiplication for equal sized mat a, mat b.

///add(): Overloaded element wise addition function: c=a+b where a,b are double, vec or mat; c is of type and dimension of higher parameter rank. Format: double/vec/mat c=add(double/vec/mat a, double/vec/mat b)
vec add(vec a, vec b){for(unsigned i=0;i<a.size();i++){a[i]+=b[i];} return a;} //Element wise addition for equal sized vec a, vec b.
vec add(vec a, double s){for(auto&d:a){d+=s;} return a;}                       //Element wise addition of double s to vec a.
vec add(double s, vec a){for(auto&d:a){d+=s;} return a;}                       //Element wise addition of double s to vec a (reverse args order).
mat add(mat a, mat b){for(unsigned i=0;i<a.size();i++){for(unsigned j=0;j<a[0].size();j++){a[i][j]+=b[i][j];}} return a;} //Element wise addition for equal sized mat a, mat b.
mat add(double s, mat m){for(auto&v:m){for(auto&d:v){d+=s;}} return m;}        //Element wise addition of double s to mat a.
mat add(mat m, double s){for(auto&v:m){for(auto&d:v){d+=s;}} return m;}        //Element wise addition of double s to mat a (reverse args order).

///div(): Overloaded element wise division of: double/vec, double/mat, mat/mat, mat/double, vec/vec, vec/double.
vec div(vec a, vec b){for(unsigned i=0;i<a.size();i++){a[i]/=b[i];} return a;} //Element wise division of vec a by vec b.
vec div(vec v, double s){for(auto &d:v){d/=s;} return v;}                      //Element wise division of vec a by double s.
vec div(double s, vec v){for(auto &d:v){d=s/d;} return v;}                     //Element wise division of double s by vec a.
mat div(mat a, mat b){for(unsigned i=0;i<a.size();i++){for(unsigned j=0;j<a[0].size();j++){a[i][j]/=b[i][j];}} return a;}                  //Element wise division of mat a by mat b.
mat div(mat m, double s){for(auto&v:m){for(auto&d:v){d/=s;}} return m;}        //Element wise division of mat a by double s.
mat div(double s, mat m){for(auto&v:m){for(auto&d:v){d=s/d;}} return m;}       //Element wise division of double s by mat a.

///dot(): Overloaded inner product of two mat/vec types, and/or their scalar product with a double. Format: mat/vec a=dot(double/vec/mat A, double/vec/mat B, double s=1) ////
double dot(vec a, vec b){double s=0; for (unsigned i=0;i<a.size();i++){s+=a[i]*b[i];} return s;} //Inner/dot product of vec a, vec b.
double dot(vec v){double s=0; for(auto d:v){s+=d*d;} return s;}                //Inner/dot product of vec v with itself.
vec dot(double s, vec v){for(auto&d:v){d*=s;} return v;}                       //Product of double s with elements of vec a.
vec dot(vec v, double s){for(auto&d:v){d*=s;} return v;}                       //Product of double s with elements of vec a (reversed args order).
mat dot(mat m, double s){for(auto&v:m){for(auto&d:v){d*=s;}} return m;}        //Product of double s with elements of mat a.
mat dot(double s, mat m){for(auto&v:m){for(auto&d:v){d*=s;}} return m;}        //Product of double s with elements of mat a.
vec dot(mat m, vec v){int r=m.size(), c=m[0].size(); vec o(r); for(int i=0;i<r;i++){for(int j=0;j<c;j++){o[i]+=m[i][j]*v[j];}} return o;}  //Inner/dot product of mat m with vec a.
vec dot(vec v, mat m){int r=m.size(), c=m[0].size(); vec o(c); for(int i=0;i<r;i++){for(int j=0;j<c;j++){o[j]+=v[i]*m[i][j];}} return o;}  //Inner/dot product of vec a with mat m.
mat dot(mat a, mat b){int o=a.size(), l=a[0].size(), n=b[0].size(); mat m=arr(o,n); for(int i=0;i<o;i++){for(int j=0;j<n;j++){for(int k=0;k<l;k++){m[i][j]+=a[i][k]*b[k][j];}}} return m;} //Inner product of mat a, mat b.

///Operator overloads for arithmetic operations with double/vec/mat objects. Notice: * operator between vec*vec, mat*vec, vec*mat, mat*mat is defined as an inner/dot product. Use mult() for element wise products.///
vec operator+(vec a, vec b){return add(a,b);}
vec operator+(vec a, double b){return add(a,b);}
vec operator+(double a, vec b){return add(a,b);}
mat operator+(mat a, mat b){return add(a,b);}
mat operator+(mat a, double b){return add(a,b);}
mat operator+(double a, mat b){return add(a,b);}

vec operator-(vec a, vec b){return add(a,dot(-1.,b));}
vec operator-(vec a, double b){return add(a,-b);}
vec operator-(double a, vec b){return add(a,dot(-1.,b));}
mat operator-(mat a, mat b){return add(a,dot(-1.,b));}
mat operator-(mat a, double b){return add(a,-b);}
mat operator-(double a, mat b){return add(a,dot(-1.,b));}

double operator*(vec a, vec b){return dot(a,b);}
vec operator*(vec a, double b){return dot(a,b);}
vec operator*(double a, vec b){return dot(a,b);}
vec operator*(vec a, mat b){return dot(a,b);}
vec operator*(mat a, vec b){return dot(a,b);}
mat operator*(mat a, mat b){return dot(a,b);}
mat operator*(mat a, double b){return dot(a,b);}
mat operator*(double a, mat b){return dot(a,b);}

vec operator/(vec a, vec b){return div(a,b);}
vec operator/(vec a, double b){return div(a,b);}
vec operator/(double a, vec b){return div(a,b);}
mat operator/(mat a, mat b){return div(a,b);}
mat operator/(mat a, double b){return div(a,b);}
mat operator/(double a, mat b){return div(a,b);}

vec operator^(vec a, int n){return pow(a,n);}
mat operator^(mat a, int n){return pow(a,n);}

///interp(): Overloaded Linear/Bilinear/Trilinear interpolation functions. Variable format: (x,y,z)= interpolation coordinates   (x0,y0,z0)=floating coordinate offsets   (di,dj,dk)=grid axis dimensions   (hi,hj,hk)= axial step sizes in units of x,y,z
//1D interpolation format: double interp(vec &f, double x, int di, double h, double x0=0)
double interp(vec &f, double x, int di, double h, double x0=0){
    double si,ii,c0,c1,c; int i;                                               //Declare index and coordinate variables.
    ii=fmin(fmax((x-x0)/h,0),(di-1)); i=ii; si=ii-i;                           //Set x ind coord (ii) constrained to range [0:di], find (int)i immediately below ii, and distance between them 'si'.
    c0=f[i]; c1=f[i+1]; c=c0*(1-si)+c1*si; return c;}                          //Extract grid points c0,c1; 1D interpolate between c0, c1, as c, return c.

//2D interpolation format: double interp(vec &f, double x, double y, int di, int dj, double hi, double hj=0, double x0=0, double y0=0)
double interp(vec &f, double x, double y, int di, int dj, double hi, double hj=0, double x0=0, double y0=0){
    double si,sj,ii,jj,c00,c01,c10,c11,c0,c1,c; int i,j;                       //Declare index and coordinate variables.
    if(fabs(hj)<1e-12){hj=hi;}                                                 //If only one index spacing is given, set every axis spacing to the same spacing.
    ii=fmin(fmax((x-x0)/hi,0),(di-1)); i=ii;  si=ii-i;                         //Set x ind coord (ii) constrained to range [0:di], find (int)i immediately below ii, and distance between them 'si'.
    jj=fmin(fmax((y-y0)/hj,0),(dj-1)); j=jj;  sj=jj-j;                         //Set y ind coord (jj) constrained to range [0:dj], find (int)j immediately below jj, and distance between them 'sj'.
    c00=f[i+(j+0)*di]   ; c10=f[i+1+(j+0)*di];                                 //Extract grid points c00, c10.
    c01=f[i+(j+1)*di]   ; c11=f[i+1+(j+1)*di];                                 //Extract grid points c01, c11
    c0=c00*(1-si)+c10*si; c1=c01*(1-si)+c11*si;                                //1D interpolation between grid neighboring lines.
    c=c0*(1-sj)+c1*sj; return c;}                                              //Compute 2D interpolated value between 2 lines and return.

//3D interpolation format:   double interp(vec &f, double x, double y, double z, int di, int dj, int dk, double hi, double hj=0, double hk=0, double x0=0, double y0=0, double z0=0)
double interp(vec &f, double x, double y, double z, int di, int dj, int dk, double hi, double hj=0, double hk=0, double x0=0, double y0=0, double z0=0){
    double si,sj,sk,ii,jj,kk,c000,c100,c010,c110,c001,c101,c011,c111,c00,c01,c10,c11,c0,c1,c; int i,j,k;//Declare index and coordinate variables.
    if(fabs(hj+hk)<1e-12){hj=hi; hk=hi;}                                       //If only one index spacing is given, set every axis spacing to the same spacing.
    ii=fmin(fmax((x-x0)/hi,0),(di-2)); i=ii;  si=ii-i;                         //Set x ind coord (ii) constrained to range [0:di], find (int)i immediately below ii, and distance between them 'si'.
    jj=fmin(fmax((y-y0)/hj,0),(dj-2)); j=jj;  sj=jj-j;                         //Set y ind coord (jj) constrained to range [0:dj], find (int)j immediately below jj, and distance between them 'sj'.
    kk=fmin(fmax((z-z0)/hk,0),(dk-2)); k=kk;  sk=kk-k;                         //Set z ind coord (kk) constrained to range [0:dk], find (int)k immediately below kk, and distance between them 'sk'.
    c000=f[i+(j+0)*di+(k+0)*di*dj]; c100=f[i+1+(j+0)*di+(k+0)*di*dj];          //Extract grid points c000, c100.
    c010=f[i+(j+1)*di+(k+0)*di*dj]; c110=f[i+1+(j+1)*di+(k+0)*di*dj];          //Extract grid points c010, c110.
    c001=f[i+(j+0)*di+(k+1)*di*dj]; c101=f[i+1+(j+0)*di+(k+1)*di*dj];          //Extract grid points c001, c101.
    c011=f[i+(j+1)*di+(k+1)*di*dj]; c111=f[i+1+(j+1)*di+(k+1)*di*dj];          //Extract grid points c011, c111.
    c00=c000*(1-si)+c100*si       ; c01=c001*(1-si)+c101*si;                   //1D interpolation between grid neighboring lines.
    c10=c010*(1-si)+c110*si       ; c11=c011*(1-si)+c111*si;                   //1D interpolation between grid neighboring lines.
    c0=c00*(1-sj)+c10*sj          ; c1=c01*(1-sj)+c11*sj;                      //2D interpolation between grid neighboring planes.
    c=c0*(1-sk)+c1*sk             ; return c;}                                 //Compute 3D interpolated value between 2 planes and return.

//Overloaded 3D interpolation format:  double interp(vec &f, vec c, ivec d, vec h, vec o={0,0,0}) vec c={x,y,z} ; ivec d=vector<int> with axial dimensions {di,dj,dk} ; vec h=vector<double> with axial grid spacing {hi,hj,hk}. If "h" is passed with 1 element, all axial step sizes are assumed isotropic.
double interp(vec &f, vec c, ivec d, vec h, vec o={0,0,0}){if(h.size()==1){h={h[0],h[0],h[0]};} return interp(f,c[0],c[1],c[2],d[0],d[1],d[2],h[0],h[1],h[2],o[0],o[1],o[2]);} //If grid spacing vector contains only 1 value, assign it to all 3 axes (isotropic grid). Return 3D interpolation using vectorized inputs.

///bounded(): Overloaded bounds checking function for types double and vec. upper and lower bounds can be defined either as double or vec. Returns true if in bounds; false if out of bounds.///
bool bounded(double x, double l, double h){if(x>=l && x<=h){return true;} return false;}                                      //Bounds checking function for doubles. format: bool bounded(double test_value,double low, double high);
bool bounded(vec v, double l, double h){for(double x:v){if(x<l || x>h){return false;}} return true;}                          //Bounds checking function for vector of doubles against single bound. If any value in test vector is outside bounds, returns false; else returns true. Format: bool bounded(vector<double> test_vector, double low, double high)
bool bounded(vec v, vec l, vec h){for(unsigned i=0; i<v.size();i++){if(v[i]<l[i] || v[i]>h[i]){return false;}} return true;}  //Bounds checking function for vector of doubles against vector bounds. If any value in test vec is outside bounds at corresponding indices of low and high, returns false; else returns true. Format: bool bounded(vector<double> test_vector, vector<double> low, vector<double> high)

///Fast multithreaded matrix m transpose product with self MTM(a)=dot(transpose(m),m); Format: mat MTM(mat a). a=rectangular matrix to be transposed and dotted with self.///
mat MTM(mat a){ int c=a[0].size(), r=a.size(); mat m(c,vec(c)); thread tr[c]; a=T(a);                                  // Query input matrix dimensions as c,r; Declare output matrix m; declare thread array; transpose input matrix as self.
    auto F=[&](int i){for(int j=i;j<c;j++){double s=0; for(int n=0;n<r;n++){s+=a[i][n]*a[j][n];} m[i][j]=m[j][i]=s;}}; // Define thread lambda function. Each thread evaluates corresponding index row of inner product of initial matrix transpose with initial matrix.
    for(int i=0;i<c;i++){tr[i]=thread(F,i);} for(auto&t:tr){t.join();} return m;}                                      // Execute 1 thread per column. Join threads, return transformed matrix.

///TIMER class for measuring time between calls .start() and .stop().  TIMER.start() starts timer.  TIMER.stop() returns elapsed time in seconds as type double. TIMER.pstop() prints elapsed time to console.
class TIMER{public: struct timespec a,b;                                                                               // Declare timespec structs for start and stop times as 'a','b'.
    void   start(void){clock_gettime(CLOCK_MONOTONIC,&a);}                                                             // Define void class method "start()", which registers current time.
    double  stop(void){clock_gettime(CLOCK_MONOTONIC,&b); return b.tv_sec-a.tv_sec+(b.tv_nsec-a.tv_nsec)/1e+9;}        // Define void class method "stop()", which registers current time, measures offset from start in seconds, returns calculated offset as double.
    void   pstop(std::string name="timer"){clock_gettime(CLOCK_MONOTONIC,&b); print(1.*(b.tv_sec-a.tv_sec+(b.tv_nsec-a.tv_nsec)/1e+9),name);} }; // Define void class method "pstop()", which registers current time, measures offset from start in seconds, prints computed offset, returns calculated offset as double.


///fit() IS UNTESTED. USE AT OWN RISK.///
/////fit(): Levenberg-Marquardt Algorithm for iterating non linear least squares fitting with self adjusting damping.///
//vec fit(auto F, auto a, vec B, vec y, int i_max=10, double trig=.99, vec s={1e-6}, double dmp=.001){
//        int mi=0, n=B.size(); vec f=F(a,B), t,b; double ssq=dot(y-f),q; mat J=arr(y.size(),n); //Initialize model output f, set start vector B=b, initial ssq, misstep mi=0, initialize jacobian.
//        if(s.size()<n){s=s[0]*fabs(b)+1e-9;}                                   //Initialize step vector s.
//        for(int i=0; i<i_max; i++){                                            //Loop i_max times.
//            for(int c=0;c<n;c++){vec r=b,l=b; r[c]+=s[c]/2; l[c]-=s[c]/2;      //Set right and left hand size step vectors for column c as 'r' and 'l'. Step elements at c by +/-s[c].
//                col(J,c,(F(a,r)-F(a,l))/s[c]);}                                //Compute Jacobian column J[:,c]=d(F[:])/d(B[c]) Set J[:,c]=(F(a,r)-F(a,l))/s[c].
//            while(mi<=i_max){                                                  //Loop while i<i_max.
//                b=B+inv(T(J)*J+dmp*diag(diag(T(J)*J)))*T(J)*(y-f);             //Solve for test b= B+inv(J'J+D)*J'*(y-f(a,b)). J=Jacobian, J'=J.transpose, D=damping matrix.
//                t=F(a,b); q=dot(y-t); if(q/ssq<trig){f=t; ssq=q; B=b; break;}  //Evaluate test values: t=F(a,b), and q=(y-f)^2.if q/ssq<trig, accept q as new ssq, b as new B, break out of loop.
//                dmp*=10; mi++;}                                                //If step was unsuccessful, increase damping by factor of 10, increment misstep index mi.
//            if(mi>i_max){break;} else if(mi>0){dmp/=10; mi--;}} return B;}     //If number of missteps==iter_max, terminate inversion; else if damping is elevated (mi>0), set: d=d/10; mi--; increment i.


///TIMER: struct timespec aa,bb; clock_gettime(CLOCK_MONOTONIC,&aa); ###CODE HERE###;  clock_gettime(CLOCK_MONOTONIC,&bb); print(bb.tv_sec-aa.tv_sec+(bb.tv_nsec-aa.tv_nsec)/1e+9,"timer");
//int main(){ ///zmath test script. To compile and execute, run: " g++ -o zmath zmath.cpp && ./zmath "
//////vec v1={1,2,3}, v2={7,8,9}, v3={4,2}, v4={3,5}; mat m1={{3,2,1},{2,3,2},{1,2,3}}, m2={{4,1,3},{2,8,6},{4,6,9},{9,8,5}}, m3={{3,2},{1,2}}, m4={{2,5},{3,7}}, m5={{2,1,2},{3,1,4}}; //Define vec,mat test objects.
//mat m=loadtxt("J_1027.txt",1027); mat a=MTM(m)+diag(vec(1027,1));
//struct timespec aa,bb; clock_gettime(CLOCK_MONOTONIC,&aa);
//mat i; for(double kk=0; kk<10.1; kk++){i=i+inv(a);}
//clock_gettime(CLOCK_MONOTONIC,&bb); print(bb.tv_sec-aa.tv_sec+(bb.tv_nsec-aa.tv_nsec)/1e+9,"timer 1");
//double S=0; for(vec v:i){S+=dot(v);} print(S);
//
//}
//print(MTM(m3),"MTM"); print(T(m3)*m3,"T(m)*m");}
//print(m3,"m"); print(det(m3),"det|m|");                                        //Determinant test.
//print(m1,"m"); print(inv(m1),"inv(m)");                                        //Matrix inverse test.
//print(m3,"m"); print(T(m3),"m^T");                                             //Matrix transpose test.
//print(add(m1,m2),"mat mat sum"   );                                            //Matrix sum test.
//print(add(v1,v2),"vec vec sum"   );                                            //Vector sum test.
//print(dot(v1,v2),"vec vec dot");                                               //Dot product test.
//print(dot(2.,v2),"scalar vec dot") ;                                           //Scalar vector product test.
//print(dot(2.,m1),"scalar mat dot") ;                                           //Scalar Matrix product test.
//print(dot(v2,m2),"mat vec dot") ;                                              //Matrix vector product test.
//print(1.2*dot(m1,m2),"mat mat dot") ;                                          //Matrix product test.
//print(v3); savetxt(v3); vec v4=loadtxt(); print(v4);                           //Vector save/load test.
//print(m1); savetxt(m1); mat m6=loadtxt("m.txt",3); print(m6);                  //Matrix save/load test.
//struct timespec aa,bb; clock_gettime(CLOCK_MONOTONIC,&aa); print("start timer"); clock_gettime(CLOCK_MONOTONIC,&bb); print(bb.tv_sec-aa.tv_sec+(bb.tv_nsec-aa.tv_nsec)/1e+9,"timer"); //Timer commands.
//ers=sqrt(diag(inv(T(J)*J)*I(m,SSQ.back()/(m-n))));                             //compute fit uncertainties ers=sqrt(diag(c)), where c=inv(JT*J)*(ssq/(n-m)*I); n=data m=params. Format: vec e=Ers()
