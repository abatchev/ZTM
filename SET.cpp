/// No non-STL dependencies. Language standard: C++14.
/// This library is released for public use as material supplementary to PhD thesis: Development of Specialized Non-Linear Inversion Algorithms, Basis Functions, Eikonal Solvers, and Their Integration for Use in Joint Seismic and Gravitational Tomographic Inversion
/// Author: Zagid Abatchev - abatchev@ucla.edu ; PI: Paul Davis. Department of Earth, Planetary, and Space Sciences; University of California, Los Angeles; Date of publication: March 2019.///
/// This library is published as is and may have unresolved bugs. Use at own discretion. Users are encouraged to adapt some or all of this code without restriction. Refer to Github page for updated source code, located at: https://github.com/abatchev/ZTM .

/// SET.cpp is a self contained C++14 Standard Earth library.///
/// NOTICE: If SET() is instantiated but lookup table "set.dat" has not been pre-generated, generation is done at runtime, and may take several hours. This is done only once. ///

#include <stdlib.h>                                                            //Include stdlib for malloc() headers.
#include <math.h>                                                              //Include STL math library.
#include <vector>                                                              //Include STL vector library.
#include <stdio.h>                                                             //Include STL file i/o library headers.

///iasp91(): Kennet & Engdahl 1991 standard earth velocity interpolation function, returns iasp91 p velocity as a function of depth, with option of linear blurring at discontinuities. ///
///Format: double v=iasp91(double d, double w=1); v=P wave velocity(km/s); d=depth(km); w=linear discontinuity blurring (km).                                                           ///
double iasp91(double d, double w=1){
    //Kennet&Engdahl 1991 Table 1: Standard earth depths and polynomial coefficients of p velocity. Form: v=c0+c1*x+c2*x^2+c3*x^3; x=(6371-d)/6371, c= coefficients of row corresponding to d>c[0] && d<c[1]
    const double c[11][6]={                                                    //Row format: c[n]=[d_min, d_max, c0, c1, c2, c3].
    { -100.0,   20.0,  5.80000,   0.00000,   0.00000,   0.00000},              //Upper crust range and coefficients.
    {   20.0,   35.0,  6.50000,   0.00000,   0.00000,   0.00000},              //Lower crust range and coefficients.
    {   35.0,  120.0,  8.78541,  -0.74953,   0.00000,   0.00000},              //Lithosphere range and coefficients.
    {  120.0,  210.0, 25.41389, -17.69722,   0.00000,   0.00000},              //Aesthenosphere range and coefficients.
    {  210.0,  410.0, 30.78765, -23.25415,   0.00000,   0.00000},              //Mid mesosphere range and coefficients.
    {  410.0,  660.0, 29.38896, -21.40656,   0.00000,   0.00000},              //Mid mesosphere range and coefficients.
    {  660.0,  760.0, 25.96984, -16.93412,   0.00000,   0.00000},              //Mid mesosphere range and coefficients.
    {  760.0, 2740.0, 25.14860, -41.15380,  51.99320, -26.60830},              //Lower mesosphere range and coefficients.
    { 2740.0, 2889.0, 14.49470,  -1.47089,   0.00000,   0.00000},              //Lower mesosphere range and coefficients.
    { 2889.0, 5153.9, 10.03904,   3.75665, -13.67046,   0.00000},              //Outer core range and coefficients.
    { 5153.9, 6371.0, 11.24094,   0.00000,  -4.09689,   0.00000}};             //Inner core range and coefficients.
    if(2*fabs(d-  20)<=w ){return  5.80+0.70*(d-  20+.5*w)/w;}                 //Check if depth is within w/2 km of 20km crustal discontinuity. If so, extend the boundary over range w and linearly interpolate.
    if(2*fabs(d-  35)<=w ){return  6.50+1.54*(d-  35+.5*w)/w;}                 //Check if depth is within w/2 km of continental Moho(35km). If so, extend the boundary over range w and linearly interpolate.
    if(2*fabs(d- 410)<=w ){return  9.03+0.33*(d- 410+.5*w)/w;}                 //Check if depth is within w/2 km of 410km discontinuity. If so, extend the boundary over range w and linearly interpolate.
    if(2*fabs(d- 660)<=w ){return 10.20+0.59*(d- 660+.5*w)/w;}                 //Check if depth is within w/2 km of upper/lower mantle discontinuity. If so, extend the boundary over range w and linearly interpolate.
    if(2*fabs(d-2889)<=w ){return 13.69-5.68*(d-2889+.5*w)/w;}                 //Check if depth is within w/2 km of core mantle boundary. If so, extend the boundary overrange w and linearly interpolate.
    if(2*fabs(d-5154)<=w ){return 10.26-0.17*(d-5154+.5*w)/w;}                 //Check if depth is within w/2 km of inner/outer core boundary. If so, extend the boundary overrange w and linearly interpolate.
    double x=(6371-d)/6371;                                                    //Compute normalized radius x=(6371-depth)/6371.
    for(int n=0; n<11; n++){                                                   //Advance through coefficient array c, find row [n] that corresponds to queried depth. Depth range of each row specified in c[][0] : c[][1]
        if(d>=c[n][0] && d<=c[n][1]){                                          //Check if d>c[n][0] and d<c[n][1] ; depth between values of first and second column.
            return c[n][2]+c[n][3]*x+c[n][4]*x*x+c[n][5]*x*x*x;}}              //Return v=c0+c1*x+c2*x^2+c3*x^3 where c0=c[n][2], c1=c[n][3], c2=c[n][4], c3=c[n][5].
    return 1;}                                                                 //If an undefined depth was queried, return v=1 to signify error.

 ///interp3d(): Trilinear interpolation function described at:'https://en.wikipedia.org/wiki/Trilinear_interpolation' ///
 ///Format:   double interp3d(std::vector<float> &f, double x, double y, double z, int di, int dj, int dk, double hi, double hj=0, double hk=0, double x0=0, double y0=0, double z0=0)
 ///Variables: (x,y,z)= interpolation coordinates   (x0,y0,z0)=floating coordinate offsets   (di,dj,dk)=grid axis dimensions   (hi,hj,hk)= axial step sizes in units of x,y,z
double interp3d(std::vector<float> &f, double x, double y, double z, int di, int dj, int dk, double hi, double hj=0, double hk=0, double x0=0, double y0=0, double z0=0){
    double si,sj,sk,ii,jj,kk,c000,c100,c010,c110,c001,c101,c011,c111,c00,c01,c10,c11,c0,c1,c; int i,j,k;//Declare index and coordinate variables.
    if(fabs(hj+hk)<.0000000001){hj=hi; hk=hi;}                                 //If only one index spacing is given, set every axis spacing to the same spacing.
    ii=fmin(fmax((x-x0)/hi,0),(di-1)); i=ii;  si=ii-i;                         //Set x index coordinate (ii) constrained to range [0:di], find (int)i immediately below ii, and distance between them 'si'.
    jj=fmin(fmax((y-y0)/hj,0),(dj-1)); j=jj;  sj=jj-j;                         //Set y index coordinate (jj) constrained to range [0:dj], find (int)j immediately below jj, and distance between them 'sj'.
    kk=fmin(fmax((z-z0)/hk,0),(dk-1)); k=kk;  sk=kk-k;                         //Set z index coordinate (kk) constrained to range [0:dk], find (int)k immediately below kk, and distance between them 'sk'.
    c000=f[i+(j+0)*di+(k+0)*di*dj]; c100=f[i+1+(j+0)*di+(k+0)*di*dj];          //Extract grid points c000, c100.
    c010=f[i+(j+1)*di+(k+0)*di*dj]; c110=f[i+1+(j+1)*di+(k+0)*di*dj];          //Extract grid points c010, c110.
    c001=f[i+(j+0)*di+(k+1)*di*dj]; c101=f[i+1+(j+0)*di+(k+1)*di*dj];          //Extract grid points c001, c101.
    c011=f[i+(j+1)*di+(k+1)*di*dj]; c111=f[i+1+(j+1)*di+(k+1)*di*dj];          //Extract grid points c011, c111.
    c00=c000*(1-si)+c100*si       ; c01=c001*(1-si)+c101*si;                   //1D interpolation between grid neighboring lines.
    c10=c010*(1-si)+c110*si       ; c11=c011*(1-si)+c111*si;                   //1D interpolation between grid neighboring lines.
    c0=c00*(1-sj)+c10*sj          ; c1=c01*(1-sj)+c11*sj;                      //2D interpolation between grid neighboring planes.
    c=c0*(1-sk)+c1*sk             ; return c;}                                 //Compute 3D interpolated value between 2 planes and return.


 ///gca(): Haversine Function for great circle arc calculation. Format: double degrees_arc=gca(double lat1, double lon1, double lat2, double lon2); ///
double gca(double lat1, double lon1, double lat2, double lon2){double rad=57.2958; return 2*asin(sqrt(pow(sin((lat2-lat1)/(2*rad)),2)+cos(lat1/rad)*cos(lat2/rad)*pow(sin((lon2-lon1)/(2*rad)),2)))*rad;}

///geo2enu(): Lat, Lon, Dep(km) to ENU x,y,z(km) coordinate transformation function. Uses lat_ref,lon_ref,d_ref as ENU origin reference. ///
///Format: std::vector<double> xyz= geo2enu(double lat_ref, double lon_ref, double d_ref, double lat, double lon, double d)              ///
std::vector<double> geo2enu(double lat_ref, double lon_ref, double d_ref, double lat, double lon, double d){
    double xi,yi,zi, xr,yr,zr, phi,theta, x,y,z, rad=57.295, R=6371;           //Initialize variables, refine radian to degree conversion, and earth radius R.
    xi=(R-d)*cos(lon/rad)*cos(lat/rad);                                        //Find ECEF Cartesian sample x coordinate using transformation: x=r*cos(lon)*cos(lat)
    yi=(R-d)*sin(lon/rad)*cos(lat/rad);                                        //Find ECEF Cartesian sample y coordinate using transformation: y=r*sin(lon)*cos(lat)
    zi=(R-d)*sin(lat/rad);                                                     //Find ECEF Cartesian sample z coordinate using transformation: z=r*sin(lat)
    xr=(R-d_ref)*cos(lon_ref/rad)*cos(lat_ref/rad);                            //Find ECEF Cartesian reference x coordinate using transformation: x=r*cos(lon)*cos(lat)
    yr=(R-d_ref)*sin(lon_ref/rad)*cos(lat_ref/rad);                            //Find ECEF Cartesian reference y coordinate using transformation: y=r*sin(lon)*cos(lat)
    zr=(R-d_ref)*sin(lat_ref/rad);                                             //Find ECEF Cartesian reference z coordinate using transformation: z=r*sin(lat)
    phi=lat_ref/rad ; theta=lon_ref/rad ;                                      //Convert reference latitude and longitude from degrees to radians, as phi,theta respectively.
    double m[9]={-sin(theta)         ,      cos(theta)     ,     0   ,         //Transformation matrix from ECEF(Earth Centered Earth Fixed) to ENU coordinates. phi=latitude, theta=longitude.
                 -sin(phi)*cos(theta), -sin(phi)*sin(theta), cos(phi),
                 cos(phi)*cos(theta) ,  cos(phi)*sin(theta), sin(phi)};
    double v[3]={xi-xr, yi-yr, zi-zr};                                         //Define difference vector 'v', specifying the ECEF difference vector between reference and sample points.
    x=m[0]*v[0]+m[1]*v[1]+m[2]*v[2];                                           //Take inner product of rotation matrix with sample-reference ECEF difference vector to give transformed x,y,z coordinates.
    y=m[3]*v[0]+m[4]*v[1]+m[5]*v[2];
    z=m[6]*v[0]+m[7]*v[1]+m[8]*v[2];
    return {x,y,z};}                                                           //Return the transformed coordinates with +x facing East, +y facing North, and +z facing Up, relative to reference point.

///enu2geo(): ENU x,y,z (km) to geodetic lat(deg),lon(deg),depth(km) coordinate transformation function. Uses lat_ref, lon_ref, d_ref as ENU origin reference.///
///Format: std::vector<double> coords= enu2geo(double lat_ref, double lon_ref, double d_ref, double x, double y, double z);                                   ///
std::vector<double> enu2geo(double lat_ref, double lon_ref, double d_ref, double x, double y, double z){
    double xi,yi,zi, xr,yr,zr, phi,theta, lat,lon,r, rad=57.295, R=6371;       //Initialize variables, define radian to degree conversion, and earth radius R.
    xr=(R-d_ref)*cos(lon_ref/rad)*cos(lat_ref/rad);                            //Find ECEF Cartesian reference x coordinate using transformation: x=r*cos(lon)*cos(lat)
    yr=(R-d_ref)*sin(lon_ref/rad)*cos(lat_ref/rad);                            //Find ECEF Cartesian reference y coordinate using transformation: y=r*sin(lon)*cos(lat)
    zr=(R-d_ref)*sin(lat_ref/rad);                                             //Find ECEF Cartesian reference z coordinate using transformation: z=r*sin(lat)
    phi=lat_ref/rad ; theta=lon_ref/rad ;                                      //Convert reference latitude and longitude from degrees to radians, as phi,theta respectively.
    double m[9]={-sin(theta), -sin(phi)*cos(theta), cos(phi)*cos(theta),       //Transformation matrix from ENU to ECEF coordinates. phi=latitude, theta=longitude.
                 cos(theta) , -sin(phi)*sin(theta), cos(phi)*sin(theta),
                       0    ,      cos(phi)       ,       sin(phi)     };
    xi=m[0]*x+m[1]*y+m[2]*z+xr;                                                //Take inner product of rotation matrix with sample-reference ECEF difference vector to give transformed x,y,z coordinates.
    yi=m[3]*x+m[4]*y+m[5]*z+yr;
    zi=m[6]*x+m[7]*y+m[8]*z+zr;
    r=sqrtf(xi*xi+yi*yi+zi*zi); lat=rad*asin(zi/r); lon=rad*atan2(yi,xi);      //Transform sample ECEF coordinates to spherical [lat,lon,radius] in [degrees,degrees,km] sample coordinates.
    return {lat, lon, R-r}; }                                                  //Return transformed latitude(deg), longitude(deg), and depth(km) of ENU sample coordinates.

///Overloaded enu2geo, geo2enu functions using coordinate vectors instead of individual doubles.///
std::vector<double> geo2enu(std::vector<double> r, std::vector<double> c){return geo2enu(r[0],r[1],r[2],c[0],c[1],c[2]);}
std::vector<double> enu2geo(std::vector<double> r, std::vector<double> c){return enu2geo(r[0],r[1],r[2],c[0],c[1],c[2]);}

///FMM2DP is a class implementation of Sethian & Popovici: 'Fast marching methods', SIAM REVIEW 1999, and a variation of Vidale, 1990: "Finite difference calculation of traveltimes in three dimensions."///
///FMM models the first arrival times from either a disturbance at defined point source or from an initialized time field. Axes r,theta correspond to indices i,j, with dimensions di,dj and spacing h,w. ///
///{r=i*h; theta=j*h}; Mapping to 1D: n[i,j]=i+j*di; to 2D: i=n%di, j=n/di. The slowness field must be defined as u[n]=1/v[i+j*di+k*di*dj] where: r[n]=(n%di)/h, theta[n]=(n/di)/w.                       ///
///Required method order: 1. event E= FMM();  2. Init();  3. March();  Method formatting described in detail above each method.                                                                           ///
class FMM2DP{   /////////////////////////           Start FMM class.           /////////////////////////
public:                                                                        //Set all elements public.
    int di, dj, N; double h, w;                                                //Declare axis dimensions 'di(radial axis)', 'dj(polar axis)', grid size 'N', spacing 'h'(radial), w(polar).
    std::vector<double> u, t; std::vector<char> s;                             //Declare slowness vector u[], time vector t[], state vector s[] where s[]={0=far, 1=near(heap), 2=frozen, 3=boundary(in heap)}.
    std::vector<int> heap, hind;                                               //Declare binary minheap 'heap' to hold indices of heap nodes, heap index array 'hind' for storing index of t[n] in 'heap'.

    ///Heap update method [Sedgewick, Wayne: Algorithms, 4th ed. p315-319]. Exchanges heap[i] and heap[j], update linked list hind. Format: void exch(int i, int j); i=heap index 1. j=heap index 2.///
    void exch(int i, int j){
        int q=heap[i]; heap[i]=heap[j]; heap[j]=q;                             //Swap values in heap[i] and heap[j].
        hind[heap[i]]=i; hind[heap[j]]=j;}                                     //Update linked list hind with new ranks of nodes in heap[i] and heap[j].

    ///Heap update method [Sedgewick, Wayne: Algorithms, 4th ed. p315-319]. Propagates heap[k] up binary minheap until heap condition is met. Format: void swim( int k ); k=heap index.              ///
    void swim(int k){while(k>1 && t[heap[k/2]]>t[heap[k]]){exch(k/2,k); k/=2;}}//While(k>1 & t[heap[k/2]]>t[heap[k]]) Exchange child/parent inds, set active ind where violating parent ind previously was.

    ///Heap update method [Sedgewick, Wayne: Algorithms, 4th ed. p315-319]. Propagates heap[k] down binary minheap until heap condition is met. Format: void sink(unsigned int k); k=heap index.     ///
    void sink(unsigned k){
       while (2*k+2<=heap.size()){                                             //Swim down heap while test index is in the upper/first half of 'heap'.
          unsigned j=2*k;                                                      //Set index for comparison at 2x test index.
          if(j<heap.size() && t[heap[j]]>t[heap[j+1]]){j++;}                   //If time at test index t[heap[j]] is more than t[heap[j+1]]], advance to second child node j+1.
          if(t[heap[j]]>=t[heap[k]]){break;}                                   //If t[heap[j]]>t[heap[k]], break.
          exch(k,j); k=j;} }                                                   //Exchange test index k with index j. Set new k at j.

    ///Heap insertion method: Inserts index 'n' to end of binary minheap and bubbles it up heap until heap condition is restored. Format: insrt(int n); n=index of t to insert.                      ///
    void add(int n){int j=heap.size(); heap.push_back(n); hind[n]=j; swim(j);} //Add n to end of heap; Restore heap structure using swim() on last element of heap.

    ///Vidale() is a Vidale plane wave Eikonal solver, for recomputing arrival time at t[n] using more points than Eikonal. Format: void Vidale(int n);  n=grid point index at which to recompute time. ///
    void Vidale(int n){
        int na,nb; double a,b,c,m,q,l;                                         //Declare upwind axial neighbor index offsets i,j, neighbor times a,b,c,  test time tt.
        na=1; if(t[n+na]>t[n-na]){na=-na;} nb=di; if(t[n+nb]>t[n-nb]){nb=-nb;} //Compute the offset index of the upwind neighbors on i and j axis.
        a=t[n+na+nb]; b=t[n+na]; c=t[n+nb];                                    //Extract upwind axial and diagonal neighbor times: a(upwind diagonal), b(upwind azimuthal), c(upwind radial).
        q=pow(w*(n%di)*h,2); l=h*h;                                            //Compute squares of axial grid spacings, with hh= radial spacing^2, qq=linearized azimuthal spacing^2.
        if(s[n+na]==2 && s[n+nb]==2 && s[n+na+nb]==2){                         //In cases where a 3 point square of upwind neighbors is frozen, travel time is computed using an anisotropic variant of (Vidale, 1990).
            m=pow((u[n%di]+u[(n+na)%di])/2,2);                                  //Compute square of the average slowness 'm'.
            t[n]=fmin(t[n], ((a-b+c)*l+(a+b-c)*q+2*sqrt(l*q*(m*(l+q)-(b-c)*(b-c))))/(l+q));}}//Solve for tt: ((c+tt-a-b)/2h)^2+((b+tt-a-c)/(w*r))^2=u^2. If tt<t[n], set as new t[n].

    ///Eikonal() computes time at t[n] using an upwind Eikonal scheme (Sethian & Popovici: 'Fast marching methods', SIAM REVIEW 1999). Algorith: Find upwind neighbor times (t[neighbor]<t[n]) on both axes; find    ///
    ///average slowness of target and upwind neighbors. If fabs(b-a)<s (2 immediate neighbors are upwind) solve for tt: (tt-b)^2/h^2+(tt-c)^2/q^2=u^2  ;  Else, solve: (tt-a)^2=u^2  ;  Set:  Set t[n]=min(tt,t[n]). ///                                                                                                                                                                      ///
    ///Format: void Eikonal(int n);   n=grid index at which to recompute time.                                                                                                                                       ///
    void Eikonal(int n){
        int na,nb; double a,b,f,q,l,tt;                                        //Declare upwind axial neighbor index offsets ii/jj/kk, upwind axial times a/b/c, test time tt, and time step 'f'.
        na=1; if(t[n+na]>t[n-na]){na=-na;} nb=di; if(t[n+nb]>t[n-nb]){nb=-nb;} //Compute the offset index of the upwind neighbor on i axis as 'na', on j axis as 'nb'.
        a=t[n+na]; b=t[n+nb]; f=(2*u[n%di]+u[(n+na)%di])/3;  q=w*(n%di)*h;     //Set slowness using the average of u[n] and u at each upwind neighbor. Assign upwind axial neighbor times on radial/polar axes as a/b.
        l=h*f; if(b<a){l=q*f;}                                                 //Set the primary propagation axis unit time (if propagation is primarily radial, s=h*f, if primarily azimuthal; s=w*f.
        if(fabs(b-a)<l){tt=(b*h*h+a*q*q+h*q*sqrt((h*h+q*q)*f*f-(a-b)*(a-b)))/(h*h+q*q);}//If upwind neighbor times differ by less than the primary axial time step, solve (tt-b)^2/h^2+(tt-c)^2/q^2=u^2 for tt.
        else {tt=fmin(a+f*h, b+f*q);}                                          //Else, set tt=a+f.
        t[n]=fmin(t[n],tt);}                                                   //If test time 'tt' is faster than current time at t[n], set test time as new t[n].

    ///.Init(): Initializes times, states, and heap near a point source. Nodes near shell of radius 'r0' centered at (x0,y0) are set to s=1(active) and added to minheap 'heap'; points inside are set to s=2(frozen).///                                                            ///
    ///Format:   void .Init(double i0, double j0, double r0) ; r0=grid distance to which to initialize. i0,j0= grid points coords of internal source using mapping: r=i*h, a=j*w.                                     ///                                                                                                   //
    void Init(double i0, double j0, double r0){                                //Internal point source initializer
        t.assign(N,10000); int i,j; double r;                                  //Initialize time grid t, with N grid points, initialized at value t_max.
        for(int n=0; n<N; n++){                                                //Iterate through every index to seed 't', 's', and 'heap' within an 'init_radius' of the source with initial values.
            i=n%di ; j=n/di ; r=h*sqrt(i0*i0+i*i-2*i0*i*cos(w*(j-j0)));        //Compute axial indices i,j at linear index n, and separation distance in unit coordinates 's'.
            if(i==0 || i==di-1 || j==0 || j==dj-1){s[n]=3;}                    //Deactivate all grid boundary nodes setting their s=3 to prevent their addition to 'heap'.
            else if(r<=r0){                                                    //Else, if n is located inside the initialization radius 'r_0'
                t[n]=r*(u[i]+u[i0]+u[(i0+i)/2])/3;  s[n]=2;                    //Compute direct travel time from source coordinates to n using average slowness, set s[n]=2 (frozen).
                if(r>r0-1.5*h){s[n]=1; add(n);}}}}                             //If n is within 1.5 grid points of the initialization shell radius and not a border point, set state s=1(active) and add n to minheap.

    ///.March(): Implements Sethian and Popovici's Fast Marching Method(FMM) and Vidale's finite difference method. Algorithm:                                           ///
    ///  Start from {t, s, heap} initialized by Init(), iterating until len(heap)=0 (all nodes are converted to "frozen"(s=2)).                                          ///
    ///  Extract index of minimal active t[n0=heap[0]] and freeze the node(s[n0]=2), apply Vidale's scheme to t[n0].                                                     ///
    ///  Update times at every non frozen(s[n]<2) bordering node 'n' by applying upwind Eikonal scheme.                                                                  ///
    ///     If n was already in the heap(s[n]==1), restore heap structure by looking up nodes position in heap in hind[n], and using swim(hind[n]).                      ///
    ///     Else if n was 'far'(s==0), add n to the end of heap and restore heap structure using swim(heap[end]).                                                        ///
    ///  Frozen node n0 is removed from heap, heap structure is restored by swapping heap[0]:heap[-1], removing heap[-1], de-indexing hind[n0], and sinking new heap[0]. ///
    ///Format: .March();                                                                                                                                                 ///
    void March(){
        while(heap.size()>0){                                                  //March active nodes downwind until 'heap' is empty (indicating that all nodes are set to their minimal time "frozen" s).
            int n, n0=heap[0], i0=n0%di, j0=n0/di; s[n0]=2; Vidale(n0);        //Extract minimal time index n0 from top of heap, its axial indices i0/j0, Recompute the travel time using Vidale(), set state s=2(frozen)
            for(int i=i0-1; i<=i0+1; i++){for(int j=j0-1; j<=j0+1; j++){       //Iterate through all points neighboring n0. (i axis neighbor iteration).
                n=i+j*di;                                                      //Calculate neighboring node's linear index n.
                if(s[n]<2){                                                    //Proceed if neighboring/test node is not frozen(inludes all finalized times and bounding nodes).
                    Eikonal(n);                                                //Execute Eikonal() to update the node's time based on 1st order axial upwind neighboring node times.
                    if(s[n]==1){swim(hind[n]);}                                //If n is 'active'(s[n]==1), update heap by running swim() with updated t[n] up from its current position in heap indexed at hind[n].
                    else{s[n]=1; add(n);}}}}                                   //Else, if neighboring node is currently "far"(s=0), set s=1 (active), add its to index to heap.
            exch(0,heap.size()-1); heap.pop_back(); sink(0);}}                 //Move last element of 'heap' to front, delete last element of 'heap', restore heap structure by sinking heap[0].

    ///FMM constructor. Create 1D binary minheap 'heap' to hold active indices n (s[n]=2), and array 'hind' to hold heap indices of n. Format: FMM( int di, int dj, int dk, float h, std::vector<float> u). ///
    ///  {di,dj,dk}={x dimension, y dimension, z dimension};  h=grid point spacing(km);  u=1D rastered slowness grid: u[n]=1/v[i+j*di+k*di*dj]; x[n]=(n%di)/h, y[n]=((n/di)%dj)/h, z[n]=(n/(di*dj))/h       ///
    FMM2DP(int _di, int _dj, double _h, std::vector<double>_u){di=_di; dj=_dj; h=_h; w=3.141593/dj; u=_u; N=di*dj; s.assign(N,0); hind.assign(N,0);}//Pass params, Initialize states 's'(0=far, 1=active(heap), 2=frozen), 'hind' to store indices of t in heap.
};/////////////////////////                   End FMM class                    /////////////////////////


///NOTICE: If SET() is instantiated but lookup table "set.dat" has not been pre-generated, generation is done at runtime, and may take several hours. This is done only once. ///
///SET: Standard earth first arrival time lookup table generation and interpolation class. When Instantiated, loads lookup table t[] from 'set.dat'. If unsuccessful, .t_generate() generates a dense 2D full earth  ///
///traveltime grid in polar coordinates (to r=6371km, arc=180 degrees), with 'bins' grid points per km radially, 100*bins grid points per degree arc on the polar axis. "blur" sets the width of the region (in km)  ///
///over which discontinuities are linearly interpolated. Velocity model is described by iasp91()(Kennet & Engdahl 1991 standard earth velocity). NOTE: After scoping out, t must be freed using free(InstanceName.t) ///
class SET{      /////////////////////////           Start SET class.           /////////////////////////
public:
    std::vector<float> t; int dr, da, ds, h;                                   //Declare iasp first arrival lookup table t[], source/arc/receiver dims {ds/da/dr}, and spacing h.

    ///.time(): Overloaded iasp91 first arrival function. Returns standard earth travel time between source/receiver using either 2 sets of {lat,lon,d}, or a source depth, great circle arc, and receiver depth. ///
    ///Format: double travel_time(s)=time( double source_depth(km), double arc_angle(degrees), double receiver_depth(km)=0).                                                                                      ///
    double time(double source_depth, double arc_angle, double receiver_depth=0){
        double x=receiver_depth, y=arc_angle, z=source_depth;                  //Set source depth as z coordinate(index k), constrain to z range of the grid; Bring outliers in to the face which they bound.
        if(x+z<19 && y<1){return sqrt((x-z)*(x-z)+12364.31*y*y)/5.8  ;}        //If combined source and receiver depths are less than 18km, and the arc is less than 1 degree, return direct ray travel time.                                                                                          //
        else             {return interp3d(t,x,y,z,dr,da,ds,h,.01*h,h);}}       //Else, interpolate times in lookup table 't' and return.
    ///Format: double travel_time(s)=time(double lat1, double lon1, double depth1, double lat2, double lon2, double depth2)                                                                                       ///
    double time(double lat1, double lon1, double d1, double lat2, double lon2, double d2){return time(d1,gca(lat1,lon1,lat2,lon2),d2);}
    ///Format: double travel_time(s)=time(std::vector<double> {lat1,lon1,depth1}, std::vector<double> {lat2,lon2,depth2})                                                                                       ///
    double time(std::vector<double> s, std::vector<double> r){return time(s[2],gca(s[0],s[1],r[0],r[1]),r[2]);}

    ///SET constructor/generator. Imports standard earth first arrival lookup table t[] from binary file 'set.dat'. If t is not imported successfully, calls .t_generate() to generate/save t[]. ///
    ///Format: SET Instance = SET(int bins=2, double blur=2); bins=grid points per km used for traveltime generation. blur=iasp91 discontinuity linear blurring (width in km).         ///
    SET(int bins=2 , double blur=2){
        dr=140; da=3600; ds=140; h=5; t.assign(dr*da*ds,0);                    //Set receiver/arcs/source dims in t, receiver/arc/source spacing (km)/(.01 deg)/(km), stride(km spacing for sources/receivers).
        FILE *o=fopen("set.dat","a+"); int l=fread(t.data(),4,t.size(),o);     //Open binary file SET.dat, read in array t[] containing travel times through IASP91.
        if(l!=dr*da*ds){ printf("generating times \n");                        //If number of elements successfully read in does not match the specified tt_len, generate a new traveltime set.
            int R=6371, b=bins, di=b*R+3, dj=18000*b; std::vector<double>u(di);//Set earth radius R, generation grid radial points/km, radial/polar dims 'di/dj', Declare radial IASP91 slowness array 'u'.
            for(int i=0; i<di; i++){u[i]=1/iasp91(R-i*1./b,blur);}             //Fill radial slowness array with 1/iasp91(d) where d ranges from 0 to center of earth in increments 1/blur (km).
            for(int s=0; s<ds; s++){ printf("d= %i \n",s*h);                   //Iterate through each source depth. Print source depth of currently generated time grid.
            FMM2DP E(di,dj,1./b,u); E.Init(b*(R-s*h),1,9); E.March();          //Instantiate event E at specified source coordinates and grid dims, initialize the fields, and march.
            for(int a=0; a<da; a++){ for(int r=0; r<dr; r++){                  //Iterate through each tabulated arc angle and receiver angle.
                t[r+a*dr+s*dr*da]=E.t[R*b-r*b*h+a*b*h*di+di];}}}               //Set travel time lookup table tt[] to the corresponding grid point on FMM2DP grid.
            fwrite(t.data(),4,t.size(),o); printf("generation complete\n");}   //Save generated travel time array t as binary set.dat. To print to text: fprintf(out,"%.4f ",tt[n])
        fclose(o);}                                                            //Close open binary set.dat
};   /////////////////////////           End SET class.           /////////////////////////

//int main(){SET S; printf("\nSET t= %.3f \n",S.time(241.55,3.138));}          //To compile:' g++ -o SET SET.cpp -std=c++11 -O3 && ./SET '    //vec dat(&S.t[0], &S.t[2*140*20]); savetxt(dat,140,"data");
