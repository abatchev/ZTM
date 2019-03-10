/// No non-STL dependencies. Language standard: C++14.
/// This library is released for public use as material supplementary to PhD thesis: Development of Specialized Non-Linear Inversion Algorithms, Basis Functions, Eikonal Solvers, and Their Integration for Use in Joint Seismic and Gravitational Tomographic Inversion
/// Author: Zagid Abatchev - abatchev@ucla.edu ; PI: Paul Davis. Department of Earth, Planetary, and Space Sciences; University of California, Los Angeles; Date of publication: March 2019.///
/// This library is published as is and may have unresolved bugs. Use at own discretion. Users are encouraged to adapt some or all of this code without restriction. Refer to Github page for updated source code, located at: https://github.com/abatchev/ZTM .

///FMM is a C++ implementation of Sethian & Popovici: 'Fast marching methods', SIAM REVIEW 1999, and a variant of Vidale, 1990: "Finite difference calculation of traveltimes in three dimensions."///
///FMM models the first arrival times from either a disturbance at defined point source or from an initialized time field.
///No non-STL dependencies. Language standard: C++14.
///Axes {x,y,z} correspond to indices {i,j,k} with dimensions {di,dj,dk} and spacing 'h'. {x=i*h; y=j*h; z=k*h}; Mapping to 1D: n[i,j,k]=i+j*di+k*di*dj; to 3D: i=(n%di), j=((n/di)%dj)/h , k=(n/(di*dj)).
///Slowness field must be defined as u[n]=1/v[i+j*di+k*di*dj] where: x[n]=(n%di)/h, y[n]=((n/di)%dj)/h, z[n]=(n/(di*dj))/h.
///Required method order: 1. FMM();  2. Init();  3. March();  Formatting described in detail above each method.
#include<math.h>                                                               //Include STL math library.
#include<vector>                                                               //Include STL vector library.

class FMM{public: /// //////////////////////          Start FMM class. All members public.          ////////////////////// ///
    int di, dj, dk, N; double h, t_max;                                        //Declare axial dimensions {di(x axis), dj(y axis), dk(z axis)}, grid size 'N', grid spacing 'h', max time 't_max'.
    std::vector<double> u, t; std::vector<char> s;                             //Declare slowness vector u, time vector t, state vector s{0=far, 1=near(heap), 2=frozen, 3=boundary}.
    std::vector<int> heap, hn;                                                 //Declare binary minheap 'heap' to hold indices of heap nodes, heap index array 'hn' for storing index of t[n] in 'heap'.

    ///Heap update method [Sedgewick, Wayne: Algorithms, 4th ed. p315-319]. Exchange heap[i] and heap[j], update linked list hn. Format: void exch(int i, int j); i, j: heap indices to swap.///
    void exch(int i, int j){int q;
        q=heap[i]; heap[i]=heap[j]; heap[j]=q; hn[heap[i]]=i; hn[heap[j]]=j;}  //Swap values in heap[i] and heap[j], Update linked list hn with new ranks of nodes in heap[i] and heap[j].

    ///Heap update method [Sedgewick, Wayne: Algorithms, 4th ed. p315-319]. Propagate heap[k] up binary minheap until heap condition is met. Format: void swim( int k ); k=heap index.         ///
    void swim(int k){while(k>1 && t[heap[k/2]]>t[heap[k]]){exch(k/2,k); k/=2;}}//Exchange child/parent inds, set active ind where violating parent ind previously was.

    ///Heap update method [Sedgewick, Wayne: Algorithms, 4th ed. p315-319]. Propagate heap[k] down binary minheap until heap condition is met. Format: void sink(unsigned int k); k=heap index.///
    void sink(unsigned k){
       while(2*k+2<=heap.size()){ unsigned j=2*k;                              //Swim down heap while test index is in the upper/first half of 'heap'. Set index for comparison at 2x test index.
          if(j<heap.size() && t[heap[j]]>t[heap[j+1]]){j++;}                   //If time at test index t[heap[j]] is more than t[heap[j+1]]], advance to second child node j+1.
          if(t[heap[j]]>=t[heap[k]]){break;}                                   //If t[heap[j]]>=t[heap[k]], break.
          exch(k,j); k=j;} }                                                   //Exchange test index k with index j. Set new k at j.

    ///Heap insertion method: Insert index 'n' to end of binary minheap, bubble it up heap until heap condition is restored. Format: insrt(int n); n=index of t to insert.               ///
    void add(int n){int j=heap.size(); heap.push_back(n); hn[n]=j; swim(j);}   //Add n to end of heap; Restore heap structure using swim() on last element of heap.

    ///Eikonal() computes time at t[n] using a 2nd order upwind Eikonal scheme (Sethian & Popovici: 'Fast marching methods', SIAM REVIEW 1999) Algorithm:                                ///
    ///  Find and sort upwind (points where t<t[n]) neighbor times on each axis from least to greatest as a<b<c; Set time step f=h*u[n] using spacing h, slowness u[n].                  ///
    ///  If a+f>c: solve (t-a)^2+(t-b)^2+(t-c)^2=f^2 for t. set t[n]=min(t,t[n]). If a+f>b: solve (t-a)^2+(t-b)^2=f^2 for t. Set t[n]=min(t,t[n]). If a+f<t[n], set t[n]=a+f;            ///
    ///Format: void Eikonal(int n);   n=grid index at which to recompute time.                                                                                                           ///
    void Eikonal(int n){
        int na,nb,nc; double a,b,c,q,f=h*u[n];                                 //Declare upwind axial neighbor index offsets na/nb/nc, upwind axial times a/b/c, temp holder 'q', set unit step time f=h*u[n].
        na=1    ; if(t[n+na]>t[n-na]){na=-na;} a=t[n+na];                      //Compute the offset index of the upwind neighbor on i axis as 'na'; Assign axis i upwind neighbor time 'a'.
        nb=di   ; if(t[n+nb]>t[n-nb]){nb=-nb;} b=t[n+nb];                      //Compute the offset index of the upwind neighbor on j axis as 'nb'; Assign axis j upwind neighbor time 'b'.
        nc=dj*di; if(t[n+nc]>t[n-nc]){nc=-nc;} c=t[n+nc];                      //Compute the offset index of the upwind neighbor on k axis as 'nc'; Assign axis k upwind neighbor time 'c'.
        if(a>b){q=a;a=b;b=q;} if(b>c){q=b;b=c;c=q;} if(a>b){q=a;a=b;b=q;}      //Sort upwind neighbor times so that a<b<c; 'a' corresponds to major axis of propagation, 'c' to minor axis.
        if(a+f>c){t[n]=fmin(t[n],(a+b+c+sqrt(3*f*f-2*(a*a+b*b+c*c-a*b-b*c-c*a)))/3);}//If (a+f>c): Solve (t-a)^2+(t-b)^2+(t-c)^2=f^2 for t. If t<t[n], set t[n]=t.
        else if(a+f>b){t[n]=fmin(t[n],(a+b+sqrt(2*f*f-(a-b)*(a-b)))/2);}       //if (a+f>b): Solve (t-a)^2+(t-b)^2=f^2 for t. If t<t[n], set t[n]=t.
        if(a+f<t[n]){t[n]=a+f;}}                                               //If a+f<t[n], set t[n]=a+f.

    ///Eikonal2() computes time at t[n] using a 2nd order upwind Eikonal scheme (Sethian & Popovici: 'Fast marching methods', SIAM REVIEW 1999) Algorithm:                                        ///
    ///  Find upwind(t[neighbor]<t[n]) neighbor t on each axis, sort so that a<b<c. Assign further coaxial upwind neighbors aa/bb/cc; Set time step f=h*u[n] using spacing h, slowness u[n].      ///
    ///  Else, set t=a+u*h. if t<t[n], set t[n]=t///                                                                                                                                              ///
    ///  If consecutive axial neighbors are frozen on 3 axes, solve: ua^2+ub^2+uc^2=f^2: ua=3*t-4*a+aa; ub=3*t-4*b+bb; uc=3*t-4*c+cc; Set t[n]=min(t,t[n]).                                       ///
    ///  If a+f>c: solve (t-a)^2+(t-b)^2+(t-c)^2=f^2. Set t[n]=min(t,t[n]).                                                                                                                       ///
    ///  If consecutive axial neighbors are frozen on 2 axes, solve: ua^2+ub^2=u^2 for t: ua=3*t-4*a+aa; ub=3*t-4*b+bb; t=(4*a-aa+4*b-bb+sqrt(8*f*f-(4*a-aa-4*b+bb)^2))/6. Set t[n]=min(t,t[n]).  ///
    ///  Else if a+f>b: solve (t-a)^2+(t-b)^2=f^2 for t. Set t[n]=min(t,t[n]); Set t[n]=min(a+f,t[n]).                                                                                            ///
    ///Format: void Eikonal(int n) ; n=grid index at which to recompute time.                                                                                                                     ///
    void Eikonal2(int n){
        int na,nb,nc,q; double a,aa,b,bb,c,cc,f=h*u[n];                        //Declare upwind axial neighbor index offsets na/nb/nc, upwind axial times a/aa, b/bb, c/cc, time step 'f', temp holder 'q'.
        na=1    ; if(t[n+na]>t[n-na]){na=-na;}                                 //Compute the offset index of the upwind neighbor on i axis as 'na';
        nb=di   ; if(t[n+nb]>t[n-nb]){nb=-nb;}                                 //Compute the offset index of the upwind neighbor on j axis as 'nb'; Assign axis j upwind neighbor time 'b', next neigbor time 'bb'
        nc=dj*di; if(t[n+nc]>t[n-nc]){nc=-nc;}                                 //Compute the offset index of the upwind neighbor on k axis as 'nc'; Assign axis k upwind neighbor time 'c', next neigbor time 'cc'
        if(t[n+na]>t[n+nb]){q=na;na=nb;nb=q;}                                  //Sort upwind neighbor times so that a<b<c; 'a' corresponds to major axis of propagation, 'c' to minor axis.
        if(t[n+nb]>t[n+nc]){q=nb;nb=nc;nc=q;}                                  //Sort upwind neighbor times so that a<b<c; 'a' corresponds to major axis of propagation, 'c' to minor axis.
        if(t[n+na]>t[n+nb]){q=na;na=nb;nb=q;}                                  //Sort upwind neighbor times so that a<b<c; 'a' corresponds to major axis of propagation, 'c' to minor axis.
        a=t[n+na];aa=t[n+2*na]; b=t[n+nb];bb=t[n+2*nb]; c=t[n+nc];cc=t[n+2*nc];//Assign sorted axial upwind neighbor times 'a/b/c', next neigbor times 'aa/bb/cc'.
        if(a+f>c && s[n+nc]==2 && s[n+2*nc]==2 && s[n+nb]==2 && s[n+2*nb]==2 && s[n+na]==2 && s[n+2*na]==2){ //If 2 consecutive axial neighbors are frozen on 3 axes:
            t[n]=fmin(t[n], (4*(a+b+c)-(aa+bb+cc)+sqrt(fabs(pow(aa+bb+cc-4*(a+b+c),2)-3*(pow(aa-4*a,2)+pow(bb-4*b,2)+pow(cc-4*c,2)-4*f*f))))/9);}//Solve 3D order 2 Eikonal: (3*t-4*a+aa)^2+(3*t-4*b+bb)^2+(3*t-4*c+cc)^2=u^2.
        if(a+f>c){t[n]=fmin(t[n],(a+b+c+sqrt(3*f*f-2*(a*a+b*b+c*c-a*b-b*c-c*a)))/3);} //Else if a+f>c Solve (t-a)^2+(t-b)^2+(t-c)^2=f^2 for t. If t<t[n], set t[n]=t.
        if(a+f>b && s[n+nb]==2 && s[n+2*nb]==2 && s[n+na]==2 && s[n+2*na]==2){t[n]=fmin(t[n],(4*a-aa+4*b-bb+sqrt(fabs(8*f*f-pow(4*a-aa-4*b+bb,2))))/6);}//Solve 2D order 2 Eikonal: ua^2+ub^2+uc^2=u^2 ; ua=3*t-4*a+aa; ub=3*t-4*b+bb; uc=(t-tc)/h;
        if(a+f>b){t[n]=fmin(t[n],(a+b+sqrt(2*f*f-(a-b)*(a-b)))/2);}            //Else if a+f>b solve (t-a)^2+(t-b)^2=f^2 for t. If t<t[n], set t[n]=t.
        if(a+f<t[n]){t[n]=a+f;}}                                               //Test t=a+f. If t<t[n], set t[n]=t.

    ///Vidale() is an Eikonal solver for computing time at t[n]. Applies eq(2) from (Vidale, 1990: Finite difference calculation of traveltimes in three dimensions.)                                         ///
    ///Find upwind (t<t[n]) neighbor index offsets on each axis, sort so that t[a]<t[b]<t[c]. If 7 required points are frozen apply eq2. Else if 6 points required for scheme B are frozen, solve eq(3).      ///
    ///Format: void Vidale(int n);  n=grid point index at which to recompute time.                                                                                                                            ///
    void Vidale(int n){
        int a,b,c,q; double f;                                                                            //Declare upwind axial neighbor index offsets a,b,c, temp placeholder q, time step 'f'.
        a=-1; if(t[n+a]>t[n-a]){a=-a;} b=-di; if(t[n+b]>t[n-b]){b=-b;} c=-di*dj; if(t[n+c]>t[n-c]){c=-c;} //Compute offset indices of upwind neighbor on i/j/k axes as a/b/c.
        if(t[n+a]>t[n+b]){q=a;a=b;b=q;} if(t[n+b]>t[n+c]){q=b;b=c;c=q;} if(t[n+a]>t[n+b]){q=a;a=b;b=q;}   //Sort upwind neighbor index offsets to enforce t[n+a]<t[n+b]<t[n+c].
        if(s[n+c]==2 && s[n+b+c]==2 && s[n+c+a]==2 && s[n+a+b]==2 && s[n+a+b+c]==2 && s[n+b]==2 && s[n+a]==2){//If 7 upwind neighbors required in "scheme A" are frozen, apply (eq2)(Vidale, 1990).
            f=h*(u[n]+u[n+a]+u[n+b]+u[n+c]+u[n+a+b]+u[n+b+c]+u[n+c+a]+u[n+a+b+c])/8;                      //Compute unit step time f=(h*u_avg)^2 where h=grid spacing, u_av= 8 point slowness average.
            t[n]=fmin(t[n], t[n+a+b+c]+sqrt(fabs(6*f*f-pow(t[n+a]-t[n+b],2)-pow(t[n+b]-t[n+c],2)-pow(t[n+c]-t[n+a],2)-pow(t[n+a+b]-t[n+b+c],2)-pow(t[n+b+c]-t[n+c+a],2)-pow(t[n+c+a]-t[n+a+b],2))/2));}//apply (eq2).
        else if(s[n+a]==2 && s[n+b]==2 && s[n+a+b]==2 && s[n+a+b+c]==2 && s[n+a+b-c]==2){                 //If 5 upwind neighbors required in "scheme B" are frozen, apply (eq3)(Vidale, 1990).
            f=h*(u[n]+u[n+a]+u[n+b]+u[n+a+b]+u[n+a+b+c]+u[n+a+b-c])/6;                                    //Compute unit step time f=(h*u_avg)^2 where h=grid spacing, u_av= 6 point slowness average.
            t[n]=fmin(t[n],t[n+a+b]+sqrt(fmax(0, 2*f*f-.5*pow(t[n+a+b+c]-t[n+a+b-c],2)-pow(t[n+a]-t[n+b],2))));}} //apply (eq3). Accept as new t[n] if smaller than current t[n].

    ///Init(): Overloaded FMM method to initialize times, states, and binary minheap around a point source, or pre-initialized times.                                                                            ///
    ///  If 4 doubles are passed, initialize point source case. Nodes near shell of radius 'r0' centered at (x0,y0,z0) are set to s=1(active) and added to minheap 'heap'; points inside are set to s=2(frozen). ///
    ///  If std::vector<double> is passed, time field 't' is set equal to passed array, times less than t_max are set to s=1 and added to 'heap'.                                                                ///
    ///Format(internal point source):   void .Init(double x0, double y0, double z0, double r0) ; r0=grid distance to which to initialize. x0,y0,z0= coords of internal source using mapping: x=i*h, y=j*h, z=k*h.///                                                                                                   //
    void Init(double x0, double y0, double z0, double r0){
        int n0=fmin(N,fmax(0,rint(x0/h)+rint(y0/h)*di+rint(z0/h)*di*dj));      //Compute the grid index n_start nearest to source coordinates, limiting it to the available coordinate set
        for(int n=0; n<N; n++){ int i=n%di, j=(n/di)%dj, k=n/(di*dj);          //Iterate through every index to seed 't', 's', 'heap' within 'r0' of source with initial values. Compute axial indices i,j,k at n.
            double r=sqrt(pow(i*h-x0,2)+pow(j*h-y0,2)+pow(k*h-z0,2));          //Compute square of the separation distance r_sq between current point n (i/j/k) and source coordinates (x_0,y_0,z_0).
            if(i==0||i==di-1 || j==0||j==dj-1 || k==0||k==dk-1){s[n]=3;}       //Deactivate all grid boundary nodes setting their s=3 to prevent their addition to 'heap'.
            else if(r<=r0){                                                    //Else, if n is located inside the initialization radius 'r_0'
                t[n]=(u[n]+u[n0])*r/2; s[n]=2; if(r>r0-1.5){s[n]=1; add(n);}}}}//Find line travel time from source to n using mean u, set s[n]=2(frozen). If n is near initialization shell and not on boundary, set s=1(active), add to heap.
    ///Format(pre-computed time field): void .Init(std::vector<double> t0) ; t0= initialization time vector with size d*dj*dk using mapping: x=i*h, y=j*h, z=k*h.                                                ///
    void Init(std::vector<double> t_init){ t=t_init;                           //Pre computed time field initializer. Set grid time array t=t_init
        for(int n=0; n<N; n++){ int i=n%di, j=(n/di)%dj, k=n/(di*dj);          //Iterate through every index to seed 't', 's', and 'heap' using a defined initial time field t. Compute axial indices i,j,k at n.
            if(i==0||i==di-1 || j==0||j==dj-1 || k==0||k==dk-1){s[n]=3;}       //Deactivate all grid boundary nodes setting their s=3 to prevent their use or addition to 'heap'.
            else if(t[n]+1<t_max){s[n]=1; add(n);}}}                           //If t[n]<t_max and n is not a border point, set state s=1(active) and add n to minheap.
    ///FAST INIT version for local events. btm sets the lowest bound of where travel times will be evaluated at. Everything below is frozen///
    void Init(double x0, double y0, double z0, double r0, int btm){
        int n0=fmin(N,fmax(0,rint(x0/h)+rint(y0/h)*di+rint(z0/h)*di*dj));      //Compute the grid index n_start nearest to source coordinates, limiting it to the available coordinate set
        for(int n=0; n<N; n++){ int i=n%di, j=(n/di)%dj, k=n/(di*dj);          //Iterate through every index to seed 't', 's', 'heap' within 'r0' of source with initial values. Compute axial indices i,j,k at n.
            double r=sqrt(pow(i*h-x0,2)+pow(j*h-y0,2)+pow(k*h-z0,2));          //Compute square of the separation distance r_sq between current point n (i/j/k) and source coordinates (x_0,y_0,z_0).
            if(i==0||i==di-1 || j==0||j==dj-1 || k<=btm||k==dk-1){s[n]=3;}     //Deactivate all grid boundary nodes setting their s=3 to prevent their addition to 'heap'.
            else if(r<=r0){                                                    //Else, if n is located inside the initialization radius 'r_0'
                t[n]=(u[n]+u[n0])*r/2; s[n]=2; if(r>r0-1.5){s[n]=1; add(n);}}}}//Find line travel time from source to n using mean u, set s[n]=2(frozen). If n is near initialization shell and not on boundary, set s=1(active), add to heap.

    ///March(): Implements Sethian and Popovici's Fast Marching Method Algorithm (FMM) and Vidale's finite difference method. Algorithm initializes {t,s,heap}, iterates until all (s=0) nodes are converted to (s=2).                ///
    ///  Extract index of minimal active t[n0=heap[0]] and freeze the node(s[n0]=2); apply Vidale's scheme to t[n0].                                                ///
    ///  Update times at every non frozen(s[n]<2) bordering node 'n' by applying upwind Eikonal scheme.                                                             ///
    ///     If n was already in the heap(s[n]==1), restore heap structure by looking up nodes position in heap in hn[n], and using swim(hn[n]).                     ///
    ///     Else if n was 'far'(s==0), add n to the end of heap and restore heap structure using swim(heap[end]).                                                   ///
    ///  Frozen node n0 is removed from heap, heap structure is restored by swapping heap[0]:heap[-1], removing heap[-1], sinking new heap[0].                      ///
    ///Format: event.March(int type=1);  type=model type. (0){1st order Eikonal}  (1){1st order Eikonal+Vidale finite difference recompute}  (2){2nd order Eikonal} ///                                                                                                                                                             //
    void March(int type=1){
        while(heap.size()>0){                                                  //March active nodes downwind until 'heap' is empty (indicating that all nodes are set to their minimal time "frozen" s).
            int n,n0=heap[0]; s[n0]=2; if(type==1){Vidale(n0);}                //Extract minimal time index n0 from top of heap, set state s=2(frozen). Recompute t[n0] using Vidale's method.
            for(int i:{-1,0,1}){for(int j:{-1,0,1}){for(int k:{-1,0,1}){       //Iterate through all points neighboring n0.
                n=n0+i+j*di+k*di*dj;                                           //Calculate neighboring node's linear index n.
                if(s[n]<2){if(type<2){Eikonal(n);} else{Eikonal2(n);}          //If neighboring/test node is not frozen(includes all finalized times and bounding nodes), Run 1st order Eikonal update.
                    if(s[n]==1){swim(hn[n]);} else{s[n]=1; add(n);} }}}}       //If n is active(s[n]==1), update heap using swim() with updated t[n] up from current position in heap indexed at hn[n]. Else, if far(s=0), set s=1, add index to heap.
            exch(0,heap.size()-1); heap.pop_back(); sink(0);}}                 //Move last element of 'heap' to front, delete last element of 'heap', restore heap structure by sinking heap[0].

    ///FMM constructor. Creates binary minheap 'heap' to hold active indices(s[n]=1), array 'hn' to hold heap indices. Format: FMM( int di, int dj, int dk, double h, std::vector<double> u).       ///
    ///  {di,dj,dk}={x dimension, y dimension, z dimension};  h=grid point spacing(km);  u=1D rastered slowness grid: u[n]=1/v[i+j*di+k*di*dj]; x[n]=(n%di)/h, y[n]=((n/di)%dj)/h, z[n]=(n/(di*dj))/h ///
    FMM(int _di, int _dj, int _dk, double _h, std::vector<double>_u){
        di=_di; dj=_dj; dk=_dk; h=_h; u=_u; N=di*dj*dk; t_max=1e+9;            //Pass constructor parameters to class members: di, dj, dk, h, u, type; compute total number of grid points N.
        s.assign(N,0); hn.assign(N,0); t.assign(N,t_max); }                    //Initialize vectors 's' to hold t states (0=far, 1=active(heap), 2=frozen), 'hn' to store indices of t in heap, and time array t.
};/// //////////////////////                   End FMM class                   ////////////////////// ///


//#include"zmath.cpp"
/////Shell test. To compile and run:' g++ -o FMM FMM.cpp -std=c++11 -O3 && ./FMM ' . use -lpthread flag to compile multithreaded functions in zmath.cpp.
//int main(){int i,j,k,di,dj,dk,e0,e1,e2,f,m; double h,x0,y0,z0,s0; vec a,t,to,u; std::string id="out.txt";
//    x0=300; y0=300; z0=30; s0=50; e0=61; e1=61; e2=35; f=1; m=2; id="tmt.txt"; //Set source x0,y0,z0 (x=i*h,y=j*h,z=k*h), init radius s0, grid dims [e0,e1,e2](i,j,k), fineness factor f, model type m, file name "id".
//    di=(e0-1)*f+1; dj=(e1-1)*f+1; dk=(e2-1)*f+1; h=10./f;                      //Set model parameters: grid dims di, dj, dk, spacing h.
//    u=loadtxt("u_111.txt");                                                    //Uncomment to import u field from a .txt.
//    //u=vec(e0*e1*e2); for(i=0; i<e0; i++){for(j=0; j<e1; j++){for(k=0; k<e2; k++){u[i+j*e0+k*e0*e1]=7-2.*exp(-(pow(i-20,2)+pow(j-25,2)+pow(k-15,2))/49.);}}} u=1/u; savetxt(u,"u_synth.txt"); //Uncomment to generate synthetic u field.
//    if(f>1){a=u; u=vec(di*dj*dk,1); for(i=0; i<di; i++){for(j=0; j<dj; j++){for(k=0; k<dk; k++){u[i+j*di+k*di*dj]=interp(a,1.*i,1.*j,1.*k,e0,e1,e2,1.*f);}}}} //If f>1, generate linearly interpolated fine grid.
//    FMM E(di,dj,dk,h,u); E.Init(x0,y0,z0,s0); E.March(m); to=E.t;              //Launch FMM as instance E, initialize at start coords
//    if(f>1){for(i=1; i<e0-1; i++){for(j=1; j<e1-1; j++){for(k=1; k<e2-1; k++){to[i+j*e0+k*e0*e1]=E.t[f*i+f*j*di+f*k*di*dj];}}}} //If fineness factor>1, populate baseline grid.
//    savetxt(to,id);}                                                          //Save file, return.
