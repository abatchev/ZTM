/// Non-STL dependencies including "zmath.cpp", "FMM.cpp", "NNI.cpp", and "SET.cpp" needed to run TFM are published together with this file. Language standard: C++14. Requires C++ multithreading flag in compilation. Multi threaded compilation on Ubuntu can be done using -lpthread flag in GCC.
/// This library is released for public use as material supplementary to PhD thesis: Development of Specialized Non-Linear Inversion Algorithms, Basis Functions, Eikonal Solvers, and Their Integration for Use in Joint Seismic and Gravitational Tomographic Inversion
/// Author: Zagid Abatchev - abatchev@ucla.edu ; PI: Paul Davis. Department of Earth, Planetary, and Space Sciences; University of California, Los Angeles; Date of publication: March 2019.///
/// This library is published as is and may have unresolved bugs. Use at own discretion. Users are encouraged to adapt some or all of this code without restriction. Refer to Github page for updated source code, located at: https://github.com/abatchev/ZTM .

/// TFM is a C++ seismic tomography package for 3D inversion of first arrival seismic data.///
/// TFM requires "TFM.cpp", "zmath.cpp", "FMM.cpp", "NNI.cpp", and "SET.cpp" to be in the same directory. Additionally, if gravity inversion is necessary, the same directory must contain gravity data and earth topo data files PGM2012B.txt and ETOPO5.bin.
/// TFM uses standard earth traveltime lookup table SET.dat which needs to be generated only once and is done automatically. This process may take several hours to complete if set.dat is not in the same directory as TFM.cpp.
/// If "grav" and "GPU" options are set to true, TFM attempts to compile GPU gravity model GFM_GPU.cu using Nvidia's NVCC compiler. If it is successful, gravity forward modeling is done on a CUDA capable GPU. If it is not successful, the inversion reverts to CPU based gravity modeling and proceeds normally.
/// Inversion parameters are set inside main().
/// Build+run inversion command: "g++ -o TFM TFM.cpp -lpthread -O3 && ./TFM"

/// Inversion input format:
///     Inversion data must be located in a subdirectory below TFM.cpp named "input". Input data must be named "E#.txt" where # is event index from 0 to last event index.
///     Inversion data is imported from a set of 2D arrays holding station information and pick times with one file per event.
///     Each row in the input file holds arrival time data and event/station information. They may additionally hold waveform data after column index 20. They are imported as space separated .txt files with each row formatted as:
///     [([0] pick time(s)),([1] pick quality(1-3)),([2] blank for TFM),([3] st lat),([4] st lon),([5] st alt(km)),([6] ev lat),([7] ev lon),([8] ev dep(km)),([9]ev offset t(s)),([10-15] ev year,month,day,hour(utc) minute,second),([16] headwave pick sample),([17] pick sample),([18] sample rate),([19] fitted model time(s)),([20]trace sample offset)];
///     a[] columns index 10 and above may be left as zeros, however input data must contain at least 20 columns.
///     a[n][20:], if defined, may contain waveform data or any other data relevant to the pick.
///     Each imported data set is filtered (with bad picks and events removed), and stacked in to full data array 'a' with a common assigned event indices for all picks corresponding to the same events.

/// Inversion output format (primary types, description of secondary types is in out()):
///     Each inversion run is output in to a separate subdirectory with a time stamp in its name. All output files carry the inversion iteration index in their name, as outtype_I.txt where outtype is the type of output information, and I is indexed from 0 upwards, with index 0 indicating initial conditions (pre-inversion).
///     a_I.txt:
///         Stacked event data 'a' is output with a 1 row header which includes:
///         [central reference lat,central reference lon,i axis indices, j axis indices, k axis indices, grid spacing, local events, teleseismic events, Moho Voronoi nodes, mantle Voronoi nodes, inverted points, scaled misfit(normalized ssq), pick time uncertainty, gravity data uncertainty,gravity inversion boolean (1 if yes, 0 if no)]};
///         All other output rows of a_I.txt are formatted as described above.
///     b_I.txt:
///         Complete parameter inversion vector, concatenating the following vectors in 1D:
///             apriori estimated local event coordinates (repeating) as template: [event lat, event lon, event depth(km), event offset time]
///             apriori estimated teleseismic event offset times (repeating): [event offset time]
///             if moho==True : [estimated initial crust velocity(km/s)]
///             if grav==True : starting estimates for: [mantle v_p, mantle crust density contrast (g/cc), mantle density perturbation scaling to velocity perturbation, linear gravity offset (mgal)]
///             if moho==True : repeating for each o2 node: Voronoi nodes defining Moho manifold: [grid i coordinate, grid j coordinate, grid k coordinate]
///             repeating for each o3 node: Voronoi nodes defining the mantle velocity field: [grid i coordinate, grid j coordinate, grid k coordinate, vp (km/s)]
///
///     J_I.txt: Inversion Jacobian.
///     dmp_I.txt: Damping matrix used in the last iteration.
///     u_I.txt: Rendered 1/v_p field as a flattened 1d vector where: field value at [i,j,k]=u[i+j*di+k*di*dj].
///     M_I.txt: Rendered Moho depth mapping as a di x dj array.
///     dG_I.txt: Rendered gravity perturbation field in mgal with same indexing as u_I.txt
///     drho_I.txt: Rendered density perturbation field in g/cc with same indexing as u_I.txt
///     sig_I.txt: Estimated uncertainties corresponding to each element in b_I.txt, of the same dimensionality and units as the corresponding b values.
///     o2_I.txt: if moho==True : 2D array with Voronoi nodes defining Moho manifold with rows defined as: [grid i coordinate, grid j coordinate, grid k coordinate]
///     o3_I.txt: 2D array with Voronoi nodes defining the mantle velocity field with rows defined as: [grid i coordinate, grid j coordinate, grid k coordinate, vp (km/s)]
///     sea.txt, topo.txt: Grid indexed 2D arrays describing sea level k coordinate, and topography k coordinate respectively.

#include "zmath.cpp"                                                           //Include extended STL compatible math library zmath.
#include "FMM.cpp"                                                             //Include FMM source code, containing first arrival Fast Marching Method forward model class FMM.
#include "NNI.cpp"                                                             //Include NNI source code, containing Discrete Sibson Interpolation/Natural Neighbor interpolation functions NNI(), NNIm().
#include "SET.cpp"                                                             //Include SET source code, containing standard earth travel times, great circle arc, and coordinate transformation functions.
struct args{mat a; vec r; int di,dj,dk; double h; int nl,nt,no2,no3; vec sig; vec tsig; vec D; bool grav; double gsig=25; int bfr=10;};//Define TFM parameter struct args.
SET ET; TIMER tmr; int ctr,ftype; string dir; vec GD,sea,topo,M0,RC; bool GPU; //Instantiate SET containing standard earth travel time function .time(), output counter ctr, path dir, forward model type ftype.

vec g2c(vec g){vec v=geo2enu(RC[0],RC[1],RC[2],g[0],g[1],g[2]); return {(v[0]+RC[3])/RC[6],(v[1]+RC[4])/RC[6],(v[2]+RC[5])/RC[6]};}          //Define shorthand for geodetic to enu conversion using globally set reference R. Takes geodetic coordinate vector g  (lat,lon,depth), returns ENU Cartesian vector.
vec c2g(vec c){return enu2geo(RC[0],RC[1],RC[2],c[0]*RC[6]-RC[3],c[1]*RC[6]-RC[4],c[2]*RC[6]-RC[5]);}                                        //Define shorthand for enu to  geodetic conversion using globally set reference R. Takes ENU Cartesian coordinate vector c, returns Geodetic vector (lat,lon,depth).
bool inside(vec c,int s=0){if(c[0]>s && c[0]<RC[7]-1-s && c[1]>s && c[1]<RC[8]-1-s && c[2]>s && c[2]<RC[9]-1-s){return true;} return false;} //Define point test function inside(). Checks if input point coordinate vector c is inside globally set boundaries. Returns true if inside, false if outside. if int s>0, test boundaries are shrunk by 2*s grid points on each axis.

///FM(): Forward Model class.///
class FM{public:
    int nl,nt,di,dj,dk,da,no2,no3,n,o,i,j,k; double h,vl,vh,mk,mm,d,v0,v1,ssq,g0; mat a; vec u,v,M,O,o2,o3,r,t,G,drho,y,f,res,G64,vv;//Declare intrinsic and object variables.
    FM(struct args S, vec B, int fg=1){                                        //Declare initializer/Forward model method.
        a=S.a; di=S.di; dj=S.dj; dk=S.dk; h=S.h; nl=S.nl; nt=S.nt; da=a.size();//Extract data 'a', reference coords 'r', E/N/U axis dims di/dj/dk, grid spacing 'h'(km)
        thread R[nl+nt]; no2=S.no2; no3=S.no3; r=S.r; t=vec(da); G=GD; g0=0;   //Extract local events nl, tele events nt, 2d nodes no2, 3d nodes no3, data points na. Set gravity field G=GD.
        o3=vec(B.end()-4*no3,B.end()); v=NNIm(o3,di,dj,dk); u=1/(v+1e-12); M=vec(di*dj,1.*dk); //Extract o3 Voronoi node set from b, render 3d 1/velocity field 'u'. convert to slowness u=1/v
        if(no2>0){o2=vec(B.end()-4*no3-3*no2,B.end()-4*no3); M=NNI(o2,di,dj); vh=B[4*nl+nt]; //if 2 layer model is used (no2 defined), render moho k coordinates M using o2, extract crustal velocity vh from B.
            for(n=0;n<di*dj*dk;n++){k=n/(di*dj); mk=M[n%(di*dj)];              //iterrate through every point in grid, compute the moho k coordinate at current i,j.
                if(k>mk){u[n]=1/vh;}                                           //if point is above moho at i,j: set u=1/vcrust;
                if(fabs(k-mk)<1){u[n]=1/(v[n]+(k-mk+1)*(vh-v[n])/2);}}}        //Blend crustal v with mantle v at the Moho, write to u.
        ///P wave arrival forward model: Iterate through each event; propagate first arrivals though a Cartesian grid using FMM(); interpolate for station times.///
        auto Q=[&](int n){int g,i,j,k,l; FMM E(di,dj,dk,h,u); vec w,f=E.t,c,e,L,v,s;//Declare lambda function Q(int n) which captures entire local scope by reference, takes argument int n setting event index; define local vars.
            if(n<nl){e={B[4*n],B[4*n+1],fabs(B[4*n+2])};}                      //If event index corresponds to a local/floating source, extract source e={lat,lon,dep} from b.
            else{for(vec&q:a){if(rint(q[2])==n){e={q[6],q[7],q[8]}; break;}}}  //If event index corresponds to a fixed/teleseismic source, extract source e={lat,lon,dep} from first row in 'a' corresponding to ev index.
            c=g2c(e); if(inside(c)){E.Init(h*c[0],h*c[1],h*c[2],5*h,fmax(0,c[2]-4));} //fmax(0,c[2]-4)  //If source is inside grid, initialize FMM with local source to a 5*h initialized radius.
            else{for(g=0; g<di*dj*dk; g++){i=g%di; j=(g/di)%dj; k=g/(di*dj);   //Initialize standard earth teleseismic travel times at grid bounds, iterating through each index n, converting to axial coords i,j,k.
                if(k<2 || i<2||i>di-3 || j<2||j>dj-3){f[g]=ET.time(e,c2g({1.*i,1.*j,1.*k}));}}//If index g is on the bottom face or any of the side faces, set E.t[h] using SET.time() from source coordinates to g
                E.Init(f);}                                                    //Initialize FMM using pre-initialized teleseismic times at bounds.
            E.March(ftype); if(n<nl){f=E.t+B[4*n+3];} else{f=E.t+B[3*nl+n];}   //Execute FMM march, set model time vector t=FMM times+ corresponding type offset time in B[].
            for(int m=0; m<da; m++){w=a[m]; c={w[3],w[4],-w[5]}; l=rint(w[2]); //Iterate through each row/data point in 'a', extract as vec v. extract data point's associated event count 'l' in v[2].
                if(l==n){t[m]=interp(f,g2c(c),{di,dj,dk},{1});}}};             //If thread event 'n' corresponds to data point event 'l', write interpolated model travel time to model time vector mt[].
        for(int n=0;n<nl+nt;n++){R[n]=thread(Q,n);}                            //Iterate through each event index 'n' in separate threads. exit FM, end class FM. To debug use: for(int n=0;n<nl+nt;n++){Q(n);}}};
        ///Gravity forward model.///
        if(S.grav){int fh=1, fv=2, ii=(di-1)*fh+1, jj=(dj-1)*fh+1, kk=(dk-1)*fv+1, N=ii*jj*kk, bfr=S.bfr,l; drho=vec(N); ivec nn; vec surf=topo+sea; //Set up moho perturbation density field drho parameters.
            double vm=B[4*nl+nt+1], ro=-B[4*nl+nt+2], dr=B[4*nl+nt+3], g0=B[4*nl+nt+4]; //Extract mantle v, offset density ro, velocity-density scaling ratio dr, and gravity offset g0
            for(int n=0;n<N;n++){double i=(n%ii)/(1.*fh), j=((n/ii)%jj)/(1.*fh), k=(n/(ii*jj))/(1.*fv), li=1./fv, m0,m1,r0,r1,rp; //Iterate through partitioned velocity grid, partitioning by factor fh horizontally, fv vertically.
                m0=interp(sea,i,j,di,dj,1)-3.5; m1=interp(M,i,j,di,dj,1); r0=0; r1=dr*(interp(v,i,j,k,di,dj,dk,1)-vm); //Use linear interpolation to set grid point moho depth, terrain altitude, and velocity scaled density
                if(k>m0+li){r0=ro;} else if(fabs(k-m0)<li){r0=r0+(k-m0+li)*(r0+ro)/(2*li);} //If k>iasp91 moho depth+ h, set r0=ro
                if(k>m1+li){r1=ro;} else if(fabs(k-m1)<li){r1=r1+(k-m1+li)*(r1+ro)/(2*li);} //If k>rendered moho depth+ h, set r0=ro
                if(fabs(r0-r1)>.001 && c2g({i,j,k})[2]>10){drho[n]=r1-r0; nn.push_back(n);}} savetxt(drho,dir+"drho.dat"); //If density perturbation is greater than .001 g/cc, add perturbed grid point index to set nn. After finishing iteration, save rendered density perturbation field.
            if(GPU){vector<float>R32, S32, G32(di*dj); FILE *fR; FILE *fS; FILE *fG; ///For GPU rendered gravity field, declare file pointers and call compiled CUDA executable GFM_GPU from shell.///
                for(auto d:drho){R32.push_back(d);} fR=fopen("/dev/shm/R.dat","w"); fwrite(R32.data(),4,R32.size(),fR); fclose(fR); //Convert to float32 and write drho to binary as R.dat to shared system memory.
                for(auto d:surf){S32.push_back(d);} fS=fopen("/dev/shm/S.dat","w"); fwrite(S32.data(),4,S32.size(),fS); fclose(fS); //Convert to float32 and write G to binary as S.dat to shared system memory.
                FILE *gpu=popen("./GFM","r"); pclose(gpu);                      //Execute compiled GPU gravity field model in shell.
                fG=fopen("/dev/shm/G.dat","r"); l=fread(G32.data(),4,di*dj,fG); fclose(fG); for(int m=0; m<di*dj; m++){G[m]=G32[m];}} //Open file containing float32 binary of rendered gravity field, read data, populate float42 (double) vector G with imported data.
            else{thread GR[dj-2*bfr];                                          ///For CPU rendered gravity field, initialize thread array GR and use lambda function GFM()///
                auto GFM=[&](int j0){for(int i0=bfr; i0<di-bfr; i0+=1){int m=i0+j0*di; double k0=surf[m]; G[m]=0;    //Define multi threaded lambda function DQ to compute surface gravity anomaly due to perturbed density field.
                    for(int n:nn){double i=(n%ii)/(1.*fh), j=((n/ii)%jj)/(1.*fh), k=(n/(ii*jj))/(1.*fv);             //Iterate through fine grid with i, j axes further partitioned by factor fh, k axis by factor fv.
                    G[m]+=(66.7*drho[n]*(k0-k)/pow((i0-i)*(i0-i)+(j0-j)*(j0-j)+(k0-k)*(k0-k),1.5))/(fh*fh*fv);}}};   //Add gravity perturbation due to every grid block using dg=G*rho*vol/dist^2 *delta_k/dist integrated over entire perturbation field
                for(int j=0;j<dj-2*bfr;j++){GR[j]=thread(GFM,j+bfr);} for(auto&w:GR){w.join();}}}//Launch dj-10 threads to compute perturbation field inside grid, with 10 point buffer. struct timespec aa,bb; clock_gettime(CLOCK_MONOTONIC,&aa); clock_gettime(CLOCK_MONOTONIC,&bb); print(bb.tv_sec-aa.tv_sec+(bb.tv_nsec-aa.tv_nsec)/1e+9," gravity timer");
        for(auto&w:R){w.join();}                                               //Join seismic forward model threads, allowing gravity to execute concurrently.
        y=col(a,0)/S.tsig; f=t/S.tsig;                                         //Normalize data and model time to estimated uncertainty.
        if(S.grav){for(int i=S.bfr;i<di-S.bfr;i++){for(int j=S.bfr;j<dj-S.bfr;j++){int nn=i+j*di; append(y,GD[nn]/S.gsig); append(f,(G[nn]+g0)/S.gsig);}}} //If gravity selected, append normalized gravity data to y, normalized gravity model data to f.
        res=y-f; ssq=dot(res);}};                                              //Compute sum of squares of residuals ssq.

///Output function for saving variables Format: void output(args S, vec b, int ind); args S: static args(tele event, station, grid data). vec b: dynamic invertible args(local evt coords, offset ts, Voronoi nodes).///
void out(args &S, vec b, int n, mat J=mat{{-1}}){
    double di=S.di,dj=S.dj,dk=S.dk,nl=S.nl,nt=S.nt,no2=S.no2,no3=S.no3,m,mm; vec v,st,w; mat a,src,sta; string I; int l;//Declare and extract parameters from S.
    I=to_string(n)+".txt"; FM F(S,b); col(S.a,19,F.t); a=S.a;                  //Set output path p, convert passed index integer to string I
    if(n>=40){for(vec &w:a){w=vec(w.begin(),w.begin()+20);}}                   //If iteration count is >100, cut down data array a to 20 columns to save space.
    m=a.size(); v={S.r[0],S.r[1],di,dj,dk,S.h,nl,nt,no2,no3,m,F.ssq,S.tsig[0],S.gsig,(double)((int)S.grav)}; //Define inversion header v, to be attached to top of S.da.
    for(l=0;l<4*nl;l+=4){append(src,g2c({b[l],b[l+1],fabs(b[l+2])}));} savetxt(src,dir+"src_"+I); //Generate local event hypocenter grid coordinate array src, save as .txt.
    for(l=0;l<m;l++){st=g2c({a[l][3],a[l][4],-a[l][5]}); mm=1; for(vec &w:sta){mm=fmin(mm,dot(w-st));} if(mm>.001){append(sta,st);}} savetxt(sta,dir+"sta_"+I);
    append(v,vec(a[0].size()-v.size())); append(v,a); savetxt(a,dir+"a_"+I);   //Extender header v to width of a, stack on top of station data S.a as output matrix a, save as a_ind.txt.
    savetxt(b,dir+"b_"+I); savetxt(S.sig,dir+"sig_"+I);                        //Save parameter vector b as b_#.txt and b estimated uncertainties sig.
    savetxt(J,dir+"J_"+I); savetxt(S.D,dir+"dmp_"+I);                          //save Jacobian as J_#.txt. Default J is {{0}}, damping vector used for inversion D as dmp.
    savetxt(F.u,dir+"u_"+I); savetxt(arr(F.o3,F.o3.size()/4,4),dir+"o3_"+I);   //Save rendered p velocity as a flattened 1d vector where: field value at [i,j,k]=u[i+j*di+k*di*dj]. Save Voronoi node arrays o3, o2 as no2 x 3  and no3 x 4 arrays with row format: [x coord, y coord, z coord, amplitude]. savetxt(arr(F.o2,no2,3),dir+"o2_"+I);
    savetxt(arr(F.o2,F.o2.size()/3,3),dir+"o2_"+I);                            //Save moho Voronoi node field
    savetxt(arr(F.M,S.di,S.dj),dir+"m_"+I);                                    //Save rendered moho depth mapping as a di x dj array with row format.
    savetxt(arr((F.G-GD),S.di,S.dj),dir+"dG_"+I); savetxt(F.drho,dir+"drho_"+I); //If gravity selected, save gravity anomaly field, density perturbation field.
    savetxt(arr(sea,S.di,S.dj),dir+"sea.txt"); savetxt(arr(topo,S.di,S.dj),dir+"topo.txt"); //Save sea level k map, topography k map.
    print(vec{(double)n,m,nl,nt,no2,no3,F.ssq},"\nout [ind,pts,nl,nt,no2,no3,ssq]");}  //Execute forward model FM using parsed args S, vec b. Compute ssq vs data 'q', print index and ssq.

///Jac(): Efficient Jacobian Function for TFM(). Identify regions corresponding to zero parameter covariance to avoid calculating null space elements.///
mat Jac(args &S, vec b, vec s){int da,df,n,nl,nt,j,e; vec d,r,f; mat J,a;        //Declare function, variables.
    f=FM(S,b).f; a=S.a; da=a.size(); df=f.size(); nl=S.nl; nt=S.nt; j=b.size(); J=arr(j,df);//Extract data 'a', reference coords 'r', E/N/U axis dims di/dj/dk, grid spacing 'h'(km), FMM init radius 's0'(km).
    s={}; for(int n=0; n<nl;n++){append(s,{.01,.01,2,1});}  for(int n=0; n<nt;n++){append(s,1.);}                     ///IMPORTANT: FIXED STEP OPTION FOR JACOBIAN
    if(S.no2>0){append(s,.05); if(S.grav){append(s,{.05,.02,.01,1});} for(int n=0; n<S.no2;n++){append(s,{1,1,.5});}} ///IMPORTANT: FIXED STEP OPTION FOR JACOBIAN
    for(int n=0; n<S.no3;n++){append(s,{1,1,1,.05});}                                                                 ///IMPORTANT: FIXED STEP OPTION FOR JACOBIAN
    for(n=0; n<3; n++){                                                        //Iterate through parameter steps in each local location spatial coordinate (x,y,z).
        r=b; for(e=0; e<nl; e++){r[4*e+n]+=s[4*e+n];} d=FM(S,r).f-f;           //Execute 3 local event location parameter steps in [lat,lon,depth]
        for(int m=0; m<da; m++){e=rint(a[m][2]);                               //Iterate through each data point.
            if(e<nl){J[4*e+n][m]=d[m]/s[4*e+n]; J[4*e+3][m]=1;}                //If source is local, set J element using location step differences in D. Offset time Jacobian is always linear(1) for all events.
            else{J[3*nl+e][m]=1;}}}                                            //Set linear offset time jacobian elements for teleseismic events.
    for(int n=4*nl+nt; n<j; n++){if(s[n]>1e-6){r=b; r[n]+=s[n]; J[n]=(FM(S,r).f-f)/s[n];}}//If uncertainty is reasonably larger than zero (1e-6), Compute Jacobian column corresponding to location parameters (dense).
    return T(J);}                                                              //Return Jacobian (transpose all elements to return true jacobian J[m][n] where m is data index(row) and n is parameter index(column).

///Creation function for o3 nodes. Compute Jacobian, add random test nodes, assigning interpolated or nearest neighbor scalar field and test estimated uncertainty of new node amplitude appending existing J. Accept best trial.///
void create(args &S, vec &B, int nn, int ntr, mat J={{-1}}){ print("\ncreate");
    int i,n,n1,n0,da; double e,d,dn,r,ds,q; mat j; vec f,fr,b,br,b0,v,y,u,E;   //Declare variables
    n1=B.size(); n0=n1-4*S.no3; da=S.a.size(); y=col(S.a,0); ds=.2;            //Extract b index of first o3 node "n0", last node "n1", a[] row count da, pick data y, Set node amplitude step size 'ds', test args instance SS.
    if(J.size()<=1){J=Jac(S,B,ds*S.sig);}                                      //If a Jacobian is not passed in to create(), generate Jacobian.
    for(int ni=0; ni<nn; ni++){ n=B.size()+4; e=1e+9; u=FM(S,B).u;             //Iterate through 'nn' sets of node creation trials. Initialize node scalar field error and benchmark error E,E0 at 3km/s.
        append(S.sig,{4,4,4,1}); S.no3++; append(B,{0,0,0,1});                 //Generate trial args struct copy SS, add new node param uncertainties to sig, increment node count no3, copy vec b as bb.
        J=T(J); for(i=0;i<4;i++){append(J,vec(da));} J=T(J);                   //Add 4 new columns to jacobian J, one for each new node parameter.
        for(int t=0; t<ntr; t++){b=B; j=J;                                     //Iterate through "trials" trials.
            b[n-4]=randint(S.di); b[n-3]=randint(S.dj); b[n-2]=randint(S.dk);  //Assign random i/j/k for trial node.
            b[n-1]=1/u[b[n-4]+b[n-3]*S.di+b[n-2]*S.di*S.dj];                   //Set node u equal to local grid point u, add node to test beta bt.
            f=FM(S,b).t; q=dot(y-f);                                           //Evaluate first arrival time grid 'f' at current "b". Calculate SSQ as 'q'.
            for(i=n-4; i<n; i++){br=b; r=ds*S.sig[i]; br[i]+=r; col(j,i,(FM(S,br).t-f)/r);}  //Update Jacobian to correspond to new test node by calculating a step vector and finite differences with respect to each test node parameter.
            d=0; do{E=diag(inv(MTM(j)+I(n,d)))*q/(da-n); d=10*(d+1e-3);} while(!(dot(E)>1e-8));  //take inverse of J'J with increasing damping factor until a non singular inverse is found. Et[i]=SSQ/(da-len(b))*inv(jj')[i][i].
            if(E[i-1]<e){B=b; J=j; e=E[i-1];}} print(vec{B[n-4],B[n-3],B[n-2],B[n-1],e},"new node: [i,j,k,v,E]"); }                                  // print(vec{B[n-4],B[n-3],B[n-2],B[n-1],e},"new node: [i,j,k,v,E]");
            out(S,B,ctr++);}                                                   //export results.

 ///LMA(): Levenberg Marquardt Algorithm.
class LMA{public: double ssq,q,qt,mse; vec y,Q,r; mat J,JJ,M,A; int n,nd,nb;   //Declare class LMA, class variables ssq, ssqt, f,ft,SSQ,db,bt,r,J,B.
void fit(args &S, vec &b, mat dm, double ds=.1, int II=4, double g=.999){      //Declare sole LMA method fit().
    q=FM(S,b).ssq; Q={q}; nd=dm.size(); nb=nd-1; print("\nlma");               //Initialize model values f at starting params b, their ssq=dot(F(S,b)-y,F(S,b)-y), ssq log SSQ, b log B.
    for(int ii=0; ii<II; ii++){mat B=arr(nd,0); thread tr[nd]; r=FM(S,b).res;  //Iterate for 'II' inversions. Initialize damping condition grid search multithreading parameters.
        J=Jac(S,b,ds*S.sig); JJ=MTM(J); M=diag(diag(JJ)); A=diag(1/(S.sig^2)); //Compute Jacobian at start b with step size=frac of each b uncertainty 'sig'/'s', f=current model times.
        auto F=[&](int n,vec d){B[n]=b+inv(JJ+d[0]*M+d[1]*A)*T(J)*r;};         //Inline damped Gauss-Newton step calculation function.
        for(n=0;n<nd;n++){tr[n]=thread(F,n,dm[n]);} for(auto&t:tr){t.join();}  //Launch and then join 'nd' parallel threads to test all damping coefficient conditions in dm.
        if(S.no2>0){for(vec&v:B){int m, db=b.size(), o3o=db-4*S.no3, o2o=o3o-3*S.no2; //If 2 layer velocity paarametrization is used, iterate through every test vector v in B to enforce parameter limits; define indices of first o2, o3 elements, last element.
            for(m=o2o;m<o3o;m+=3){v[m+2]=fmin(fmax(v[m+2],S.dk-11),S.dk-4);}   //Enforce moho depth limits.
            for(m=o3o;m<db;m+=4){v[m+3]=fmin(fmax(v[m+3],7.5),9.5);} }}        //Enforce mantle node velocity limits
        for(n=0;n<nd;n++){qt=FM(S,B[n]).ssq; if(qt<q){q=qt; b=B[n]; nb=n;}}    //Iterate through every damping condition tested. If SSQT[]<ssq: set ssq=SSQT[], b=tested b BT[], f=tested f(FT[]); alt: ///for(vec v:B){qt=dot(y-FM(S,v).t); if(qt<q){q=qt; b=v;}}  print(mat{dm[n],{q,qt}},"damping, qt, t");
        S.D=diag(dm[nb][0]*M+dm[nb][1]*A); print(dm[nb],"best damping");       //Compute the damping matrix diagonal from the damping conditions used in iteration.
        Q.push_back(q); out(S,b,ctr++,J); if(Q[ii+1]/Q[ii]>g){break;}}}};      //Add best ssq to SSQ. If ratio of current ssq/previous ssq >g, fit is no longer improving; break iteration loop.

///Annihilation function for o3 nodes. Iteratively find node in o3 with highest uncertainty of scalar field amplitude, delete node, update J.///
void annihilate(args &S, vec &b, mat &J, int rmv){print("\nannihilate");
    for(int ind=rmv; ind>0; ind--){ if(S.no3==0){print("no nodes"); return;}   //Iterate through "rmv" loops to remove one node per loop.
        int db=b.size(), da=S.a.size(), m=db-4; double d=0; vec V;             //Extract side of b, rows in a, new b length m (after each 4 parameter node is removed), damping factor d, scaled error vector V. Generate identity damping matrix.
        do{V=diag(inv(MTM(J)+I(db,d))); d=10*(d+1e-6);} while(!(dot(V)>1e-8)); //take inverse of J'J with increasing damping factor until a non singular inverse is found.
        for(int n=db-4*S.no3; n<db; n+=4){if(V[n+3]>V[m+3]){m=n;}}             //Iterate through o3 parameters in b corresponding to amplitude, to find the one with the largest relative uncertainty E[m].
        print(vec(b.begin()+m,b.begin()+m+4),"worst {i,j,k,v}");               //Print parameters of worst resolved o3 node before annihilating it.
        b.erase(b.begin()+m,b.begin()+m+4); S.no3--;                           //Delete node m with largest uncertainty from b, decrement node count no3.
        S.sig.erase(S.sig.begin()+m,S.sig.begin()+m+4);                        //Delete apriori uncertainty estimates from S.sig corresponding to node m.
        J=T(J); J.erase(J.begin()+m,J.begin()+m+4); J=T(J);} out(S,b,ctr++);}  //Delete Jacobian columns corresponding to node m parameters.

///trim: Goodness of fit maximization through data reduction:
void trim(struct args &S, vec b, double R){vec r=fabs(FM(S,b).t-T(S.a)[0]); print("\ntrim"); //Compute absolute residual vector r using current S and b.
    for(int i=S.a.size()-1;i>=0;i--){if(r[i]>R){S.a.erase(S.a.begin()+i);}}}   //Iterate through each row in a, backwards. If absolute residual/misfit exceeds R, elimiate row from a.

///Tomographic Inversion Algorithm. Build+run inversion command: "g++ -o TFM TFM.cpp -lpthread -O3 && ./TFM"  ///
int main(){
        int runs=20; DRE.seed(seed); for(int run=0; run<runs;run++){           //Set number of monte carlo runs in "runs", seed random number generator, iterate for runs specified.
        int di,dj,dk,II,wa,nlm,ntr,n,nl,nt,evs,no2,no3,l,bfr; double h,trm,i,j,k,nc,mm,q,gsig,stsig,airy,vmantle,vcrust,drhomc,drhov; string root,name; mat a,e,dm,BG,dmp; LMA lma; vec t,b,sig,c,s,y,nn,r,o3sig,o2sig,srcsig; args S; vector<string>fl; FILE *F; bool grav,moho; //Declare variables.
        name="v66MC"; dir=mkdir("output/"+name+"/");                           //create new directory 'dir' in root to read data and populate with inversion results (optionally pass str with dir name). Return full path 'dir'.
        ftype=1;                                                               //Set seismic forward model type (described in FMM).
        di=61; dj=61; dk=33; h=10;                                             //Set grid dimensions [di,dj,dk]=[E,N,U], grid spacing h(km);
        evs=132;                                                               //Set total number of event files to import.
        r={-15.25,-72.5};                                                      //Set central reference {lat,lon} coordinate vector 'r',
        wa=2000;                                                               //Set inputdata array column rank (can be 20+).
        airy=5; vmantle=8, vcrust=6;                                           //Set airy compensation ratio for Moho initialization, initialization/background mantle v, crustal v,
        stsig=.3; srcsig={.4,.4,40,100}; o2sig={2,2,2}; o3sig={2,2,2,.5};      //Set apriori uncertainties of: arrival pick times "stsig", local source coordinates 'srcsig', moho mapping Voronoi nodes 'o2', Mantle velocity mapping Voronoi nodes 'o3'.
        moho=true;                                                             //Specify whether to use 2 layer model with Moho, or perturbed halfspace.
        grav=true; GPU=true;                                                   //Specify whether to use gravity data and inversion as "grav", whether attempt to use GPU gravity modeling "GPU".
        gsig=20; bfr=10;                                                       //Set gravity data uncertainty values in milligal as 'gsig', gravity grid boundary buffer "bfr" at which to set null perturbation;
        drhomc=.4; drhov=.3;                                                   //Set starting mantle-crust density contrast "drhomc", perturbed velocity scaling factor for density perturbation "drhov"
        dmp={{0,.3,1,3,10,100},{10000.,100000.}};                              //Set inversion damping coefficient array dmp.
        trm=1.2; II=1; nlm=5; ntr=20;                                          //Set outlier pick trim threshold 'trm', number of meta inversion iterations to run "II", number of damped gauss newton inversions per II cycle "nlm", number of MCMC creation trials per node ntr.
        nn={1,15};                                                             //Set number of new nodes and destroyed nodes per II cycle as nn={created,annihilated}.
        if(runs>1){airy=rand(4.,6); vmantle=rand(7.7,8.3); vcrust=rand(5.9,6.2); drhomc=rand(.3,.45); drhov=rand(.1,.4); gsig=rand(10,40); //If multiple runs used, set ranges for random starting parameters.
            print(vec{airy,vmantle,vcrust,drhomc,drhov,gsig},"airy,vmantle,vcrust,drhomc,drhov,gsig");} //Print randomly generated starting parameters

        ctr=1; nc=-1; nl=nt=no2=no3=0;                                         //Initialize counters
        if(GPU && grav){F=popen("nvcc -o GFM GFM_GPU.cu -O3  2>&1","r"); string gs; char ch[64]; if(F){while(!feof(F)){if(fgets(ch,64,F)!=NULL){gs.append(ch);}}} pclose(F);  //If GPU option is selected, compile GPU gravity forward model in local file GFM_GPU.cu using nvcc.
            if(string::npos==gs.find("error")){print("\nGPU gravity model compiled");} else{GPU=false; print("GPU error:\n\n"+gs+"\nUsing CPU gravity model");}}
        for(auto dlm:dmp[0]){for(auto dba:dmp[1]){append(dm,{dlm,dba});}}       //Iterate through 2 sets of damping coefficients in to create damping matrix dm.
        for(int n=0;n<evs;n++){fl.push_back("input/E"+to_string(n)+".txt");}    //Populate string vector fl with input file paths. Data name format: "EN.txt". Iterate through evs indices.
        RC={r[0],r[1],(dk-3)*h/2-5,(di-1)*h/2,(dj-1)*h/2,(dk-1)*h/2,h,di*1.,dj*1.,dk*1.};  //Set offset/central coordinate vector RC.
        ///Event data import and format.///
        for(auto E:fl){e=loadtxt(E,wa); s=e[0];                                ///Local: Extract event data e from each file in "files" as 20 column matrix. Compute ENU coordinates of the event as c={E,N,U}={x,y,z}.
            if(inside(g2c({s[6],s[7],s[8]}))){nl++; nc++;                      //If the ENU coordinates of event source are within the grid limits, add event data as a local floating source.
                append(b,{s[6],s[7],s[8],s[9]}); append(sig,srcsig);           //Append pick file source coordinates in to parameter vec b, and to uncertainty vector sig {.02",.02",2km,1s} .
                for(vec p:e){if(inside(g2c({p[3],p[4],p[5]}),1)){p[2]=nc; append(a,p);}}}}//Iterate through each station data point "vec p" as rows of mat e; Transform station geodetic coordinates to ENU coordinates c. If station inside grid, append local evt station data point to 'a' as a row.
        for(auto E:fl){e=loadtxt(E,wa); s=e[0];                                ///Teleseismic: Extract event data e from each file in "files" as 20 column matrix. Compute ENU coordinates of the event as c={E,N,U}={x,y,z}.
            if(gca(s[6],s[7],r[0],r[1])>20){nt++; nc++;                        //If the event source is >10 degrees from reference, add event data as a fixed teleseismic source. (previous cutoff was 7 degrees)
                append(b,s[9]); append(sig,srcsig[3]);                         //Insert pick file source time in to parameter vec b, insert 1s uncertainty in to sig.
                for(vec p:e){if(inside(g2c({p[3],p[4],p[5]}),1)){p[2]=nc; append(a,p);}}}}//Iterate through each station data point "vec p" as rows of mat e; Transform station geodetic coordinates to ENU coordinates c. If station inside grid, append local evt station data point to 'a' as a row.
        ///Import global Digital Elevation Map and Bougher anomaly map. Use either isostatic equilibrium or bougher anomaly to generate initial moho M0. ///
        sea=topo=GD=M0=vec(di*dj); BG=loadtxt("PGM2012B.txt",3);               //Load cropped WGM2012 Bouguer anomaly data as a 3 column matrix as [lat,lon,dg]
        vector<int16_t>H(9331200); F=fopen("ETOPO5.bin","r"); l=fread(H.data(),2,H.size(),F); fclose(F);//Import binary digital elevation map as H. File format: int 16, global map, 5 min spacing starting from 90N 0E (dim= 12*12*180*360).
        for(int n=0;n<di*dj;n++){c=c2g({n%di+0.,n/di+0.,dk-3.}); c[2]=0; sea[n]=g2c(c)[2];//iterate through each i,j, coordinate, transform semiterminal top sheet i,j,k to geodetic coordinates.
            topo[n]=H[(int)((90-c[0])*12.)*360*12+(int)(std::fmod((360+c[1]),360.)*12.)]/(1000*h); topo[n]=fmax(topo[n],-.5); M0[n]=sea[n]-35/h-airy*topo[n]; //Use local elevation to calculate isostatic moho depth using crustal thickening of 5km:1km altitude. //savetxt(moho_init,dir+"moho_init.txt");
            q=mm=1e-15; for(vec&v:BG){if(pow(c[0]-v[1],2)+pow(c[1]-v[0],2)<.04){q+=v[2]; mm++;}} GD[n]=q/mm;} //average bouguer anomalies from all defined points within .2 degrees summing them as q, with count mm.
        for(int i=0; i<di; i++){for(int j=0; j<dj; j++){if(i<bfr || i>=di-bfr || j<bfr || j>=dj-bfr ){GD[i+j*di]=0;}}} //Zero gravity data field everywhere within bfr points of boundaries.
        if(moho){append(b,vcrust); append(sig,.2);                             //Add Seismic model crustal velocity,
            if(grav){append(b,{vmantle,drhomc,drhov,0}); append(sig,{.1,.1,.05,100});}//If gravity inversion selected, add velocity:density contrast factor, and gravity background velocity parameters for generating drho to b, corresponding uncertainties to sig.
            ///Generate o2/2d moho depth parameters, add to 'b'.///
            //l=10; for(i=l;i<di-l;i+=40){for(j=l;j<dj-l;j+=40){append(b,{i,j,M0[(int)(.1+i+j*di)]}); append(sig,{1e-7,1e-7,1e-7}); no2++;}} //Add grid of o2 Voronoi nodes to b (frozen moho corner nodes set at isostatic equilibrium depth.)
            for(vec w:a){s=g2c({w[3],w[4],-w[5]}); mm=100; q=4; s=vec{s[0],s[1]}+rand(-q,q,2); s[0]=fmin(fmax(s[0],2),di-3); s[1]=fmin(fmax(s[1],2),dj-3); //Create test o3 Voronoi nodes to b by taking random perturbations of each initial event location.
                for(l=b.size()-3*no2; l<b.size(); l+=3){mm=fmin(mm,dot(vec{b[l],b[l+1]}-vec{s[0],s[1]}));} //Find distance to nearest accepted o3 node.
                if(mm>16){append(b,{s[0],s[1],M0[(int)s[0]+di*(int)s[1]]}); append(sig,o2sig); no2++;}} } //Accept node if its more than dmin points from nearest o3 node.
        else{for(vec w:a){s=g2c({w[3],w[4],-w[5]}); mm=100; q=4; s=s+rand(-q,q,3); s[0]=fmin(fmax(s[0],2),di-3); s[1]=fmin(fmax(s[1],2),dj-3); s[2]=fmin(fmax(s[2],2),dk-3); //If single space model is used, generate conditional velocity nodes near stations.
                for(l=b.size()-4*no3; l<b.size(); l+=4){mm=fmin(mm,dot(vec{b[l],b[l+1],b[l+2]}-s));} //Iterate through generated node set and find closest accepted node with respect to conditional node. Record nearest neighbor distance^2 as mm.
                if(mm>16){append(b,{s[0],s[1],s[2],vmantle}); append(sig,o3sig); no3++;}}} //If nearest neighbor is more than 4 grid points away, accept conditional node, add to b.
        ///Generate o3/3d velocity field parameters, add to 'b'.///
        l=12; for(i=l;i<di-l;i+=l){for(j=l;j<dj-l;j+=l){for(double k:{7,14,20,25}){append(b,{i,j,k,vmantle}); append(sig,o3sig); no3++;}}}//Add grid of o3 Voronoi nodes to b, their uncertainties to sig, update S.no3. {6,13,19,24,27}
        for(n=0; n<4*nl; n+=4){s=g2c({b[n],b[n+1],b[n+2]}); mm=100; q=4; s=s+rand(-q,q,3); s[0]=fmin(fmax(s[0],2),di-3); s[1]=fmin(fmax(s[1],2),dj-3); s[2]=fmin(fmax(s[2],2),dk-6); //Create test o3 Voronoi nodes to b by taking random perturbations of each initial event location.
            for(l=b.size()-4*no3; l<b.size(); l+=4){mm=fmin(mm,dot(vec{b[l],b[l+1],b[l+2]}-s));} //Iterate through generated node set and find closest accepted node with respect to conditional node. Record nearest neighbor distance^2 as mm.
            if(mm>16){append(b,{s[0],s[1],s[2],vmantle}); append(sig,o3sig); no3++;}}  //If nearest neighbor is more than 4 grid points away, accept conditional node, add to b.
        S={a,r,di,dj,dk,h,nl,nt,no2,no3,sig,stsig/col(a,1),b*0,grav,gsig,bfr}; //Pack TFM arguments in to struct S
        out(S,b,0);                                                            // output initial state.
        ///Run alternating series of LMA and MCMC toggling between Levenberg Marquardt gradient based inversion and Markov Chain Monte Carlo reassignment of worst fitting nodes.///
        for(int I=1;I<=II;I++){                                                //Loop creation/migration/annihilation cycle for a total II iterations.
            lma.fit(S,b,dm,.25,nlm); out(S,b,I+1000);                          //Execute a multi iteration LMA gradient based inversion, export results.
            trim(S,b,trm); annihilate(S,b,lma.J,nn[1]); lma.fit(S,b,dm,.25,nlm/2);  //Trim data set, annihilate nn[1] of the highest error nodes, execute 1 iteration of LMA, export results.
            for(l=0;l<2;l++){create(S,b,nn[0],ntr,lma.J); lma.fit(S,b,dm,.25,nlm/2);}}}}//Create new nodes, execute 1 iteration of LMA. Loop 3 times.

///main() synthetic time generation///
//args SS={a,r,di,dj,dk,h,nl,nt,no2,0,sig}; vec sb=1*b;                      //SYNTHETIC Define synthetic copies of args SS and parameters sb
//for(int n=0; n<SS.no2; n++){ i=randint(di); j=randint(dj);                 //SYNTHETIC Iterate through no2 nodes, randomly assign i,j.
//  double m=1, c=0, i0=(i+m*j-m*c)/(1+m*m), j0=m*i0+c, d;                   //SYNTHETIC Define synthetic perturbation line slope m and intercept c.
//	d=35*exp(-(pow(i-i0,2)+pow(j-j0,2))/400); append(sb,{i,j,d});}           //SYNTHETIC Set synthetic node depth as gaussian of distance to ridge line.
//output(SS,sb,111); col(a,0,TFM(SS,sb)+nrand(0,.15,a.size()));              //SYNTHETIC Generate synthetic data set, set pick times as noisy synthetic model time set.

///error vector err
//vec err=sqrt(q/(S.a.size()-b.size())*diag(inv(T(J)*J)));                   //Compute error vector err.
