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
