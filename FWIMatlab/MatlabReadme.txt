FWIMatlab: Based on SimpleFWI code by TristanvanLeeuwen   (https://github.com/TristanvanLeeuwen/SimpleFWI)

This code performs Full Waveform Inversion in the frequency domain, solving the problem: 

 m_w = \min_{m} \sum_{i} \frac{1}{2}||R u_i - d_i||_2^2 + \frac{\alpha}{2}||Lm||_2^2, 

where 

 u_i=A(m)^{-1}q_i , 

which involves the sum over sources (s) and frequencies (omega) (i = s & omega),

A(m) is the discretized Helmholtz operator \omega^2 m + \nabla^2 with impedance boundary conditions, 

R is the sampling operator (linear interpolant), q_i are the source vectors, 

d_i are the observed data and L is the discretized \nabla.


Functions: 
Data.m:          Returns synthetic data d (computed on a finer grid to avoid inverse crime)
Misfit.m:        Input the medium parameters m and data d_i 
                 Outputs the misfit value & its gradient (& function handle to compute the action of the Gauss-Newton hessian).
BBiter.m:        Barzilai-Borwein optimisation algorithm  
FWIfull.m:       Calls above functions to compute optimal model mw
Continuation. m: Frequency Continuation FWI - input type, number of groups and overlap between groups 


How to run an example:

1). Run startup.m

2). Run Example1.m. 

Different examples may be produced by editing the following variables in Example1 :
   - Grid: Spacing & number of nodes can be edited by changing h & n respectively 
   - Subsurface: Model m edited by changing backgound wavespeed (v0) and obstacle wavespeed (dv) 
   - Frequency: Frequencies over which the misfit function is summed stored in f (more frequencies = better reconstruction)
   - Receiver locations: xr and zr (x and z coordinate)
   - Source locations: xs and zs (x and z coordinate)
   - FWI parameters: fwitol (FWI converges when ||gradient||_2 < fwitol) & maxit (maximum FWI iterations)
   - Regularisation parameter: alpha
   - Noise added to synthetic data by uncommenting 'Add Noise' section
   - Frequency Continuation: Comment out Original FWI and uncomment Frequency Continuation
                             -> type=1 for sequential, type=2 for progressive    