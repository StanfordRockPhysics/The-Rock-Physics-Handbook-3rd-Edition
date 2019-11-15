%RPHtools 
%Matlab files for The Rock Physics Handbook 2nd edition
%Mavko, Mukerji, and Dvorkin
%Cambridge University Press, 2009
%
% Anisotropy
%   bkusc      - Backus average for thin layered TI anisotropy
%   bkuslog    - Computes Backus Average from log data
%   CSiso      - 6x6 stiffness and compliance tensors for isotropic rock
%   c2anis     - Anisotropy parameters (Thomsen) from stiffnesses Cij
%   c2sti      - Convert between stiffness and compliance for TI symmetry
%   c2vti      - Stiffness Cij to velocity for TI media
%   cti2v      - Stiffness Cij to fast/slow velocities and Thomsen parameters
%   ezbond     - Coordinate rotation for stiffness tensor using Bond transformation
%   Rorsym     - reflectivity in the symmetry plane for interfaces between
%                2 orthorhombic media
%   Rruger     - Reflectivity AVOZ in weakly anisotropic media 
%   v2cti      - Velocities to stiffnesses Cij for TI media
%
% AVO
%   avopp      - Reflectivity vs. angle, P-to-P for single interface
%   avops      - Reflectivity vs. angle, P-to-S for single interface
%   avo_abe    - AVO intercept and gradient attributes for P-P and P-S
%   eimp       - Far-offset elastic impedance for P-P and P-S (reflection angle)
%   eimp2      - Far-offset elastic impedance for P-P and P-S (incidence angle)
% 
%
% Effective Medium
%   berrysc    - Self-Consistent (Berryman) effective elastic moduli
%   berryscm   - Self-Consistent (Berryman) multi-component (> 2)
%                effective elastic moduli
%   berryscp   - Effective elastic moduli vs. pressure for multi-component
%                composite using Berryman's Self-Consistent method
%   bkus       - Backus average anisotropic stiffnesses for thin layers
%   bound      - Upper and Lower bounds on elastic moduli of aggregates
%   critpor    - Velocity, density, moduli at critical porosity
%   dem        - Differential Effective Medium elastic moduli
%   dem1       - Differential Effective Medium elastic moduli at a single
%                volume fraction of constituents
%   echeng     - Eshelby-Cheng TI stiffnesses Cijkl for cracked rock
%   hash       - Hashin-Shtrikman upper and lower bound plot for moduli
%   hashv      - Hashin-Shtrikman upper and lower bound plot for velocities
%   hudson     - Hudson anisotropic Cijkl for cracked rock
%   hudson3    - Hudson anisotropic Cijkl for cracked rock with 3 crack sets
%   hudsoncone - Hudson anisotropic Cijkl for conical crack distributions
%   hudsonF    - Hudson anisotropic Cijkl for Fisher crack distributions
%   
%
% Fluid Effects, Velocity Dispersion, & Attenuation
%   biot       - Velocity dispersion and attenuation from Biot theory
%   biothf     - High frequency limiting velocity from Biot theory
%   biothfb    - Approximate high frequency limiting velocity from Biot theory
%   bkti       - Brown & Korringa anisotropic (TI) saturated compliances
%   BKc2c      - Brown & Korringa (stiffness) arbritrary anisotropy and fluids
%   BKd2s      - Brown & Korringa arbritrary anisotropy dry to sat.
%   BKs2d      - Brown & Korringa arbritrary anisotropy sat. to dry
%   BKs2s      - Brown & Korringa (compliance) arbritrary anisotropy and fluids
%   co2prop    - CO2 properties versus temperature and pressure
%   flprop     - Reservoir fluid properties from Batzle-Wang relations
%   flpropui   - Reservoir fluid properties from Batzle-Wang relations
%              - with interactive input/output dialog boxes.
%   gassmnk    - Gassmann fluid substitution (bulk moduli)
%   gassmnv    - Gassmann fluid substitution (velocities)
%   mmti       - Unrelaxed wet frame anisotropic (TI) compliances from dry
%   patchw     - White's patchy model with Dutta-Ode correction
%   squirt     - Mavko squirt model for high frequency saturated velocities
%   stdlin     - Standard linear viscoelastic solid
%
%  
% Granular Media
%   Cem        - Dvorkin cementation model
%   hertzmind  - Bulk and shear moduli of sphere packs using Hertz-Mindlin model
%   hertzmindv - Vp, Vs of sphere packs using Hertz-Mindlin model
%   Johnson    - Norris-Johnson model for stress-induced anisotropy in 
%                sphere packs
%   John_Makse - Johnson-Makse model for stress-induced anisotropy in sphere packs
%   Unconsol   - Modified Hashin-Shtrikman lower bound with Hertz-Mindlin
%                for unconsolidated sediments
%   walton     - Bulk and shear moduli of sphere packs using Walton's model
%   waltonv    - Velocities of sphere packs using Walton's model
%
%
% Moduli and Velocities
%   ku2v       - Velocities from bulk and shear moduli and density
%   lm2v       - Velocities from Lame's constants and density
%   v2ku       - Bulk and shear moduli and Poisson's ratio from velocities
%   v2lm       - Lame's constants and Poisson's ratio from velocities 
%
% Permeability
%   PermMenu   - Menu tool to select different permeability models
%   BernabeE   - Bernabe permeability model
%   Bloch      - Bloch empirical permeability model
%   CoatDum    - Coates-Dumanoir permeability model
%   Coates     - Coates permeability model
%   FredrichE  - Fredrich et al. permeability model
%   KozCarmE   - Kozeny-Carman
%   ModKozCarm - Modified Kozeny-Carman with percolation porosity
%   Owolabi    - Owolabi et al. permeability model
%   PandaLake  - Panda-Lake permeability model
%   PandaLakeKCE - Panda-Lake Kozeny-Carman model
%   RevilE     - Revil et al. shaly permeability model
%   Timur      - Timur permeability model
%   Tixier     - Tixier permeability model
%   WylGregE   - Wyllie Gregory Kozeny-Carman model
%
% Seismic and Signal Processing
%   fftplot    - Amplitude and phase spectra plot
%   iatrib     - Instantaneous attributes (amplitude, phase, frequency)   
%
%
% Statistics
%   bayesclass - Bayes classification based on pdf
%   hist2d     - Bivariate histogram of data
%   monte      - Monte-Carlo draws from non-parametric marginal cdf 
%                                  (followed by linear regression)
%   monteccdf  - Monte-Carlo draws from non-parametric conditional cdfs 
%   pdfbayes   - Non-parametric pdf estimation, Bayes' error & Information
%
% Synthetic Seismograms and Wave Propagation
%   ezseis     - Quick normal incidence seismic section from velocity and
%                density
%   kenfdisp   - Velocity dispersion in 1-D layered media
%   kenfrtt    - Exact traveltimes in 1-D layered media for normally
%                incident plane waves
%   kennet     - Invariant imbedding wave propagation in 1-D layered media for
%                plane waves at normal incidence
%   pgator     - Propagator matrix wave propagation in 1-D layered media for
%                plane waves at normal incidence
%   sourcewvlt - Default source wavelet 
%
%   
% Utilities
%   blockav    - Block average of logs or signals
%   loadlas    - Load LAS formatted file
