% run_PHT3D_1Dprofile_quickstart.m
% (starting point: run_PHT3D_1Dprofile_quickstart3_ForNextTime_test.m)
%
% 11/3/15
%
% - uses PHREEQC to first equilibrate ic w/ minerals and charge balance
%
% Sets up model by custom-creating the following input files:
%   1) *nam: namelist file 
%   2) *_ph.dat: PHREEQC interface package file (to link PHT3D to PHREEQC)
%   3) postfix: additional PHREEQC optiocecns 
% (the following files 4)-8) are MT3DMS input file formats)
%   4) *btn: basic transport package file
%   5) *ssm: source/sink mixing package file (not geochem reactiions, 
%      corresponds to stress packages from MODFLOW)
%   6) *dsp: dispersion package file
% (the following files 7)-8) are set up the same each time in this script)
%   7) *adv:advection package file
%   8) *gcg: GCG solver package file (iterative scheme to address stability
%      problems due to time step)
%
% _Orgcsource (3/31/16): adds another immobile component that is source
%   that continuously replenishes Orgcsed.  For steady-state:
% _dike (4/1/16): No Orgcsed in dike
%
% Assumes pht3d_databas.dat (geochem database) file already exists; specify
% in 'use_file_databas'.
%       
%
% Uses the following functions:
%   - SC_InitCond_chem_2()
%       Specifies observed conc, equilibrates with PHREEQC
%   - Alk2DIC()
%       Called by SC_InitCond_chem(), converts alkalinity to total C(4)
%   - generate_ic_PHREEQC_f_062215a()
%       Called by SC_InitCond_chem_2(), generates equilibrated ic using
%          PHREEQC
%   - go_PHREEQC2_4()
%       Called by generate_ic_PHREEQC_f_062215a(), runs PHREEQC
%   - PHT3D_ph_f_general4()
%       Creates pht3d_ph.dat
%   - PHT3D_btn_f_051715a()
%       Creates btn, dsp, ssm files using the specified inputs
%

clear all, close all; fclose all;

% ============================
% SECTION 1: INPUTS
% ============================

% 1A) DIRECTORY AND FILE NAMES

fl_gcng = 0;  % 1: Crystal, 0: Patrick

% ********** CUSTOMIZE FILE NAMES TO YOUR COMPUTER!! **********************
% Make sure to set the following: 
%   matlab_dir, PHT3D_exe, phrq_exe, sim_dir, flo_file, use_file_databas
if fl_gcng
    % - matlab_dir: directory with matlab functions 
    % matlab_dir = '.';  % use '.' if all your functions and scripts in the same directory
    matlab_dir = '/home/gcng/workspace/matlab_files/my_toolbox/PHT3D_functions';  

    % - choose one of these
    slashstr = '/'; % for Linux or Windows/Cygwin
    % slashstr = '\'; % for Windows/MS-DOS 

    % - Full path PHT3D executable
    PHT3D_exe = '/home/gcng/workspace/Models/PHT3D/src/pht3dv211';

    % - Full path PHREEQC executable
    phrq_exe = '/home/gcng/workspace/Models/PHREEQC/phreeqc-3.1.7-9213/bin/phreeqc';

    % - sim_dir: Directory with input files and where you run simulations
    sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/Minntac/test1/';

    % - flo_file: MODFLOW flo simulation file (full path)
%     flo_file = 'C:\Hydro_Modeling\MINNTAC_MATLAB_FILES\test2_1D_3\test.flo'; % 2D, no recharge
    flo_file = '/home/gcng/workspace/ModelRuns_scratch/MODFLOW_projects/Minntac/test3/test2_1D_3/test_MW12.flo';
    fl_rech = 1;  % 1: for recharge, 0 for no recharge in the MODFLOW simulation

    % - use_file_databas: geochem database, this is copied into 'pht3d_datab.dat' in sim_dir for simulation
%     use_file_databas = '/home/gcng/Documents/Teaching/ESCI5980_HydModeling/Fall2015/Lectures/L16_PHT3D_2/PHT3D_files/pht3d_datab.dat';
    use_file_databas = '/home/gcng/workspace/ProjectFiles/DNR_sulfate/PHT3D_files/database_files/pht3d_datab.dat_160331';
    % ************ (end of computer-specific file specifications) *************
else    
    % - matlab_dir: directory with matlab functions 
    % matlab_dir = '.';  % use '.' if all your functions and scripts in the same directory
    matlab_dir = 'C:\Hydro_Modeling\PHT3D_files\';  

    % - choose one of these
    slashstr = '/'; % for Linux or Windows/Cygwin
    % slashstr = '\'; % for Windows/MS-DOS 

    % - Full path PHT3D executable
    PHT3D_exe = '/cygdrive/c/Hydro_Modeling/pht3dv210/bin/pht3dv210.exe';

    % - Full path PHREEQC executable
    phrq_exe = 'C:\Hydro_Modeling\phreeqc';

    % - sim_dir: Directory with input files and where you run simulations
    sim_dir = 'C:\Hydro_Modeling\GW003_pht3d_dir\';

    % - flo_file: MODFLOW flo simulation file (full path)
    % flo_file = '/home/gcng/workspace/ModelRuns_scratch/MODFLOW_projects/ESCI5980/test2_1D_3/test.flo';
    flo_file = 'C:\Hydro_Modeling\MINNTAC_MATLAB_FILES\GW003_MODFLOW_dir\test.flo'; % 2D, no recharge
    % flo_file = '/home/gcng/workspace/ModelRuns_scratch/MODFLOW_projects/ESCI5980/Asst5_3/test.flo'; % 2D, w/ recharge
    fl_rech = 1;  % 1: for recharge, 0 for no recharge in the MODFLOW simulation

    % - use_file_databas: geochem database, this is copied into 'pht3d_datab.dat' in sim_dir for simulation
    use_file_databas = 'C:\Hydro_Modeling\pht3d_dir\pht3d_datab.dat';
    % ************ (end of computer-specific file specifications) *************
end

% - Name file will contain names of all input files, including MODFLOW flo file, 
%   Many files depend on resolution; 'suffix' is for names of resolution-
%   dependent files (btn, dsp, ssm)
nam_fil = 'pht3d_2D.nam';
suffix = '2D';  % suffix to file names


% 1B) (*DO NOT CHANGE*: SET UP A TIMER AND SOME OTHER FILE ADMINISTRATION)
a = clock;
if a(5) < 10
    fprintf('\n\nStart time: (%d/%d/%d) %d:0%d\n', a(2), a(3), a(1), a(4), a(5));
else
    fprintf('\n\nStart time: (%d/%d/%d) %d:%d\n', a(2), a(3), a(1), a(4), a(5));
end
addpath(matlab_dir);


% 1C) SET TIME PARAMETERS
% timprs = [0:30:1800]; % print out times [Start Time:Increment:End Time] [d]
% timprs = [0,1,365:365:365*6]; % print out times [d]
timprs = [0:365/10:365*12]; % print out times [d]
nstp = 20;  % number of time steps, incr for better numerical performance, decr for faster simulations


% 1D) PHYSICAL INPUTS AND PARAMETERS 

% -- Domain info 
nlay = 75;
ncol = 150;   
nrow = 1;
domain_bot_elev = -27.2; % m
domain_top_elev = 0; % top of domain must be at least this elev (include extra space for WT mov't)
domain_len = 224; % [m]

y_scale = 200/nlay; %ratio set by initial harcoded discretization of 200 rows by 400 columns
x_scale = 400/ncol; %ratio set by initial harcoded discretization of 200 rows by 400 columns


% -- General properties
tempC = 9.14; % water temperature, geochem rxs are temp-sensitive
por = 0.35;   % porosity



% -- Dispersion 

% - Calculate molecular diffusion coefficient

% Tstar: tortuosity 
% typical Tstar from Fitts (p. 530):
%   0.1 (clays) to 0.7 (sands) [Marsily 1986] 
%   0.56 to 0.8 (granular soils) [Bear 1972]
Tstar = 0.5;

% D: molecular diffusion coefficient
% Typical D from Fitts (Table 11.5, p. 531), increases with temp
%   D for SO42-: 1.1e-5 to 2.1e-5 cm2/s at 20degC [Li and Gregory, 1974]
D = 1e-5 / 1e4 * 3600 * 24;  % cm2/s -> m2/d

% Dstar: molecular diffusion coefficient in porous media
% Example Dstar values: Bemidji 3e-10 m2/s (2.6e-5 m2/d)
Dstar = Tstar*D;


% - Calculate mechanical dispersion coefficient

% alphaL: longitudinal dispersivity
% Typical alphaL: 
%   Bemidji: alphaL = 1 m
%   Cape Cod: alphaL = 0.96 m
%   Borden: alphaL = 0.36 m
%   Zheng textbook: Fig. 11.3 p. 303: alphaL vs. scale, mostly alphaL: 1e-2 
%       to 1e2 m, alphaL = 1e-2 m is for 1 m scale (this is similar to lab
%       column values)
long_disp = 1;  % m
hor2longdisp = 0.018; % ratio horiz transverse disp / long dispersivity (Garabedian: 0.018 m)
vert2longdisp = 0.0015; % ratio vertical transverse disp / long dispersivity (Garabedian: 0.0015 m)


% 1E) SET UP ORGANIC CARBON INPUTS, KINETIC MODEL THAT CONTROLS REDOX RXNS

% Set parameters for degradable organic carbon (Orgcsed)
% ** Modify to control kinetic rate of redox reactions
% Current model: 1st order decay relative to Orgcsed concentration; faster rate for aerobic degradation

% - set Orgsed ic conc [mol/g]
Orgcsed_conc_ic = 0.00066;  % [mol/Lw] * por / rho_b = [mol/g] * current value: arbitrary!!!

% - kinetic parameters (we're assuming 1st order decay)
% (increase logK for faster decay)
% Orgcsed_log10K = log10([1e5, 0.075e-6]);  % rate-limiting, 1. aerobic, 2. anaerobic (orig)
%Orgcsed_log10K = log10([0.075e-6, 0.075e-7]);  % rate-limiting, 1. aerobic, 2. anaerobic 
Orgcsed_log10K = log10([0.1e-6, 0.075e-6]);  % rate-limiting, 1. aerobic, 2. anaerobic 

% - set Orgcsource to replenish Orgcsed (Orgcsed is what actually degrades)
fl_Orgcsource = 1; % 0: ignore all parts related to Orgcsource
if fl_Orgcsource
    % - set Orgcsource ic conc [mol/g] 
    Orgcsource_conc_ic = Orgcsed_conc_ic * 10000;  % an amount that won't run out

    % - kinetic parameters (we're assuming 1st order decay), use same params as
    % Orcsed for steady-state replenishment
    % (increase logK for faster decay)
    Orgcsource_log10K = Orgcsed_log10K;  % rate-limiting, 1. aerobic, 2. anaerobic 
end

% 1F) SET UP GEOCHEMICAL INITIAL AND BOUNDARY CONDITIONS
% ***NOTE: Most of i.c. and b.c. should be set in function: Minntac_InitCond_chem.m ***

% -- Generate equilibrated and charge-balanced solutions based on observations:
% (nrow,ncol,nlay,n_comp)
[mob_eq_comp, mob_eq_ic_z, mob_eq_extra_z, min_eq_comp, min_eq_ic_z, catex_comp, catex_ic_z, ...
    surf_comp, surf_ic_z, surf_par, surf_cpl, surf_calc_type] = ...
    GW003_Minntac_InitCond_chem_red_wFeS(sim_dir, phrq_exe, use_file_databas, por, tempC);
if isempty(mob_eq_comp)
    fprintf('Minntac_InitCond_chem did not return valid results, exiting... \n')
    return
end

% -- Most likely no need to change this block: 
% Initialize arrays for initial conditions, boundary conditions will 
% automatically be set to the initial conditions in the corresponding grid cells
n_mob_eq = length(mob_eq_comp);
n_min_eq = length(min_eq_comp);
n_catex = length(catex_comp);
n_surf = length(surf_comp);
mob_eq_ic = zeros(nrow,ncol,nlay,n_mob_eq); % mobile equil components
mob_eq_extra = cell(n_mob_eq,1); % to specify charge-balance component
min_eq_ic = zeros(nrow,ncol,nlay,n_min_eq); % mineral equil components 
fl_mob_eq_const_rech = zeros(n_mob_eq,1); % default: apply constant concentration in recharge over space 
mob_eq_const_rech = zeros(n_mob_eq,1); % to specify space-constant concentration in recharge   
mob_eq_distr_rech = zeros(nrow,ncol,n_mob_eq); % to specify distributed concentration in recharge over space
catex_ic = zeros(nrow,ncol,nlay,n_catex); 
surf_ic = zeros(nrow,ncol,nlay,n_surf); 

% -- *** Set initial conditions here according the the various equilibrated
% and charge-balanced solutions (different "zones") from
% Minntac_InitCond_chem() function ***
for ii = 1: n_mob_eq
    % Rest of domain
    mob_eq_ic(:,:,1:end,ii) = mob_eq_ic_z(2,ii); 
    mob_eq_extra(ii) = mob_eq_extra_z(2,ii); 
    
    % Cell 2 boundary
    mob_eq_ic(:,1,:,ii) = mob_eq_ic_z(1,ii); % cell 2 
    
    % recharge 
    mob_eq_const_rech(ii) = mob_eq_ic_z(3,ii); % spatially constant
    mob_eq_distr_rech(:,:,ii) = mob_eq_ic_z(3,ii); % spatially distributed
    mob_eq_distr_rech(1,1:(400/x_scale),ii) = mob_eq_ic_z(4,ii); % recharge concentration thru perimeter dike
end
for ii = 1: n_min_eq
    % Rest of domain
    min_eq_ic(:,:,1:end,ii) = min_eq_ic_z(2,ii);     

    % Cell 2
    min_eq_ic(:,1,:,ii) = min_eq_ic_z(1,ii); % cell 2 
end
for ii = 1: n_catex
    % Rest of domain
    catex_ic(:,:,:,ii) = catex_ic_z(2,ii); 
end
for ii = 1: n_surf
    % Rest of domain
    surf_ic(:,:,:,ii) = surf_ic_z(2,ii); 
end

%==================================================
% RECHARGE CONCENTRATIONS: SULFATE & CHLORIDE (mg/L)
%==================================================
ind_sulfate = find(strcmp(mob_eq_comp, 'S(6)'));
ind_DO = find(strcmp(mob_eq_comp, 'O(0)'));
% mob_eq_distr_rech(1,1:(216/x_scale),ind_sulfate) = 10 % recharge concentration of Sulfate from perimeter dike
% mob_eq_distr_rech(1,1:(216/x_scale),ind_sulfate) = mob_eq_ic_z(3,ind_sulfate); % recharge concentration of Sulfate from perimeter dike
%mob_eq_distr_rech(1,(216/x_scale):(400/x_scale),11) = 0 % recharge concetrations of Sulfate from natural land surface
%mob_eq_distr_rech(1,1:(216/x_scale),3) = 0 % recharge concentrations of Cl from perimeter dike
%mob_eq_distr_rech(1,(216/x_scale):(400/x_scale),3) = 0 % recharge concentrations of Cl from natural land surface

% force anoxic
mob_eq_ic(:,2:end,1:end,ii) = 0; % rest of domain (not Cell water)
mob_eq_distr_rech(1,1:(240/x_scale),ind_DO) = 0; 


% 1G) ADDITIONAL OUTPUT VARIABLES FOR 'SELECT' FILE, SET IN POSTFIX FILE
% Default outputs will include most information, below block allows additional more detailed outputs 
% -- additional .sel output variables (other than input components)
addl_sel_outlist.total = {'C(-4)'}; % total concentration [mol/Lw] of component of specified oxidation state
addl_sel_outlist.mol = {}; % concentration [mol/Lw] of specified compound (e.g. HCO3- instead of total C(4))
addl_sel_outlist.equilphase = {}; % concentration [mol/Lw] of specified equilibrium mineral or gas phase
addl_sel_outlist.si = {}; % saturation index of specified mineral or gas phase
addl_sel_outlist.gas = {}; % concentration [mol/Lw] of specified equilibrium gas phase


% 1H) MISCELLANEOUS INPUTS
fl_setup_only = 1;  % 0: to directly run simulation from this matlab script



%% === In general, do not change below... =================================

% ============================
% SECTION 2: MODEL SET-UP
% ============================

% 2A) DIRECTORY AND FILE SETTINGS

% - file names
file_databas = [sim_dir, 'pht3d_datab.dat'];  % cannot change this name
btn_file = [sim_dir, 'pht3dbtn_', suffix, '.dat'];
ssm_file = [sim_dir, 'pht3dssm_', suffix, '.dat'];
dsp_file = [sim_dir, 'pht3ddsp_', suffix, '.dat'];
ph_file = [sim_dir, 'pht3d_ph.dat'];

% in this script, these files do not change
adv_file = [sim_dir, slashstr, 'pht3dadv.dat'];
gcg_file = [sim_dir, slashstr, 'pht3dgcg.dat'];

% specify general output file name
out_file = [sim_dir, slashstr, 'pht3d.out'];

% specify additional "select" output file, allows for extra detailed ouput information
sel_file = 'out_X.sel';


% - go to simulation dir (check to make sure not over-writing, no current programs running)
if ~exist(sim_dir, 'dir')
    mkdir(sim_dir);
end
%if exist([sim_dir, slashstr, sel_file], 'file')    
%    fprintf('sim_dir has .sel file(s)!  Could be job running or unsaved job there.  Exiting...\n');
%    fprintf('(sim_dir %s) \n', sim_dir);
%    return
%end
curr_dir = pwd;
cd(sim_dir);
addpath(curr_dir);

% set database file (must be in sim_dir with file name 'pht3d_datab.dat')
if ~strcmp(use_file_databas, file_databas)
    copyfile(use_file_databas, file_databas);
end


% 2B) DOMAIN PARAMETERS
htop = domain_top_elev;
DELC = domain_len / ncol; % thickness of column 
DELR = DELC; % thickness of row (doesn't matter in x-section)
dz = repmat((domain_top_elev - domain_bot_elev)/nlay, nlay, 1);


n_par_max = 10; % max number kinetic parameters allowed


% 2C) SET REMAINING VARIABLES FOR GEOCHEMISTRY INITIAL AND BOUNDARY CONDITIONS
% (Units -- aq: mol/L_w, user-defined immob (e.g. bacteria, napl): mol/L_w, 
% minerals (and gases?): mol/L_v, exchangers and surfaces: mol/L_v)

% - mobile kinetic components 
n_mob_kin_max = 10;
mob_kin_comp = cell(n_mob_kin_max,1);
mob_kin_ic = zeros(nrow,ncol,nlay,n_mob_kin_max);
mob_kin_par = nan(n_par_max, n_mob_kin_max);
mob_kin_formula = cell(n_mob_kin_max,1);
ii = 0;
n_mob_kin = ii;
mob_kin_comp = mob_kin_comp(1:n_mob_kin);
mob_kin_ic = mob_kin_ic(:,:,:,1:n_mob_kin); 

% - mobile equil components
% (already set in INPUT section)

% - immobile kinetic componentscomponents
n_imob_kin_max = 10;
imob_kin_comp = cell(n_imob_kin_max,1);
imob_kin_ic = zeros(nrow,ncol,nlay,n_imob_kin_max); 
imob_kin_par = nan(n_par_max, n_mob_kin_max);
imob_kin_formula = cell(n_mob_kin_max,1);
ii = 0;
% No Orgcsed in dike!!
ii = ii+1; imob_kin_comp{ii} = 'Orgcsed'; % ********   
imob_kin_par(1:length(Orgcsed_log10K),ii) = 10.^Orgcsed_log10K;
imob_kin_ic(:,:,:,ii) = Orgcsed_conc_ic; 
imob_kin_ic(:,1:round(160/x_scale),1:round(66/y_scale),ii) = 0;
imob_kin_ic(:,round(160/x_scale):round(240/x_scale),1:round(88/y_scale),ii) = 0;
imob_kin_ic(:,round(240/x_scale):round(400/x_scale),1:round(66/y_scale),ii) = 0;
% imob_kin_formula{ii} = 'Orgcsed -1.0 Orgc 1.0 ';
imob_kin_formula{ii} = 'Orgcsed -1.0 CH2O 1.0 ';
if fl_Orgcsource
    ii = ii+1; imob_kin_comp{ii} = 'Orgcsource'; % ********   
    imob_kin_par(1:length(Orgcsource_log10K),ii) = 10.^Orgcsource_log10K;
    imob_kin_ic(:,:,:,ii) = Orgcsource_conc_ic; 
    imob_kin_ic(:,1:round(160/x_scale),1:round(66/y_scale),ii) = 0;
    imob_kin_ic(:,round(160/x_scale):round(240/x_scale),1:round(88/y_scale),ii) = 0;
    imob_kin_ic(:,round(240/x_scale):round(400/x_scale),1:round(66/y_scale),ii) = 0;
    imob_kin_formula{ii} = 'Orgcsource -1.0 Orgcsed 1.0 ';
end
n_imob_kin = ii;
imob_kin_comp = imob_kin_comp(1:n_imob_kin);
imob_kin_ic = imob_kin_ic(:,:,:,1:n_imob_kin); 

% - mineral eq components (ic specified below)
% (already set in INPUT section)

% - catex components (catex_exch is exchanger site)
% (already set in INPUT section)

% - surface complexation components
% (already set in INPUT section)

% - mineral kinetic components
n_min_kin_max = 10;
min_kin_comp = cell(n_min_kin_max,1);
min_kin_ic = zeros(nrow,ncol,nlay,n_min_kin_max); 
min_kin_par = zeros(1,n_min_kin_max); 
ii = 0;
n_min_kin = ii;
min_kin_comp = min_kin_comp(1:n_min_kin);
min_kin_ic = min_kin_ic(:,:,:,1:n_min_kin); 
min_kin_par = min_kin_par(1,1:n_min_kin); 


% 2D) CREATE INPUT FILES

% -- create btn, ssm, and disp files 
eff_por = ones(nrow,ncol,nlay) * por;
PHT3D_btn_f_051715a(btn_file, ssm_file, dsp_file, ...
    nrow, ncol, nlay, DELR, DELC, htop, dz, eff_por, ...
    long_disp, hor2longdisp, vert2longdisp, Dstar, ...
    mob_kin_comp, mob_kin_ic, ...
    mob_eq_comp, mob_eq_ic, fl_mob_eq_const_rech, mob_eq_const_rech, mob_eq_distr_rech, ...
    imob_kin_comp, imob_kin_ic, ...
    min_eq_comp, min_eq_ic, ...
    catex_comp, catex_ic, ...
    surf_comp, surf_ic, ...
    timprs, nstp, fl_rech);    


% -- Create nam file
fid = fopen(nam_fil, 'wt');
fprintf(fid, 'List   7    %s\n', out_file);
fprintf(fid, 'FTL    66  %s\n', flo_file);
fprintf(fid, 'BTN    31    %s\n', btn_file);
fprintf(fid, 'ADV    32    %s\n', adv_file);
fprintf(fid, 'DSP    33    %s\n', dsp_file);
fprintf(fid, 'SSM    34    %s\n', ssm_file);
fprintf(fid, 'GCG    35    %s\n', gcg_file);
fprintf(fid, 'PHC    64    %s\n', ph_file);
fclose(fid);


% -- Create files that don't change (gcg and adv)
% - Sets numerical method inputs for advection term
MIXELM = 2; % 0: FD, 1: MOC, 2: MMOC, 3: HMOC
fid = fopen(adv_file, 'wt');
fprintf(fid, '         %d       .75      5000         0\n', MIXELM);
fprintf(fid, '         3        .5\n');
fprintf(fid, '         1         0        15\n');
fclose(fid);

fid = fopen(gcg_file, 'wt');
fprintf(fid, '        10       500         1         0\n');
fprintf(fid, '         1    .00001         0\n');
fclose(fid);


% -- create _ph file (this added to handle alkalinity)
fl_catex_toequil = 1;
PHT3D_ph_f_general4(ph_file, ...
    tempC, ...
    mob_kin_comp, mob_kin_par, mob_kin_formula, ...
    mob_eq_comp, mob_eq_extra, ...
    imob_kin_comp, imob_kin_par, imob_kin_formula, ...
    min_eq_comp, ...
    catex_comp, fl_catex_toequil, ...
    surf_comp, surf_par, surf_cpl, surf_calc_type, ...
    min_kin_comp, min_kin_par)


% -- create postfix file (what to output to .sel file)
% - select output list
sel_outlist.total = [mob_kin_comp; mob_eq_comp(~strncmp('p',mob_eq_comp,1)); imob_kin_comp];
sel_outlist.mol = catex_comp;  % better to put desired surf species in add_sel_outlist.mol
% sel_outlist.mol = [catex_comp; surf_comp];
sel_outlist.equilphase = min_eq_comp;
sel_outlist.si = {'Siderite', 'FeS(ppt)'}; 
% sel_outlist.si = []; 
sel_outlist.gas = []; 
% sel_outlist.si = {'CH4(g)', 'CO2(g)', 'Ntwo(g)', 'O2(g)'}; 
% sel_outlist.gas = {'CH4(g)', 'CO2(g)', 'Ntwo(g)', 'O2(g)'}; 

% include any additional components
sel_outlist.total = [sel_outlist.total; addl_sel_outlist.total(:)];
sel_outlist.mol = [sel_outlist.mol; addl_sel_outlist.mol(:)];
sel_outlist.equilphase = [sel_outlist.equilphase; addl_sel_outlist.equilphase(:)];
sel_outlist.si = [sel_outlist.si; addl_sel_outlist.si(:)]; 
sel_outlist.gas = [sel_outlist.gas; addl_sel_outlist.gas(:)]; 

postfix_fil = [sim_dir, 'postfix.phrq'];
fid = fopen(postfix_fil, 'wt');
fprintf(fid, 'SELECTED_OUTPUT \n');
fprintf(fid, '    -file %s \n', sel_file);
fprintf(fid, '    -reset false \n');
fprintf(fid, '    -time  true \n');
fprintf(fid, '    -ph \n');
fprintf(fid, '    -pe \n');
if ~isempty(sel_outlist.total)
    fprintf(fid, '    -total ');
    for ii = 1: length(sel_outlist.total)
        fprintf(fid, ' %s', sel_outlist.total{ii});
    end
    fprintf(fid, '\n');
end
if ~isempty(sel_outlist.mol)
    fprintf(fid, '    -molalities ');
    for ii = 1: length(sel_outlist.mol)
        fprintf(fid, ' %s', sel_outlist.mol{ii});
    end
    fprintf(fid, '\n');
end
if ~isempty(sel_outlist.equilphase)
    fprintf(fid, '    -equilibrium_phases ');
    for ii = 1: length(sel_outlist.equilphase)
        fprintf(fid, ' %s', sel_outlist.equilphase{ii});
    end
    fprintf(fid, '\n');
end
if ~isempty(sel_outlist.si)
    fprintf(fid, '    -saturation_indices ');
    for ii = 1: length(sel_outlist.si)
        fprintf(fid, ' %s', sel_outlist.si{ii});
    end
    fprintf(fid, '\n');
end
if ~isempty(sel_outlist.gas)
    fprintf(fid, '    -gases ');
    for ii = 1: length(sel_outlist.gas)
        fprintf(fid, ' %s', sel_outlist.gas{ii});
    end
    fprintf(fid, '\n');
end
fprintf(fid, 'END\n');

%% ------------------------------------------------------------------------
% run simulation

if isunix
    command = [PHT3D_exe, ' ', nam_fil, ' &> out.txt'];
else
    command = [PHT3D_exe, ' ', nam_fil];
end    
if fl_setup_only
    fprintf('Ready for manual execution in %s!\n', sim_dir);
    fprintf('(%s)\n', command)
else
    system(command);
end
fprintf('\n');


cd(curr_dir);

rmpath(matlab_dir);


