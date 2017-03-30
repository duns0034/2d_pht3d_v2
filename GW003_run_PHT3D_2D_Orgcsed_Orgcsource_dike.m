% run_PHT3D_1Dprofile_quickstart.m
%
% 11/3/15
%
% - uses PHREEQC to first equilibrate ic w/ minerals and charge balance
%
% Sets up model by custom-creating the following input files:
%   1) *nam: namelist file 
%   2) *_ph.dat: PHREEQC interface package file (to link PHT3D to PHREEQC)
%   3) postfix: additional PHREEQC options 
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
% Assumes pht3d_databas.dat (geochem database) file already exists; specify
% in 'use_file_databas'.
%       
%
% Uses the following functions:
%   - SC_InitCond_chem_2()
%       SpecifiesOr observed conc, equilibrates with PHREEQC
%   - Alk2DIC()
%       Called by SC_InitCond_chem(), converts alkalinity to total C(4)
%   - generate_ic_PHREEQC_f_062215a()
%       Called by SC_InitCond_chem_2(), generates equilibrated ic using
%          PHREEQC
%   - go_PHREEQC2_4()
%       Called by generate_ic_PHREEQC_f_062215a(), runs PHREEQC
%   - PHT3D_ph_f_general5()
%       Creates pht3d_ph.dat
%   - PHT3D_btn_f_051715a()
%       Creates btn, dsp, ssm files using the specified inputs
%
% run_PHT3D_1Dprofile_quickstart_v2: includes As as redox species
% run_PHT3D_1Dprofile_quickstart_v3: includes As sorption
% run_PHT3D_1Dprofile_quickstart_v4: includes temperature-controlled Orgc degradation


clear all, close all; fclose all;

a = clock;
if a(5) < 10
    fprintf('\n\nStart time: (%d/%d/%d) %d:0%d\n', a(2), a(3), a(1), a(4), a(5));
else
    fprintf('\n\nStart time: (%d/%d/%d) %d:%d\n', a(2), a(3), a(1), a(4), a(5));
end


% **************** CUSTOMIZE TO YOUR COMPUTER!! ***************************
% Make sure to set the following: 
%   matlab_dir, PHT3D_exe, phrq_exe, sim_dir, flo_file, use_file_databas

fl_gcng = 0;  % 1: Crystal, 0: Aubrey

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
    sim_dir = '/home/gcng/workspace/ModelRuns_scratch/PHT3D_projects/CapeCod/Arsenic/test1_c2/';

    % - flo_file: MODFLOW flo simulation file (full path)
%     flo_file = 'C:\Hydro_Modeling\MINNTAC_MATLAB_FILES\test2_1D_3\test.flo'; % 2D, no recharge
%     flo_file = '/home/gcng/workspace/ModelRuns_scratch/MODFLOW_projects/ESCI5980/Asst6_oppdir/test.flo';
    flo_file = '/home/gcng/workspace/ModelRuns_scratch/MODFLOW_projects/CapeCod/Arsenic/test_1D_1/test_AJD_AsMN_50m_200lay_downwardflux.flo';    
%     fl_rech = 1;  % 1: for recharge, 0 for no recharge in the MODFLOW simulation

    % - use_file_databas: geochem database, this is copied into 'pht3d_datab.dat' in sim_dir for simulation
%     use_file_databas = '/home/gcng/Documents/Teaching/ESCI5980_HydModeling/Fall2015/Lectures/L16_PHT3D_2/PHT3D_files/pht3d_datab.dat';
    use_file_databas = '/home/gcng/workspace/ProjectFiles/CapeCod/arsenic/database_files/pht3d_datab_As_160606_AJDedits_gcng.txt';
    % ************ (end of computer-specific file specifications) *************
else
        % - choose one of these
    %slashstr = '/'; % for Linux or Windows/Cygwin
    slashstr = '\'; % for Windows/MS-DOS 
    
    % - matlab_dir: directory with matlab functions 
    matlab_dir = pwd;  % current directory

    % - Full path PHT3D executable
    PHT3D_exe = 'C:\pht3dv210\bin\pht3dv210';

    % - Full path PHREEQC executable
    phrq_exe = 'C:\USGS\phreeqc-3.3.7-11094-x64\bin\phreeqc';

    % - sim_dir: Directory with input files and where you run simulations
    sim_dir = 'C:\Users\Owner\Desktop\Arsenic_MN\Arsenic_MN_Project\OTT3_Model\OTT3_Test_Models\Test_Model_03142017\';

    % - flo_file: MODFLOW flo simulation file (full path)
    flo_file = 'C:\Users\Owner\Desktop\Arsenic_MN\Arsenic_MN_Project\OTT3_Model\OTT3_Test_Models\Test_Model_03142017\test_100m_200lay_noflow.flo';

    % - use_file_databas: geochem database, this is ultimately copied into 
    %   'pht3d_datab.dat' in sim_dir for simulation
    use_file_databas = 'C:\Users\Owner\Desktop\Arsenic_MN\Arsenic_MN_Project\OTT3_Model\OTT3_Test_Models\Test_Model_03142017\pht3d_datab_As_03062017_OTT3.dat';
    % ************* (end to CUSTOMIZE TO YOUR COMPUTER!!) *********************
end

% - This will contain names of all input files, including MODFLOW flo file, 
%   Many files depend on resolution; 'suffix' is for names of resolution-
%   dependent files (btn, dsp, ssm)
nam_fil = 'pht3d_1D.nam';
suffix = '1D';  % suffix to file names

%fl_NoChargeBal = 1; % will always do charge balance for first initialization

% -- other entries for _ph, btn files

ctr=0;
fl_ReDo = 0;

% end_ctr = 15; % any number, can loop thru tempC_v
end_ctr = 1; % any number, can loop thru tempC_v

OperSplit = 3; % 2 for reaction step every user-specified time step, 
%       3 for reaction every system-determined transport step (all runs
%       previous to 160920 were OperSplit=2, hard-coded in
%       PHT3D_ph_f_general5

 while(1)
    ctr=ctr+1;
    if fl_ReDo
        % redo with different charge balance, don't advance ctr
        ctr=ctr-1;
    end

% %Two year binary temperature model.    
tempC_kin_v = [10 10 10 10];  % for kinetics
% tempC_kin_v = [25 25 2 25];  % for kinetics
% tempC_kin_v = [25];  % for kinetics
tempC_eq_v = [10];  % for equil

%     tempC=tempC_v(ctr);

%One year month-dependent temperature model
%     tempC_v = [2; 4; 5; 7; 10; 15; 20; 25; 23; 19; 13; 8];
% 
%     tempC=tempC_v(ctr);

    tempC_eq=tempC_eq_v(mod(ctr-1,length(tempC_eq_v))+1);
    tempC_kin=tempC_kin_v(mod(ctr-1,length(tempC_kin_v))+1);
    fprintf('ctr=%d, tempC_eq=%g, tempC_kin=%g\n', ctr, tempC_eq, tempC_kin);

    units = 'mol/kgw';

    if ctr == 1
%         timprs=[0 1 30 60 90 120 150 180 210 240 270 300 330 365*[100:2000:20000]];
%         timprs=[0 1 30 60 90 120 150 180 210 240 270 300 330 365*[100:1500:15000]];
%         timprs=[0 1 30 60 90 120 150 180 210 240 270 300 330 365*[100:1000:10000]];
%         timprs=[0 1 30 60 90 120 150 180 210 240 270 300 330 365*[100:100:2000]];
%         timprs=[0 1 30 60 90 120 150 180 210 240 270 300 330 365*[100:100:1000]];
        timprs=[0 1 365*1 365*[100:100:1000]];
%         timprs=[0 1 30 60 90 120 150 180 210 240 270 300 330 365*[1:100]];
%         timprs=[0 1 30 60 90 120 150 180 210 240 270 300 330 365*[1:5]];
%         timprs=[0 1 30 60 90 120 150 180 210 240 270 300 330 365];
%         timprs=[1 91.25 182.5];
%          timprs=[0 1 30 60 90];
    else
%         timprs=[0 30 60 90 120 150 180];
        timprs=[91.25 182.5];
    end
        
    %timprs_v=[0:179; 180:359; 360:539; 540:719];
    %timprs_v=[0:29; 30:59; 60:89; 90:119; 120:149; 150:179; 180:209; 210:239; 240:269; 270:299; 300:329; 330:359];
    %timprs_v=[0 29; 30 59; 60 89; 90 119; 120 149; 150 179; 180 209; 210 239; 240 269; 270 299; 300 329; 330 359];
    
%         timprs=timprs_v(ctr,1:end); 

%timprs = [0 1 90 180 365 365*[2:10]]; % print out times [d]
%timprs = [0: 1 : 30]; % print out times [d]
% nstp = round(2/180*timprs(end));  % number of steps = steps per 182.5 days times total number of days in run
nstp = 50;
% nstp = round(10/180*timprs(end));  % often used 50 for longer runs, 10 is faster
por = 0.20;   % porosity, 0.39 assumed in Smith et al. 2013
rho_b = 2650; % g/L (bulk density), 1864g/L from Smith et al. 2013
pe_value = 12;

fl_Yexch = 1; % 1 to include Y cation exchanger

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
 %long_disp = 1e-2;  % m
long_disp = 1;  % m
hor2longdisp = 0.018; % ratio horiz transverse disp / long dispersivity (Garabedian: 0.018 m)
vert2longdisp = 0.0015; % ratio vertical transverse disp / long dispersivity (Garabedian: 0.0015 m)
hivert2longdisp = 0.015; % ratio vertical transverse disp / long dispersivity (Garabedian: 0.0015 m), used to create impression of increased flux in aquifer


% Domain info ------------------------------------------------------------
% - Domain parameters (for ba6, dis)
nlay = 200;
% nlay = 100;
% nlay = 50;
ncol = 1;   
nrow = 1;
% domain_len = 50; % meters
domain_bot_elev = -100; % meters
domain_top_elev = 0; % top of domain must be at least this elev (include extra space for WT mov't)
% DELR = domain_len / ncol;  DELC = DELR; % discretization for 2D profile
% dz = repmat((domain_top_elev - domain_bot_elev)/nlay, nlay, 1);
DELR = 100;  DELC = 100; % arbitrary for 1D profile
dz = repmat((domain_top_elev - domain_bot_elev)/nlay, nlay, 1);


% -- various flags
fl_setup_only = 0;  % 0: to directly run simulation from this matlab script


%% Orgcsed-----------------------------------------------------------------

% for i=length(ctr)
%     if ctr==1

% -- set high limit for Orgsed ic conc (zone 2) [mol/g]
% hi_Orgcsed = 0.3*por/rho_b;  % [mol/Lw] * por / rho_b = [mol/g]
%hi_Orgcsed = 0.1*por/rho_b;  % [mol/Lw] * por / rho_b = [mol/g]
% hi_Orgcsed = 0.03*por/rho_b;  % [mol/Lw] * por / rho_b = [mol/g]
hi_Orgcsed = 0.01*por/rho_b;  % [mol/Lw] * por / rho_b = [mol/g]
% - kinetic parameters
Orgcsed_log10K = log10([1e3, 0.075e-8]);  % rate-limiting, 1. aerobic, 2. anaerobic 

%% Orgc--------------------------------------------------------------------

% -- set high limit for Orgsed ic conc (zone 2) [mol/g]
% hi_Orgc = 0.3;  % [mol/Lw]
% hi_Orgc = 0.1;  % [mol/Lw]
hi_Orgc = 0.03;  % [mol/Lw]
% hi_Orgc = 0;  % [mol/Lw]
% - kinetic parameters
Orgc_log10K = log10([5e-4, 4e-9]);  % rate-limiting, 1. aerobic, 2. anaerobic (starting point is Cape Cod reoxy model)
% - temperature-dependency parameters for: lnk = m*(1/T_K) + b, (T in K, k in s-1)
% (Arrhenius model, based on conversation with Doug 6/1/16)
% Set these params:
Pt1_T_k_Aer = [9+273.15; 0.016/3600]; % T [K], 1st order k_Aer [1/s] (point 1)
dT_mult_Aer = [10; 2];  % dT [K], multiply rate [-] (to get aerobic point 2)
An_Aer_ratio_9degC = 1e-6;  % k_An / k_Aer at 9degC (to get anaerobic point 1)
dT_mult_An = [10; 5];  % dT [degK], multiply rate [-] (to get anaerobic point 2)

% do not change below:
Pt2_T_k_Aer = [Pt1_T_k_Aer(1)+dT_mult_Aer(1); Pt1_T_k_Aer(2)*dT_mult_Aer(2)]; % T [K], 1st order k_Aer [1/s] (point 2)
Pt1_T_k_An = [Pt1_T_k_Aer(1); Pt1_T_k_Aer(2)*An_Aer_ratio_9degC]; % T [K], 1st order k_An [1/s] (point 1)
Pt2_T_k_An = [Pt1_T_k_An(1)+dT_mult_An(1); Pt1_T_k_An(2)*dT_mult_An(2)]; % T [K], 1st order k_An [1/s] (point 2)
for ii = 1:2 % loop thru Aerobic, Anaerobic
    if ii == 1
        Pt1_T_k = Pt1_T_k_Aer;
        Pt2_T_k = Pt2_T_k_Aer;
    else
        Pt1_T_k = Pt1_T_k_An;
        Pt2_T_k = Pt2_T_k_An;
    end
    Pt1 = [1/Pt1_T_k(1), log10(Pt1_T_k(2))];
    Pt2 = [1/Pt2_T_k(1), log10(Pt2_T_k(2))];
    m = (Pt2(2)-Pt1(2)) / (Pt2(1)-Pt1(1));  % slope for: lnk = m*(1/T_K) + b
    b = -m*Pt1(1) + Pt1(2);
    if ii == 1
        Orgc_Aer_m = m;
        Orgc_Aer_b = b;
    else
        Orgc_An_m = m;
        Orgc_An_b = b;
    end
end
Orgc_log10K = [Orgc_Aer_m*(1/(273.15+tempC_kin))+Orgc_Aer_b, Orgc_An_m*(1/(273.15+tempC_kin))+Orgc_An_b];  % rate-limiting, 1. aerobic, 2. anaerobic (starting point is Cape Cod reoxy model)

%     end
% end

%% Pyrite------------------------------------------------------------------
% 
% %Example of KINETICS data block for pyrite rate:
% %	KINETICS 1
% %	Pyrite  
% %		-tol    1e-8
% %		-m0 	5.e-4
% %		-m	5.e-4
% %		-parms -5.0 	0.1 	.5 	-0.11
% 
% Pyrite_log10K = log10([1e-5, 1.2589, 3.1623, 0.7762]);  
%% ------------------------------------------------------------------------
% -- additional .sel output variables (other than input components)
addl_sel_outlist.total = {'C(-4)'};
% addl_sel_outlist.mol = {'Hfo_wH2AsO4', 'Hfo_wHAsO4-', 'Hfo_wOHAsO4-3', 'Hfo_wH2AsO3', 'Hfo_sOFe+', 'Hfo_wOFe+', 'Hfo_wOFeOH', 'FeY2'};
addl_sel_outlist.mol = {'Hfo_wH2AsO4', 'Hfo_wHAsO4-', 'Hfo_wOHAsO4-3', 'Hfo_wH2AsO3', 'Hfo_wH2PO4', 'Hfo_wHPO4-', 'Hfo_wPO4-2'};
addl_sel_outlist.equilphase = {};
addl_sel_outlist.si = {}; 
addl_sel_outlist.gas = {}; 

sel_file = 'out_X.sel';

%% In general, do not change below...


%% -- Directory settings:
% - Enter simulation files here:
phrq_sim_dir = sim_dir;


%% file names

file_databas = [sim_dir, 'pht3d_datab.dat'];  % cannot change this name
btn_file = [sim_dir, 'pht3dbtn_', suffix, '.dat'];
ssm_file = [sim_dir, 'pht3dssm_', suffix, '.dat'];
dsp_file = [sim_dir, 'pht3ddsp_', suffix, '.dat'];
ph_file = [sim_dir, 'pht3d_ph.dat'];

% in this script, these files do not change
adv_file = [sim_dir, slashstr, 'pht3dadv.dat'];
gcg_file = [sim_dir, slashstr, 'pht3dgcg.dat'];

% to specify output file name
out_file = [sim_dir, slashstr, 'pht3d.out'];


%% ----------------------------------------------------------------------
% Generally don't need to change below here:

addpath(matlab_dir);

% - go to simulation dir (check to make sure not over-writing, no current
% programs running)
if ~exist(sim_dir, 'dir')
    mkdir(sim_dir);
end
% if exist([sim_dir, slashstr, sel_file], 'file')    
%     fprintf('sim_dir has .sel file(s)!  Could be job running or unsaved job there.  Exiting...\n');
%     fprintf('(sim_dir %s) \n', sim_dir);
%     return
% end
curr_dir = pwd;
cd(sim_dir);
addpath(curr_dir);

if ~strcmp(use_file_databas, file_databas)
    copyfile(use_file_databas, file_databas);
end

htop = domain_top_elev;

n_par_max = 10; % max number kinetic parameters allowed

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

if ctr==1

% -- SET INITIAL MODEL CHEMISTRY FOR BTN FILE
% (Units -- aq: mol/L_w, user-defined immob (e.g. bacteria, napl): mol/L_w, 
% minerals (and gases?): mol/L_v, exchangers and surfaces: mol/L_v)

% - mobile kinetic components
n_mob_kin_max = 10;
mob_kin_comp = cell(n_mob_kin_max,1);
mob_kin_ic = zeros(nrow,ncol,nlay,n_mob_kin_max);
mob_kin_par = nan(n_par_max, n_mob_kin_max);
mob_kin_formula = cell(n_mob_kin_max,1);
ii = 0;
% ii=ii+1; mob_kin_comp{ii} = 'Orgc'; % *************************************
% mob_kin_par(1:length(Orgc_log10K),ii) = 10.^Orgc_log10K;
% mob_kin_ic(:,:,1:50,ii) = 0;
% mob_kin_ic(:,:,51:end,ii) = 0;
% 
% mob_kin_formula{ii} = 'Orgc -1.0 CH2O 1.0 ';
n_mob_kin = ii;
mob_kin_comp = mob_kin_comp(1:n_mob_kin);
mob_kin_ic = mob_kin_ic(:,:,:,1:n_mob_kin);  

% - mobile equil components
n_mob_eq_max = 50;
mob_eq_comp = cell(n_mob_eq_max,1);
% - (ic from SC_InitCond_chem.m)    
mob_eq_ic = zeros(nrow,ncol,nlay,n_mob_eq_max);
% - (constant conc bc, not recharge conc bc)
fl_mob_eq_const_rech = ones(n_mob_eq_max,1); % default cont rech
mob_eq_const_rech = zeros(n_mob_eq_max,1);    
mob_eq_distr_rech = zeros(nrow,ncol,n_mob_eq_max); 
mob_eq_extra = cell(n_mob_eq_max,1);
ii = 0;
% ii=ii+1; mob_eq_comp{ii} = 'Al'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
% ii=ii+1; mob_eq_comp{ii} = 'B'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
%ii=ii+1; mob_eq_comp{ii} = 'H(1)'; % ***********************************
ii=ii+1; mob_eq_comp{ii} = 'As(3)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
ii=ii+1; mob_eq_comp{ii} = 'As(5)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
ii=ii+1; mob_eq_comp{ii} = 'C(4)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
%ii=ii+1; mob_eq_comp{ii} = 'C(-4)'; % ***********************************
ii=ii+1; mob_eq_comp{ii} = 'Ca'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 21; % mg/L
ii=ii+1; mob_eq_comp{ii} = 'Cl'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 2.42; % mg/L
%            mob_eq_extra{ii} = 'charge';
ii=ii+1; mob_eq_comp{ii} = 'Fe(2)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 8.9; % mg/L
ii=ii+1; mob_eq_comp{ii} = 'Fe(3)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 8.9; % mg/L
ii=ii+1; mob_eq_comp{ii} = 'K'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 2.89; % mg/L
ii=ii+1; mob_eq_comp{ii} = 'Mg'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 21.1; % mg/L
ii=ii+1; mob_eq_comp{ii} = 'Mn(2)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 3.29; % mg/L
% ii=ii+1; mob_eq_comp{ii} = 'N'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 1.35; % mg/L
% ii=ii+1; mob_eq_comp{ii} = 'N(5)'; % ***********************************
% ii=ii+1; mob_eq_comp{ii} = 'N(0)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
ii=ii+1; mob_eq_comp{ii} = 'Amm'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
ii=ii+1; mob_eq_comp{ii} = 'Na'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 6.89; % mg/L
%if ~fl_NoChargeBal
%          mob_eq_extra{ii} = 'charge';
%end
ii=ii+1; mob_eq_comp{ii} = 'O(0)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
ii=ii+1; mob_eq_comp{ii} = 'P'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0.2; % mg/L
% ii=ii+1; mob_eq_comp{ii} = 'S(-2)'; % ***********************************
% % mob_eq_ic(:,:,:,ii) = 0.067; % mg/L
ii=ii+1; mob_eq_comp{ii} = 'S(6)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 26.6; % mg/L
ii=ii+1; mob_eq_comp{ii} = 'S(-2)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 26.6; % mg/L
ii=ii+1; mob_eq_comp{ii} = 'Si'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 46.3; % mg/L
ii=ii+1; mob_eq_comp{ii} = 'N(5)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
ii=ii+1; mob_eq_comp{ii} = 'N(3)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
ii=ii+1; mob_eq_comp{ii} = 'N(0)'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
ii=ii+1; mob_eq_comp{ii} = 'pH'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 6.92;
ii=ii+1; mob_eq_comp{ii} = 'pe'; % ***********************************
% mob_eq_ic(:,:,:,ii) = 0;
% ii=ii+1; mob_eq_comp{ii} = 'Orgc'; % ***********************************
n_mob_eq = ii;
mob_eq_comp = mob_eq_comp(1:n_mob_eq);
mob_eq_ic = mob_eq_ic(:,:,:,1:n_mob_eq);
mob_eq_const_rech = mob_eq_const_rech(1:n_mob_eq);
mob_eq_distr_rech = mob_eq_distr_rech(:,:,1:n_mob_eq);
mob_eq_extra = mob_eq_extra(1:n_mob_eq);

% - immobile kinetic components
n_imob_kin_max = 10;
imob_kin_comp = cell(n_imob_kin_max,1);
imob_kin_ic = zeros(nrow,ncol,nlay,n_imob_kin_max); 
imob_kin_par = nan(n_par_max, n_mob_kin_max);
imob_kin_formula = cell(n_mob_kin_max,1);
ii = 0;
ii = ii+1; imob_kin_comp{ii} = 'Orgcsed'; % ********     
imob_kin_par(1:length(Orgcsed_log10K),ii) = 10.^Orgcsed_log10K;
imob_kin_ic(:,:,1:100,ii) = hi_Orgcsed;
imob_kin_ic(:,:,101:200,ii) = 0;

% imob_kin_formula{ii} = 'Orgcsed -1.0 Orgc 1.0 ';  %For ever 1 mol
% %orgcsed, get 1 mol Orgc
imob_kin_formula{ii} = 'Orgcsed -1.0 CH2O 1.0 ';

% ii = ii+1; imob_kin_comp{ii} = 'Arsenopyrite'; % ********     
% imob_kin_par(1:length(Orgcsed_log10K),ii) = 10.^Orgcsed_log10K;
% imob_kin_ic(:,:,1:50,ii) = 1e-6*rho_b;
% imob_kin_ic(:,:,51:100,ii) = 0;

% imob_kin_formula{ii} = 'Arsenopyrite -1.0 Fe+2 1.0 As+3 0.008 SO4-2 1.0 ';

% ii = ii+1; imob_kin_comp{ii} = 'Pyrite'; % ********     
% imob_kin_par(1:length(Pyrite_log10K),ii) = 10.^Pyrite_log10K;
% imob_kin_ic(:,:,1:50,ii) = 1e-6;
% imob_kin_ic(:,:,51:100,ii) = 0;
% imob_kin_formula{ii} = 'Pyrite -1.0 Fe+2 1.0 SO4-2 2.0 ';

n_imob_kin = ii;
imob_kin_comp = imob_kin_comp(1:n_imob_kin);
imob_kin_ic = imob_kin_ic(:,:,:,1:n_imob_kin); 

% - mineral eq components (ic specified below)
n_min_eq_max = 10;
min_eq_comp = cell(n_min_eq_max,1);
min_eq_ic = zeros(nrow,ncol,nlay,n_min_eq_max); 
ii = 0;
ii=ii+1; min_eq_comp{ii} = 'Orpiment'; % **********************************
    min_eq_ic(:,:,1:100,ii) = 1.03e-20; % (mol/Lv)
ii=ii+1; min_eq_comp{ii} = 'Arsenopyrite'; % ******************************
    min_eq_ic(:,:,1:100,ii) = 6.23e-20; % (mol/Lv)
ii=ii+1; min_eq_comp{ii} = 'Fe(OH)3(a)'; % ********************************
    min_eq_ic(:,:,101:200,ii) = 6.23e-5; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Gypsum'; % ************************************
% min_eq_ic(:,:,51:100,ii) = 1e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Pyrite'; % ************************************
% min_eq_ic(:,:,1:50,ii) = 1e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Calite'; % ************************************
% min_eq_ic(:,:,51:100,ii) = 1e-20*rho_b; % (mol/Lv)  <-- DO NOT add
% calcite into AsMn model.  Causes model to crash due to charge imbalance.  
% ii=ii+1; min_eq_comp{ii} = 'Gibbsite'; % **********************************
% min_eq_ic(:,:,:,ii) = 30e-6*rho_b; % (mol/Lv)
 %ii=ii+1; min_eq_comp{ii} = 'Goethite'; % **********************************
 %min_eq_ic(:,:,:,ii) = 50e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Pyrolusite'; % ********************************
% min_eq_ic(:,:,:,ii) = .25e-6*rho_b;  % 
% ii=ii+1; min_eq_comp{ii} = 'Siderite'; % **********************************
% min_eq_ic(:,:,:,ii) = 1e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Rhodochrosite'; % *****************************
% min_eq_ic(:,:,:,ii) = 1e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'FeS(ppt)'; % **********************************
% min_eq_ic(:,:,1:50,ii) = 1e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Mackinawite'; % *******************************
% min_eq_ic(:,:,1:50,ii) = 1e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Hydroxyapatite'; % ****************************
% min_eq_ic(:,:,51:100,ii) = 1e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Hematite'; % **********************************
% min_eq_ic(:,:,51:100,ii) = 1e-6*rho_b; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Vivianite'; % *********************************
% min_eq_ic(:,:,51:100,ii) = 1e-6*rho_b; % (mol/Lv)
n_min_eq = ii;
min_eq_comp = min_eq_comp(1:n_min_eq);
min_eq_ic = min_eq_ic(:,:,:,1:n_min_eq); 

% - catex components (catex_exch is exchanger site)
n_catex_exch_max = 10;
catex_exch_comp = cell(n_catex_exch_max,1);
catex_exch_ic = zeros(nrow,ncol,nlay,n_catex_exch_max); 
ii = 0;
if fl_Yexch
    ii=ii+1; catex_exch_comp{ii} = 'Y'; % ***********************************
    catex_exch_ic(:,:,:,ii) = 0.020725*por;  % Doug's xls "Sorption Notes" mol/L_w, email 4/12/13
end
n_catex_exch = ii;
catex_exch_comp = catex_exch_comp(1:n_catex_exch);
catex_exch_ic = catex_exch_ic(:,:,:,1:n_catex_exch); 

% - catex components (catex is exchanger sorbed compound)
% (ic from PHREEQC batch run)
n_catex_max = 10;
catex_comp = cell(n_catex_max,1);
catex_ic = zeros(nrow,ncol,nlay,n_catex_max); 
ii = 0;
if fl_Yexch
    ii=ii+1; catex_comp{ii} = 'HY'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'NaY'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'KY'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'AmmHY'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'CaY2'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'MgY2'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'FeY2'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'MnY2'; % ***********************************
%     ii=ii+1; catex_comp{ii} = 'SrY2'; % ***********************************
%     ii=ii+1; catex_comp{ii} = 'BaY2'; % ***********************************
end
n_catex = ii;
catex_comp = catex_comp(1:n_catex);
catex_ic = catex_ic(:,:,:,1:n_catex);

% - surface complexation components
% **** WARNING: not set up for surfaces w/
%               multiple sites (would need to enter only 1 surf_calc_type
%               per surface in _ph file)
n_surf_max = 10;
surf_comp = cell(n_surf_max,1);
% surf_ic: number of sites [mol/Lv] if not coupled to phase or reactant,
%          (number of sites per mol x porosity) [mol/mol] if coupled
surf_ic = zeros(nrow,ncol,nlay,n_surf_max); 
% surf_par(1,ii): SurfArea ([m2/g] if not coupled to phase or reactant,
%              [m2/mol] if coupled to phase or reactant)
% surf_par(2,ii): mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_par = zeros(2,n_surf_max); 
% surf_cpl{1,ii}: Empty if not coupled to phase or reactant
% surf_cpl{1,ii}: Name of pure phase or kin reactant coupled to
% surf_cpl{2,ii}: 'equilibrium_phase' to couple to pure phase, or
%                    'kinetic_reactant' to couple to kinetic reactant
surf_cpl = cell(2,n_surf_max); 
surf_calc_type = ''; % must have same type for all surfaces, 3 options: '', '-no_edl', 'diffuse_layer'
ii = 0;
% 6/6/16: 
% Arsenic_Lakes_parameter_calcs_20160528_dbk_corrected_gcngnotes.xlsx
% ("Sorption" tab)
ii=ii+1; surf_comp{ii} = 'Hfo_wOH'; % ***********************************
surf_ic(:,:,:,ii) = 0.2*por;  % number of sites, [mol/mol]*por if coupled
surf_par(1,ii) = 1243.5/0.004665;  % SurfArea ([m2/mol] if coupled to phase or reactant)
surf_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_cpl{1,ii} = 'Fe(OH)3(a)';
surf_cpl{2,ii} = 'equilibrium_phase';
ii=ii+1; surf_comp{ii} = 'Hfo_sOH'; % ***********************************
surf_ic(:,:,:,ii) = 0.005*por;  % number of sites, [mol/mol]*por if coupled
surf_par(1,ii) = 1243.5/0.00012;  % SurfArea ([m2/mol] if coupled to phase or reactant)
surf_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_cpl{1,ii} = 'Fe(OH)3(a)';
surf_cpl{2,ii} = 'equilibrium_phase';  

n_surf = ii;
surf_comp = surf_comp(1:n_surf);
surf_ic = surf_ic(:,:,:,1:n_surf); 
surf_par = surf_par(:,1:n_surf); 

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

n_comp_all = [n_mob_kin; n_mob_eq; n_imob_kin; n_min_eq; n_catex; n_surf];

% -- (done setting model chemical components)


% - select output list
phrq_sel_outlist.total = [mob_kin_comp; mob_eq_comp(~strncmp('p',mob_eq_comp,1)); imob_kin_comp];
phrq_sel_outlist.mol = [catex_comp; surf_comp];
phrq_sel_outlist.equilphase = min_eq_comp;
phrq_sel_outlist.si = []; phrq_sel_outlist.gas = [];    

sel_outlist.total = [mob_kin_comp; mob_eq_comp(~strncmp('p',mob_eq_comp,1)); imob_kin_comp];
sel_outlist.mol = catex_comp;  % better to put desired surf species in add_sel_outlist.mol
% sel_outlist.mol = [catex_comp; surf_comp];
sel_outlist.equilphase = min_eq_comp;
% sel_outlist.si = {'FeS(ppt)'}; 
sel_outlist.si = {'Fe(OH)3(a)'};
% sel_outlist.si = {'Pyrite'};
% sel_outlist.si = []; 
sel_outlist.gas = []; 
% sel_outlist.si = {'CH4(g)', 'CO2(g)', 'Ntwo(g)', 'O2(g)'}; 
% sel_outlist.gas = {'CH4(g)', 'CO2(g)', 'Ntwo(g)', 'O2(g)'}; 

% iunclude any additional components
sel_outlist.total = [sel_outlist.total; addl_sel_outlist.total(:)];
sel_outlist.mol = [sel_outlist.mol; addl_sel_outlist.mol(:)];
sel_outlist.equilphase = [sel_outlist.equilphase; addl_sel_outlist.equilphase(:)];
sel_outlist.si = [sel_outlist.si; addl_sel_outlist.si(:)]; 
sel_outlist.gas = [sel_outlist.gas; addl_sel_outlist.gas(:)]; 


% -- Get initial conditions: equilibrated and charge-balanced solutions
% based on observations:

%         [mob_eq_comp_ic, mob_eq_sw_ic, mob_eq_pw_ic, catex_comp_ic, catex_pw_ic, surf_sp_comp_ic, surf_sp_pw_ic] = ...
%             As_InitCond_chem_3_gcng_AJD_combined_8N14E_Lake(sim_dir, phrq_exe, tempC, fl_Yexch);
%         [mob_eq_comp_ic, mob_eq_sw_ic, mob_eq_pw_ic, catex_comp_ic, catex_pw_ic, surf_sp_comp_ic, surf_sp_pw_ic] = ...
%             As_InitCond_chem_3_gcng_AJD_combined(sim_dir, phrq_exe, tempC_eq, fl_Yexch);
        [mob_eq_comp_ic, mob_eq_sw_ic, mob_eq_pw_ic, catex_comp_ic, catex_pw_ic, surf_sp_comp_ic, surf_sp_pw_ic] = ...
           As_InitCond_chem_4_testmodel_hi_lo_As_test4(sim_dir, phrq_exe, tempC_eq, fl_Yexch);

        for ii = 1: n_mob_eq
            ind = find(strcmp(mob_eq_comp(ii), mob_eq_comp_ic));
        %     mob_eq_const_rech(ii) = mob_eq_sw_ic(ind);
            mob_eq_ic(:,:,1:100,ii) = mob_eq_sw_ic(ind);  % BC
                    % mob_eq_ic(:,:,1:38,ii) = mob_eq_sw_ic(ind); % BC
                    % mob_eq_ic(:,:,43:44,ii) = mob_eq_sw_ic(ind);
                    % mob_eq_ic(:,:,63:64,ii) = mob_eq_sw_ic(ind);
                    % mob_eq_ic(:,:,67:80,ii) = mob_eq_sw_ic(ind);
                    % mob_eq_ic(:,:,87:120,ii) = mob_eq_sw_ic(ind);
        %     mob_eq_ic(:,:,2:end,ii) = mob_eq_sw_ic(ind); % **** CHANGE THIS BACK TO PW!!! ******
        %     mob_eq_ic(:,:,1,ii) = mob_eq_pw_ic(ind); % **** CHANGE THIS BACK TO SW!!! ******
            mob_eq_ic(:,:,101:end,ii) = mob_eq_pw_ic(ind);    
                    % mob_eq_ic(:,:,39:42,ii) = mob_eq_pw_ic(ind);
                    % mob_eq_ic(:,:,45:62,ii) = mob_eq_pw_ic(ind);
                    % mob_eq_ic(:,:,65:66,ii) = mob_eq_pw_ic(ind);
                    % mob_eq_ic(:,:,81:86,ii) = mob_eq_pw_ic(ind);
        end
        for ii = 1: n_catex  % cation exchange components (sorption)
            % - Surface boundary: no cation exchangers in sw
            % - Middle bulk of domain
            ind = find(strcmp(catex_comp(ii), catex_comp_ic));
            catex_ic(:,:,2:end,ii) = catex_pw_ic(ind); 
        %     % - Bottom boundary
        %     if ~isempty(ind_pwgw)
        %         catex_ic(:,:,ind_pwgw:end,ii) = catex_ic_z(3,ii);                  
        end
else
    % for ctr > 1, only change Orgc parameters
    ii = strcmp(mob_kin_comp, 'Orgc');
    mob_kin_par(1:length(Orgc_log10K), ii) = 10.^Orgc_log10K;
    
end
             
% create _ph file (this added to handle alkalinity)
fl_catex_toequil = zeros(n_catex); % 0 b/c already equilibrated exchanger
PHT3D_ph_f_general6(ph_file, ...
    OperSplit, tempC_eq, ...
    mob_kin_comp, mob_kin_par, mob_kin_formula, ...
    mob_eq_comp, mob_eq_extra, ...
    imob_kin_comp, imob_kin_par, imob_kin_formula, ...
    min_eq_comp, ...
    catex_comp, fl_catex_toequil, ...
    surf_comp, surf_par, surf_cpl, surf_calc_type, ...
    min_kin_comp, min_kin_par)

% create btn, ssm, and disp files 
fl_rech = 0;
eff_por = ones(nrow,ncol,nlay) * por;
PHT3D_btn_f_160608a(btn_file, ssm_file, dsp_file, ...
    nrow, ncol, nlay, DELR, DELC, htop, dz, eff_por, ...
    long_disp, hor2longdisp, vert2longdisp, Dstar, ...
    mob_kin_comp, mob_kin_ic, ...
    mob_eq_comp, mob_eq_ic, fl_mob_eq_const_rech, mob_eq_const_rech, mob_eq_distr_rech, ...
    imob_kin_comp, imob_kin_ic, ...
    min_eq_comp, min_eq_ic, ...
    catex_comp, catex_ic, ...
    surf_comp, surf_ic, ...
    timprs, nstp, fl_rech);    

% create postfix file (what to output to .sel file)
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

% -- Check if model crashed!!
phout_fil = [sim_dir, slashstr, 'phout.dat'];
fid = fopen(phout_fil);
while(~feof(fid))
    str = fgets(fid);
end
fclose(fid);
if strncmp(str, 'Stopping.', 9);
    % model did not complete
    fprintf('Model did not complete, redo with different charge balance \n');
    fl_ReDo = 1;
    Na_ind = find(strcmp(mob_eq_comp, 'Na'),1);
    Cl_ind = find(strcmp(mob_eq_comp, 'Cl'),1);
    ind = find(strcmp(mob_eq_extra, 'charge'));
    if ind == Na_ind
        mob_eq_extra{Na_ind} = '';
        mob_eq_extra{Cl_ind} = 'charge';
    elseif ind == Cl_ind
        mob_eq_extra{Cl_ind} = '';
        mob_eq_extra{Na_ind} = 'charge';        
    else
        fprintf('  ...but charge component is neither Cl nor Na. Stopping...\n');
        return
    end
else
    
    % did not crash, continue...
    fl_ReDo = 0;


           %**********OBTAINING NEW INIT COND***************            

    % n_comp_all = [n_mob_kin; n_mob_eq; n_imob_kin; n_min_eq; n_catex]; % old version

    comp_ctr=0;

        for cc = 1: 6
            dist_ic = zeros(nrow,ncol,nlay,n_comp_all(cc));
            switch cc
                case 1
                    comp_name_i = mob_kin_comp;
                case 2
                    comp_name_i = mob_eq_comp;
                case 3
                    comp_name_i = imob_kin_comp;
                case 4
                    comp_name_i = min_eq_comp;
                case 5
                    comp_name_i = catex_comp;
                case 6 
                    comp_name_i = surf_comp;

            end 

                    for ss = 1: n_comp_all(cc)
                comp_ctr = comp_ctr + 1;

                if comp_ctr >= 100
                    numstr = num2str(comp_ctr);
                elseif comp_ctr >= 10
                    numstr = [ '0', num2str(comp_ctr)];
                else
                    numstr = [ '00', num2str(comp_ctr)];
                end
                fil = ['PHT3D', numstr, '.ACN'];
    %             [status, result] = system(['dir C:\Users\Owner\Desktop\Arsenic_MN\Test_Models\As_into_Model_20160611\', fil]); 
    %             if status ~= 0
    %                 fprintf('No .ACN file found for %s!  Exiting... \n', comp_name_i{ss});
    %                 return
    %             end
                x = load(fil);
                tmax = length(x)/(ncol*nrow*nlay);
                data = reshape(x, ncol, nrow, nlay, tmax);  % (ncol,nrow,nlay,ntime);
                %take final time
                data = data(:,:,:,end);            

                dist_ic(:,:,:,ss) = permute(data, [2 1 3]); % (nrow,ncol,nlay,n_comp);
    %             dist_ic(:,:,1,ss) = dist_ic(:,:,2,ss);
                    end

            switch cc
                case 1
                    mob_kin_ic = dist_ic;
                case 2
                    mob_eq_ic = dist_ic;
                case 3
                    imob_kin_ic = dist_ic;
                case 4
                    min_eq_ic = dist_ic;
                case 5
                    catex_ic = dist_ic;
                case 6
                    surf_ic = dist_ic;
            end
        end

    newname=['out_X',num2str(ctr),'.sel'];

    copyfile('out_X.sel', newname)

        newname1 = ['pht3dbtn_1D_', num2str(ctr), '.dat'];

        copyfile('pht3dbtn_1D.dat', newname1)

            newname2=['pht3d_ph_', num2str(ctr), '.dat'];

            copyfile('pht3d_ph.dat', newname2)
    if fl_gcng, system(['cp -p out.txt out_', num2str(ctr), '.txt']); end


        %*************CONCATENATING .SEL FILES******************    

    fclose('all');    

    %  % - read in data
    outfile = ['out_X',num2str(ctr),'.sel'];

    % formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';

    ndatalines = length(find(timprs>0)) * (nlay-1); % does not write time=0 nor constant conc BC cells

    if ctr==1

    fid = fopen(outfile);    
    %     line0 = fscanf(fid, '%c', [1 462]); %reads in column titles
        line0 = fgets(fid); %reads in first line with column titles
        col_ti = textscan(line0, '%s'); 
        col_ti = col_ti{1};
        n_col_ti = length(col_ti);
        formatSpec = repmat('%f ', 1, n_col_ti);
        data1 = cell2mat(textscan(fid, formatSpec, ndatalines)); %reads in data
    fclose(fid);

    fid1 = fopen('out_Xfinal.sel', 'w');  
        fprintf(fid1, '%s', line0);    
    %     fprintf(fid1, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n\r', data1');
        fmtstr = repmat('%12g ', 1, n_col_ti);
        fprintf(fid1, [fmtstr, '\n'], data1');
        fclose(fid1);


    elseif ctr >1


       %Array used to make time sequential for 50 layer model run for 180 days
    %    tim = ones(294,1);
    %    ntim = 15552000.*tim.*(ctr-1);

       %Array used to make time sequential for 200 layer model run for 180 days
       tim = ones(ndatalines,1);
    %    ntim = 15552000.*tim.*(ctr-1);   
       ntim = (timprs(end)*3600*24).*tim.*(ctr-1);   

    fid = fopen(outfile);
       data2 = cell2mat(textscan(fid, formatSpec, ndatalines, 'headerLines', 1)); %reads in data
       data2 = [data2(1:end,1)+ntim, data2(1:end,2:end)]; %Reformats matrix with the correct time
    fclose(fid);

    fid1 = fopen('out_Xfinal.sel', 'a'); 
    %     fprintf(fid1, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n\r', data2');
        fprintf(fid1, [fmtstr, '\n'], data2');
        fclose(fid1);   
    end
end            
if ctr==end_ctr
    return
end 
end
 
fclose('all');

cd(curr_dir);

rmpath(matlab_dir);
