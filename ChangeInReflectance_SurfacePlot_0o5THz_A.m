% For un-polarized light:

% Housekeeping
clear
clc

% Principal equation:
% rtotal = (r1 + r2*e^((-4*1i*pi*n1*d1)/wavelength)) / (1 + r1*r2*e^((-4*1i*pi*n1*d1)/wavelength)

% List of refractive index values:

% ABS (Acrylonitrile butadiene styrene) - 1.57
% Bendlay - 1.532
% HDPE (High density polyethylene) - 1.532
% Nylon - 1.72
% Polystyrene - 1.56
% PP (Polypropylene) - 1.495
% PLA (Polylactic acid) - 1.89
% HRFZ-Si

% Source of RI values: "Busch, S.F., Weidenbach, M., Fey, M. et al. Optical Properties of 3D Printable Plastics in the THz Regime and their Application for 3D Printed THz Optics. J Infrared Milli Terahz Waves 35, 993–997 (2014). https://doi.org/10.1007/s10762-014-0113-9"

list={'ABS','Bendlay','HDPE','Nylon','Polystyrene','PP','PLA','HRFZ-Si'}; % Defines list of avaialable imaging window materials
did_select=0;
while ~did_select   % Creates loop to ensure that user selects an imaging window material
    [material,did_select] = listdlg('PromptString', ... % Prompts user with dialogue box
    {'Select an imaging window material:'}, ... % Title of dialogue box
    'SelectionMode','single', ...   % Allows only one single selection from given options
    'ListString',list); % Sets list of options displayed
end

switch material % Looks to match index value of material chosen
    case 1  % Value corresponds to order of material options given
        x='ABS';    % Defines variable to be used in title of plot
        n1=1.57;    % Defines imaging window refractive index based on material selected
    case 2
        x='Bendlay';
        n1=1.532;
    case 3
        x='HDPE';
        n1=1.532;
    case 4
        x='Nylon';
        n1=1.72;
    case 5
        x='Polystyrene';
        n1=1.56;
    case 6
        x='PP';
        n1=1.495;
    case 7
        x='PLA';
        n1=1.89;
    case 8
        x='HRFZ-Si';
        n1=1.54;
end

% Input to define constant indexes of refraction
n0 = 1; %(Consistant value of n0=1)
clc

% Find values used in thickness array
dmin = 0.004; %(Defines minimum value plotted)
dmax = 0.005; %(Defines maximum value plotted)
dr = dmax-dmin;  %(Determines range of thickness being computed)
di = dr/1000; %(Determines increments used in thickness range)

% Find values used in refractive index array
n21min = 1.5; %(Defines minimum value plotted)
n21max = 3.0; %(Defines maximum value plotted)
nr = n21max-n21min;  %(Determines range of thickness being computed)
ni = nr/1000; %(Determines increments used in thickness range)

% Input to define wavelength of emitted radiation
wl = 0.00059958; %(1THz frequency used equates to 0.00029979m wavelength)
clc

% Input to define angle of incidence
Oi1 = 0;    % Assumes normal incidence

% Creates array of thickness and refractive index values values
darray = dmin:di:dmax;  %(Creates array of all plotted thickness values)
narray = n21min:ni:n21max;  %(Creates array of all plotted thickness values)

% Preallocates empty matrix for change in reflectance values
DeltaRmatrix = zeros(length(darray):length(narray));

% Begings loop of calculations for all values of thickness and refractive index being being tested
for j = 1:length(narray)   

    for i = 1:length(darray)

% Determines current thickness
d = darray(i);

% Determines current skin refractive indices
n21 = narray(j);
n22 = n21+0.01;

% Calculate angle of transmittance
Ot1 = asin(n0*sin(Oi1)/n1);
Oi2 = Ot1;
Ot21 = asin(n1*sin(Oi2)/n21);
Ot22 = asin(n1*sin(Oi2)/n22);

% Determine reflection of each interface for s-polarization
rs1 = (n0*cos(Oi1)-n1*cos(Ot1))/(n0*cos(Oi1)+n1*cos(Ot1));
rs21 = (n1*cos(Oi2)-n21*cos(Ot21))/(n1*cos(Oi2)+n21*cos(Ot21));
rs22 = (n1*cos(Oi2)-n22*cos(Ot22))/(n1*cos(Oi2)+n22*cos(Ot22));

% Determine reflection of each interface for p-polarization
rp1 = (n0*cos(Ot1)-n1*cos(Oi1))/(n0*cos(Ot1)+n1*cos(Oi1));
rp21 = (n1*cos(Ot21)-n21*cos(Oi2))/(n1*cos(Ot21)+n21*cos(Oi2));
rp22 = (n1*cos(Ot22)-n22*cos(Oi2))/(n1*cos(Ot22)+n22*cos(Oi2));

% Calculate overall reflection of s-polarization
rst1 = (rs1+rs21*exp((-4*1i*pi*n1*d*cos(Oi1))/wl))/(1+rs1*rs21*exp((-4*1i*pi*n1*d*cos(Oi1))/wl));
rst2 = (rs1+rs22*exp((-4*1i*pi*n1*d*cos(Oi2))/wl))/(1+rs1*rs22*exp((-4*1i*pi*n1*d*cos(Oi2))/wl));

% Calculate overall reflection of p-polarization
rpt1 = (rp1+rp21*exp((-4*1i*pi*n1*d*cos(Oi2))/wl))/(1+rp1*rp21*exp((-4*1i*pi*n1*d*cos(Oi2))/wl));
rpt2 = (rp1+rp22*exp((-4*1i*pi*n1*d*cos(Oi2))/wl))/(1+rp1*rp22*exp((-4*1i*pi*n1*d*cos(Oi2))/wl));

% Calculate overall total reflectance of system
Rs1 = rst1^2;
Rs2 = rst2^2;
Rp1 = rpt1^2;
Rp2 = rpt2^2;
R1 = (Rs1+Rp1)/2;
R2 = (Rs2+Rp2)/2;
DeltaR = (abs(R1-R2));

DeltaRmatrix(j,i) = DeltaR; % Creates matrix of change in reflectance from thickness of imaging window and refractive index of skin

    end
end

surf(darray,narray,DeltaRmatrix,EdgeColor="none")   % Generates surface plot of change in reflectance in terms of both thickness and skin refractive index
xlabel('Thickness')
ylabel('Refractive Index')
zlabel('Change in Reflectance')
title("ΔR, " + x + "-IW @0.5THz (Cons.)")   % Inserts name of imaging window material into title of plot
alpha(.8)   % Sets level of transparency of plot
colorbar    % Includes legend for color of plot and correlating change in reflectance




% FUTURE NOTES:

% DeltaRmatrix(end,:)

% Likely will want to create separate script for each material analyzed

% Creating universal programmatic method will likely be time consuming
% (and maybe a little too difficult for my current skill level)

% Select single tunnel
% Find maximum DeltaR across RI and thickness
% Find allowable
% Check using logic matrix
% Output values of 2-d darray by plugging in logic magtrix
% Plugging logic matrix into original DeltaRmatrix should give matrix of only acceptable deltaR's
% Can use this to create surface plot and find range of viable thicknesses
% May have issues with dimensions of matrices, but can confer with Dr. Marn at that point
