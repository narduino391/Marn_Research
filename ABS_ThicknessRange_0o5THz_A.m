% For un-polarized light:

% Housekeeping
clear
clc

% Principal equation:
% rtotal = (r1 + r2*e^((-4*1i*pi*n1*d1)/wavelength)) / (1 + r1*r2*e^((-4*1i*pi*n1*d1)/wavelength)

% ABS (Acrylonitrile butadiene styrene) - 1.57
% Source of RI values: "Busch, S.F., Weidenbach, M., Fey, M. et al. Optical Properties of 3D Printable Plastics in the THz Regime and their Application for 3D Printed THz Optics. J Infrared Milli Terahz Waves 35, 993–997 (2014). https://doi.org/10.1007/s10762-014-0113-9"

x = 'ABS';
n1 = 1.57;

% Input to define constant indexes of refraction
n0 = 1; %(Consistant value of n0=1)
clc

% Find values used in thickness array
dmin = 0.00335; %(Defines minimum value plotted)
dmax = 0.003525; %(Defines maximum value plotted)
dr = dmax-dmin;  %(Determines range of thickness being computed)
di = dr/1000; %(Determines increments used in thickness range)

% Find values used in refractive index array
n21min = 1.6; %(Defines minimum value plotted)
n21max = 3.1; %(Defines maximum value plotted)
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

%surf(darray,narray,DeltaRmatrix,EdgeColor="none")   % Generates surface plot of change in reflectance in terms of both thickness and skin refractive index
%xlabel('Thickness')
%ylabel('Refractive Index')
%zlabel('Change in Reflectance')
%title("ΔR, " + x + "-IW @0.5THz (Cons.)")   % Inserts name of imaging window material into title of plot
%alpha(.8)   % Sets level of transparency of plot
%colorbar    % Includes legend for color of plot and correlating change in reflectance

[maxDeltaRmatrix,maxDeltaR_thickness_index] = max(DeltaRmatrix,[],2);    % Finds maximum change in reflectance values at each skin refractive index value
maxDeltaR_thickness_index = mean(maxDeltaR_thickness_index);    % Finds thickness index of peak change in reflectance (finds mean of identical values iterated in maxDeltaR_thickness_index)
maxDeltaR_thickness = darray(maxDeltaR_thickness_index);    % Finds thickness value at given thickness index associated with peak change in reflectance
requiredDeltaRmatrix = 0.9*maxDeltaRmatrix; % At each skin refractive index, finds minimum required change in reflectance value in order to stay within 10% of peak 
tf_acceptable_DeltaR = DeltaRmatrix > requiredDeltaRmatrix;  % Generates logic matrix that checks change in reflectance at all thicknesses for each reractive indices to see if it meets the required value
darray_matrix = meshgrid(darray,narray);    % Generates matrix with the same dimensions as DeltaRmatrix, containing repeating rows of the thickness values of darray
acceptable_thickness = darray_matrix .* tf_acceptable_DeltaR;   % Creates matrix containing thicknesses at which change in reflectance is acceptable and displays zeros at thicknesses where it is not
acceptable_DeltaRmatrix = DeltaRmatrix .* tf_acceptable_DeltaR; % Creates matrix containing change in reflectance values that meet the required value and zeros when it does not

tmp = acceptable_thickness; % Creates temporary copy of acceptable_thickness matrix for slight manipulation
tmp(tmp==0) = Inf;  % Sets all 0 values to infinity
minimum_thickness = min(tmp,[],2);  % Finds minimum thickness value of each row (skin refractive index) of the tmp(acceptable_thickness) matrix

maximum_thickness = max(acceptable_thickness,[],2); % Finds maximum thickness value of each row (skin refractive index) of the acceptable_thickness matrix)

thickness_range = maximum_thickness-minimum_thickness;  % Finds range of acceptable thicknesses from each row (skin refractive index)
minimum_thickness_range = min(thickness_range); % Finds minimum range of thickness for all skin refractive indices
minimum_thickness_variation = 0.5*minimum_thickness_range;  % Cuts range in half to find accepted variance from peak

disp(x) % Displays material name for verification
fprintf('For peak change in reflectance, the target thickness is %.10fm.\n',maxDeltaR_thickness)    % Displays target thickness for peak change in reflectance
fprintf('The lowest acceptable variation from the target thickness is +/-%.10fm.\n',minimum_thickness_range)   % Displays lowest acceptable thickness variation


% Likely will want to create separate script for each material analyzed

% Creating universal programmatic method will likely be time consuming
% (and maybe a little too difficult for my current skill level)

% Select single tunnel (Check)
% Find maximum DeltaR across RI and thickness (Check)
% Find allowable (Check)
% Check using logic matrix (Check)
% Output values of 2-d darray by plugging in logic magtrix (Check)
% Plugging logic matrix into original DeltaRmatrix should give matrix of only acceptable deltaR's
% Can use this to create surface plot and find range of viable thicknesses
% May have issues with dimensions of matrices, but can confer with Dr. Marn at that point