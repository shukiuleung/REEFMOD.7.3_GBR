% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Feb 2026.
%
% Converts the predict CoTS density per reef grid into number of CoTS per tow.
% This reverses the calibration model of Moran & De'ath (1992) AUST MAR FRESHW
% implemented in f_convert_CoTS_tow_2_density.
% CoTS density is the summed density across all subadults and adults, after application of the detectability vector 
% (lower detectability for smaller CoTS).
% Major difference here is that the model is deterministic
% Model and data described in COTS_MTC_calibration.R available in:
%
% https://github.com/ymbozec/acanthaster_GBR_abundance/
% -------------------------------------------------------------------------

function MTC_per_tow = f_convert_CoTS_density_2_tow(COTS_density, total_area_cm2)

% The original model predicts SCUBA swim counts (SSC) performed on an area equivalent to a tow sample (~200m*12m = 2400m2)
% from manta tow counts (MTC) expressed in a per tow basis (ie, as averaged across a reef). The number of individual tows
% performed around the perimeter of a reef was used as weight in the regression model, to infer stochastic predictions.
% We simply reverse the relationship but do not apply stochastic predictions as this was not part of the original model
% (ie, SSC as response, MTC as predictor) 

% First adjust density to the original area, which was equivalent to a tow sample (~200m*12m = 2400m2)
SSC_per_tow = 2400.*COTS_density./(total_area_cm2/10000);

% Then apply the reversed regression model
% Note the model was build after a cubic root transformation of the predictor (MTC) and response (SSC)
MTC_per_tow = ((1/1.2008).*((SSC_per_tow.^(1/3)) - 0.8071)).^3 ;

% MTC_per_tow(MTC_per_tow<0)=0; % Force to 0 the negatives
MTC_per_tow(MTC_per_tow<0.001)=0; % Force to 0 values below 0.001 for optimisation