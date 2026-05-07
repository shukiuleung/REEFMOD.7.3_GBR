% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Feb 2026.
%
% Converts the number of CoTS per tow into a density estimate per reef grid.
% This operationalises the regression model of Moran & De'ath (1992) AUST MAR FRESHW
% of SCUBA swim counts (SSC) vs manta tow counts (MTC) and adds uncertainty around 
% CoTS per tow predictions (consistent with the model and data) based on the number of tows performed.
% Model and data described in COTS_MTC_calibration.R available in:
%
% https://github.com/ymbozec/acanthaster_GBR_abundance/
% -------------------------------------------------------------------------

function COTS_density = f_convert_CoTS_tow_2_density(MTC_per_tow, NB_tows, total_area_cm2)

% The model predicts SCUBA swim counts (SSC) performed on an area equivalent to a tow sample (~200m*12m = 2400m2)
% from manta tow counts (MTC) expressed in a per tow basis (ie, as averaged across a reef). The number of individual tows
% performed around the perimeter of a reef was used as weight in the regression model, so is used here to infer stochastic
% predictions. If NB_tows = 0, the conversion is deterministic.
% Note the model was build after a cubic root transformation of the predictor (MTC) and response (SSC).

se_fit = 0.12 ; %standard error for the mean prediction at given predictor values
% Here it is fixed to an average 0.12, although it varies non linearly between 0.10 and 0.15 in the range 0-20 CoTS per tow,
% (up to 0.28 for 100 CoTS per tow).
res_scale = 1.223346 ; % model's estimate of the residual standard deviation

% MTC_per_tow = 0:100; % for testing
% NB_tows = 20; % for testing

SSC_per_tow_cubic_root = 1.2008.*MTC_per_tow.^(1/3) + 0.8071; % deterministic prediction

% If stochastic prediction (NB_tows must be non-null)
if NB_tows > 0
    SSC_per_tow_cubic_root = normrnd(SSC_per_tow_cubic_root, sqrt(se_fit.^2 + res_scale.^2)./sqrt(NB_tows)); % add noise dependent on number of tows
end

% Output is a random estimate of absolute density of noncryptic COTS (cubic root) for a 2400m2 equivalent area.
% Needs to back transformed and scaled to the area of a reef grid
COTS_density = (total_area_cm2/10000).*(SSC_per_tow_cubic_root.^3)/2400; % (the first term is to express grid area in m2)
