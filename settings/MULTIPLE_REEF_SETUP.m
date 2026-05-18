%__________________________________________________________________________
%
% REEFMOD-GBR - populate reef parameters and initial state 
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 09/2019
%__________________________________________________________________________


%% Populate REEF parameters
for n = 1:length(META.reef_ID)
        
    REEF(n).initial_coral_cover = init_coral_cover(n,:);
    REEF(n).initial_algal_cover = init_algal_cover(n,:) ;
    REEF(n).nongrazable_substratum = init_sand_cover(n,1) ;
    REEF(n).initial_rubble_pct = 100*init_rubble_cover(n,1) ;
    
    % Default values (allows for reef-specific modifications in runmodel)
    REEF(n).juv_whole_mortality_rate = CORAL.juv_whole_mortality_rate;
    REEF(n).adol_whole_mortality_rate = CORAL.adol_whole_mortality_rate;
    % REEF(n).adult_whole_mortality_rate = CORAL.adult_whole_mortality_rate; % doesn't need this anymore as mortality
    % decreases from subadult mortality with decay rate of 3
    
    % Default values for herbivory
    REEF(n).diadema = REEF(1).diadema ;
    REEF(n).herbivory = REEF(1).herbivory ;
    REEF(n).dictyota_declines_seasonally = REEF(1).dictyota_declines_seasonally ;
    
    % Store the bleaching scenario
    if META.doing_genetics == 0
        
        REEF(n).predicted_DHWs = squeeze(DHW(n,1:META.nb_time_steps))+ single(META.doing_cooling)*META.cooling_factor*12;
        REEF(n).Topt_baseline =[];
        REEF(n).predicted_SST = [];
        
    else
        REEF(n).predicted_DHWs = squeeze(DHW(n,1:META.nb_time_steps,:))'+ single(META.doing_cooling)*META.cooling_factor*12; % with all the Topt stacks
        REEF(n).Topt_baseline = MMM_CoRTAD(META.reef_ID(n),1) - META.Topt_offset ;
        REEF(n).predicted_SST = SST_GBR(n,:) + single(META.doing_cooling)*META.cooling_factor/2;  % Decreases SST for only 3 months in summer

        for s=1:META.nb_coral_types
            
            % select genotypes from the bank and add random random variation (creates genetic diversity among species)
            REEF(n).coral(s).QTL_pool_IN = squeeze(QTL1000(META.reef_ID(n),1:META.genetics_pop_size,:,:)) + ...
                    normrnd(META.initial_push_QTL_mu(s), META.initial_push_QTL_sd(s),META.genetics_pop_size, META.genetics.nb_loci,2);            
        end
    end
    
    % Allocation to record the pool of ongoing/incoming HT phenotypes (same principle as QTL pool in genetic model)
    REEF(n).HT_pool_OUT_TMP = zeros(META.HT_pop_size,META.nb_coral_types,'single'); % temporary structure to store HT as reefs are processed iteratively
    REEF(n).HT_pool_OUT = REEF(n).HT_pool_OUT_TMP;
    REEF(n).HT_pool_IN = REEF(n).HT_pool_OUT_TMP;
    
    % Store the scenario of cyclones
    REEF(n).hurricane_chronology = CYCLONE_CAT(n,1:META.nb_time_steps);
    
    if META.randomize_hurricane_strike == 1 % randomise the incidence of a predicted cylone (based on orbital velocity thresholds):  
        REEF(n).hurricane_strike_incidence = STRIKE_INCIDENCE(n); %if 1, a cyclone predicted on a reef will always strike
    else
        REEF(n).hurricane_strike_incidence = META.hurricane_strike_incidence; % otherwise, use the default value defined in PARAMETERS
    end
    
    % Assign to every reef the defaut background density of CoTS
    if META.doing_COTS == 1
        REEF(n).COTS_background_density = META.COTS_background_density;
        
        % TINA: Assign specific mortality rate based on reef protection status
        if GBR_REEFS.GreenZone(n) == 1 %GreenZone
            REEF(n).COTS_mortality = META.COTS_mortality;
            
        else %Not GreenZone
            REEF(n).COTS_mortality = META.COTS_mortality_spec;
        end
    end
end
