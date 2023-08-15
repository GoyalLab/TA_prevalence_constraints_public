% original version written by Lea Schuh and Yogesh Goyal, 2020
% Schuh et al., Cell Systems, 2020

%function to produce Latin hypercube sampled parameter sets in all seven
%independent parameters (r_prod, r_deg, r_on, n, r_add, r_off, d) and
%create the dependent parameter k (where k = 0.95*r_prod/r_deg*d)
%
%INPUT:
%
%nruns:         number of parameter sets sampled
%n_species:     number of parameter multiples (eg. if 3 node network and
%               one requires a different parameter set per node: n_species = 3; 
%               here: n_species = 1, as we assume the parameter sets to be 
%               the same for all nodes in the network)
%type:          normal - normal parameter space; constrained - constrained
%               parameter space
%
%OUTPUT:
%
%D:             nruns * n_species*8 - matrix of nruns different LHS sampled
%               parameter sets for n_species

%%
function D = LHSsampling(nruns, n_species, type)

rng('shuffle');                             %receive different values upon new start

if isequal(type,'normal') == 1
    %for general Data1000
    min_range = [repmat(0.01,1,1),...       %production rate
        repmat(0.001,1,1),...               %degradation rate
        repmat(0.001,1,1),...               %basal on-rate of burst
        repmat(0.1,1,1),...                 %Hill coefficient n (Hill function)
        repmat(0.1,1,1),...                 %additional on-rate due to dependency to other node
        repmat(0.01,1,1),...                %off rate of burst
        repmat(2,1,1)];                     %proddiff
    
    max_range = [repmat(1,1,1),...      %production rate
        repmat(0.1,1,1),...             %degradation rate
        repmat(0.1,1,1),...             %basal on-rate of burst
        repmat(10,1,1),...              %Hill coefficient n (Hill function)
        repmat(1,1,1),...               %additional on-rate due to dependency to other node
        repmat(0.1,1,1),...             %off rate of burst
        repmat(100,1,1)];               %proddiff
    
elseif isequal(type,'constrained') == 1
    %for constrained parameter space
    min_range = [repmat(0.01,1,1),...       %production rate
        repmat(0.001,1,1),...               %degradation rate
        repmat(0.001,1,1),...               %basal on-rate of burst
        repmat(0.1,1,1),...                 %Hill coefficient n (Hill function)
        repmat(0.15,1,1),...                %additional on-rate due to dependency to other node
        repmat(0.06,1,1),...                %off rate of burst
        repmat(2,1,1)];                     %proddiff
    
    max_range = [repmat(1,1,1),...      %production rate
        repmat(0.1,1,1),...             %degradation rate
        repmat(0.025,1,1),...           %basal on-rate of burst
        repmat(10,1,1),...              %Hill coefficient n (Hill function)
        repmat(0.36,1,1),...            %additional on-rate due to dependency to other node
        repmat(0.1,1,1),...             %off rate of burst
        repmat(100,1,1)];               %proddiff
end

latinhyp = lhsdesign_modified(nruns, min_range, max_range);

latinhyp_prod = repmat(latinhyp(:,1),1,n_species);
latinhyp_deg = repmat(latinhyp(:,2),1,n_species);
latinhyp_onbasal = repmat(latinhyp(:,3),1,n_species);
latinhyp_n = repmat(latinhyp(:,4),1,n_species);
latinhyp_ondep = repmat(latinhyp(:,5),1,n_species);
latinhyp_off = repmat(latinhyp(:,6),1,n_species);
latinhyp_proddiff = repmat(latinhyp(:,7),1,n_species);

k = repmat(0.95*latinhyp(:,1)./latinhyp(:,2).*latinhyp(:,7),1,n_species); %0.95 of stationary on-state of system

D = [latinhyp_prod, latinhyp_deg, latinhyp_onbasal, latinhyp_n, latinhyp_ondep, latinhyp_off, k, latinhyp_proddiff];

end