function dydt = odefun_mutmut_B4state(t,y,rates)
% diploid, mut/mut ODE, manually excludes paralog and nonsense

r_prodbasal_A1 = rates(1);
r_prodbasal_Anonsense1 = rates(2);
r_prodbasal_Aprime1 = rates(3);
r_prodbasal_Aprime2 = rates(4);
r_prodbasal_B1 = rates(5);

r_prodon_A1 = rates(6);
r_prodon_Anonsense1 = rates(7);
r_prodon_Aprime1 = rates(8);
r_prodon_Aprime2 = rates(9);
r_prodon_B1 = rates(10);

d_Aprime1_B1 = rates(11);
d_Aprime2_B1 = rates(12);

r_deg_A1 = rates(13);
r_deg_Anonsense1 = rates(14);
r_deg_Aprime1 = rates(15);
r_deg_Aprime2 = rates(16);
r_deg_B1 = rates(17);

r_onbasal_A1 = rates(18);
r_onbasal_Anonsense1 = rates(19);
r_onbasal_Aprime1 = rates(20);
r_onbasal_Aprime2 = rates(21);
r_onbasal_B1 = rates(22);

r_nitc_byAnonsense1_A1 = rates(23);
r_nitc_byAnonsense1_Anonsense1 = rates(24);
r_nitc_byAnonsense1_Aprime1 = rates(25);
r_nitc_byAnonsense1_Aprime2 = rates(26);
r_addon_byA1_B1 = rates(27);
r_addon_byAprime1_B1 = rates(28);
r_addon_byAprime2_B1 = rates(29);

r_off_A1 = rates(30);
r_off_Anonsense1 = rates(31);
r_off_Aprime1 = rates(32);
r_off_Aprime2 = rates(33);
r_offorig_B1 = rates(34);
r_offpara1_B1 = rates(35);
r_offpara2_B1 = rates(36);

k_A1  = rates(37);
k_Anonsense1  = rates(38);
k_Aprime1  = rates(39);
k_Aprime2  = rates(40);
k_B1  = rates(41);

n_A1  = rates(42);
n_Anonsense1  = rates(43);
n_Aprime1  = rates(44);
n_Aprime2  = rates(45);
n_B1  = rates(46);

% vector: y: 
% y(1) = ona_1 (wt)
% y(2) = offa_1 (wt)
% y(3) = ona_2 (wt)
% y(4) = offa_2 (wt)
% y(5) = a 
% y(6) = anons
% y(7) = onap1_1 
% y(8) = offap1_1 
% y(9) = onap1_2 
% y(10) = offap1_2 
% y(11) = onap2_1 
% y(12) = offap2_1 
% y(13) = onap2_2 
% y(14) = offap2_2 
% y(15) = aprim1
% y(16) = aprim2
% y(17) = onborig_1
% y(18) = onbpara1_1
% y(19) = onbpara2_1
% y(20) = offb_1
% y(21) = onborig_2
% y(22) = onbpara1_2
% y(23) = onbpara2_2
% y(24) = offb_2
% y(25) = b

dydt = zeros(25,1);
% y(1) = ona_1 (mut)
dydt(1) = y(2)*(r_onbasal_Anonsense1 + r_nitc_byAnonsense1_A1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    - r_off_Anonsense1*y(1);
% y(2) = offa_1 (mut)
dydt(2) = -y(2)*((r_onbasal_Anonsense1 + r_nitc_byAnonsense1_A1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) ...
    + r_off_Anonsense1*y(1);
% y(3) = ona_2 (mut)
dydt(3) = y(4)*(r_onbasal_Anonsense1 + r_nitc_byAnonsense1_Anonsense1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    - r_off_Anonsense1*y(3);
% y(4) = offa_2 (mut)
dydt(4) = -y(4)*(r_onbasal_Anonsense1 + r_nitc_byAnonsense1_Anonsense1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    + r_off_Anonsense1*y(3);
% y(5) = a 
dydt(5) = -r_deg_A1*y(5);
% y(6) = anons
dydt(6) = r_prodon_Anonsense1*y(1) + r_prodon_Anonsense1*y(3) - r_deg_Anonsense1*y(6);
% y(7) = onap1_1 
dydt(7) = y(8)*(r_onbasal_Aprime1 + r_nitc_byAnonsense1_Aprime1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    - r_off_Aprime1*y(7);
% y(8) = offap1_1 
dydt(8) = -y(8)*((r_onbasal_Aprime1 + r_nitc_byAnonsense1_Aprime1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) ...
    + r_off_Aprime1*y(7);
% y(9) = onap1_2 
dydt(9) = y(10)*(r_onbasal_Aprime1 + r_nitc_byAnonsense1_Aprime1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    - r_off_Aprime1*y(9);
% y(10) = offap1_2 
dydt(10) = -y(10)*((r_onbasal_Aprime1 + r_nitc_byAnonsense1_Aprime1*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) ...
    + r_off_Aprime1*y(9);
% y(11) = onap2_1 
dydt(11) = y(12)*(r_onbasal_Aprime2 + r_nitc_byAnonsense1_Aprime2*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    - r_off_Aprime2*y(11);
% y(12) = offap2_1 
dydt(12) = -y(12)*((r_onbasal_Aprime2 + r_nitc_byAnonsense1_Aprime2*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) ...
    + r_off_Aprime2*y(11);
% y(13) = onap2_2 
dydt(13) = y(14)*(r_onbasal_Aprime2 + r_nitc_byAnonsense1_Aprime2*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1)) ...
    - r_off_Aprime2*y(13);
% y(14) = offap2_2 
dydt(14) = -y(14)*((r_onbasal_Aprime2 + r_nitc_byAnonsense1_Aprime2*y(6)^n_Anonsense1/(k_Anonsense1^n_Anonsense1 + y(6)^n_Anonsense1))) ...
    + r_off_Aprime2*y(13);
% y(15) = aprim1 
dydt(15) = r_prodon_Aprime1*y(7) + r_prodon_Aprime1*y(9) - r_deg_Aprime1*y(15);
% y(16) = aprim1 
dydt(16) = r_prodon_Aprime1*y(11) + r_prodon_Aprime1*y(13) - r_deg_Aprime1*y(16);
% y(17) = onborig_1
dydt(17) = (r_onbasal_B1 + r_addon_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(20) - r_offorig_B1*y(17);
% y(18) = onbpara1_1
dydt(18) = (r_onbasal_B1 + r_addon_byAprime1_B1*y(15)^n_Aprime1/(k_Aprime1^n_Aprime1 + y(15)^n_Aprime1))*y(20) - r_offpara1_B1*y(18);
% y(19) = onbpara2_1
dydt(19) = (r_onbasal_B1 + r_addon_byAprime2_B1*y(16)^n_Aprime2/(k_Aprime2^n_Aprime2 + y(16)^n_Aprime2))*y(20) - r_offpara2_B1*y(19);
% y(20) = offb_1
dydt(20) = -y(20)*((r_onbasal_B1 + r_addon_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1)) + (r_onbasal_B1 + r_addon_byAprime1_B1*y(15)^n_Aprime1/(k_Aprime1^n_Aprime1 + y(15)^n_Aprime1)) + (r_onbasal_B1 + r_addon_byAprime2_B1*y(16)^n_Aprime2/(k_Aprime2^n_Aprime2 + y(16)^n_Aprime2))) ...
    + r_offorig_B1*y(17) + r_offpara1_B1*y(18) + r_offpara2_B1*y(19);
% y(21) = onborig_2
dydt(21) = (r_onbasal_B1 + r_addon_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1))*y(24) - r_offorig_B1*y(21);
% y(22) = onbpara1_2
dydt(22) = (r_onbasal_B1 + r_addon_byAprime1_B1*y(15)^n_Aprime1/(k_Aprime1^n_Aprime1 + y(15)^n_Aprime1))*y(24) - r_offpara1_B1*y(22);
% y(23) = onbpara2_2
dydt(23) = (r_onbasal_B1 + r_addon_byAprime2_B1*y(16)^n_Aprime2/(k_Aprime2^n_Aprime2 + y(16)^n_Aprime2))*y(24) - r_offpara2_B1*y(23);
% y(24) = offb_2
dydt(24) = -y(24)*((r_onbasal_B1 + r_addon_byA1_B1*y(5)^n_A1/(k_A1^n_A1 + y(5)^n_A1)) + (r_onbasal_B1 + r_addon_byAprime1_B1*y(15)^n_Aprime1/(k_Aprime1^n_Aprime1 + y(15)^n_Aprime1)) + (r_onbasal_B1 + r_addon_byAprime2_B1*y(16)^n_Aprime2/(k_Aprime2^n_Aprime2 + y(16)^n_Aprime2))) ...
    + r_offorig_B1*y(21) + r_offpara1_B1*y(22) + r_offpara2_B1*y(23);
% y(25) = b
dydt(25) = r_prodon_B1*y(17) + (d_Aprime1_B1)*r_prodon_B1*y(18) + + (d_Aprime2_B1)*r_prodon_B1*y(19) + ...
    r_prodon_B1*y(21) + (d_Aprime1_B1)*r_prodon_B1*y(22) + + (d_Aprime2_B1)*r_prodon_B1*y(23) - r_deg_B1*y(25);
end