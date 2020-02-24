% [1] source : Poulin P., Theil F. Prediction of pharmacokinetics prior to 
% in vivo studies. II. Generic physiologically based pharmacokinetic models 
% of drug disposition. 
% https://doi.org/10.1002/jps.10128
% [2] source : Poulin P., Theil F. Prediction of pharmacokinetics prior to 
% in vivo studies. I. Mechanism?based prediction of volume of distribution
% https://doi.org/10.1002/jps.10005
% [3] source : Igari Y., Sugiyama Y., Sawada Y., Iga T., Hanano M.
% Prediction of diazepam disposition in the rat and man by a 
% physiologically based pharmacokinetic model. 
% https://doi.org/10.1007/BF01059058

% PARS = [1@Weight_rat 2@Blood-Plasma_ratio 3@Dow 4@Pow 5@fup 6@CLint 7@Dose]
PARS = [0.25 1.04 2.07 2.99 0.14 83 1.2];
[t,y] = ode15s(@(t,y)PBPK_iv_model(t,y,PARS),[0 30],zeros(13,1));

cc = parula(13);
figure(1)
titles = {'Arterial blood';'Venous blood';'Lung';'Liver';'Adipose';...
    'Bone';'Brain';'Heart';'Kidney';'Muscle';'Spleen';'Gut';'Skin'};
for j=1:13
    subplot(4,4,j)
    plot(t,y(:,j),'-','LineWidth',2,'Color',cc(j,:))
    title(titles{j},'Color','w')
    set(gca,'Color','k','XColor','w','YColor','w','Box','on','FontName','GothamPro','FontSize',12)
end
set(gcf,'Color','k')


%% PBPK iv model
function dydt = PBPK_iv_model(t,y,param)
dydt = zeros(13,1);

Wrat = param(1);

% % Drug-specific parameters

Qc = 0.235*(Wrat^0.75); % cardiac output
BP = param(2); % blood-plasma ratio
Dow = param(3); % log olive oil:water partition coefficient of both the non-ionized and ionized species at pH=7.4
Pow = param(4); % log n-octanol-buffer partition coefficient of the non-ionized species at pH=7.4
fup = param(5); % plasma protein binding; unbound fraction in plasma
fut = 1/(1 + (((1 - fup)/fup)*0.5)); % unbound fraction in tissue
CLint = param(6); % intrinsic clearance
Dose = param(7); % dose

% % Concentration in :

Cab = y(1); % arterial blood
Cvb = y(2); % venous blood
Clung = y(3); % lung
Cliver = y(4); % liver
Cadipose = y(5); % adipose
Cbone = y(6); % bone
Cbrain = y(7); % brain
Cheart = y(8); % heart
Ckidney = y(9); % kidney
Cmuscle = y(10); % muscle
Cspleen = y(11); % spleen
Cgut = y(12); % gut
Cskin = y(13); % skin

% % Species-specific parameters

% blood flow rates
Qliver = 0.175*Qc;
Qadipose = 0.07*Qc;
Qbone = 0.122*Qc;
Qbrain = 0.02*Qc;
Qheart = 0.049*Qc;
Qkidney = 0.141*Qc;
Qmuscle = 0.278*Qc;
Qspleen = 0.02*Qc;
Qgut = 0.131*Qc;
Qskin = 0.058*Qc;
% volumes
Vliver = 0.0366*Wrat; 
Vadipose = 0.07*Wrat; % 0.0761 ? source [2]
Vbone = 0.04148*Wrat; % 0.041476 ? source [2]
Vbrain = 0.0057*Wrat; 
Vheart = 0.0033*Wrat;
Vkidney = 0.0073*Wrat;
Vmuscle = 0.404*Wrat;
Vspleen = 0.002*Wrat;
Vgut = 0.027*Wrat;
Vskin = 0.19*Wrat;
Vlung = 0.005*Wrat;
Vab = 0.0272*Wrat; % 0.0449 ? source [2]
Vvb = 0.00544*Wrat;
% water
VWliver = 0.705*Vliver;
VWadipose = 0.12*Vadipose;
VWbone = 0.446*Vbone;
VWbrain = 0.788*Vbrain;
VWheart = 0.779*Vheart;
VWkidney = 0.771*Vkidney;
VWmuscle = 0.756*Vmuscle;
VWspleen = 0.771*Vspleen;
VWgut = 0.749*Vgut;
VWskin = 0.651*Vskin;
VWlung = 0.79*Vlung;
VWab = 0.96*Vab;
% neutral lipids
VNLliver = 0.0138*Vliver;
VNLadipose = 0.853*Vadipose;
VNLbone = 0.0273*Vbone;
VNLbrain = 0.0392*Vbrain;
VNLheart = 0.014*Vheart;
VNLkidney = 0.0123*Vkidney;
VNLmuscle = 0.01*Vmuscle;
VNLspleen = 0.0077*Vspleen;
VNLgut = 0.0292*Vgut;
VNLskin = 0.0239*Vskin;
VNLlung = 0.0219*Vlung;
VNLab = 0.00147*Vab;
% phospholipids
VPHliver = 0.0303*Vliver;
VPHadipose = 0.002*Vadipose;
VPHbone = 0.0027*Vbone;
VPHbrain = 0.0533*Vbrain;
VPHheart = 0.0118*Vheart;
VPHkidney = 0.0284*Vkidney;
VPHmuscle = 0.009*Vmuscle;
VPHspleen = 0.0136*Vspleen;
VPHgut = 0.0138*Vgut;
VPHskin = 0.018*Vskin;
VPHlung = 0.014*Vlung;
VPHab = 0.00083*Vab;

% % Tissue:Plasma ratios

Pliver = (((Pow*(VNLliver + 0.3*VPHliver)) + (VWliver + 0.7*VPHliver))/...
    ((Pow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/fut);
Padipose = (((Dow*(VNLadipose + 0.3*VPHadipose)) + (VWadipose + 0.7*VPHadipose))/...
    ((Dow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/1);
Pbone = (((Pow*(VNLbone + 0.3*VPHbone)) + (VWbone + 0.7*VPHbone))/...
    ((Pow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/fut);
Pbrain = (((Pow*(VNLbrain + 0.3*VPHbrain)) + (VWbrain + 0.7*VPHbrain))/...
    ((Pow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/fut);
Pheart = (((Pow*(VNLheart + 0.3*VPHheart)) + (VWheart + 0.7*VPHheart))/...
    ((Pow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/fut);
Pkidney = (((Pow*(VNLkidney + 0.3*VPHkidney)) + (VWkidney + 0.7*VPHkidney))/...
    ((Pow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/fut);
Pmuscle = (((Pow*(VNLmuscle + 0.3*VPHmuscle)) + (VWmuscle + 0.7*VPHmuscle))/...
    ((Pow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/fut);
Pspleen = (((Pow*(VNLspleen + 0.3*VPHspleen)) + (VWspleen + 0.7*VPHspleen))/...
    ((Pow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/fut);
Pskin = (((Pow*(VNLskin + 0.3*VPHskin)) + (VWskin + 0.7*VPHskin))/...
    ((Pow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/fut);
Pgut = (((Pow*(VNLgut + 0.3*VPHgut)) + (VWgut + 0.7*VPHgut))/...
    ((Pow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/fut);
Plung = (((Pow*(VNLlung + 0.3*VPHlung)) + (VWlung + 0.7*VPHlung))/...
    ((Pow*(VNLab + 0.3*VPHab)) + (VWab + 0.7*VPHab)))*(fup/fut);

% % Venous blood leaving tissue

Clung_L = Clung/(Plung/BP); % lung
Cliver_L = Cliver/(Pliver/BP); % liver
Cadipose_L = Cadipose/(Padipose/BP); % adipose
Cbone_L = Cbone/(Pbone/BP); % bone
Cbrain_L = Cbrain/(Pbrain/BP); % brain
Cheart_L = Cheart/(Pheart/BP); % heart
Ckidney_L = Ckidney/(Pkidney/BP); % kidney
Cmuscle_L = Cmuscle/(Pmuscle/BP); % muscle
Cspleen_L = Cspleen/(Pspleen/BP); % spleen
Cgut_L = Cgut/(Pgut/BP); % gut
Cskin_L = Cskin/(Pskin/BP); % skin

% % Equations

Eliver = CLint/(CLint + Qliver); % hepatic extraction ratio

% lung
dydt(3) = (Qc/Vlung)*(Cvb - Clung_L);
% arterial blood
dydt(1) = (Qc/Vab)*(Clung_L - Cab);
% venous blood
dydt(2) = (Qadipose*Cadipose_L +...
    Qbone*Cbone_L + ...
    Qbrain*Cbrain_L + ...
    Qheart*Cheart_L + ...
    Qkidney*Ckidney_L + ...
    Qmuscle*Cmuscle_L + ...
    Qliver*Cliver_L + ...
    Qskin*Cskin_L - ...
    Qc*Cvb + infusion(t,Dose))/Vvb;
% liver
dydt(4) = (Cab*(Qliver - Qgut - Qspleen) + ...
    (Qgut*Cgut_L + Qspleen*Cspleen_L - Qliver*Cliver_L))/Vliver -...
    (Cab*(Qliver - Qgut - Qspleen) +...
    (Qgut*Cgut_L + Qspleen*Cspleen_L))*Eliver/Vliver;
% adipose
dydt(5) = (Qadipose/Vadipose)*(Cab - Cadipose_L);
% bone
dydt(6) = (Qbone/Vbone)*(Cab - Cbone_L);
% brain
dydt(7) = (Qbrain/Vbrain)*(Cab - Cbrain_L);
% heart
dydt(8) = (Qheart/Vheart)*(Cab - Cheart_L);
% kidney
dydt(9) = (Qkidney/Vkidney)*(Cab - Ckidney_L);
% muscle
dydt(10) = (Qmuscle/Vmuscle)*(Cab - Cmuscle_L);
% spleen
dydt(11) = (Qspleen/Vspleen)*(Cab - Cspleen_L); % (1+Eliver)*(Qspleen/Vspleen)*(Cab - Cspleen_L);
% gut
dydt(12) = (Qgut/Vgut)*(Cab - Cgut_L);
% skin
dydt(13) = (Qskin/Vskin)*(Cab - Cskin_L);
end

function g = infusion(t,D)
% theta = 20/60;
% g = D*theta*((theta*t)^2)*((1 - theta*t)^2); ? source [3]
g = D*exp(-100*t); % 100 ? temporary constant
end