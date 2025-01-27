%% STANDARD STRATEGY

rp_i = a_i * (1 - e_i);
ra_f = a_f * (1 + e_f);
rp_t = rp_i;
ra_t = ra_f;
a_t = 0.5 * (rp_t + ra_t);
e_t = (ra_t - rp_t) / (ra_t + rp_t);

% Change orbit plane
[delta_v_pc, om_pc, th_pc] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
delta_t_plane = TOF(a_i, e_i, th_i, th_pc, mu);

% Change of pericenter argument maneuver
[delta_v_omc, th_i_omc, th_f_omc] = changePericenterArg(a_i, e_i, om_pc, om_f, mu);
delta_t_1 = TOF(a_i, e_i, th_f_omc(1), 2*pi, mu); % from first intersection to pericenter

% Plotting Vectors
vet_a = [a_i;a_i;a_i;a_t;a_f];
vet_e = [e_i;e_i;e_i;e_t;e_f];
vet_i = [i_i;i_f;i_f;i_f;i_f];
vet_OM = [OM_i;OM_f;OM_f;OM_f;OM_f];
vet_om = [om_i;om_pc;om_f;om_f;om_f];
vet_th_i = [th_i;th_pc;0;pi;pi];
vet_th_f = [th_pc;th_i_omc(1);th_f_omc(1);2*pi;-th_f];

% Vectors Chart 2
vet2_th_i = [0;0;0;0;0];
vet2_th_f = [2*pi;2*pi;2*pi;2*pi;2*pi];

% Chart 2 Options
graphopts.linewidth = [0.125 0.125 0.125 0.125 0.125];
graphopts.color = [0 0.4470 0.7410; 0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
graphopts.linestyle = '-';

% Plotting 
plotOrbitMul(vet_a, vet_e, vet_i, vet_OM, vet_om, vet_th_i, vet_th_f, 200, mu);
hold on
plotOrbitMul_GraphOpts(vet_a, vet_e, vet_i, vet_OM, vet_om, vet2_th_i, vet2_th_f, 200, mu, graphopts)
legend('','Initial Orbit','Plane Changed Orbit','Changed Pericenter Anomaly Orbit',...
    'Bitangent Elliptical Transfer Orbit','Final Orbit');
title('Standard Transfer');


%% Bielliptic transfer for elliptic orbits and change orbital plane (ALTERNATIVE STRATEGY) 

% Calculation of auxiliary values
ra_t = 1.731830*1e4;
rp_i = a_i * (1 - e_i);
rp_f = a_f * (1 - e_f);

rp_t1 = rp_i;
ra_t1 = ra_t;
a_t1 = 0.5 * (rp_t1 + ra_t1);
e_t1 = (ra_t1 - rp_t1) / (ra_t1 + rp_t1);

rp_t2 = rp_f;
ra_t2 = ra_t;
a_t2 = 0.5 * (rp_t2 + ra_t2);
e_t2 = (ra_t2 - rp_t2) / (ra_t2 + rp_t2);

% Plotting Vectors
vet_a = [a_i; a_t1; a_t2; a_t2; a_t2];
vet_e = [e_i; e_t1; e_t2; e_t2; e_t2];
vet_i = [i_i; i_i; i_i; i_f; i_f];
vet_OM = [OM_i; OM_i; OM_i; OM_f; OM_f];
vet_om = [om_i; om_i; om_i; omf_pc; om_f];
vet_th_i = [0; 0; pi; th_pc_2; th2(2)];
vet_th_f = [2*pi; pi; th_pc_2; th1(2); 2*pi+th_f];

% Chart 2 Options
graphopts2.linewidth = 1;
graphopts2.color = [0 0 0];
graphopts2.linestyle = '--';

% Vectors Chart 3
vet3_th_i = [0;0;0;0;0];
vet3_th_f = [2*pi;2*pi;2*pi;2*pi;2*pi];

% Chart 3 Options
graphopts3.linewidth = [0.125 0.125 0.125 0.125 0.125];
graphopts3.color = [0 0.4470 0.7410; 0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
graphopts3.linestyle = ':';

% Plotting
plotOrbitMul(vet_a, vet_e, vet_i, vet_OM, vet_om, vet_th_i, vet_th_f, 200, mu);
hold on
plotOrbitMul_GraphOpts(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 200, mu, graphopts2);
hold on
plotOrbitMul_GraphOpts(vet_a, vet_e, vet_i, vet_OM, vet_om, vet3_th_i, vet3_th_f, 200, mu, graphopts3);
legend('','Initial Orbit','Bitangent Orbit 1','Bitangent Orbit 2','Plane Changed Orbit',...
    'Changed Pericenter Anomaly Orbit','Final Orbit');


%% Circularization and Hohmann maneuver (ALTERNATIVE STRATEGY) 

% Calculation of auxiliary values
ra_i_1 = a_i * (1 + e_i);
rp_f_1 = a_i * (1 - 0);
rpT_1 = ra_i_1;
raT_1 = rp_f_1;
aT_1 = 0.5 * (rpT_1 + raT_1);
eT_1 = (rpT_1 - raT_1) / (raT_1 + rpT_1);

rp_i_2 = a_i * (1 - 0);
ra_f_2 = a_f * (1 + 0);
rpT_2 = rp_i_2;
raT_2 = ra_f_2;
aT_2 = 0.5 * (rpT_2 + raT_2);
eT_2 = (raT_2 - rpT_2) / (raT_2 + rpT_2);

ra_i_3 = a_f * (1 + 0);
rp_f_3 = a_f * (1 - e_f);
rpT_3 = ra_i_3;
raT_3 = rp_f_3;
aT_3 = 0.5 * (rpT_3 + raT_3);
eT_3 = (rpT_3 - raT_3) / (raT_3 + rpT_3);

% Plotting Vectors
vet_a = [a_i;aT_1;aT_2;a_f;a_f;a_f;aT_3];
vet_e = [e_i;eT_1;eT_2;0;0;0;eT_3];
vet_i = [i_i;i_i;i_i;i_i;i_f;i_f;i_f];
vet_OM = [OM_i;OM_i;OM_i;OM_i;OM_f;OM_f;OM_f];
vet_om = [om_i;om_i;om_i;om_i;om_pc;om_f;om_f];
vet_th_i = [th_i;pi;0;pi;th_pc;th2(2);pi];
vet_th_f = [pi;2*pi;pi;th_pc;th1(2);pi;2*pi];

% Chart 2 Options
graphopts2.linewidth = 2;
graphopts2.color = [1 0 1];
graphopts2.linestyle = '-';

% Vectors Chart 3
vet3_th_i = [0;0;0;0;0;0;0];
vet3_th_f = [2*pi;2*pi;2*pi;2*pi;2*pi;2*pi;2*pi];

% Chart 3 Options
graphopts3.linewidth = [0.125 0.125 0.125 0.125 0.125 0.125 0.125];
graphopts3.color = [0 0.4470 0.7410; 0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];
graphopts3.linestyle = ':';

% Chart 4 Options
graphopts4.linewidth = 0.125;
graphopts4.color = [1 0 1];
graphopts4.linestyle = ':';

% Plotting
plotOrbitMul(vet_a, vet_e, vet_i, vet_OM, vet_om, vet_th_i, vet_th_f, 200, mu);
hold on
plotOrbitMul_GraphOpts(a_f, e_f, i_f, OM_f, om_f, 0, th_f, 200, mu, graphopts2);
hold on
plotOrbitMul_GraphOpts(vet_a, vet_e, vet_i, vet_OM, vet_om, vet3_th_i, vet3_th_f, 200, mu, graphopts3);
hold on
plotOrbitMul_GraphOpts(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 200, mu, graphopts4);
legend('','Initial Orbit','Circularized orbit 1','Hohmann maneuver','Waiting orbit to plane change point',...
    'Plane Changed Orbit','Changed Pericenter Anomaly Orbit','Bitangent Orbit','Final Orbit');




