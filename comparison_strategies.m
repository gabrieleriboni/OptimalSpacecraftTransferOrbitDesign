mu = 398600;
th0 = 0;
thf = 2*pi;
dth = 1e-5;

% Initial point data, [km] [km/s]
x_i = 4643.9563;
y_i = 5721.0035;
z_i = 2608.9369;
vx_i = -5.6160;
vy_i = 2.6040;
vz_i = 4.2930;

rr_i = [x_i; y_i; z_i];
vv_i = [vx_i; vy_i; vz_i];

[a_i, e_i, i_i, OM_i, om_i, th_i] = car2par(rr_i, vv_i, mu);

% Final point data, [km] [rad]
a_f = 14160.0;
e_f = 0.2088;
i_f = 1.4870;
OM_f = 1.6100;
om_f = 2.3450;
th_f = 1.5110;

[rr_f, vv_f] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu);

%% STANDARD STRATEGY PA
[delta_v1, om_pc, th_pc] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
delta_t1 = TOF(a_i, e_i, th_i, th_pc, mu);
[delta_v2, th_i_omc, th_f_omc] = changePericenterArg(a_i, e_i, om_pc, om_f, mu);
delta_t2 = TOF(a_i, e_i, th_f_omc(1), 2*pi, mu); % from first intersection to pericenter
[delta_v1_pa, delta_v2_pa, delta_t3] = bitangentTransfer(a_i, e_i, a_f, e_f, 'pa', mu);
delta_t4 = TOF(a_i, e_i,th_pc , th_i_omc(1), mu);
delta_t5 = TOF(a_f, e_f, pi, th_f, mu);
delta_v3 = delta_v1_pa + delta_v2_pa;

deltaT_PA = delta_t1 + delta_t2 + delta_t3 + delta_t4 + delta_t5; % 29110.0717 s
deltaV_PA = abs(delta_v1) + abs(delta_v2) + abs(delta_v3); % 9.8756 km/s

%% STANDARD STRATEGY AP
[delta_v1, om_pc, th_pc] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
delta_t1 = TOF(a_i, e_i, th_i, th_pc, mu);
[delta_v2, th_i_omc, th_f_omc] = changePericenterArg(a_i, e_i, om_pc, om_f, mu);
delta_t2 = TOF(a_i, e_i, th_f_omc(2), pi, mu); % from second intersection to apocenter
[delta_v1_ap, delta_v2_ap, delta_t3] = bitangentTransfer(a_i, e_i, a_f, e_f, 'ap', mu);
delta_t4 = TOF(a_i, e_i,th_pc , th_i_omc(2), mu);
delta_t5 = TOF(a_f, e_f, 0, th_f, mu);
delta_v3 = delta_v1_ap + delta_v2_ap;

deltaT_AP = delta_t1 + delta_t2 + delta_t3 + delta_t4 + delta_t5; % 16146.5137 s
deltaV_AP = abs(delta_v1) + abs(delta_v2) + abs(delta_v3); % 9.9444 km/s

%% ALTERNATIVE 1
[~, ~, ~, ~, delta_v1, ~, ~, ~, delta_t1, ~, omf_pc] = biellipticChangePlane(a_i, e_i, a_f, e_f, 1.731830*1e4, i_i, OM_i, om_i, i_f, OM_f, th_i, mu);
[delta_v2, th1, th2] = changePericenterArg(a_f, e_f, omf_pc, om_f, mu);
delta_t2 = TOF(a_f, e_f, 0, th1(2), mu);
delta_t3 = TOF(a_f, e_f, th2(2), th_f, mu);

deltaV_ALT1 = delta_v1 + abs(delta_v2); % 8.7061 km/s
deltaT_ALT1 = delta_t1 + delta_t2 + delta_t3; % 52066.4027 s

%% ALTERNATIVE 2
ra_i = a_i * (1 + e_i);
ra_f = a_f * (1 + e_f);

deltaT0 = TOF(a_i, e_i, th_i, pi, mu);
[delta_v1_ap, delta_v2_ap, deltaT1] = bitangentTransfer(a_i, e_i, ra_i, 0, 'ap', mu);
deltaV1 = abs(delta_v1_ap) + abs(delta_v2_ap);
[delta_v1_h, delta_v2_h, deltaT2] = bitangentTransfer(ra_i, 0, ra_f, 0, 'pa', mu);
deltaVH = abs(delta_v1_h) + abs(delta_v2_h);
[delta_v_pc, om_pc, th_pc] = changeOrbitalPlane(ra_f, 0, i_i, OM_i, om_i, i_f, OM_f, mu);
deltaT3 = TOF(ra_f, 0, pi, th_pc, mu);
[delta_v_PerArg, th1, th2] = changePericenterArg(ra_f, 0, om_pc, om_f, mu);
deltaT4 = TOF(ra_f, 0, th_pc, th1(2), mu); % th1(2) = 249.2707°
deltaT5 = TOF(ra_f, 0, th2(2), pi, mu); % th2(2) = 110.7293°
[delta_v1, delta_v2, deltaT6] = bitangentTransfer(ra_f, 0, a_f, e_f, 'ap', mu);
deltaV2 = abs(delta_v1) + abs(delta_v2);
deltaT7 = TOF(a_f, e_f, 0, th_f, mu);
deltaV_ALT2 = deltaV1 + deltaVH + delta_v_pc + delta_v_PerArg + deltaV2; % 7.9673 km/s
deltaT_ALT2 = deltaT0 + deltaT1 + deltaT2 + deltaT3 + deltaT4 + deltaT5 + deltaT6 + deltaT7; % 36616.8834 s

%% ALTERNATIVE 3 VM
p_i = a_i * (1 - e_i^2);
p_f = a_f * (1 - e_f^2);
[delta_v_pc, om_pc, th_pc] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
delta_t1 = TOF(a_i, e_i, th_i, th_pc, mu); % from th_i to th_pc (initial orbit)
delta_t2 = TOF(a_i, e_i, th_pc, 0, mu); % from th_pc to pericenter (initial orbit)
th_2 = 0:0.01:th_f;
deltaT_vect = [];
deltaV_vect = [];
for i = 1 : length(th_2)
    th_t2 = th_2(i) - om_pc + om_f;
    k = p_i / (1 + e_i);
    e_t = (p_f - k - k * e_f * cos(th_2(i))) / (k + k * e_f * cos(th_2(i)) - p_f * cos(th_t2));
    p_t = p_i * (1 + e_t) / (1 + e_i);
    a_t = p_t / (1 - e_t^2);
    
    [delta_v1, ~, ~] = bitangentTransfer(a_i, e_i, a_t, e_t, 'pa', mu);
    [delta_v_sec] = secantTransfer(a_t, e_t, a_f, e_f, th_t2, th_2(i), mu);
    delta_t3 = TOF(a_t, e_t, 0, th_t2, mu); % from pericenter to th_t (transfer orbit)

    deltaV_ALT3 = abs(delta_v_pc) + abs(delta_v1) + delta_v_sec;
    deltaT_ALT3 = delta_t1 + delta_t2 + delta_t3;

    deltaT_vect = [deltaT_vect, deltaT_ALT3];
    deltaV_vect = [deltaV_vect, deltaV_ALT3];

end

delta_v_min = deltaV_vect(41);
th_2_vm = th_2(41);
th_t2_vm = th_2_vm - om_pc + om_f;
k = p_i / (1 + e_i);
e_tvm = (p_f - k - k * e_f * cos(th_2_vm)) / (k + k * e_f * cos(th_2_vm) - p_f * cos(th_t2_vm)); % 0.1905
p_tvm = p_i * (1 + e_tvm) / (1 + e_i);
a_tvm = p_tvm / (1 - e_tvm^2);

[delta_v1, ~, ~] = bitangentTransfer(a_i, e_i, a_tvm, e_tvm, 'pa', mu);
[delta_v_sec_vm] = secantTransfer(a_tvm, e_tvm, a_f, e_f, th_t2_vm, th_2_vm, mu);
delta_t3_vm = TOF(a_tvm, e_tvm, 0, th_t2_vm, mu); % from pericenter to th_t (transfer orbit)
delta_t4_vm = TOF(a_f, e_f, th_2_vm, th_f, mu);

deltaV_ALT3_vm = abs(delta_v_pc) + abs(delta_v1) + delta_v_sec_vm; % 8.5063 km/s
deltaT_ALT3_vm = delta_t1 + delta_t2 + delta_t3_vm + delta_t4_vm; % 16991.9979 s
%% ALTERNATIVE 3 TM
p_i = a_i * (1 - e_i^2);
p_f = a_f * (1 - e_f^2);
[delta_v_pc, om_pc, th_pc] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
delta_t1 = TOF(a_i, e_i, th_i, th_pc, mu); % 5330.5965 s   from th_i to th_pc (initial orbit)
delta_t2 = TOF(a_i, e_i, th_pc, 0, mu); % from th_pc to pericenter (initial orbit)

th_2 = 0:0.01:th_f;
deltaT_vect = [];
deltaV_vect = [];
for i = 1 : length(th_2)
    th_t2 = th_2(i) - om_pc + om_f;
    k = p_i / (1 + e_i);
    e_t = (p_f - k - k * e_f * cos(th_2(i))) / (k + k * e_f * cos(th_2(i)) - p_f * cos(th_t2));
    p_t = p_i * (1 + e_t) / (1 + e_i);
    a_t = p_t / (1 - e_t^2);
    
    [delta_v1, ~, ~] = bitangentTransfer(a_i, e_i, a_t, e_t, 'pa', mu);
    [delta_v_sec] = secantTransfer(a_t, e_t, a_f, e_f, th_t2, th_2(i), mu);
    delta_t3 = TOF(a_t, e_t, 0, th_t2, mu); % from pericenter to th_t (transfer orbit)

    deltaV_ALT3 = abs(delta_v_pc) + abs(delta_v1) + delta_v_sec;
    deltaT_ALT3 = delta_t1 + delta_t2 + delta_t3;

    deltaT_vect = [deltaT_vect, deltaT_ALT3];
    deltaV_vect = [deltaV_vect, deltaV_ALT3];

end

t_min = deltaT_vect(1);
th_2_tm = th_2(1);
th_t2_tm = th_2_tm - om_pc + om_f;
e_ttm = (p_f - k - k * e_f * cos(th_2_tm)) / (k + k * e_f * cos(th_2_tm) - p_f * cos(th_t2_tm));
p_ttm = p_i * (1 + e_ttm) / (1 + e_i);
a_ttm = p_ttm / (1 - e_ttm^2);
[delta_v1, ~, ~] = bitangentTransfer(a_i, e_i, a_ttm, e_ttm, 'pa', mu);
[delta_v_sec_tm] = secantTransfer(a_ttm, e_ttm, a_f, e_f, th_t2_tm, th_2_tm, mu);
delta_t3_tm = TOF(a_ttm, e_ttm, 0, th_t2_tm, mu); % from pericenter to th_t (transfer orbit)
delta_t4_tm = TOF(a_f, e_f, th_2_tm, th_f, mu);

deltaV_ALT3_tm = abs(delta_v_pc) + abs(delta_v1) + delta_v_sec_tm; % 8.8634 km/s
deltaT_ALT3_tm = delta_t1 + delta_t2 + delta_t3_tm + delta_t4_tm; % 14443.8075 s
%% ALTERNATIVE 4
% You have to plot it from 'progettoFinaleC12.m'
deltaV_ALT4_vm = dv_min;
deltaT_ALT4_vm = dt_vmin;

deltaV_ALT4_tm = dv_TOT;
deltaT_ALT4_tm = dt_min;

%% Comparison of strategies (standard and alternatives)
cost_vect = [deltaV_PA, deltaV_AP, deltaV_ALT1, deltaV_ALT2, deltaV_ALT3_vm, deltaV_ALT3_tm, deltaV_ALT4_vm, deltaV_ALT4_tm]; % speed vector [ km/s ];
time_vect_sec = [deltaT_PA, deltaT_AP, deltaT_ALT1, deltaT_ALT2, deltaT_ALT3_vm, deltaT_ALT3_tm, deltaT_ALT4_vm, deltaT_ALT4_tm]; % time vector [ s ];
time_vect_hour = time_vect_sec / 3600; % time vector, [ h ]

k = cost_vect .* time_vect_hour;
kk = [cost_vect; time_vect_hour]';
vv = 0:0.0001:100;

y1 = @(v) k(1) ./ v;
y2 = @(v) k(2) ./ v;
y3 = @(v) k(3) ./ v;
y4 = @(v) k(4) ./ v;
y5 = @(v) k(5) ./ v;
y6 = @(v) k(6) ./ v;
y7 = @(v) k(7) ./ v;
y8 = @(v) k(8) ./ v;


x = [1, 2, 3, 4, 5, 6, 7, 8];
w1 = 0.5;
w2 = 0.25;
label_text = sprintf('DeltaV [km/s]\nDeltaT [h]');

figure(1)
bar(x, cost_vect, w1)
hold on
bar(x, time_vect_hour, w2)
ylabel(label_text);
grid on
legend('DeltaV', 'DeltaT', 'Location','northwest')
ax = gca;
ax.XTick = [1 2 3 4 5 6 7 8]; 
ax.XTickLabels = {'PA','AP','ALT1','ALT2', 'ALT3 vm', 'ALT3 tm', 'ALT4 vm', 'ALT4 tm'};
hold off

figure(2)
plot(vv, y1(vv))
hold on
plot(vv, y2(vv))
plot(vv, y3(vv))
plot(vv, y4(vv))
plot(vv, y5(vv))
plot(vv, y6(vv))
plot(vv, y7(vv))
plot(vv, y8(vv))
plot(kk(:,1), kk(:,2), 'k.', 'MarkerSize',7)
ylim([0,20])
xlim([0,70])
xlabel('Delta V [km/s]');
ylabel('Delta T [h]');
grid on
legend('PA', 'AP', 'ALT1','ALT2', 'ALT3 vm', 'ALT3 tm','ALT4 vm', 'ALT4 tm','')

