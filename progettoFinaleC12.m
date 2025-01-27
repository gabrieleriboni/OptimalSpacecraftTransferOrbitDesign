clear
close all
clc

%%
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

planet = 'Earth Cloudy';
opts.Units = 'km';

% Definition of plan parameters
a = 0;
b = 0;
c = 0;

% Creation of the x and y co-ordinate grid
[x, y] = meshgrid(-1.5e4:3e4:2e4); % It varies from -10 to 10 with 0.5 step
% Calculation of z values corresponding to the plane
z = a * x + b * y + c;
 
% figure
% planet3D(planet, opts);
% light('Position',[1,-1,0]);
% hold on;
% plotOrbit(a_i, e_i , i_i, OM_i, om_i, th0, thf, dth, mu)
% plotOrbit(a_f, e_f , i_f, OM_f, om_f, th0, thf, dth, mu)
% surf(x, y, z,'FaceColor', [0.204, 0.229, 0.255], 'FaceAlpha',0.4);
% xlabel('X [km]');
% ylabel('Y [km]');
% zlabel('Z [km]');
% grid on
% legend('', 'Initial orbit', 'Final orbit')

%% STANDARD STRATEGY
% Change orbit plane
[delta_v_pc, om_pc, th_pc] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
delta_t_plane = TOF(a_i, e_i, th_i, th_pc, mu);

% Change of pericenter argument maneuver
[delta_v_omc, th_i_omc, th_f_omc] = changePericenterArg(a_i, e_i, om_pc, om_f, mu);
delta_t_1 = TOF(a_i, e_i, th_f_omc(1), 2*pi, mu); % from first intersection to pericenter
delta_t_2 = TOF(a_i, e_i, th_f_omc(2), pi, mu); % from second intersection to apocenter

% Bitangent transfer for elliptic orbits
[delta_v1_pa, delta_v2_pa, delta_t_pa] = bitangentTransfer(a_i, e_i, a_f, e_f, 'pa', mu);
[delta_v1_pp, delta_v2_pp, delta_t_pp] = bitangentTransfer(a_i, e_i, a_f, e_f, 'pp', mu);
[delta_v1_ap, delta_v2_ap, delta_t_ap] = bitangentTransfer(a_i, e_i, a_f, e_f, 'ap', mu);
[delta_v1_aa, delta_v2_aa, delta_t_aa] = bitangentTransfer(a_i, e_i, a_f, e_f, 'aa', mu);

% Second intersection and 'ap'
delta_t_pc_omc = TOF(a_i, e_i,th_pc , th_i_omc(2), mu);
delta_t_thf = TOF(a_f, e_f, 0, th_f, mu);
delta_t_TOT_ap = delta_t_plane + delta_t_pc_omc + delta_t_2 + delta_t_ap + delta_t_thf; % 16146.5137 s
delta_v_ap = delta_v1_ap + delta_v2_ap;
delta_v_TOT_ap = abs(delta_v_pc) + abs(delta_v_omc) + abs(delta_v_ap); % 9.9444 km/s

% First intersection and 'pa' (cheaper cost)
delta_t_pc_omc = TOF(a_i, e_i,th_pc , th_i_omc(1), mu);
delta_t_thf = TOF(a_f, e_f, pi, th_f, mu);
delta_t_TOT_pa = delta_t_plane + delta_t_pc_omc + delta_t_1 + delta_t_pa + delta_t_thf; % 29110.0717 s
delta_v_pa = delta_v1_pa + delta_v2_pa;
delta_v_TOT_pa = abs(delta_v_pc) + abs(delta_v_omc) + abs(delta_v_pa); % 9.8756 km/s

% Second intersection and 'aa'
delta_t_pc_omc = TOF(a_i, e_i,th_pc , th_i_omc(2), mu);
[delta_v_omc_aa, th1omc, th2omc] = changePericenterArg(a_f, e_f, om_f + pi, om_f, mu); % 2nd change of pericenter argument maneuver
delta_t_pc_omc2 = TOF(a_f, e_f, pi, th1omc(1), mu);
delta_t_thf = TOF(a_f, e_f, th2omc(1), th_f, mu);
delta_t_TOT_aa = delta_t_plane + delta_t_pc_omc + delta_t_2 + delta_t_aa + delta_t_pc_omc2 + delta_t_thf;
delta_v_aa = delta_v1_aa + delta_v2_aa;
delta_v_TOT_aa = abs(delta_v_pc) + abs(delta_v_omc) + abs(delta_v_aa) + abs(delta_v_omc_aa);

% First intersection and 'pp'
delta_t_pc_omc = TOF(a_i, e_i,th_pc , th_i_omc(1), mu);
[delta_v_omc_pp, th1omc, th2omc] = changePericenterArg(a_f, e_f, om_f + pi, om_f, mu); % 2nd change of pericenter argument maneuver
delta_t_pc_omc2 = TOF(a_f, e_f, 0, th1omc(2), mu);
delta_t_thf = TOF(a_f, e_f, th2omc(2), th_f, mu);
delta_t_TOT_pp = delta_t_plane + delta_t_pc_omc + delta_t_1 + delta_t_pp + delta_t_pc_omc2 + delta_t_thf;
delta_v_pp = delta_v1_pp + delta_v2_pp;
delta_v_TOT_pp = abs(delta_v_pc) + abs(delta_v_omc) + abs(delta_v_pp) + abs(delta_v_omc_pp);

deltaV_vect = [delta_v_TOT_pa, delta_v_TOT_ap, delta_v_TOT_aa, delta_v_TOT_pp];
deltaT_vect = [delta_t_TOT_pa, delta_t_TOT_ap, delta_t_TOT_aa, delta_t_TOT_pp];

x = [1, 2, 3, 4];
w1 = 0.5;
w2 = 0.25;
label_text = sprintf('DeltaV [km/s]\nDeltaT [h]');

figure(2)
bar(x, deltaV_vect, w1)
hold on
bar(x, deltaT_vect / 3600, w2)
ylabel(label_text);
grid on
legend('DeltaV', 'DeltaT', 'Location','northwest')
ax = gca;
ax.XTick = [1 2 3 4]; 
ax.XTickLabels = {'P-A','A-P','A-A','P-P'};
hold off

deltaV_PA = delta_v_TOT_pa; % cost of the standard strategy
deltaV_AP = delta_v_TOT_ap; % cost of the 'ap' strategy instead of 'pa'

deltaT_PA = delta_t_TOT_pa; % time for the standard strategy
deltaT_AP = delta_t_TOT_ap; % time for the 'ap' strategy instead of 'pa'

%% Bielliptic transfer for elliptic orbits and change orbital plane (1st ALTERNATIVE STRATEGY)
ra_t = 6000:0.25:80000;
delta_v = [];
delta_t = [];
for i = 1 : length(ra_t)
    [~, ~, ~, ~, delta_v_TOT, ~, ~, ~, delta_t_TOT, ~, ~] = biellipticChangePlane(a_i, e_i, a_f, e_f, ra_t(i), i_i, OM_i, om_i, i_f, OM_f, th_i, mu);
    delta_v = [delta_v; delta_v_TOT];
    delta_t = [delta_t; delta_t_TOT];
end

label_text = sprintf('DeltaV [km/s]\nDeltaT [h]');

figure(3)
plot(ra_t, delta_v, 'LineWidth',1.5)
hold on
plot(ra_t, delta_t / 3600, 'LineWidth',1.5)
plot(1.731830*1e4, 6.587182,'kx','LineWidth',1.5)
xline(1.731830*1e4)
yline(6.587182)
xlabel('Transfer orbits apocenter distance [km]');
ylabel(label_text);
grid on
legend('DeltaV', 'DeltaT', 'Location','northwest')

% Change orbit plane and Change of pericenter argument maneuver
[~, ~, ~, ~, delta_v_TOT, ~, ~, ~, delta_t_TOT, th_pc, omf_pc] = biellipticChangePlane(a_i, e_i, a_f, e_f, 1.731830*1e4, i_i, OM_i, om_i, i_f, OM_f, th_i, mu);
[delta_v, th1, th2] = changePericenterArg(a_f, e_f, omf_pc, om_f, mu);

delta_v_TOT = delta_v_TOT + abs(delta_v); % 8.7061 km/s

delta_t3 = TOF(a_f, e_f, 0, th1(2), mu);
delta_t4 = TOF(a_f, e_f, th2(2), th_f, mu);
delta_t_TOT = delta_t_TOT + delta_t3 + delta_t4; % 52066.4027 s

deltaV_ALT1 = delta_v_TOT; % cost of the 2nd alternative strategy
deltaT_ALT1 = delta_t_TOT;
%% Circularization and Hohmann maneuver (2nd ALTERNATIVE STRATEGY)
ra_i = a_i * (1 + e_i);
ra_f = a_f * (1 + e_f);
% From th_i to apocenter
deltaT0 = TOF(a_i, e_i, th_i, pi, mu);

% Circularization of the first orbit at its apocenter
[delta_v1_ap, delta_v2_ap, deltaT1] = bitangentTransfer(a_i, e_i, ra_i, 0, 'ap', mu);
deltaV1 = abs(delta_v1_ap) + abs(delta_v2_ap);

% Hohmann maneuver
[delta_v1_h, delta_v2_h, deltaT2] = bitangentTransfer(ra_i, 0, ra_f, 0, 'pa', mu);
deltaVH = abs(delta_v1_h) + abs(delta_v2_h);

% Change orbit plane
[delta_v_pc, om_pc, th_pc] = changeOrbitalPlane(ra_f, 0, i_i, OM_i, om_i, i_f, OM_f, mu);
deltaT3 = TOF(ra_f, 0, pi, th_pc, mu);

% Change of pericenter argument maneuver
[delta_v_PerArg, th1, th2] = changePericenterArg(ra_f, 0, om_pc, om_f, mu);
deltaT4 = TOF(ra_f, 0, th_pc, th1(2), mu);
deltaT5 = TOF(ra_f, 0, th2(2), pi, mu);

% Decircularization of the second orbit
[delta_v1, delta_v2, deltaT6] = bitangentTransfer(ra_f, 0, a_f, e_f, 'ap', mu);
deltaV2 = abs(delta_v1) + abs(delta_v2);
deltaT7 = TOF(a_f, e_f, 0, th_f, mu);
 
deltaV_ALT2 = deltaV1 + deltaVH + delta_v_pc + delta_v_PerArg + deltaV2; % 7.9673 km/s
deltaT_ALT2 = deltaT0 + deltaT1 + deltaT2 + deltaT3 + deltaT4 + deltaT5 + deltaT6 + deltaT7; % 36616.8834 s

%% Change of Plane, Tangent and Secant (3rd ALTERNATIVE STRATEGY)
th_f = 1.5110;
om_f = 2.3450;
om_i = 0.527014267568125;
p_i = a_i * (1 - e_i^2);
p_f = a_f * (1 - e_f^2);
[delta_v_pc, om_pc, th_pc] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
delta_t1 = TOF(a_i, e_i, th_i, th_pc, mu); % from th_i to th_pc (initial orbit)
delta_t2 = TOF(a_i, e_i, th_pc, 0, mu); % from th_pc to pericenter (initial orbit)

[rr_pi, ~] = par2car(a_i, e_i, i_f, OM_f, om_pc, 0, mu);

th_2 = 0:0.01:th_f;
deltaT_vect = [];
deltaV_vect = [];

planet = 'Earth Cloudy';
opts.Units = 'km';
figure(1)
planet3D(planet, opts);
light('Position',[1,-1,0]);
hold on;
plotOrbit(a_i, e_i , i_f, OM_f, om_pc, th0, thf, dth, mu)
plotOrbit(a_f, e_f , i_f, OM_f, om_f, th0, thf, dth, mu)
plot3(rr_pi(1),rr_pi(2),rr_pi(3),'.',MarkerSize=10,Color=[0 0 0])

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on
j=0;
for i = 1 : length(th_2)
    j=j+1;
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

    if mod(j,11)==0
    plotOrbit(a_t, e_t , i_f, OM_f, om_pc, th0, th_t2, dth, mu,'g')
    end
end
legend('', 'Initial orbit', 'Final orbit','Initial point','Auxiliary orbits')
hold off

delta_v_min = deltaV_vect(41);
th_2_vm = th_2(41);
th_t2_vm = th_2_vm - om_pc + om_f;
k = p_i / (1 + e_i);
e_tvm = (p_f - k - k * e_f * cos(th_2_vm)) / (k + k * e_f * cos(th_2_vm) - p_f * cos(th_t2_vm));
p_tvm = p_i * (1 + e_tvm) / (1 + e_i);
a_tvm = p_tvm / (1 - e_tvm^2);

[delta_v1, ~, ~] = bitangentTransfer(a_i, e_i, a_tvm, e_tvm, 'pa', mu);
[delta_v_sec_vm] = secantTransfer(a_tvm, e_tvm, a_f, e_f, th_t2_vm, th_2_vm, mu);
delta_t3_vm = TOF(a_tvm, e_tvm, 0, th_t2_vm, mu); % from pericenter to th_t (transfer orbit)
delta_t4_vm = TOF(a_f, e_f, th_2_vm, th_f, mu);

deltaV_ALT3_vm = abs(delta_v_pc) + abs(delta_v1) + delta_v_sec_vm; % 8.5063 km/s
deltaT_ALT3_vm = delta_t1 + delta_t2 + delta_t3_vm + delta_t4_vm; % 16991.9979 s

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

[rr_vm, ~] = par2car(a_f, e_f, i_f, OM_f, om_f, th_2_vm, mu);
[rr_tm, ~] = par2car(a_f, e_f, i_f, OM_f, om_f, th_2_tm, mu);


planet = 'Earth Cloudy';
opts.Units = 'km';
figure(2)
planet3D(planet, opts);
light('Position',[1,-1,0]);
hold on;
plotOrbit(a_i, e_i , i_f, OM_f, om_pc, th0, thf, dth, mu)
plotOrbit(a_f, e_f , i_f, OM_f, om_f, th0, thf, dth, mu)
plotOrbit(a_tvm, e_tvm, i_f, OM_f, om_pc, th0, th_t2_vm, dth, mu, 'm')
plotOrbit(a_ttm, e_ttm , i_f, OM_f, om_pc, th0, th_t2_tm, dth, mu,'c')
plot3(rr_pi(1),rr_pi(2),rr_pi(3),'.',MarkerSize=10,Color=[0 0 0])
plot3(rr_vm(1), rr_vm(2), rr_vm(3), 'mx' , 'LineWidth',1)
plot3(rr_tm(1), rr_tm(2), rr_tm(3), 'cx' , 'LineWidth',1)
plot3(rr_f(1), rr_f(2), rr_f(3), 'k.' , 'MarkerSize',10)
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on
legend('','Initial orbit on the final plane', 'Final orbit', 'dV min orbit', 'dT min orbit', 'Initial point', 'Intersection point dV min', 'Intersection point dT min', 'Final point' )
%% Direct transfer (4th ALTERNATIVE STRATEGY)

% Searching for an orbit that goes directly from the first point to the last one
hht_direction = cross(rr_i, rr_f);

it = acos(hht_direction(3) / norm(hht_direction));
di=i_f-i_i;

N = cross([0; 0; 1], hht_direction);
N = N / norm(N);

OMt = 2*pi - acos(N(1));

[dVplane, om1f, theta_changePlane] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, it, OMt, mu);

% rr_i * rr_f = r_i * r_f * cos(dtheta)

r_i = norm(rr_i);
r_f = norm(rr_f);
dtheta = acos(dot(rr_i, rr_f) / (r_i * r_f));

tht_vect = linspace(0, 2*pi, 360);

th1t_vect = [];
th2t_vect = [];
dt_VECT = [];
dv1_vect = [];
dv2_vect = [];
dt_min = 1e6;
dv_min = 1e3;
k = 0;

planet3D(planet, opts);
light('Position',[1,-1,0]);
hold on;
plotOrbit(a_i, e_i , i_i, OM_i, om_i, th0, thf, dth, mu,'b')
plotOrbit(a_f, e_f , i_f, OM_f, om_f, th0, thf, dth, mu,'r')
plot3(rr_i(1),rr_i(2),rr_i(3),'.',MarkerSize=10,Color=[0 0 0])
plot3(rr_f(1),rr_f(2),rr_f(3),'.',MarkerSize=10,Color=[0 0 0])

for th1_t = tht_vect

    th2_t = th1_t + dtheta;
    
    if th2_t > 2*pi
        th2_t = th2_t - 2*pi;
    end

    omt = om1f + th_i - th1_t;
    if omt < 0
        omt = omt + 2*pi;
    end

    et = (r_f - r_i) / ((r_i * cos(th1_t)) - (r_f * cos(th2_t)));

    if et >= 0 && et < 1
        pt = r_i + r_i * et * cos(th1_t);
        at = pt / (1 - et^2);
        rpt = pt / (1 + et);
        
        
        if rpt >= (6378 + 100)
            k = k + 1;

            [rt1, vt1] = par2car(at, et, it, OMt, omt, th1_t, mu);
            [rt2, vt2] = par2car(at, et, it, OMt, omt, th2_t, mu);
            th1t_vect = [th1t_vect, th1_t];
            th2t_vect = [th2t_vect, th2_t];
            dv1 = abs(norm(vt1 - vv_i));
            dv2 = abs(norm(vv_f - vt2));
            dv1_vect = [dv1_vect, dv1];
            dv2_vect = [dv2_vect, dv2];
            dv_direct = dv1 + dv2;
            dt_direct = TOF(at, et, th1_t, th2_t, mu);
            dt_VECT=[dt_VECT;dt_direct];

            if dv_direct < dv_min
                dv_min = dv_direct; % 23.8082 km/s
                k_dvmin = k;
            end
            if dt_direct < dt_min
                dt_min = dt_direct; % 2199.6887 s
                k_dtmin = k;
            end
            if mod(k, 13) == 0
            plotOrbit(at, et, it, OMt, omt, th1_t, th2_t, 1e-3, mu,'g')
            end
        end
    end
end

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on
legend('', 'Initial orbit', 'Final orbit','Initial point','Final Point','Auxiliary orbits')


th1_t_v=th1t_vect(k_dvmin);
th2_t_v=th2t_vect(k_dvmin);
etvm = (r_f - r_i) / ((r_i * cos(th1_t_v)) - (r_f * cos(th2_t_v)));
ptvm = r_i + r_i * et * cos(th1_t_v);
atvm = ptvm / (1 - etvm^2);
dt_vmin = TOF(atvm, etvm, th1_t_v, th2_t_v, mu); % 5149.5967 s

omtvm = om1f + th_i - th1_t_v;
if omtvm < 0
    omtvm = omtvm + 2*pi;
end

figure
planet3D(planet, opts);
light('Position',[1,-1,0]);
hold on
plotOrbit(a_i, e_i , i_i, OM_i, om_i, th0, thf, dth, mu,'b')
plotOrbit(a_f, e_f , i_f, OM_f, om_f, th0, thf, dth, mu,'r')
plot3(rr_i(1),rr_i(2),rr_i(3),'.',MarkerSize=10,Color=[0 0 0])
plot3(rr_f(1),rr_f(2),rr_f(3),'.',MarkerSize=10,Color=[0 0 0])
plotOrbit(atvm, etvm, it, OMt, omtvm, th1_t_v, th2_t_v, 1e-3, mu,'g')

th1_t_t=th1t_vect(k_dtmin);
th2_t_t=th2t_vect(k_dtmin);

omttm = om1f + th_i - th1_t_t;

if omttm < 0
    omttm = omttm + 2*pi;
end

ettm = (r_f - r_i) / ((r_i * cos(th1_t_t)) - (r_f * cos(th2_t_t)));
pttm = r_i + r_i * ettm * cos(th1_t_t);
attm = pttm / (1 - ettm^2);
rptm = pttm / (1 + ettm);
[rt1_t, vt1_t] = par2car(attm, ettm, it, OMt, omttm, th1_t_t, mu);
[rt2_t, vt2_t] = par2car(attm, ettm, it, OMt, omttm, th2_t_t, mu);
dv1 = abs(norm(vt1_t - vv_i));
dv2 = abs(norm(vv_f - vt2_t));
dv_TOT = dv1 + dv2; % 27.1790 km/s

plotOrbit(attm, ettm, it, OMt, omttm, th1_t_t, th2_t_t, 1e-3, mu,'m')
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on
legend('', 'Initial orbit', 'Final orbit','Initial point','Final Point','Delta v min orbit', 'Delta t min orbit')