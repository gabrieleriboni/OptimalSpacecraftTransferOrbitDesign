clear
close all
clc

%% SECOND ORBIT CHARACTERIZATION 
mu = 398600;
th0 = 0;
thf = 2*pi;
dth = 1e-5;

a_f = 14160.0;
e_f = 0.2088;
i_f = 1.4870;
OM_f = 1.6100;
om_f = 2.3450;
th_f = 1.5110;

[rr_f, vv_f] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu);

p_f = a_f * (1 - e_f^2);
rp_f = a_f * (1 - e_f);
ra_f = a_f * (1 + e_f);
T_f = 2 * pi *sqrt((a_f^3) / mu);
E = -mu / 2 / a_f;

r_f = norm(rr_f);
hh_f = cross(rr_f,vv_f);
ee_f = (cross(vv_f,hh_f) / mu) - rr_f/r_f;
[rr_p, vv_p] = par2car(a_f, e_f, i_f, OM_f, om_f, 0, mu);
[rr_a, vv_a] = par2car(a_f, e_f, i_f, OM_f, om_f, pi, mu);

planet = 'Earth Cloudy';
opts.Units = 'km';

figure
planet3D(planet, opts);
light('Position',[1,-1,0]);
hold on;
plotOrbit(a_f, e_f , i_f, OM_f, om_f, th0, thf, dth, mu)
plot3(rr_f(1), rr_f(2), rr_f(3), 'o', 'Color',[1, 0.5, 0] , 'LineWidth',3)
plot3(rr_p(1), rr_p(2), rr_p(3), 'gx', 'LineWidth',1)
text(rr_p(1)+1.5e3, rr_p(2)+1.5e3, rr_p(3)+1.5e3, 'P', 'FontSize', 12, 'FontWeight', 'bold')
plot3(rr_a(1), rr_a(2), rr_a(3), 'gx', 'LineWidth',1)
text(rr_a(1)-2e3, rr_a(2)+1.5e3, rr_a(3)-1.5e3, 'A', 'FontSize', 12, 'FontWeight', 'bold')
quiver3(0, 0, 0, hh_f(1)*0.25, hh_f(2)*0.25, hh_f(3)*0.25, 'm', 'LineWidth',2)
quiver3(0, 0, 0, ee_f(1)*0.75e5, ee_f(2)*0.75e5, ee_f(3)*0.75e5, 'b', 'LineWidth',2)
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on
legend('', 'Final orbit', 'Final point','Pericenter', 'Apocenter','Specific angular momentum vector','Eccentricity vector')
leg = legend;
leg.FontSize = 8;
leg.Position = [0.7, 0.7, 0.2, 0.2];