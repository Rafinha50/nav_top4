clc;
clear all;
load trabalho_missao_3.mat;

# Coordenadas
lat_dec = geodesica_ini(1);
long_dec = geodesica_ini(2);
h_dec = geodesica_ini(3);

# Conversao Geodesica -> Cartesiano
a = 6378137;
e = 0.0818;
N0 = a/(sqrt(1-e^2*(sind(lat_dec))^2));
x0(1) = (N0+h_dec)*cosd(lat_dec)*cosd(long_dec);
y0(1) = (N0+h_dec)*cosd(lat_dec)*sind(long_dec);
z0(1) = (N0*(1-e^2)+h_dec)*sind(lat_dec);

rx(1) = 0;
ry(1) = 0;
rz(1) = 0;

#Determinacao das velocidades iniciais
m = length(ax_ms2);
vx = zeros(1,m);
vy = zeros(1,m);
vz = zeros(1,m);

vx(1) = 0;
vy(1) = 0;
vz(1) = 0;

time = zeros(1,m);

# Descobrir a rotacao da atitude
m = length(ax_ms2);
phi = zeros(1,m);
theta = zeros(1,m);
psi = zeros(1,m);

# Determinacao da atitude inicial
phi(1) = atitude_ini(1);
theta(1) = atitude_ini(2);
psi(1) = atitude_ini(3);

for i = 1 : m
  #0 Determinacao da matriz de rotacao do girometro
  Rb2n = [ 1 sind(phi(i))*tand(theta(i)) cosd(phi(i))*tand(theta(i));
           0        cosd(phi(i))                -sind(phi(i));
           0 sind(phi(i))/cosd(theta(i)) cosd(phi(i))/cosd(theta(i))];

  #1 Determinacao das derivadas dos angulos de atitude
  d_ang = Rb2n * [wx_degs(i);
                  wy_degs(i);
                  wz_degs(i)];

  #2 Determinacao dos angulos phi e theta segundo os girometros
  phi(i+1) = phi(i) + d_ang(1) * dt(1);
  theta(i+1) = theta(i) + d_ang(2) * dt(1);

  #3 Determinacao do psi com o mag
  Cb2n = [cosd(theta(i)) sind(phi(i))*sind(theta(i)) cosd(phi(i))*sind(theta(i));
          0                       cosd(phi(i))             -sind(phi(i));
          -sind(theta(i)) sind(phi(i))*cosd(theta(i)) cosd(phi(i))*cosd(theta(i))]; #matriz de rotacao para os girometros com psi=0

  #4 Determinacao dos ms para encontrar o psi
  m = Cb2n * [mx_uT(i);
              my_uT(i);
              mz_uT(i)];

  #5 Determinacao do psi verdadeiro com a declinacao
  psi(i+1) = atan2d(-m(2) , m(1)) + declina(i);

  #6 Determinacao das atitudes finais com os accs
  Cn2b = [cosd(psi(i))*cosd(theta(i)) cosd(psi(i))*sind(theta(i))*sind(phi(i))-sind(psi(i))*cosd(phi(i)) cosd(psi(i))*sind(theta(i))*cosd(phi(i))+sind(psi(i))*sind(phi(i));
          sind(psi(i))*cosd(theta(i)) sind(psi(i))*sind(theta(i))*sind(phi(i))+cosd(psi(i))*cosd(phi(i)) sind(psi(i))*sind(theta(i))*cosd(phi(i))-cosd(psi(i))*sind(phi(i));
          -sind(theta(i))             cosd(theta(i))*sind(phi(i))                                        cosd(theta(i))*cosd(phi(i))]; #matriz de rotacao para os aceler?metros

  #7 Determinacao das aceleracaes
  acc_n = Cn2b * [ax_ms2(i);
              ay_ms2(i);
              az_ms2(i)];

  #8 Retirando a gravidade do az
  acc_n(3) = acc_n(3) - g(i);

  #9 Determinacao das velocidades
  vx(i+1) = vx(i) + acc_n(1) * dt(i);
  vy(i+1) = vy(i) + acc_n(2) * dt(i);
  vz(i+1) = vz(i) + acc_n(3) * dt(1);

  #10 Determinacao das posicaes
  rx(i+1) = rx(i) + vx(i+1) * dt(i);
  ry(i+1) = ry(i) + vy(i+1) * dt(i);
  rz(i+1) = rz(i) + vz(i+1) * dt(i);

  time(i+1) = time(i) + dt(i);
endfor

ep2 = e^2/(1-e^2);
b = a*sqrt(1-e^2);

R1 = [cosd(geodesica_ini(1)) , sind(geodesica_ini(1)) , 0 ;
     -sind(geodesica_ini(1)) , cosd(geodesica_ini(1)) , 0 ;
          0      ,     0      , 1 ];

R2 = [-sind(geodesica_ini(2)) ,   0   ,  cosd(geodesica_ini(2)) ;
                0             ,   1   ,           0             ;
      -cosd(geodesica_ini(2)) ,   0   ,  -sind(geodesica_ini(2))];

Re2n = R2*R1;
ECEF_0 = [x0 ; y0 ; z0];
ENU = [ ry(end) ; rx(end) ; -rz(end) ];

ECEF_fin = Re2n'*ENU + ECEF_0;
x = ECEF_fin(1);
y = ECEF_fin(2);
z = ECEF_fin(3);

long_final = atan2d(y,x);
P = sqrt(x^2 + y^2);
theta_geof = atan2d(z*a,P*b);
lat_final = atan2d((z + (ep2*b*(sind(theta_geof))^3)) , (P - (e^2*a*(cosd(theta_geof))^3)));
N_final = a/sqrt(1 - (e^2*sind(lat_final))^2);

h_final = P/cosd(lat_final) - N_final;

aeroporto_final = [ lat_final ; long_final ; h_final ];

# Plot posicao cartesiana
figure;
subplot(3,1,1);
plot(time,rx);
ylabel('R_x (m)');
subplot(3,1,2);
plot(time,ry);
ylabel('R_y (m)');
subplot(3,1,3);
plot(time,rz);
ylabel('R_z (m)');






