%% FUNCIONES

% Función incline


%% DEFINICION TRANSMISORES Y RECEPTORES

%posición Ti
x_i=3;
y_i=0.5;
z_i=4.5;

%posición Rj1 
x_j=3;
y_j=1;
z_j=1.8;

%posición del Rj2     
x_j2=3;
y_j2=1.5;
z_j2=1.8;

%posición del Rj3    
x_j3=4;
y_j3=2;
z_j3=1.8;

%posición del Rj4    
x_j4=2.2;
y_j4=2.5;
z_j4=1.8;

%posición del Rj5    
x_j5=1;
y_j5=2.5;
z_j5=1.6;

%% PARAMETROS DE  SIMULACION

%posición del área reflectiva (pared) usar para puntos específicos de elementos reflectantes
% % first reflection_point 
x_w1=2;
y_w1=2;
z_w1=3;
% % second reflection_point
x_w2=1.5;
y_w2=2.3;
z_w2=4.0;
% % third reflection_point
x_w=5;
y_w=2.9;
z_w=3.5;

Aw=1; % área del elemento reflectante
pw=0.6; % coeficiente de reflexión del área reflectiva
ang_rad = 60; % semi-ángulo de mitad de potencia del LED
m = -log(2)/log(abs(cos(ang_rad*pi/180))); % número Lambertiano
Ap = 0.0001; % área física del receptor (1 cm^2)
eta = 1.5; % índice de refracción del PD
fov = 70; % campo de visión

% parámetros de los ángulos de inclinación y rotación
beta_i = 45; % ángulo de inclinación del LED con respecto al eje z
alpha_i = 45; % ángulo de rotación del LED con respecto al eje x

beta_j = 45; % ángulo de inclinación del PD con respecto al eje z
alpha_j = 45; % ángulo de rotación del PD con respecto al eje x

% wall rotation angles para ángulo específico
alpha_w = 20;
beta_w = 70;

% distribución uniforme de probabilidad
gv = [20 60];
fv = [0 5];
W = 6; % ancho del obstáculo
H = 5; % altura del obstáculo
X = 6; % longitud de la mina
Y = 3; % altura de la mina
es = 5; % epsilon 
t = 5*10^-9; % valor del tiempo
gymma = 0.017; % coeficiente de reflexión
g = 0.72; % responsividad
f = 0.5; % índice de refracción
kr = [0.1 0.01]; % distribución uniforme de probabilidad
km = [0 10]; % distribución uniforme de probabilidad
ks = kr+km; % distribución uniforme de probabilidad
N = 70; % número de dispersores
p = 0.1; % parámetro utilizado en el cálculo de Gn (ecuación 3.21)

c = 3*10^8; % velocidad de la luz

Sampling_time = 0.25e-9; % tiempo de muestreo (0.25 nanosegundos)
time = 0:Sampling_time:35e-9; % vector de tiempo para observar el CIR
time = round(time*1e12)*1e-12; % redondeo del tiempo a 12 cifras significativas
t_rise = 0.5e-9; % tiempo de subida del PD
t_fall = 1e-9; % tiempo de bajada del PD
h_led = 10*(exp(-time/t_fall) - exp(-time/(t_rise))); % respuesta al impulso del PD


%% Cálculo de ángulos de incidencia e irradiancia entre los enlaces 

% Ángulo de incidencia del enlace entre LED y PD para 5 posiciones del PD
% 1st position of receiver
incidencia_radian1 = incline(x_i, y_i, z_i, x_j, y_j, z_j, alpha_j, beta_j);
incidencia1 = rad2deg(incidencia_radian1); % conversión de radianes a grados
% 2nd position of receiver
incidencia_radian2 = incline(x_i, y_i, z_i, x_j2, y_j2, z_j2, alpha_j, beta_j);
incidencia2 = rad2deg(incidencia_radian2);
% 3rd position of receiver
incidencia_radian3 = incline(x_i, y_i, z_i, x_j3, y_j3, z_j3, alpha_j, beta_j);
incidencia3 = rad2deg(incidencia_radian3);
% 4th position of receiver
incidencia_radian4 = incline(x_i, y_i, z_i, x_j4, y_j4, z_j4, alpha_j, beta_j);
incidencia4 = rad2deg(incidencia_radian4);
% 5th position of receiver
incidencia_radian5 = incline(x_i, y_i, z_i, x_j5, y_j5, z_j5, alpha_j, beta_j);
incidencia5 = rad2deg(incidencia_radian5);

% Ángulo de irradiancia del enlace entre LED y PD para 5 posiciones del PD
% 1st position of receiver
irradiancia_radian1 = rotacion(x_j, y_j, z_j, x_i, y_i, z_i, alpha_i, beta_i);
irradiancia1 = rad2deg(irradiancia_radian1); % conversión de radianes a grados
% 2nd position of receiver
irradiancia_radian2 = rotacion(x_j2, y_j2, z_j2, x_i, y_i, z_i, alpha_i, beta_i);
irradiancia2 = rad2deg(irradiancia_radian2);
% 3rd position of receiver
irradiancia_radian3 = rotacion(x_j3, y_j3, z_j3, x_i, y_i, z_i, alpha_i, beta_i);
irradiancia3 = rad2deg(irradiancia_radian3);
% 4th position of receiver
irradiancia_radian4 = rotacion(x_j4, y_j4, z_j4, x_i, y_i, z_i, alpha_i, beta_i);
irradiancia4 = rad2deg(irradiancia_radian4);
% 5th position of receiver
irradiancia_radian5 = rotacion(x_j5, y_j5, z_j5, x_i, y_i, z_i, alpha_i, beta_i);
irradiancia5 = rad2deg(irradiancia_radian5);

% Ángulo de irradiancia del enlace entre LED y pared
irradiancialw_radian = rotacion(x_w, y_w, z_w, x_i, y_i, z_i, alpha_i, beta_i);
irradiancialw = rad2deg(irradiancialw_radian); % conversión de radianes a grados

% Ángulo de incidencia del enlace entre LED y pared
incidencialw_radian = incline(x_i, y_i, z_i, x_w, y_w, z_w, alpha_w, beta_w);
incidencialw = rad2deg(incidencialw_radian); % conversión de radianes a grados

% Ángulo de irradiancia del enlace entre pared y PD
irradianciaw_radian = rotacion(x_j, y_j, z_j, x_w, y_w, z_w, alpha_w, beta_w);
irradianciaw = rad2deg(irradianciaw_radian); % conversión de radianes a grados

% Ángulo de incidencia del enlace de pared a PD
% 1st position of receiver
incidenciaw_radian11 = incline(x_w1, y_w1, z_w1, x_j, y_j, z_j, alpha_j, beta_j);
incidenciaw11 = rad2deg(incidenciaw_radian11); % conversión de radianes a grados

incidenciaw_radian21 = incline(x_w2, y_w2, z_w2, x_j, y_j, z_j, alpha_j, beta_j);
incidenciaw21 = rad2deg(incidenciaw_radian21);

incidenciaw_radian31 = incline(x_w, y_w, z_w, x_j, y_j, z_j, alpha_j, beta_j);
incidenciaw31 = rad2deg(incidenciaw_radian31);

% 2nd position of receiver
incidenciaw_radian12 = incline(x_w1, y_w1, z_w1, x_j2, y_j2, z_j2, alpha_j, beta_j);
incidenciaw12 = rad2deg(incidenciaw_radian12); % conversión de radianes a grados

incidenciaw_radian22 = incline(x_w2, y_w2, z_w2, x_j2, y_j2, z_j2, alpha_j, beta_j);
incidenciaw22 = rad2deg(incidenciaw_radian22);

incidenciaw_radian32 = incline(x_w, y_w, z_w, x_j2, y_j2, z_j2, alpha_j, beta_j);
incidenciaw32 = rad2deg(incidenciaw_radian32);

% 3rd position of receiver
incidenciaw_radian13 = incline(x_w1, y_w1, z_w1, x_j3, y_j3, z_j3, alpha_j, beta_j);
incidenciaw13 = rad2deg(incidenciaw_radian13); % conversión de radianes a grados

incidenciaw_radian23 = incline(x_w2, y_w2, z_w2, x_j3, y_j3, z_j3, alpha_j, beta_j);
incidenciaw23 = rad2deg(incidenciaw_radian23);

incidenciaw_radian33 = incline(x_w, y_w, z_w, x_j3, y_j3, z_j3, alpha_j, beta_j);
incidenciaw33 = rad2deg(incidenciaw_radian33);

% 4th position of receiver
incidenciaw_radian14 = incline(x_w1, y_w1, z_w1, x_j4, y_j4, z_j4, alpha_j, beta_j);
incidenciaw14 = rad2deg(incidenciaw_radian14); % conversión de radianes a grados

incidenciaw_radian24 = incline(x_w2, y_w2, z_w2, x_j4, y_j4, z_j4, alpha_j, beta_j);
incidenciaw24 = rad2deg(incidenciaw_radian24);

incidenciaw_radian34 = incline(x_w, y_w, z_w, x_j4, y_j4, z_j4, alpha_j, beta_j);
incidenciaw34 = rad2deg(incidenciaw_radian34);

% 5th position of receiver
incidenciaw_radian15 = incline(x_w1, y_w1, z_w1, x_j5, y_j5, z_j5, alpha_j, beta_j);
incidenciaw15 = rad2deg(incidenciaw_radian15); % conversión de radianes a grados

incidenciaw_radian25 = incline(x_w2, y_w2, z_w2, x_j5, y_j5, z_j5, alpha_j, beta_j);
incidenciaw25 = rad2deg(incidenciaw_radian25);

incidenciaw_radian35 = incline(x_w, y_w, z_w, x_j5, y_j5, z_j5, alpha_j, beta_j);
incidenciaw35 = rad2deg(incidenciaw_radian35);

% número de valores phi y theta generados en función de la cantidad de dispersores y el FOV
numero_valores = N*fov; % tamaño del vector de valores phi y theta

%% Calcular CIR directo

% distancia directa LED a PD
D_ij1 = sqrt((x_j - x_i)^2 + (y_j - y_i)^2 + (z_j - z_i)^2);
D_ij2 = sqrt((x_j2 - x_i)^2 + (y_j2 - y_i)^2 + (z_j2 - z_i)^2);
D_ij3 = sqrt((x_j3 - x_i)^2 + (y_j3 - y_i)^2 + (z_j3 - z_i)^2);
D_ij4 = sqrt((x_j4 - x_i)^2 + (y_j4 - y_i)^2 + (z_j4 - z_i)^2);
D_ij5 = sqrt((x_j5 - x_i)^2 + (y_j5 - y_i)^2 + (z_j5 - z_i)^2);

% función de transferencia directa entre LED y PD

% 1st position of receiver
H_direct1 = ((m+1)*Ap)/(2*pi*D_ij1^2)*(cos(irradiancia1*pi/180)^m)*T(incidencia1, fov, eta, g)*(cos(incidencia1*pi/180));
H_direct1_total = H_direct1.*h_led;

% 2nd position of receiver
H_direct2 = ((m+1)*Ap)/(2*pi*D_ij2^2)*(cos(irradiancia2*pi/180)^m)*T(incidencia2, fov, eta, g)*(cos(incidencia2*pi/180));
H_direct2_total = H_direct2.*h_led;

% 3rd position of receiver
H_direct3 = ((m+1)*Ap)/(2*pi*D_ij3^2)*(cos(irradiancia3*pi/180)^m)*T(incidencia3, fov, eta, g)*(cos(incidencia3*pi/180));
H_direct3_total = H_direct3.*h_led;

% 4th position of receiver
H_direct4 = ((m+1)*Ap)/(2*pi*D_ij4^2)*(cos(irradiancia4*pi/180)^m)*T(incidencia4, fov, eta, g)*(cos(incidencia4*pi/180));
H_direct4_total = H_direct4.*h_led;

% 5th position of receiver
H_direct5 = ((m+1)*Ap)/(2*pi*D_ij5^2)*(cos(irradiancia5*pi/180)^m)*T(incidencia5, fov, eta, g)*(cos(incidencia5*pi/180));
H_direct5_total = H_direct5.*h_led;

% resuelve el enlace directo para 5 posiciones del PD
H_direct = [H_direct1_total; H_direct2_total; H_direct3_total; H_direct4_total; H_direct5_total];

%% CIR indirecto (reflexión)

% distancia indirecta LED -> pared -> PD

% 1st position of receiver
D_iw1 = sqrt((x_i - x_w1)^2 + (y_i - y_w1)^2 + (z_i - z_w1)^2);
D_w1j = sqrt((x_w1 - x_j)^2 + (y_w1 - y_j)^2 + (z_w1 - z_j)^2);
D_ij_indirect1 = D_iw1 + D_w1j;

D_iw2 = sqrt((x_i - x_w2)^2 + (y_i - y_w2)^2 + (z_i - z_w2)^2);
D_w2j = sqrt((x_w2 - x_j)^2 + (y_w2 - y_j)^2 + (z_w2 - z_j)^2);
D_ij_indirect2 = D_iw2 + D_w2j;

D_iw3 = sqrt((x_i - x_w)^2 + (y_i - y_w)^2 + (z_i - z_w)^2);
D_w3j = sqrt((x_w - x_j)^2 + (y_w - y_j)^2 + (z_w - z_j)^2);
D_ij_indirect3 = D_iw3 + D_w3j;

% 2nd position of receiver
D_w1j2 = sqrt((x_w1 - x_j2)^2 + (y_w1 - y_j2)^2 + (z_w1 - z_j2)^2);
D_ij_indirect4 = D_iw1 + D_w1j2;

D_w2j2 = sqrt((x_w2 - x_j2)^2 + (y_w2 - y_j2)^2 + (z_w2 - z_j2)^2);
D_ij_indirect5 = D_iw2 + D_w2j2;

D_w3j2 = sqrt((x_w - x_j2)^2 + (y_w - y_j2)^2 + (z_w - z_j2)^2);
D_ij_indirect6 = D_iw3 + D_w3j2;

% 3rd position of receiver
D_w1j3 = sqrt((x_w1 - x_j3)^2 + (y_w1 - y_j3)^2 + (z_w1 - z_j3)^2);
D_ij_indirect7 = D_iw1 + D_w1j3;

D_w2j3 = sqrt((x_w2 - x_j3)^2 + (y_w2 - y_j3)^2 + (z_w2 - z_j3)^2);
D_ij_indirect8 = D_iw2 + D_w2j3;

D_w3j3 = sqrt((x_w - x_j3)^2 + (y_w - y_j3)^2 + (z_w - z_j3)^2);
D_ij_indirect9 = D_iw3 + D_w3j3;

% 4th position of receiver
D_w1j4 = sqrt((x_w1 - x_j4)^2 + (y_w1 - y_j4)^2 + (z_w1 - z_j4)^2);
D_ij_indirect10 = D_iw1 + D_w1j4;

D_w2j4 = sqrt((x_w2 - x_j4)^2 + (y_w2 - y_j4)^2 + (z_w2 - z_j4)^2);
D_ij_indirect11 = D_iw2 + D_w2j4;

D_w3j4 = sqrt((x_w - x_j4)^2 + (y_w - y_j4)^2 + (z_w - z_j4)^2);
D_ij_indirect12 = D_iw3 + D_w3j4;

% 5th position of receiver
D_w1j5 = sqrt((x_w1 - x_j5)^2 + (y_w1 - y_j5)^2 + (z_w1 - z_j5)^2);
D_ij_indirect13 = D_iw1 + D_w1j5;

D_w2j5 = sqrt((x_w2 - x_j5)^2 + (y_w2 - y_j5)^2 + (z_w2 - z_j5)^2);
D_ij_indirect14 = D_iw2 + D_w2j5;

D_w3j5 = sqrt((x_w - x_j5)^2 + (y_w - y_j5)^2 + (z_w - z_j5)^2);
D_ij_indirect15 = D_iw3 + D_w3j5;

% retardo de los tiempos
retardo = [D_ij1/c; D_ij2/c; D_ij3/c; D_ij4/c; D_ij5/c; D_ij_indirect1/c; D_ij_indirect2/c; D_ij_indirect3/c; D_ij_indirect4/c; D_ij_indirect5/c; D_ij_indirect6/c; D_ij_indirect7/c; D_ij_indirect8/c; D_ij_indirect9/c; D_ij_indirect10/c; D_ij_indirect11/c; D_ij_indirect12/c; D_ij_indirect13/c; D_ij_indirect14/c; D_ij_indirect15/c];

% función de transferencia indirecta entre LED y PD
rho = 0.8; % coeficiente de reflexión

% 1st position of receiver
H_indirect1 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect1^2)*cos(irradianciaw1*pi/180)^m*T(incidenciaw1, fov, eta, g)*cos(incidenciaw1*pi/180)*cos(irradianciaw1*pi/180);
H_indirect1_total = H_indirect1.*h_led;

H_indirect2 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect2^2)*cos(irradianciaw2*pi/180)^m*T(incidenciaw2, fov, eta, g)*cos(incidenciaw2*pi/180)*cos(irradianciaw2*pi/180);
H_indirect2_total = H_indirect2.*h_led;

H_indirect3 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect3^2)*cos(irradianciaw3*pi/180)^m*T(incidenciaw3, fov, eta, g)*cos(incidenciaw3*pi/180)*cos(irradianciaw3*pi/180);
H_indirect3_total = H_indirect3.*h_led;

% 2nd position of receiver
H_indirect4 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect4^2)*cos(irradianciaw12*pi/180)^m*T(incidenciaw12, fov, eta, g)*cos(incidenciaw12*pi/180)*cos(irradianciaw12*pi/180);
H_indirect4_total = H_indirect4.*h_led;

H_indirect5 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect5^2)*cos(irradianciaw22*pi/180)^m*T(incidenciaw22, fov, eta, g)*cos(incidenciaw22*pi/180)*cos(irradianciaw22*pi/180);
H_indirect5_total = H_indirect5.*h_led;

H_indirect6 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect6^2)*cos(irradianciaw32*pi/180)^m*T(incidenciaw32, fov, eta, g)*cos(incidenciaw32*pi/180)*cos(irradianciaw32*pi/180);
H_indirect6_total = H_indirect6.*h_led;

% 3rd position of receiver
H_indirect7 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect7^2)*cos(irradianciaw13*pi/180)^m*T(incidenciaw13, fov, eta, g)*cos(incidenciaw13*pi/180)*cos(irradianciaw13*pi/180);
H_indirect7_total = H_indirect7.*h_led;

H_indirect8 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect8^2)*cos(irradianciaw23*pi/180)^m*T(incidenciaw23, fov, eta, g)*cos(incidenciaw23*pi/180)*cos(irradianciaw23*pi/180);
H_indirect8_total = H_indirect8.*h_led;

H_indirect9 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect9^2)*cos(irradianciaw33*pi/180)^m*T(incidenciaw33, fov, eta, g)*cos(incidenciaw33*pi/180)*cos(irradianciaw33*pi/180);
H_indirect9_total = H_indirect9.*h_led;

% 4th position of receiver
H_indirect10 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect10^2)*cos(irradianciaw14*pi/180)^m*T(incidenciaw14, fov, eta, g)*cos(incidenciaw14*pi/180)*cos(irradianciaw14*pi/180);
H_indirect10_total = H_indirect10.*h_led;

H_indirect11 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect11^2)*cos(irradianciaw24*pi/180)^m*T(incidenciaw24, fov, eta, g)*cos(incidenciaw24*pi/180)*cos(irradianciaw24*pi/180);
H_indirect11_total = H_indirect11.*h_led;

H_indirect12 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect12^2)*cos(irradianciaw34*pi/180)^m*T(incidenciaw34, fov, eta, g)*cos(incidenciaw34*pi/180)*cos(irradianciaw34*pi/180);
H_indirect12_total = H_indirect12.*h_led;

% 5th position of receiver
H_indirect13 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect13^2)*cos(irradianciaw15*pi/180)^m*T(incidenciaw15, fov, eta, g)*cos(incidenciaw15*pi/180)*cos(irradianciaw15*pi/180);
H_indirect13_total = H_indirect13.*h_led;

H_indirect14 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect14^2)*cos(irradianciaw25*pi/180)^m*T(incidenciaw25, fov, eta, g)*cos(incidenciaw25*pi/180)*cos(irradianciaw25*pi/180);
H_indirect14_total = H_indirect14.*h_led;

H_indirect15 = ((m+1)*Ap*rho*Aw)/(2*pi*D_ij_indirect15^2)*cos(irradianciaw35*pi/180)^m*T(incidenciaw35, fov, eta, g)*cos(incidenciaw35*pi/180)*cos(irradianciaw35*pi/180);
H_indirect15_total = H_indirect15.*h_led;

% resuelve el enlace indirecto para 5 posiciones del PD
H_indirect = [H_indirect1_total; H_indirect2_total; H_indirect3_total; H_indirect4_total; H_indirect5_total; H_indirect6_total; H_indirect7_total; H_indirect8_total; H_indirect9_total; H_indirect10_total; H_indirect11_total; H_indirect12_total; H_indirect13_total; H_indirect14_total; H_indirect15_total];

% respuesta impulso total entre LED y PD
h_total = [H_direct; H_indirect];
time = [0:retardo(1):5e-9]; % tiempo del retardo

% función de transferencia del canal óptico
h = conv(h_total, rectpuls(time));

% plot respuesta impulso
plot(time, h);
title('Respuesta Impulso del Canal');
xlabel('Tiempo (s)');
ylabel('Amplitud');

% Cálculo de la tasa de error de bit (BER)
R = 1e9; % tasa de datos en bits por segundo
N0 = 1e-22; % densidad espectral de ruido
SNR = (sum(h_total.^2)*R)/(N0*B); % relación señal a ruido

BER = qfunc(sqrt(2*SNR)); % función Q para calcular la BER
disp(['La tasa de error de bit es: ', num2str(BER)]);

% Fin del script



