




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
time = round(time); % redondeo del tiempo a 12 cifras significativas
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
incidenciaw12 = rad2deg(incidenciaw_radian12);

incidenciaw_radian22 = incline(x_w2, y_w2, z_w2, x_j2, y_j2, z_j2, alpha_j, beta_j);
incidenciaw22 = rad2deg(incidenciaw_radian22);

incidenciaw_radian32 = incline(x_w, y_w, z_w, x_j2, y_j2, z_j2, alpha_j, beta_j);
incidenciaw32 = rad2deg(incidenciaw_radian32);

% 3rd position of receiver
incidenciaw_radian13 = incline(x_w1, y_w1, z_w1, x_j3, y_j3, z_j3, alpha_j, beta_j);
incidenciaw13 = rad2deg(incidenciaw_radian13);

incidenciaw_radian23 = incline(x_w2, y_w2, z_w2, x_j3, y_j3, z_j3, alpha_j, beta_j);
incidenciaw23 = rad2deg(incidenciaw_radian23);

incidenciaw_radian33 = incline(x_w, y_w, z_w, x_j3, y_j3, z_j3, alpha_j, beta_j);
incidenciaw33 = rad2deg(incidenciaw_radian33);

% 4th position of receiver
incidenciaw_radian14 = incline(x_w1, y_w1, z_w1, x_j4, y_j4, z_j4, alpha_j, beta_j);
incidenciaw14 = rad2deg(incidenciaw_radian14);

incidenciaw_radian24 = incline(x_w2, y_w2, z_w2, x_j4, y_j4, z_j4, alpha_j, beta_j);
incidenciaw24 = rad2deg(incidenciaw_radian24);

incidenciaw_radian34 = incline(x_w, y_w, z_w, x_j4, y_j4, z_j4, alpha_j, beta_j);
incidenciaw34 = rad2deg(incidenciaw_radian34);

% 5th position of receiver
incidenciaw_radian15 = incline(x_w1, y_w1, z_w1, x_j5, y_j5, z_j5, alpha_j, beta_j);
incidenciaw15 = rad2deg(incidenciaw_radian15);

incidenciaw_radian25 = incline(x_w2, y_w2, z_w2, x_j5, y_j5, z_j5, alpha_j, beta_j);
incidenciaw25 = rad2deg(incidenciaw_radian25);

incidenciaw_radian35 = incline(x_w, y_w, z_w, x_j5, y_j5, z_j5, alpha_j, beta_j);
incidenciaw35 = rad2deg(incidenciaw_radian35);



%% Calculo Respuesta al Impulso

% Llamada a las funciones para calcular las contribuciones LoS y nLoS de la respuesta al impulso
[h1,t1]=los(time,x_i,y_i,z_i,x_j,y_j,z_j,m,Ap,eta,fov,t,h_led);
[h2,t2]=los(time,x_i,y_i,z_i,x_j2,y_j2,z_j2,m,Ap,eta,fov,t,h_led);
[h3,t3]=los(time,x_i,y_i,z_i,x_j3,y_j3,z_j3,m,Ap,eta,fov,t,h_led);
[h4,t4]=los(time,x_i,y_i,z_i,x_j4,y_j4,z_j4,m,Ap,eta,fov,t,h_led);
[h5,t5]=los(time,x_i,y_i,z_i,x_j5,y_j5,z_j5,m,Ap,eta,fov,t,h_led);

[t_los_r1,h_los_r1] = los(time,x_i,y_i,z_i,x_j,y_j,z_j,m,Ap,eta,fov,t,h_led);
[t_los_r2,h_los_r2] = los(time,x_i,y_i,z_i,x_j2,y_j2,z_j2,m,Ap,eta,fov,t,h_led);
[t_los_r3,h_los_r3] = los(time,x_i,y_i,z_i,x_j3,y_j3,z_j3,m,Ap,eta,fov,t,h_led);
[t_los_r4,h_los_r4] = los(time,x_i,y_i,z_i,x_j4,y_j4,z_j4,m,Ap,eta,fov,t,h_led);
[t_los_r5,h_los_r5] = los(time,x_i,y_i,z_i,x_j5,y_j5,z_j5,m,Ap,eta,fov,t,h_led);

[t_nlos_r1,h_nlos_r1] = nlos(time,x_i,y_i,z_i,x_w1,y_w1,z_w1,x_w,y_w,z_w,x_w2,y_w2,z_w2,x_j,y_j,z_j,m,Ap,eta,fov,t,h_led,Aw,pw,alpha_w,beta_w);
[t_nlos_r2,h_nlos_r2] = nlos(time,x_i,y_i,z_i,x_w1,y_w1,z_w1,x_w,y_w,z_w,x_w2,y_w2,z_w2,x_j2,y_j2,z_j2,m,Ap,eta,fov,t,h_led,Aw,pw,alpha_w,beta_w);
[t_nlos_r3,h_nlos_r3] = nlos(time,x_i,y_i,z_i,x_w1,y_w1,z_w1,x_w,y_w,z_w,x_w2,y_w2,z_w2,x_j3,y_j3,z_j3,m,Ap,eta,fov,t,h_led,Aw,pw,alpha_w,beta_w);
[t_nlos_r4,h_nlos_r4] = nlos(time,x_i,y_i,z_i,x_w1,y_w1,z_w1,x_w,y_w,z_w,x_w2,y_w2,z_w2,x_j4,y_j4,z_j4,m,Ap,eta,fov,t,h_led,Aw,pw,alpha_w,beta_w);
[t_nlos_r5,h_nlos_r5] = nlos(time,x_i,y_i,z_i,x_w1,y_w1,z_w1,x_w,y_w,z_w,x_w2,y_w2,z_w2,x_j5,y_j5,z_j5,m,Ap,eta,fov,t,h_led,Aw,pw,alpha_w,beta_w);

% Contribución total de la respuesta al impulso (LoS + nLoS)
h_total_r1 = h_los_r1 + h_nlos_r1;
h_total_r2 = h_los_r2 + h_nlos_r2;
h_total_r3 = h_los_r3 + h_nlos_r3;
h_total_r4 = h_los_r4 + h_nlos_r4;
h_total_r5 = h_los_r5 + h_nlos_r5;


%% GRAFICOS

% Gráfica de las respuestas al impulso para cada receptor
figure;

subplot(2,3,1);
plot(t1,h1,'r');
hold on;
plot(t_nlos_r1,h_nlos_r1,'g');
plot(time,h_total_r1,'b');
title('PD1');
xlabel('Time (s)');
ylabel('h(t)');
legend('LoS','NLoS','Total');

subplot(2,3,2);
plot(t2,h2,'r');
hold on;
plot(t_nlos_r2,h_nlos_r2,'g');
plot(time,h_total_r2,'b');
title('PD2');
xlabel('Time (s)');
ylabel('h(t)');
legend('LoS','NLoS','Total');

subplot(2,3,3);
plot(t3,h3,'r');
hold on;
plot(t_nlos_r3,h_nlos_r3,'g');
plot(time,h_total_r3,'b');
title('PD3');
xlabel('Time (s)');
ylabel('h(t)');
legend('LoS','NLoS','Total');

subplot(2,3,4);
plot(t4,h4,'r');
hold on;
plot(t_nlos_r4,h_nlos_r4,'g');
plot(time,h_total_r4,'b');
title('PD4');
xlabel('Time (s)');
ylabel('h(t)');
legend('LoS','NLoS','Total');

subplot(2,3,5);
plot(t5,h5,'r');
hold on;
plot(t_nlos_r5,h_nlos_r5,'g');
plot(time,h_total_r5,'b');
title('PD5');
xlabel('Time (s)');
ylabel('h(t)');
legend('LoS','NLoS','Total');



%% FUNCIONES

% Función incline

function incidencia_radian = incline(x_i, y_i, z_i, x_j, y_j, z_j, alpha_j, beta_j)
    % Vector desde el transmisor al receptor
    dx = x_j - x_i;
    dy = y_j - y_i;
    dz = z_j - z_i;
    
    % Vector normal del receptor
    nx = sind(beta_j) * cosd(alpha_j);
    ny = sind(beta_j) * sind(alpha_j);
    nz = cosd(beta_j);
    
    % Vector unitario desde el transmisor al receptor
    mag = sqrt(dx^2 + dy^2 + dz^2);
    ux = dx / mag;
    uy = dy / mag;
    uz = dz / mag;
    
    % Ángulo de incidencia (en radianes)
    incidencia_radian = acosd(nx*ux + ny*uy + nz*uz);
end

% Funcion rotacion


function irradiancia_radian = rotacion(x_j, y_j, z_j, x_i, y_i, z_i, alpha_i, beta_i)
    % Vector desde el receptor al transmisor
    dx = x_i - x_j;
    dy = y_i - y_j;
    dz = z_i - z_j;
    
    % Vector normal del transmisor
    nx = sind(beta_i) * cosd(alpha_i);
    ny = sind(beta_i) * sind(alpha_i);
    nz = cosd(beta_i);
    
    % Vector unitario desde el receptor al transmisor
    mag = sqrt(dx^2 + dy^2 + dz^2);
    ux = dx / mag;
    uy = dy / mag;
    uz = dz / mag;
    
    % Ángulo de irradiancia (en radianes)
    irradiancia_radian = acosd(nx*ux + ny*uy + nz*uz);
end


% Funcion los

function [h_los, t_los] = los(time, x_i, y_i, z_i, x_j, y_j, z_j, m, Ap, eta, fov, t, h_led)
    % Distancia entre el transmisor y el receptor
    d = sqrt((x_j - x_i)^2 + (y_j - y_i)^2 + (z_j - z_i)^2);
    
    % Ángulo de incidencia
    theta = acosd((z_j - z_i) / d);
    
    % Ángulo de irradiancia
    phi = theta;
    
    % Verificar si está dentro del campo de visión
    if theta <= fov
        % Respuesta al impulso de LoS
        h_los = (m + 1) * Ap / (2 * pi * d^2) * cosd(phi)^m * cosd(theta);
    else
        h_los = 0;
    end
    
    % Convolución con la respuesta al impulso del LED
    h_los = h_los * conv(h_led, exp(-time/t));
    t_los = time;
end


% Funcion nlos


function [h_nlos, t_nlos] = nlos(time, x_i, y_i, z_i, x_w1, y_w1, z_w1, x_w, y_w, z_w, x_w2, y_w2, z_w2, x_j, y_j, z_j, m, Ap, eta, fov, t, h_led, Aw, pw, alpha_w, beta_w)
    % Distancia entre el transmisor y el punto reflectante
    d1 = sqrt((x_w - x_i)^2 + (y_w - y_i)^2 + (z_w - z_i)^2);
    
    % Distancia entre el punto reflectante y el receptor
    d2 = sqrt((x_j - x_w)^2 + (y_j - y_w)^2 + (z_j - z_w)^2);
    
    % Ángulo de incidencia en la pared
    theta_w = acosd((z_w - z_i) / d1);
    
    % Ángulo de irradiancia desde la pared al receptor
    phi_w = acosd((z_j - z_w) / d2);
    
    % Verificar si está dentro del campo de visión
    if phi_w <= fov
        % Respuesta al impulso de nLoS
        h_nlos = (m + 1) * Ap * Aw * pw / (2 * pi * d1^2 * d2^2) * cosd(theta_w)^m * cosd(phi_w);
    else
        h_nlos = 0;
    end
    
    % Convolución con la respuesta al impulso del LED
    h_nlos = h_nlos * conv(h_led, exp(-time/t));
    t_nlos = time;
end





