function [H_mat,Final_response] = get_channel(receptores_mimo, trans, parameters)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
alpha_i=parameters.alpha_i;                                                  % Angulo en eje X 
ang_rad=parameters.ang_rad;                                               %semi-ángulo de mitad de potencia del LED
m=parameters.m;                   %número Lambertiano

beta_i = parameters.beta_i; % Angulo Eje Z
% Propiedades Receptores

Ap=parameters.Ap;                                     % Area del PD
eta=parameters.eta;                                  % Refraccion
fov=parameters.fov;                                   % FOV


% Propiedades Paredes
Aw=parameters.Aw;                                     %área del elemento reflectante
pw=parameters.pw;                                    %coeficiente de reflexión del área reflectiva

% Propiedades Obstaculo
W=parameters.W;                                        % width of obstacle
H=parameters.H;                                        % hight of obstacle
gv=parameters.gv;                                 %uniform distribution probability  
fv=parameters.fv;                                   %uniform distribution probability   
es=parameters.es;                                       % epsilon 
g=0.72;                                     %Responsivity
gymma=0.017;                                %Reflection coefficient
p=0.1;                                      % parameter is used in calculation of Gn below equation 3.21 
N=70;                                       %Number de scatters

% Propiedades Tunel
X=parameters.X;                                        %Lenght of mine
Y=parameters.Y;                                        %height of mine

% Propiedades Scattering

kr=parameters.kr;                              %uniform distribution probability
km=parameters.km;                                  %uniform distribution probability
ks=parameters.ks;                                   %uniform distribution probability
N=parameters.N;                                       %Number de scatters

% Propiedades Simulacion
t=parameters.t;                                  % value of time
c=parameters.c;                                  % speed of light
Sampling_time=parameters.Sampling_time;                      %(0.25 nano segundo de muestreo, se disminuye para mayor numero de muestras)
time=parameters.time;                  
t_rise=parameters.t_rise;                              %time rise subida del PD
t_fall=parameters.t_fall;                                %time bajada del PD
h_led=parameters.h_led; %respuesta al impulso del PD

f = 0.5; % Refractive Index


% Angulo de Incidencia
incidencia_radian = incline(trans(1),trans(2),trans(3),receptores_mimo(1),receptores_mimo(2), ...
receptores_mimo(3),receptores_mimo(4),receptores_mimo(5));
incidencia = rad2deg(incidencia_radian);

% Angulo de Irradiancia
irradiancia_radian = rotacion(trans(1),trans(2),trans(3),receptores_mimo(1),receptores_mimo(2), ...
receptores_mimo(3),alpha_i, beta_i);
irradiancia=rad2deg(irradiancia_radian);

% LOS
h_vector=zeros(1,length(time));
tlos = 0;
[m_HLoS1,dl1]=HLoS_direct(trans(1),trans(2),trans(3),receptores_mimo(1), ...
receptores_mimo(2),receptores_mimo(3),Ap,eta,alpha_i,receptores_mimo(4),beta_i,receptores_mimo(5), ...
incidencia_radian,incidencia,m,fov,gv,fv,W,H,X,Y,t,es,c);
ind = round(dl1,9);
tlos = dl1; 
index = find(ind==time);
h_vector(index) = h_vector(index)+m_HLoS1;
Final_responsel=conv(h_led,h_vector);


% WALL 1
a=1; %indice para ir guardando componentes nLoS, delay, angulo de incidencia, angulo de inclinación y angulo de rotación, 
lx=6; ly=3; lz=3.5;  %Area del escenario Minero
Nx=lx*3; Ny=ly*3; Nz=round(lz*3); %numero de grid de cada superficie
dA=lz*ly/(Ny*Nz); % area de cada grid
x=0:lx/Nx:lx; %vector que cubre todo el grid en x
y=0:ly/Ny:ly;  %vector que cubre todo el grid en y
z=0:lz/Nz:lz;  %vector que cubre todo el grid en z

HnLos=zeros(1,3*length(time));  %vector donde se guarda la CIR del nLoS
tem=zeros(1,3*length(time));    %vector donde se guarda el delay del nLoS
incide=zeros(1,3*length(time)); %vector donde se el angulo de indicencia CIR del nLoS
Piw = zeros(1,3*length(time));
Pwj = zeros(1,3*length(time));

for kk=1:Nx+1  
for ll=1:Nz+1
WP1=[x(kk) 0.2 z(ll)];    % ubicación de la pared fija en y se mueven en x y z
r(a) = randi([0 90],1,1); % angulo de rotacion del elemento reflectante aleatorio uniforme
s(a) = randi([0 90],1,1); % angulo de inclinación del elemento reflectante aleatorio uniforme

% Si no se cumple esta condicion, el impacto del NLoS es infimo
if abs(x(kk)-receptores_mimo(1))<=1.5 && z(ll)>=1.5
    Pwj(a) = Shadowing(receptores_mimo(1),receptores_mimo(2),receptores_mimo(3),x(kk),3,z(ll),gv,fv,W,H,X,Y,t,es);
    Piw(a) = Shadowing(trans(1),trans(2),trans(3),x(kk),3,z(ll),gv,fv,W,H,X,Y,t,es);
    incidenciaw_pru=incline(x(kk),0.2,z(ll),receptores_mimo(1),receptores_mimo(2),receptores_mimo(3),receptores_mimo(4),receptores_mimo(5)); %angulo de radiancia del LED al elemento reflectante
 incidenciawpru=rad2deg(incidenciaw_pru);%conversion de radianes a grados

   %Función que calcula el CIR y el delay del nLoS
    [m_HnLoSp,dnp] = HnLos_calculation_opt(trans(1),trans(2),trans(3),receptores_mimo(1), ...
    receptores_mimo(2),receptores_mimo(3),x(kk),0.2,z(ll),dA/70,pw,alpha_i,receptores_mimo(4),r(a),beta_i,receptores_mimo(5), ...
    s(a),Ap,incidenciawpru,incidenciaw_pru,eta,m,fov,gv,fv,W,H,X,Y,t,es,c,Piw(a),Pwj(a));

    incide(a)=incidenciawpru; %guarda angulo de incidencia
    HnLos(a)=m_HnLoSp; %guarda CIR del nLoS
    tem(a)= dnp; %guarda delay del nLoS

    
end
a=a+1;
%a
end
end

h_vector3=zeros(1,length(time)); %vector donde se guarda el CIR del nLoS dependiendo del delay
tem=round(tem,9); %redondeo del vector de delays
indno=zeros(1,length(time));   %vector donde se guardara el indice del delay del nLoS = tiempo
indexno=1;

for l=1:size(HnLos,1)
var=tem(l);  %guarda en una variable el delay
var=round(var,9); %redonde el delay a 9 cifras significativas
indexno=find(time==var); %encuentra el indice cuando el delay es igual al tiempo
if ~isempty(indexno)
    indno(l)=indexno;   %guarda ese indice en el vector de indices
    h_vector3(indexno)=HnLos(l)+h_vector3(indexno); % coloca el CIR del nlos en el indice correspondiente
end
end


Final_responsen=conv(h_led,h_vector3);

% WALL 2
b=1; %indice para ir guardando componentes nLoS, delay, angulo de incidencia, angulo de inclinación y angulo de rotación, 
lx=6; ly=3; lz=3.5;  %Area del escenario Minero
Nx=lx*3; Ny=ly*3; Nz=round(lz*3); %numero de grid de cada superficie
dA=lz*ly/(Ny*Nz); % area de cada grid
x=0:lx/Nx:lx; %vector que cubre todo el grid en x
y=0:ly/Ny:ly;  %vector que cubre todo el grid en y
z=0:lz/Nz:lz;  %vector que cubre todo el grid en z

HnLos1=zeros(1,length(time)*3);  %vector donde se guarda la CIR del nLoS
tem1=zeros(1,length(time)*3);    %vector donde se guarda el delay del nLoS
incide1=zeros(1,length(time)*3); %vector donde se el angulo de indicencia CIR del nLoS

for kk=1:Nx+1  
for ll=1:Nz+1
WP1=[x(kk) 0.2 z(ll)];    % ubicación de la pared fija en y se mueven en x y z
r(b) = randi([0 90],1,1); % angulo de rotacion del elemento reflectante aleatorio uniforme
s(b) = randi([0 90],1,1); % angulo de inclinación del elemento reflectante aleatorio uniforme

% Si no se cumple esta condicion, el impacto del NLoS es infimo
if abs(x(kk)-receptores_mimo(1))<=1.5 && z(ll)>=1.5
    Pwj(b) = Shadowing(receptores_mimo(1),receptores_mimo(2),receptores_mimo(3),x(kk),3,z(ll),gv,fv,W,H,X,Y,t,es);
    Piw(b) = Shadowing(trans(1),trans(2),trans(3),x(kk),3,z(ll),gv,fv,W,H,X,Y,t,es);
    incidenciaw_pru=incline(x(kk),0.2,z(ll),receptores_mimo(1),receptores_mimo(2),receptores_mimo(3),receptores_mimo(4),receptores_mimo(5)); %angulo de radiancia del LED al elemento reflectante
 incidenciawpru=rad2deg(incidenciaw_pru);%conversion de radianes a grados

   %Función que calcula el CIR y el delay del nLoS
    [m_HnLoSp,dnp] = HnLos_calculation_opt(trans(1),trans(2),trans(3),receptores_mimo(1), ...
    receptores_mimo(2),receptores_mimo(3),x(kk),0.2,z(ll),dA/70,pw,alpha_i,receptores_mimo(4),r(b),beta_i,receptores_mimo(5), ...
    s(b),Ap,incidenciawpru,incidenciaw_pru,eta,m,fov,gv,fv,W,H,X,Y,t,es,c,Piw(b),Pwj(b));

    incide1(b)=incidenciawpru; %guarda angulo de incidencia
    HnLos1(b)=m_HnLoSp; %guarda CIR del nLoS
    tem1(b)= dnp; %guarda delay del nLoS

    
end
b=b+1;
%a
end
end

h_vector1=zeros(1,length(time)); %vector donde se guarda el CIR del nLoS dependiendo del delay
tem1=round(tem1,9); %redondeo del vector de delays
indno1=zeros(1,length(time));   %vector donde se guardara el indice del delay del nLoS = tiempo
indexno=1;

for l=1:size(HnLos1,1)
    var=tem1(l);  %guarda en una variable el delay
    var=round(var,9); %redonde el delay a 9 cifras significativas
    indexno=find(time==var); %encuentra el indice cuando el delay es igual al tiempo
if ~isempty(indexno)
    indno1(l)=indexno;   %guarda ese indice en el vect
    h_vector1(indexno)=h_vector1(indexno)+HnLos1(l); % coloca el CIR del nlos en el indice correspondiente
end
end

Final_responsen1=conv(h_led,h_vector1);

% SCATTERING
Hsca=zeros(1,length(time));
h_vector2=zeros(1,length(time));
temsca=zeros(1,length(time));

[m_Hscat1,ds1]=H_scater(trans(1),trans(2),trans(3),receptores_mimo(1) ...
,receptores_mimo(2),receptores_mimo(3),Ap,m,f,g,gymma,kr,km,ks,p,N,irradiancia,c,alpha_i,beta_i)   ;
inds=zeros(1,length(time));
indexs=1;
for l=1:length(m_Hscat1)
 Hsca(l)=m_Hscat1(l);
temsca(l)= ds1(l);   
end
for l=1:length(time)
    vars=temsca(l);
    vars=round(vars,9);
    indexs=find(time==vars);
    if ~isempty(indexs)
    inds(l)=indexs;
    h_vector2(indexs)=h_vector2(indexs)+Hsca(l);
    end
end
Final_responses = conv(h_led,h_vector2);


Final_response=Final_responsel+Final_responsen+Final_responsen1+Final_responses;

H_mat = get_dc_gain(Final_response, 1/Sampling_time);

%H_mat1 = squeeze(sum(Final_responsel,4))
%H_mat2 = squeeze(sum(Final_responsen,4))
%H_mat3 = squeeze(sum(Final_responsen1,4))
%H_mat4 = squeeze(sum(Final_responses,4));




end


function DC_mat = get_dc_gain(matrix, fs)
    spectrum = abs(fftshift(fft(matrix)));
    precision = fs/length(matrix);
    ffunction = linspace(-fs/2+precision/2, fs/2-precision/2, length(matrix));
    DC_mat = squeeze(spectrum(abs(ffunction)<1));
end


% function to plot for one impulse
function multi_pulse_plot= graph_draw(dm_total,m_total_Hscat,t1,D1,N)

for i=1:N
    y1=dirac(t1-dm_total(i));
    multi_pulse_plot= plottt(y1,t1,m_total_Hscat(i),D1);
end

end



function ang_inc=incline(x1,y1,z1,x2,y2,z2,alpha,beta)
v=[x1-x2,y1-y2,z1-z2];
v2 = -v/sqrt(sum(v.^2));
beta = 90-beta;
Ntilt=[cosd(alpha)*sind(beta),sind(alpha)*sind(beta),cosd(beta)];
vec_center=[0.1*cosd(alpha)*sind(beta),0.1*sind(alpha)*sind(beta),0.1*cosd(beta)];


d_p=dot_product(v,Ntilt);
d=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
ang_inc= acos(d_p/d);
%alphav = (-(x1-x2)/abs(x1-x2)+1)/2*180+atan(v(2)/v(1))*180/pi;
%betav = (-v(3)/abs(v(3))+1)/2*180+atan(sqrt(v(1)^2+v(2)^2)/v(3))*180/pi;
%ang_inc_og = ang_inc;
discriminator = 2*(vec_center*v2')/(v2*v2');

if discriminator>0
    ang_inc = pi;
end

end
%función para calcular el angulo de irradiancia%
function ang_incidencia=rotacion(x1,y1,z1,x2,y2,z2,alpha,beta)

v=[x1-x2,y1-y2,z1-z2];
Ntilt=[cosd(alpha)*sind(beta),sind(alpha)*sind(beta),-cosd(beta)];
d_p=dot_product(v,Ntilt);
d=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
ang_incidencia=acos(d_p/d);

end


%funcion to calculate the HnLos for three reflection points
function [HnLoS_Total,delta_t] =HnLos_calculation_total(x_i,y_i,z_i,x_j,y_j,z_j,x_w1,y_w1,z_w1,x_w2,y_w2,z_w2,x_w3,y_w3,z_w3,Aw,pw,alpha_i,alpha_j,alpha_w1,alpha_w2,alpha_w3,beta_i,beta_j,beta_w1,beta_w2,beta_w3,Ap,inc1,inc_r1,inc2,inc_r2,inc3,inc_r3,eta,m,fov,gv,fv,W,H,X,Y,t,es,c);       
         
[m_HnLoS1,dm1]=HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x_w3,y_w3,z_w3,Aw,pw,alpha_i,alpha_j,alpha_w3,beta_i,beta_j,beta_w3,Ap,inc1,inc_r1,eta,m,fov,gv,fv,W,H,X,Y,t,es,c);
[m_HnLoS2,dm2]=HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x_w1,y_w1,z_w1,Aw,pw,alpha_i,alpha_j,alpha_w1,beta_i,beta_j,beta_w1,Ap,inc2,inc_r2,eta,m,fov,gv,fv,W,H,X,Y,t,es,c);
[m_HnLoS3,dm3]=HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x_w2,y_w2,z_w2,Aw,pw,alpha_i,alpha_j,alpha_w2,beta_i,beta_j,beta_w2,Ap,inc3,inc_r3,eta,m,fov,gv,fv,W,H,X,Y,t,es,c);
HnLoS_Total=[m_HnLoS1,m_HnLoS2,m_HnLoS3];
delta_t=[dm1,dm2,dm3];
      
end

function [Psh] = Shadowing(x_i,y_i,z_i,x_w,y_w,z_w,gv,fv,W,H,X,Y,t,es)
dv_iw= dv(x_i,y_i,x_w,y_w,fv);
sv_iw= sv(x_i,y_i,z_i,x_w,y_w,z_w,fv);
Psh=P_expt(gv,fv,W,H,X,Y,t,es, dv_iw,sv_iw);
end
%funcion que calcula el HNLoS
function[m_HnLoS,dm]=HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x_w,y_w,z_w,Aw,pw,alpha_i,alpha_j,alpha_w,beta_i,beta_j,beta_w,Ap,inc,inc_r,eta,m,fov,gv,fv,W,H,X,Y,t,es,c)

dv_iw= dv(x_i,y_i,x_w,y_w,fv);
sv_iw= sv(x_i,y_i,z_i,x_w,y_w,z_w,fv);
Piw=P_expt(gv,fv,W,H,X,Y,t,es, dv_iw,sv_iw);

dv_wj= dv(x_w,y_w,x_j,y_j,fv);
sv_wj= sv(x_w,y_w,z_w,x_j,y_j,z_j,fv);
Pwj=P_expt(gv,fv,W,H,X,Y,t,es, dv_wj,sv_wj);


g=gain(eta,inc,inc_r,fov);
   
[v1,d1]=point_to_vector(x_i,y_i,z_i,x_w,y_w,z_w);
Nnorm1=norm_vec_trans(alpha_i,beta_i);
p1=dot_product(v1,Nnorm1);
[v2,d2]=point_to_vector(x_w,y_w,z_w,x_i,y_i,z_i);
Nnorm2=norm_vec_receiver(alpha_w,beta_w);
p2=dot_product(v2,Nnorm2);
[v3,d3]=point_to_vector(x_w,y_w,z_w,x_j,y_j,z_j);
Nnorm3=norm_vec_receiver(alpha_w,beta_w);
p3=dot_product(v3,Nnorm3);
[v4,d4]=point_to_vector(x_j,y_j,z_j,x_w,y_w,z_w);
Nnorm4=norm_vec_receiver(alpha_j,beta_j);
p4=dot_product(v4,Nnorm4);


digits(2);
dm=((d1+d3)/c);
dm=vpa(dm);
dm=double(subs(dm));
m_HnLoS= abs(((m+1)*Ap*Aw*pw*p1*p2*p3*p4*g*Piw*Pwj)/((d1^2)*(d3^2)*d1*d2*d3*d4))
m_HnLoS=vpa(m_HnLoS);
m_HnLoS=double(subs(m_HnLoS));

end

function[m_HnLoS,dm]=HnLos_calculation_opt(x_i,y_i,z_i,x_j,y_j,z_j,x_w,y_w,z_w,Aw,pw,alpha_i,alpha_j,alpha_w,beta_i,beta_j,beta_w,Ap,inc,inc_r,eta,m,fov,gv,fv,W,H,X,Y,t,es,c, Piw, Pwj)
g=gain(eta,inc,inc_r,fov);

[v1,d1]=point_to_vector(x_i,y_i,z_i,x_w,y_w,z_w);
Nnorm1=norm_vec_trans(alpha_i,beta_i);
p1=dot_product(v1,Nnorm1);
[v2,d2]=point_to_vector(x_w,y_w,z_w,x_i,y_i,z_i);
Nnorm2=norm_vec_receiver_wall(alpha_w,beta_w);
p2=dot_product(v2,Nnorm2);
[v3,d3]=point_to_vector(x_w,y_w,z_w,x_j,y_j,z_j);
Nnorm3=norm_vec_receiver_wall(alpha_w,beta_w);
p3=dot_product(v3,Nnorm3);
[v4,d4]=point_to_vector(x_j,y_j,z_j,x_w,y_w,z_w);
Nnorm4=norm_vec_receiver_pd(alpha_j,beta_j);
p4=dot_product(v4,Nnorm4);
digits(2);
dm=((d1+d3)/c);
dm=vpa(dm);
dm=double(subs(dm));
m_HnLoS= abs(((m+1)*Ap*Aw*pw*p1*p2*p3*p4*g*Piw*Pwj)/((d1^2)*(d3^2)*d1*d2*d3*d4));
m_HnLoS=vpa(m_HnLoS);
m_HnLoS=double(subs(m_HnLoS));
%if x_w>2 && x_w<4 && z_w>1.7
%fprintf('z_w: %f, x_w: %f \n', z_w, x_w)
%fprintf('p1: %f, p2: %f, p3: %f, p4: %f, tot: %f \n', p1, p2, p3, p4, p1*p2*p3*p4*1e8)
%fprintf('H_NLos: %f \n', m_HnLoS*1e8) 
%end
end
function[m_HnLoS,dm,debug]=HnLos_calculation_opt_debug(x_i,y_i,z_i,x_j,y_j,z_j,x_w,y_w,z_w,Aw,pw,alpha_i,alpha_j,alpha_w,beta_i,beta_j,beta_w,Ap,inc,inc_r,eta,m,fov,gv,fv,W,H,X,Y,t,es,c, Piw, Pwj)
g=gain(eta,inc,inc_r,fov);
[v1,d1]=point_to_vector(x_i,y_i,z_i,x_w,y_w,z_w);
Nnorm1=norm_vec_trans(alpha_i,beta_i);
p1=dot_product(v1,Nnorm1);
[v2,d2]=point_to_vector(x_w,y_w,z_w,x_i,y_i,z_i);
Nnorm2=norm_vec_receiver_wall(alpha_w,beta_w);
p2=dot_product(v2,Nnorm2);
[v3,d3]=point_to_vector(x_w,y_w,z_w,x_j,y_j,z_j);
Nnorm3=norm_vec_receiver_wall(alpha_w,beta_w);
p3=dot_product(v3,Nnorm3);
[v4,d4]=point_to_vector(x_j,y_j,z_j,x_w,y_w,z_w);
Nnorm4=norm_vec_receiver_pd(alpha_j,beta_j);
p4=dot_product(v4,Nnorm4);
digits(2);
dm=((d1+d3)/c);
dm=vpa(dm);
dm=double(subs(dm));
m_HnLoS= abs(((m+1)*Ap*Aw*pw*p1*p2*p3*p4*g*Piw*Pwj)/((d1^2)*(d3^2)*d1*d2*d3*d4));
m_HnLoS=vpa(m_HnLoS);
m_HnLoS=double(subs(m_HnLoS));
debug = [p1,p2,p3,p4,m_HnLoS,z_w,x_w];
end
%función para calcular el HLoS%
function [m_HLoS,dm]=HLoS_direct(x_i,y_i,z_i,x_j,y_j,z_j,Ap,eta,alpha_i,alpha_j,beta_i,beta_j,incidencia,incidencia_r,m,fov,gv,fv,W,H,X,Y,t,es,c)
   dv_ij= dv(x_i,y_i,x_j,y_j,fv);
   sv_ij= sv(x_i,y_i,z_i,x_j,y_j,z_j,fv);
   Pij=P_expt(gv,fv,W,H,X,Y,t,es, dv_ij,sv_ij);
   [v1,d1]=point_to_vector(x_i,y_i,z_i,x_j,y_j,z_j);
   Nnorm1=norm_vec_trans(alpha_i,beta_i);
   p1=dot_product(v1,Nnorm1);
   [v2,d2]=point_to_vector(x_j,y_j,z_j,x_i,y_i,z_i);
   Nnorm2=norm_vec_receiver_pd(alpha_j,beta_j);
   p2=dot_product(v2,Nnorm2);
   g=gain(eta,incidencia,incidencia_r,fov);
   digits(2);
   dm= d1/c;
   dm=vpa(dm);
   dm=double(subs(dm));
   if (incidencia>=0) && (incidencia<=2*fov)
      m_HLoS=abs(((m+1)*Ap/(2*3.1416*d1^2))*(p1^m/d1)*(p2/d2)*g* Pij);
      m_HLoS=vpa(m_HLoS);
      m_HLoS=double(subs(m_HLoS));
   else
      m_HLoS=0;
   end
end


% Function to calculate the  Hscat for 40 scatering point
function [m_total_Hscat,dm_total]=H_scater(x_i,y_i,z_i,x_j,y_j,z_j,Ap,m,f,g,gymma,kr,km,ks,p,N,theta_ij,c,alpha_i,beta_i)
         [v1,dij]=point_to_vector(x_i,y_i,z_i,x_j,y_j,z_j);
         
         dm_total=zeros(1,(N+4));
         m_total_Hscat=zeros(1,N);
         for i =1:N
            
             Rr=0.5;
             rn=Rr*rand(1,1);
             theta_sn_j=randi([-180 180]);
             B_i_sn=Bisn(theta_sn_j,theta_ij);
             
             
             xs=rn*cosd(B_i_sn);
             ys=rn*sind(B_i_sn);
             zs=rn*cosd(theta_sn_j);
             phi_i_sn_radian=phi_scater(xs,ys,zs,x_i,y_i,z_i,alpha_i,beta_i);
             phi_i_sn=rad2deg(phi_i_sn_radian);
             
             di_sn=sqrt(rn^2+dij^2-2*rn*dij*cosd(B_i_sn));
             Di_j=di_sn+rn;
             Gn=Gain_n(f,g,gymma,phi_i_sn,kr,km,ks,p,N);
               
             digits(2);
             dm=(Di_j/c);
             dm=vpa(dm);
             dm=double(subs(dm));
             dm_total(i)=+dm;
             
             if (theta_sn_j>=-180) && (theta_sn_j<=180)
                 Hscat=abs(((m+1)*Ap*Gn/(2*3.1416*Di_j^2))*(cosd(phi_i_sn))^m*cosd(theta_sn_j));
                 Hscat=vpa(Hscat);
                 Hscat=double(subs(Hscat));
             else
                Hscat=0;
             end
             m_total_Hscat(i)=+Hscat;
        
         
         end
end






%Function to calculate the  Pij for shadowing model
function Pij=P_expt(gv,fv,W,H,X,Y,t,es,d_v,s_v)
         syms w h x y p E
         
         if(gv(1)>=2*d_v) &&(gv(2)>=s_v)
             w_int=int(gv(1),w,0,W);
             h_int=int(gv(2),h,0,H);
             A=[w_int h_int];
             x_int=int(fv(1),x,0,X);
             y_int=int(fv(2),y,0,Y);
             B=[x_int y_int];
             exp_value=dot(A,B);
             f=p*t;
             est=-es*exp_value;
             d=vpa(subs(f,p,est),8);
             f=exp(E);
             Pij=vpa(subs(f,E,d),4);
         else
             Pij=0;
         end
         
         
end

%function to calculate the d(xv,yv)  
function d_v= dv(x1,y1,x2,y2,fv)
          d_v= abs((y1-y2)*fv(1)-(x1-x2)*fv(2)-x2*y1+x1*y2)/sqrt((y1-y2)^2+(x1-x2)^2);
end


%function to calculate the s(xv,yv) 
function s_v= sv(x1,y1,z1,x2,y2,z2,fv)
         if(z1<=z2)
             s_v=((y1-y2)^2+(x1-x2)^2+(fv(1)-x1)^2+(fv(2)-y1)^2-[(fv(1)-x2)^2+(fv(2)-y2)^2]/2*sqrt((y1-y2)^2+(x1-x2)^2))+z1;
         else
             s_v=((y1-y2)^2+(x2-x1)^2+(fv(1)-x2)^2+(fv(2)-y2)^2-[(fv(1)-x1)^2+(fv(2)-y1)^2]/2*sqrt((y1-y2)^2+(x2-x1)^2))+z2;
             
         end
end


%funcion para calcular la ganancia
function g=gain(eta,incide,incide_r,fov)
    fov_radian = pi*fov/180;
  if  (incide_r<=fov) && (incide_r>=0)
      g=(eta^2)/(sin(fov_radian)^2);
  else
      g=0;
  end
end


% funcion para calcular el producto punto
function f=dot_product(a,b)
      f=dot(a,b);
end

% Calculo del vector normal para receptores
function m=norm_vec_receiver_wall(alpha,beta)
      m=[cosd(alpha)*sind(beta),sind(alpha)*sind(beta),cosd(beta)];
end

function m=norm_vec_receiver_pd(alpha,beta)
      m=[cosd(alpha)*sind(beta),sind(alpha)*sind(beta),cosd(beta)];
end

% Calculo del vector normal para transmisores
function n=norm_vec_trans(alpha,beta)
      n=[cosd(alpha)*sind(beta),sind(alpha)*sind(beta),-cosd(beta)];
end


%función para calcular vector y distancia entre dos puntos
function [Vec,lenght]= point_to_vector(x1,y1,z1,x2,y2,z2)
      Vec=[x2-x1,y2-y1,z2-z1];
      lenght=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
      
end


%calculation of Gn gain
function Gn=Gain_n(f,g,gymma,phi_i_sn,kr,km,ks,p,N)
        p_mie=pmie(f,g,phi_i_sn);
        p_ray=pray(gymma,phi_i_sn);
        p_total=(kr/ks)*p_ray+(km/ks)*p_mie;
        f_scat=p_total*sind(phi_i_sn);
        Gn=p*f_scat/N;
end

%calculation of pmie
function p_mie=pmie(f,g,phi_i_sn)
         p_mie=(1-g^2/4*pi)*(1/(1+g^2-2*g*cosd(phi_i_sn))^1.5+ f*3*(cosd(phi_i_sn))^2-1/2*(1+g^2)^1.5);
end


%calculation of pray
function p_ray=pray(gymma,phi_i_sn)
         p_ray=3*[1+3*gymma+(1-gymma)*(cosd(phi_i_sn))^2]/(16*pi*(1+2*gymma));
end

%scatering angles phi and theta
function phi_i_sn=phi_scater(x1,y1,z1,x2,y2,z2,alpha,beta)
v=[x1-x2,y1-y2,z1-z2];
Ntilt=[cosd(alpha)*sind(beta),sind(alpha)*sind(beta),-cosd(beta)];
d_p=dot_product(v,Ntilt);
d=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
phi_i_sn=acos(d_p/d);
end

%scatering angle beta 
function B_i_sn=Bisn(theta_sn_j,theta_ij)
if(theta_ij<theta_sn_j)
    B_i_sn=theta_sn_j-theta_ij;
else
   B_i_sn= theta_ij-theta_sn_j;
end
end

% function to plot for one impulse
function plott= plottt(y1,t,dij,D)
         d = y1;
         idx = d == Inf; % find Inf
         d(idx) = dij;
         plott=plot3(t,D*ones(size(t)),d);
end

% Distance function
function D= DIST(x1,y1,z1,x2,y2,z2)
      
      D=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
      
end