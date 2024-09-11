%% MAIN
% Propiedades 

%Propiedades Transmisores
% [X,Y,Z]
trans = [4,0.5,3];  
% Angulo en eje Z
parameters.beta_i=45;                                     
% Angulo en eje X 
parameters.alpha_i=90;                                    
%semi-ángulo de mitad de potencia del LED
parameters.ang_rad = 60;          
%número Lambertiano
parameters.m= -log(2)/log(abs(cos(parameters.ang_rad*pi/180))); 

% Propiedades Receptores

% Area del PD
parameters.Ap = 0.0001;        
% Refraccion
parameters.eta = 1.5;
% FOV
parameters.fov = 60;                  
% Radio Esfera
r = 0.05;                  
% Angulo Esfera Eje X
angle = 30;                 
% Angulo Esfera en Radianes
ele = 60*pi/180; 
% Angulo Esfera en Grados
eled = 60;                    

% Propiedades Paredes

%Area del elemento reflectante
parameters.Aw=1;        
%Coeficiente de reflexión del área reflectiva
parameters.pw=0.6;                                   

% Propiedades Obstaculo
% Ancho (Width)
parameters.W=6;     
% Altura (Height)
parameters.H=2.5;    
% Distribucion Uniforme 
parameters.gv=[20 60];       
% Distribucion Uniforme 
parameters.fv=[0 5];                    
% Epsilon
parameters.es=5;                                        

% Propiedades Tunel
% Largo
parameters.X=6;               
% Altura
parameters.Y=3;                                        

% Propiedades Scattering
% Distribucion Uniforme
parameters.kr=[0.1 0.01];        
% Distribucion Uniforme
parameters.km=[0 10];                    
% Distribucion Uniforme
parameters.ks=parameters.kr+parameters.km;  
% Numero de Scatters
parameters.N=70;    

% Propiedades Simulacion
% Tiempo
parameters.t=5*10^-9;      
% Velocidad de la Luz
parameters.c=3*10^8 ;    
%(0.25 nano segundo de muestreo, se disminuye para mayor numero de muestras)
parameters.Sampling_time=0.25e-9;                      
% Vector de Tiempo
parameters.time=0:parameters.Sampling_time:35e-9;    
% Redondeo
parameters.time=round(parameters.time,12);   
% Time Rise
parameters.t_rise=0.5e-9;      
% Time Low
parameters.t_fall=1e-9;                                
% Impulso PD
parameters.h_led=10*(exp(-parameters.time/parameters.t_fall)-exp(-parameters.time/(parameters.t_rise))); 


%% Calculo de H y h
cX = 3;
cY = 1;
receiver_center = [cX, cY, 1.8];
PD = PD_Position(receiver_center, r, angle, ele, eled);
[Hmat, Final_response] = get_channel(PD, trans, parameters);

%%
function PDs = PD_Position(center, r, angle, ele, eled)


    PDs = [center(1)+r*cosd(angle(1))*sin(ele(1)), ...
        center(2)+r*sind(angle(1))*sin(ele(1)), ...
        center(3)+r*cos(ele(1)),angle(1),90-eled(1)];

end


