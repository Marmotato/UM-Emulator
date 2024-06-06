import numpy as np

# posición Ti
x_i = 3
y_i = 0.5
z_i = 4.5

# posición Rj1    
x_j = 3
y_j = 1
z_j = 1.8

# posición del Rj2     
x_j2 = 3
y_j2 = 1.5
z_j2 = 1.8

# posición del Rj3    
x_j3 = 4
y_j3 = 2
z_j3 = 1.8

# posición del Rj4    
x_j4 = 2.2
y_j4 = 2.5
z_j4 = 1.8

# posición del Rj5    
x_j5 = 1
y_j5 = 2.5
z_j5 = 1.6

# posición del área reflectiva (pared) usar para puntos específicos de elementos reflectantes

# first reflection_point 
x_w1 = 2
y_w1 = 2
z_w1 = 3

# second reflection_point
x_w2 = 1.5
y_w2 = 2.3
z_w2 = 4.0

# third reflection_point
x_w = 5
y_w = 2.9
z_w = 3.5

# position of scattering particle
# x_s = 4
# y_s = 4.5
# z_s = 2.5

# parámetros de la simulación

Aw = 1                                        # área del elemento reflectante
pw = 0.6                                      # coeficiente de reflexión del área reflectiva
ang_rad = 60                                  # semi-ángulo de mitad de potencia del LED
m = int(-np.log(2) / np.log(abs(np.cos(ang_rad * np.pi / 180))));   # número Lambertiano
Ap = 0.0001                                   # área física del receptor (1 cm^2)
eta = 1.5                                     # índice de refracción del PD
fov = 70                                      # field of view

# parámetros de los ángulos de inclinación y rotación

beta_i = 45                                   # ángulo de inclinación del LED con respecto al eje z
alpha_i = 45                                  # ángulo de rotación del LED con respecto al eje x

beta_j = 45                                   # ángulo de inclinación del PD con respecto al eje z
alpha_j = 45                                  # ángulo de rotación del PD con respecto al eje x

# wall rotation angles para ángulo específico                   
alpha_w = 20
beta_w = 70

gv = [20, 60]                                 # distribución uniforme de probabilidad  
fv = [0, 5]                                   # distribución uniforme de probabilidad   
W = 6                                         # ancho del obstáculo
H = 5                                         # altura del obstáculo
X = 6                                         # Longitud de la mina
Y = 3                                         # altura de la mina
es = 5                                        # epsilon 
t = 5 * 10**-9                                # valor del tiempo
gymma = 0.017                                 # coeficiente de reflexión
g = 0.72                                      # responsividad
f = 0.5                                       # índice de refracción
kr = [0.1, 0.01]                              # distribución uniforme de probabilidad
km = [0, 10]                                  # distribución uniforme de probabilidad
ks = np.add(kr, km)                           # distribución uniforme de probabilidad
N = 70                                        # Número de scatterers
# rn=sqrt((x_j-x_s)^2+(y_j-y_s)^2+(z_j-z_s)^2);
p = 0.1                                       # el parámetro se utiliza en el cálculo de Gn debajo de la ecuación 3.21 

c = 3 * 10**8                                 # velocidad de la luz

Sampling_time = 0.25e-9                       # 0.25 nano segundos de muestreo, se disminuye para mayor número de muestras
time = np.arange(0, 35e-9, Sampling_time)    # vector de tiempo donde se observará el CIR
time = np.round(time)                         # redondeo del tiempo 
t_rise = 0.5e-9                               # subida del tiempo de subida del PD
t_fall = 1e-9                                 # tiempo de caída del PD
h_led = 10 * (np.exp(-time / t_fall) - np.exp(-time / t_rise))  # respuesta al impulso del PD

