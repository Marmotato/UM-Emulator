## Copyright (C) 2024 anl
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} los (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: anl <anl@anl-HP-Pavilion-Gaming-Laptop-16-a0xxx>
## Created: 2024-06-05
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
endfunction