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
## @deftypefn {} {@var{retval} =} nlos (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: anl <anl@anl-HP-Pavilion-Gaming-Laptop-16-a0xxx>
## Created: 2024-06-05
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
endfunction