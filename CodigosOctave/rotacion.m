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
## @deftypefn {} {@var{retval} =} rotacion (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: anl <anl@anl-HP-Pavilion-Gaming-Laptop-16-a0xxx>
## Created: 2024-06-05

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
endfunction
