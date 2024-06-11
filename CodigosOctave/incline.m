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
## @deftypefn {} {@var{retval} =} incline (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: anl <anl@anl-HP-Pavilion-Gaming-Laptop-16-a0xxx>
## Created: 2024-06-05

function incidencia_radian = incline (x_i, y_i, z_i, x_j, y_j, z_j, alpha_j, beta_j)
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
endfunction
