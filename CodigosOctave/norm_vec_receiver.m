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
## @deftypefn {} {@var{retval} =} norm_vec_receiver (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: anl <anl@anl-HP-Pavilion-Gaming-Laptop-16-a0xxx>
## Created: 2024-06-06

function m=norm_vec_receiver(alpha,beta)
      m=[cosd(alpha)*sind(beta),sind(alpha)*sind(beta),cosd(beta)];
endfunction
