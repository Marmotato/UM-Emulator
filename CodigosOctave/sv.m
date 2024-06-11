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
## @deftypefn {} {@var{retval} =} sv (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: anl <anl@anl-HP-Pavilion-Gaming-Laptop-16-a0xxx>
## Created: 2024-06-06

%function to calculate the s(xv,yv) 
function s_v= sv(x1,y1,z1,x2,y2,z2,fv)
         if(z1<=z2)
             s_v=((y1-y2)^2+(x1-x2)^2+(fv(1)-x1)^2+(fv(2)-y1)^2-[(fv(1)-x2)^2+(fv(2)-y2)^2]/2*sqrt((y1-y2)^2+(x1-x2)^2))+z1;
         else
             s_v=((y1-y2)^2+(x2-x1)^2+(fv(1)-x2)^2+(fv(2)-y2)^2-[(fv(1)-x1)^2+(fv(2)-y1)^2]/2*sqrt((y1-y2)^2+(x2-x1)^2))+z2;
             
         end
endfunction
