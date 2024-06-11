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
## @deftypefn {} {@var{retval} =} P_expt (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: anl <anl@anl-HP-Pavilion-Gaming-Laptop-16-a0xxx>
## Created: 2024-06-06
%Function to calculate the  Pij for shadowing model
function Pij=P_expt(gv,fv,W,H,X,Y,t,es,d_v,s_v)
  
    p=0.1;
         
    if (gv(1) >= 2 * d_v) && (gv(2) >= s_v)
        w_int = gv(1) * W / 2;
        h_int = gv(2) * H / 2;
        A = [w_int, h_int];
        x_int = fv(1) * X / 2;
        y_int = fv(2) * Y / 2;
        B = [x_int, y_int];
        exp_value = dot(A, B);
        f = p * t;
        est = -es * exp_value;
        d = exp(est);
        Pij = d;
    else
        Pij = 0;
    end
         
endfunction
