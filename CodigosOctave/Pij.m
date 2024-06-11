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
## @deftypefn {} {@var{retval} =} Pij (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: anl <anl@anl-HP-Pavilion-Gaming-Laptop-16-a0xxx>
## Created: 2024-06-06
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
endfunction
