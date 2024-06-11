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
## @deftypefn {} {@var{retval} =} HLoS_direct (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: anl <anl@anl-HP-Pavilion-Gaming-Laptop-16-a0xxx>
## Created: 2024-06-06


function [m_HLoS,dm]=HLoS_direct(x_i,y_i,z_i,x_j,y_j,z_j,Ap,eta,alpha_i,alpha_j,beta_i,beta_j,incidencia,incidencia_r,m,fov,gv,fv,W,H,X,Y,t,es,c)
   dv_ij= dv(x_i,y_i,x_j,y_j,fv);
   sv_ij= sv(x_i,y_i,z_i,x_j,y_j,z_j,fv);
   Pij=P_expt(gv,fv,W,H,X,Y,t,es, dv_ij,sv_ij);
   [v1,d1]=point_to_vector(x_i,y_i,z_i,x_j,y_j,z_j);
   Nnorm1=norm_vec_trans(alpha_i,beta_i);
   p1=dot_product(v1,Nnorm1);
   [v2,d2]=point_to_vector(x_j,y_j,z_j,x_i,y_i,z_i);
   Nnorm2=norm_vec_receiver(alpha_j,beta_j);
   p2=dot_product(v2,Nnorm2);
   g=gain(eta,incidencia,incidencia_r,fov);
   digits(2);
   dm= d1/c;
   dm=vpa(dm);
   dm=double(subs(dm));

   if (incidencia>=0) && (incidencia<=2*fov)
      m_HLoS=abs(((m+1)*Ap/(2*3.1416*d1^2))*(p1^m/d1)*(p2/d2)*g* Pij);
      m_HLoS=vpa(m_HLoS);
      m_HLoS=double(subs(m_HLoS));
   else
      m_HLoS=0;
   end
endfunction

