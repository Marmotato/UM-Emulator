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
## @deftypefn {} {@var{retval} =} Hnlos_calculation (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: anl <anl@anl-HP-Pavilion-Gaming-Laptop-16-a0xxx>
## Created: 2024-06-06

%funcion que calcula el HNLoS
function[m_HnLoS,dm]=HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x_w,y_w,z_w,Aw,pw,alpha_i,alpha_j,alpha_w,beta_i,beta_j,beta_w,Ap,inc,inc_r,eta,m,fov,gv,fv,W,H,X,Y,t,es,c)
   
dv_iw= dv(x_i,y_i,x_w,y_w,fv);
sv_iw= sv(x_i,y_i,z_i,x_w,y_w,z_w,fv);
Piw=P_expt(gv,fv,W,H,X,Y,t,es, dv_iw,sv_iw);
   
dv_wj= dv(x_w,y_w,x_j,y_j,fv);
sv_wj= sv(x_w,y_w,z_w,x_j,y_j,z_j,fv);
Pwj=P_expt(gv,fv,W,H,X,Y,t,es, dv_wj,sv_wj);

g=gain(eta,inc,inc_r,fov);
   
[v1,d1]=point_to_vector(x_i,y_i,z_i,x_w,y_w,z_w);
Nnorm1=norm_vec_trans(alpha_i,beta_i);
p1=dot_product(v1,Nnorm1);
[v2,d2]=point_to_vector(x_w,y_w,z_w,x_i,y_i,z_i);
Nnorm2=norm_vec_receiver(alpha_w,beta_w);
p2=dot_product(v2,Nnorm2);
[v3,d3]=point_to_vector(x_w,y_w,z_w,x_j,y_j,z_j);
Nnorm3=norm_vec_receiver(alpha_w,beta_w);
p3=dot_product(v3,Nnorm3);
[v4,d4]=point_to_vector(x_j,y_j,z_j,x_w,y_w,z_w);
Nnorm4=norm_vec_receiver(alpha_j,beta_j);
p4=dot_product(v4,Nnorm4);
    
digits(2);
dm=((d1+d3)/c);
dm=vpa(dm);
dm=double(subs(dm));
m_HnLoS= abs(((m+1)*Ap*Aw*pw*p1*p2*p3*p4*g*Piw*Pwj)/((d1^2)*(d3^2)*d1*d2*d3*d4));
m_HnLoS=vpa(m_HnLoS);
m_HnLoS=double(subs(m_HnLoS));
   
endfunction