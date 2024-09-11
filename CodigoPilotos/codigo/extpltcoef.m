% Remove pilot coefficients -- We actually do not use this function.
function [ch chp]=extpltcoef(ct,N,nb)
m=1;
ch=[];
for i=1:N:length(ct)
chp(m)=ct(i); %extracting only pilot symbols' coeeficient
ch=[ch ct(i+1:i+N-1)];
m=m+1;
end
end