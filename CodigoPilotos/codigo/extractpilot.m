%Extraxting and removing pilot from transmitted signal
function [rx rxp]=extractpilot(t,N,nb)
m=1;
rx=[];
for i=1:N:length(t)
rxp(m)=t(i); %extracting only pilot symbols
rx=[rx t(i+1:i+N-1)]; %extracting data symbols
m=m+1;
end
end