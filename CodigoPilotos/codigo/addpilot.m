function [txp]=addpilot(tx,nb,pilot,N) %Number of symbols per frame (Interpolation Factor)
m=1;k=0;
for i=1:N-1:nb
txp(m)=pilot;
txp(m+1:m+N-1)=tx(i:i+N-2);
m=m+N;
k=k+1;
end
end