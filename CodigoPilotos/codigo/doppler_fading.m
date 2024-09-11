function ct = doppler_fading(len, n_ref, vel, freq)
vel = vel/3.6;
lambda = 3e8/freq;
fmax=vel/lambda; %Max doppler shift
A=1; %amplitude
f=10000; %sampling frequency
t=0:1/f:((len/10000)-(1/f)); %sampling time
ct=zeros(1,len);
ph=2*pi* rand(1,n_ref);
theta=2*pi*rand(1,n_ref);
fd=fmax*cos(theta); %doppler shift
for n=1:len
    for i=1:n_ref
        ct(n)=ct(n)+(A*exp(j*(2*pi*fd(i)*t(n)+ph(i))));
    end
end
ct=ct/sqrt(n_ref); %channel coefficient
end



% function ct = doppler_fading(len, n_ref, vel, freq)
% lambda = 3e8/freq; %longitud de onda
% vel = vel/3.6;
% fmax = vel/lambda; %(fmax = v/lambda)
% t = linspace(0,1,len);
% an = ones(1,n_ref)*sqrt(1/n_ref); %potencia normalizada de los paths
% thetan = 2*pi*rand(1,n_ref); %fase aleatoria del path
% fDn = fmax * cos( 2*pi*rand(1,n_ref) ); %Doppler aleatorio
% ct = zeros(1,len);
% 
% for n = 1:n_ref
% ct = ct + an(n)*exp( j*( thetan(n) - 2*pi*fDn(n)*t ) );
% end
% 
% end