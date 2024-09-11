% definicion de 16QAM
A1 = sqrt(1/10);
A2 = 3*sqrt(1/10);
pilot=A1+1j*A1; % pilot symbol
inphase = [A2 A2 A2 A2 A1 A1 A1 A1];
quadr = [-A2 -A1 A1 A2 -A2 -A1 A1 A2];
inphase = [inphase;-inphase]; inphase = inphase(:);
quadr = [quadr;quadr]; quadr = quadr(:);
const = inphase + j*quadr;

% grafico
re_const = real(const);
im_const = imag(const);
figure
scatter(re_const, im_const, 'filled', 'MarkerFaceColor',[0.5 0 0.7])
xlabel("Parte Real", "FontSize", 13)
ylabel("Parte Imaginaria", "FontSize",13)
title("Constelación para 16QAM", "FontSize", 14)
grid('on')
hold on;
% Eje x (Parte Real)
plot([-1 1], [0 0], 'k--', 'LineWidth', 1);
% Eje y (Parte Imaginaria)
plot([0 0], [-1 1], 'k--', 'LineWidth', 1);
hold off;

% definicion para QPSK
A = sqrt(1/2);
pilot = A+1j*A;
inphase = [A A -A -A];
quadr = [A -A A -A];
const = inphase + j*quadr;

% grafico
re_const = real(const);
im_const = imag(const);
figure
scatter(re_const, im_const, 'filled', 'MarkerFaceColor',[0.5 0 0.7])
xlabel("Parte Real", "FontSize", 13)
ylabel("Parte Imaginaria", "FontSize", 13)
title("Constelación para QPSK", "FontSize", 14)
grid('on')
hold on
% Eje x (Parte Real)
plot([-1 1], [0 0], 'k--', 'LineWidth', 1);
% Eje y (Parte Imaginaria)
plot([0 0], [-1 1], 'k--', 'LineWidth', 1);
hold off;
