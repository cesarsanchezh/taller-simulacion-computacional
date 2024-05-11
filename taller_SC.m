clear;
clc;

Cm = 1e-6; 
gNa = 120e-3; 
gK = 36e-3; 
gL = 0.3e-3; 
ENa = 50e-3; 
EK = -77e-3; 
EL = -54.4e-3; 

tSpan = [0, 100]; 
t = linspace(tSpan(1), tSpan(2), 1000); 
dt = t(2) - t(1); 
N = length(t); 

Iext = zeros(N, 1); 
Iext(100:300) = 10; 

u = -65 * ones(N, 1); 
m = 0.05 * ones(N, 1); 
h = 0.6 * ones(N, 1); 
n = 0.32 * ones(N, 1); 

alpha_m = @(u) (2.5 - 0.1 * (u + 65)) ./ (exp(2.5 - 0.1 * (u + 65)) - 1);
beta_m = @(u) 4 * exp(-(u + 65) / 18);
alpha_h = @(u) 0.07 * exp(-(u + 65) / 20);
beta_h = @(u) 1 ./ (exp((30 - (u + 65)) / 10) + 1);
alpha_n = @(u) (0.01 * (u + 65 - 10)) ./ (exp((10 - (u + 65)) / 10) - 1);
beta_n = @(u) 0.125 * exp(-(u + 65) / 80);

for i = 1:N - 1
    INa = gNa * (m(i)^3) * h(i) * (u(i) - ENa);
    IK = gK * (n(i)^4) * (u(i) - EK);
    IL = gL * (u(i) - EL);
    
    dudt = (Iext(i) - INa - IK - IL) / Cm;
    
    u(i + 1) = u(i) + dudt * dt;
    
    m(i + 1) = m(i) + (alpha_m(u(i)) * (1 - m(i)) - beta_m(m(i))) * dt;
    h(i + 1) = h(i) + (alpha_h(u(i)) * (1 - h(i)) - beta_h(h(i))) * dt;
    n(i + 1) = n(i) + (alpha_n(u(i)) * (1 - n(i)) - beta_n(n(i))) * dt;
end

figure;
subplot(2, 1, 1);
plot(t, u);
xlabel('Tiempo (ms)');
ylabel('Potencial de Membrana (mV)');
title('Potencial de Membrana en el Modelo de Hodgkin-Huxley');

subplot(2, 1, 2);
plot(t, Iext);
xlabel('Tiempo (ms)');
ylabel('Corriente Externa');
title('Corriente Externa Aplicada');
