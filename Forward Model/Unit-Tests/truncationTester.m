addpath('..');
clear all, close all, clc

rTest = 0.7 * 0.3830;   % 0.7 * R -- a genuinely interior point, scaled to THIS R
thetaTest = pi/4;

eta_M50 = deflection(rTest, thetaTest, 13.85, 4.6985e-5, 1.4548e-3, 0.3830, 0.3, 50, 10, 10);
eta_M70 = deflection(rTest, thetaTest, 13.85, 4.6985e-5, 1.4548e-3, 0.3830, 0.3, 70, 15, 15);

fprintf('M=50,P=10,N=10: %.6f%+.6fi\n', real(eta_M50), imag(eta_M50));
fprintf('M=70,P=15,N=15: %.6f%+.6fi\n', real(eta_M70), imag(eta_M70));
fprintf('diff = %.3e\n', abs(eta_M50 - eta_M70));