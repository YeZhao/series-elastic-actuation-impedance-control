function F = criticaldamp_both_zeta_1(x,omega_1,omega_2,zeta1,zeta2)
% parameters;
utsea_v2_OL;

x = exp(x);
Kx = x(1);
Bx = x(2);
Ktau = x(3);
Btau = x(4);

% Left Hand Side
LHS1 = (IL * bM + IM * bL + IL * beta1 * Btau * k)/(IM * IL);
X = IL * (1 + beta1 * Ktau) + IM + beta1 * Btau * (bL + Bx);
LHS2 = (X * k + bL * bM)/(IM * IL);
Y = bM + beta1 * Btau * Kx;
LHS3 = ((bL + Bx) * k * (1 + beta1 * Ktau) + Y * k)/(IM * IL);
                            
                            
LHS4 = (1 + beta1 * Ktau) * k * Kx/(IM * IL);

f1 = LHS1 - 2 * (zeta1 * omega_1 + zeta2 * omega_2);
f2 = LHS2 - omega_1^2 - omega_2^2 - 4 * zeta1 * zeta2 * omega_1 * omega_2;
f3 = LHS3 - 2 * (zeta2 * omega_1 + zeta1 * omega_2) * omega_1 * omega_2;
f4 = LHS4 - (omega_1 * omega_2)^2;
F=[f1; f2; f3; f4];