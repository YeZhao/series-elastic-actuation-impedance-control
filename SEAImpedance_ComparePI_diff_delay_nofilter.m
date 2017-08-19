clear all;
clc
% clf

% parameters;
utsea_v2_OL;

%%%%%%%%%%%%%%%%%%%%%%
% tuning parameters
%%%%%%%%%%%%%%%%%%%%%%
zeta1 = 1;
zeta2 = 1;
offset = 0;

fnarray = 30;%stable

s = tf('s');
Tfv = 0;
Qv = 1/(Tfv * s + 1);

Tftau = 0;
Qtau = 1/(Tftau * s + 1);

delta = 1;

TqsArray = [0.001, 0.005, 0.01];
TqdArray = [0.001, 0.002, 0.003];
TtauArray = [0.001, 0.001, 0.001];

for i = 1 : 4
    for j = 1 : 1
        
        if i == 1
            Tqs = 0;
            Tqd = 0;
            Ttau = 0;
        else
            Tqs = TqsArray(i-1);
            Tqd = TqdArray(i-1);
            Ttau = TtauArray(i-1);
        end
        
        % fsolve function
        x0 = log([10000; 1000; 1; 0.1]);
        
        omega1 = 2 * pi * fnarray;
        omega2 = delta * omega1;
        
        options = optimset('MaxFunEvals', 200000, 'MaxIter', 100000, 'TolX', 1e-4, 'TolFun', 1e-4, 'Display','iter');
        [x,fval] = fsolve(@(x) criticaldamp_both_zeta_1(x,omega1,omega2,zeta1,zeta2), x0, options);
        x = exp(x);
        
        Kq = x(1);
        Bq = x(2);
        Ktau = x(3);
        Btau = x(4);
        
        KqArray(i) = Kq;
        BqArray(i) = Bq;
        KtauArray(i) = Ktau;
        BtauArray(i) = Btau;

        H1 = tf([omega1^2],[1, 2*zeta1*omega1, omega1^2]);
        H2 = tf([omega2^2],[1, 2*zeta2*omega2, omega2^2]);
        H_ideal(i) = H1 * H2;
                
        fourth_order_coeff = IM * IL/k;
        third_order_coeff = (IL * bM + IM * bL)/k + ...
            IL * beta1 * Btau * Qtau * ss(exp(-Ttau * s));
        second_order_coeff = IL * (1 + ss(exp(-Ttau * s)) * beta1 * Ktau)...
            + bL * beta1 * Btau * Qtau * ss(exp(-Ttau * s)) ...
            + beta1 * Btau * Bq * ss(exp(-Tqd * s)) * Qv * Qtau ...
            + IM + bL * bM/k;
        first_order_coeff = bL * (1 + ss(exp(-Ttau * s)) * beta1 * Ktau) ...
            + ss(exp(-Tqd * s))  * (1 + beta1 * Ktau) * Bq * Qv ...
            + beta1 * Btau * Kq * ss(exp(-Tqs * s)) * Qtau + bM;
        const_coeff = ss(exp(-Tqs * s)) * (1 + beta1 * Ktau) * Kq;
        
        num = (1 + beta1 * Ktau + beta1 * Btau * s) * (Bq * s + Kq);
        den = fourth_order_coeff * s^4 + third_order_coeff * s^3 + ...
            second_order_coeff * s^2 + first_order_coeff * s + const_coeff;
        
        num1 = (1 + beta1 * Ktau) * Kq;
        H_CL(i,j) = num1/den;
        
        C = Ktau + Btau * Qtau * s;
        Omega = ss(exp(-Tqd * s)) * Bq * Qv * s + ss(exp(-Tqs * s)) * Kq;
        Gamma = C * ss(exp(-Ttau * s)) + beta1^(-1);
        
        num_OL_second = beta1 * Btau * Bq * ss(exp(-Tqd * s)) * Qv * Qtau;
        num_OL_first = ss(exp(-Tqd * s)) * (1 + beta1 * Ktau) * Bq * Qv + beta1 * Btau * Qtau * Kq * ss(exp(-Tqs * s));
        num_OL_zero =  ss(exp(-Tqs * s)) * (1 + beta1 * Ktau) * Kq;
        num_OL = num_OL_second * s^2 + num_OL_first * s + num_OL_zero;
        
        den_OL = den - num_OL;
        H_OL(i,j) = num_OL/den_OL;
        
        % SEA Impedance TF
        N4 = IM * Tftau * Tfv * beta1 * k;
        N3 = beta1 * k * (IM * (Tftau + Tfv) + Tftau * Tfv * bM);
        N2 = IM * beta1 * k + beta1 * k * bM * (Tftau + Tfv) + k * ktau * ...
             (Tftau + beta1 * (Btau + Ktau * Tftau)) * (Bq * ss(exp(-Tqd * s)) + Kq * Tfv * ss(exp(-Tqs * s)));
        N1 = bM * beta1 * k + k * ktau * (Bq + Bq * Ktau * beta1) * ss(exp(-Tqd * s)) ...
             + Kq * k * ktau * (Tfv + Tftau + beta1 * (Btau + Ktau * (Tftau + Tfv))) * ss(exp(-Tqs * s));
        N0 = Kq * k * ktau * (1 + Ktau * beta1) * ss(exp(-Tqs * s));
        
        D5 = IM * Tftau * Tfv * beta1;
        D4 = IM * beta1 * (Tfv + Tftau) + Tfv * Tftau * bM * beta1;
        D3 = beta1 * IM + beta1 * bM * (Tfv + Tftau) + Tfv * beta1 * k * (Tftau + ktau * ...
            (Btau + Ktau * Tftau) * ss(exp(-Ttau * s)));
        D2 = beta1 * (bM + Tftau * k + k * ktau * (Btau + Ktau * Tftau) * ss(exp(-Ttau * s))) ...
            + Tfv * beta1 * k * (1 + Ktau * ktau * ss(exp(-Ttau * s)));
        D1 = beta1 * k * (1 + Ktau * ktau * ss(exp(-Ttau * s)));
        D0 = 0;
        
        N4Array(i) = N4;
        N3Array(i) = N3;
        N2Array(i) = N2;
        N1Array(i) = N1;
        N0Array(i) = N0;
        
        D5Array(i) = D5;
        D4Array(i) = D4;
        D3Array(i) = D3;
        D2Array(i) = D2;
        D1Array(i) = D1;
        
        IMTF(i,j) = (N4 * s^4 + N3 * s^3 + N2 * s^2 + N1 * s + N0)/(D5 * s^5 + D4 * s^4 + D3 * s^3 + D2 * s^2 + D1 * s);
        % convert from linear to rotary
        IMTF(i,j) = r_joint^2 * IMTF(i,j);
        
        % Vallery PI torque control + inner motor velocity feedback
        
        Ki = 0.002;
        Ii = 0.0005;
        Ko = Ktau;
        Io = Ktau/50;
        P = Kq;
        
        a = Ii * Ko + Io * Ki;
                
        num_PI = k * (IM * s^4 + Ki * s^3 + (Ki * Ko * P + Ii) * s^2 + a * P * s + Ii * Io * P);
        den_PI = s * (IM * s^4 + Ki * s^3 + (Ki * Ko * k + k + Ii) * s^2 + a * k * s + Ii * Io * k);
        IMTF_PI(i,j) = num_PI/den_PI;
        
        r = (IL * s^2 + bL * s)/(IL * s^2 + bL * s + k);
        P_X = 1/(IL * s^2 + bL * s);
        P_CF = (1 + beta1 * Ktau)/(1 + beta1 * Ktau * ss(exp(-Ttau * s)));
        
        P_F = beta1 * r * k/(IM * s^2 + bM * s + r * k);
        
        [Gm,Pm,Wgm,Wpm] = margin(H_OL(i));
        PM(i,j) = Pm;
    end
end


figure(7)
% rotary spring stiffness
SpringStiffness = r_joint^2 * k/s;
VirtualStiffness = r_joint^2 * ((Kq * ktau * (1 + Ktau * beta1))/(beta1 * (1 + Ktau * ktau)))/s;%
hold on;
h = bodeplot(SpringStiffness, 'y--', VirtualStiffness, 'b--', IMTF(1,1), 'r--',IMTF(2,1), 'k-',IMTF(3,1), 'm-',IMTF(4,1), 'g-');
grid
p = getoptions(h);
p.FreqUnits = 'Hz';
p.Xlabel.FontSize = 14;
p.Xlabel.FontWeight = 'bold';
p.Ylabel.FontSize = 14;
p.Ylabel.FontWeight = 'bold';
p.Title.FontSize = 11;
p.TickLabel.FontSize = 13;
p.TickLabel.FontWeight = 'bold';
% for rotary stable
p.XLim = [0.4 100];
p.YLim = {[-20 25]; [-100 0]};
h_legend = legend('Physical Stiffness k/j\omega', 'Virtual Stiffness K_{vir}/j\omega', 'No Delay', ...
    'T_{qs} = T_{qd} = T_{\tau} = 1 ms', 'T_{qs} = 5 ms, T_{qd} = 2 ms, T_{\tau} = 1 ms',...
    'T_{qs} = 10 ms, T_{qd} = 3 ms, T_{\tau} = 1 ms');
set(h_legend,'FontSize',9,'Location','NorthEast','FontWeight', 'bold');
setoptions(h,p)
h = findobj(gcf,'type','line');
set(h,'linewidth',2);
hold on;
title('SEA Impedance with different delays','FontSize',20,'FontWeight', 'bold');
% print -depsc Ideal_SEA_Impedance
% print -depsc Practical_SEA_Impedance

KqArray = KqArray * r_joint^2
BqArray = BqArray * r_joint^2
KtauArray = KtauArray/r_joint
BtauArray = BtauArray/r_joint
