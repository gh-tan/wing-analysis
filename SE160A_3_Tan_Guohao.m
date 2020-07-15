clear, close, clc

format shortG
% format longG

%% load data
data = xlsread('SE160A_Project3_InputFile');

% Stringer Definition (use loop)
% for i = 1:6 (input all values as one)
%     for j = 1:12
%         sg(i,j) = data(i,j);
%     end
% end

% Preallocation for faster computing speed
yi = zeros(1:6);
zi = zeros(1:6);
E = zeros(1:6);
Sy = zeros(1:6);
Su = zeros(1:6);
FSyST = zeros(1:6);
FSuST = zeros(1:6);
A = zeros(1:6);
Iyy = zeros(1:6);
Izz = zeros(1:6);
Iyz = zeros(1:6);
t = zeros(10:15);
G = zeros(10:15);
tau_y = zeros(10:15);
tau_u = zeros(10:15);
FSySK = zeros(10:15);
FSuSK = zeros(10:15);

for i = 1:6
    yi(i) = data(i,2); % (in)
    zi(i) = data(i,3); % (in)
    E(i) = data(i,4); % (Msi)
    Sy(i) = data(i,5); % (Ksi)
    Su(i) = data(i,6); % (Ksi)
    FSyST(i) = data(i,7); % stringer FSy
    FSuST(i) = data(i,8); % stringer FSu
    A(i) = data(i,9); %(in^2)
    Iyy(i) = data(i,10); % (in^4)
    Izz(i) = data(i,11); % (in^4)
    Iyz(i) = data(i,12); % (in^4)
end

% Skin/Spar Definition
% for u = 10:15 (input all values as one)
%     for v = 1:7
%         sk = data(u,v);
%     end
% end
% sk

for j = 10:15
    t(j) = data(j,2); % (in)
    G(j) = data(j,3); % (Msi)
    tau_y(j) = data(j,4); % (Ksi)
    tau_u(j) = data(j,5); % (Ksi)
    FSySK(j) = data(j,6); % Skin FSy
    FSuSK(j) = data(j,7); % Skin FSu
end

% Aerodynamic Properties
L = data(18,5); % (in)
c = data(19,5); % (in)
n = data(20,5);
w = data(21,5); % (lb/in)
pzo = data(22,5); % (ln/in)
pz2 = data(23,5); % (ln/in)
pz4 = data(24,5); % (ln/in)
pyo = data(25,5); % (ln/in)
py2 = data(26,5); % (ln/in)
mxo = data(27,5); % (ln-in/in)

%% Section 1: Wing Cross-Section Properties

% find EA
for i = 1:6
    E_A(i) = E(i)*A(i);
end
EA = sum(E_A); % (Msi-in)
EA_kip = EA*1000; % (kip) , here 1 Msi-in^2 = 1e6 lb = 1e3 kip

% find yc zc
for i = 1:6
    y_c(i) = E(i)*A(i)*yi(i); % (Msi-in)
    z_c(i) = E(i)*A(i)*zi(i); % (Msi-in)
end
yc = sum(y_c)/EA; %(in),in result, need to add -yi(1) in front b/c this yc is the dist from centroid
zc = sum(z_c)/EA; %(in)

% find EIyy,EIzz,EIyz
for i = 1:6
    EIyy(i) = E(i)*Iyy(i) + E(i)*A(i)*(zi(i) - zc)^2; % (Msi-in^4)
    EIzz(i) = E(i)*Izz(i) + E(i)*A(i)*(yi(i) - yc)^2; % (Msi-in^4)
    EIyz(i) = E(i)*Iyz(i) + E(i)*A(i)*(yi(i) - yc)*(zi(i) - zc); % (Msi-in^4)
end
EIyy = sum(EIyy)*1e3; % (kip-in^2)
EIzz = sum(EIzz)*1e3; % (kip-in^2)
EIyz = sum(EIyz)*1e3; % (kip-in^2)

% find GJ
b = (zi(6) - zi(2))/2; % cross section height / 2
a = yi(2)-yi(1); % distance from stringer 1 to the center of ellipse
A_total = (pi*a*b)/2 + (yi(3)-yi(2))*b*2 + 2*0.5*b*(yi(4)-yi(3)); %cross section area
S = 2*1/4*(pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)))) ... 
    + 2*(yi(3)-yi(2)) + 2*sqrt((yi(4)-yi(3))^2+b^2);
GJ = (4*A_total^2*(G(10)*1e3)*t(10))/S;

disp('Modulus Weighted Section Properties')
disp('------------------------------------------------------')
disp('The Centroid is Measured Relative to the Leading Edge')
fprintf('yc:\t%12f\tin\n',yc-yi(1))
fprintf('zc:\t%12f\tin\n',zc)
disp('------------------------------------------------------')
disp('The Inertial properties are w.r.t. the Centroidal Axes')
fprintf('EA:\t%12f\tkip-in^2\n',EA_kip)
fprintf('EIyy:\t%12f\tkip-in^2\n',EIyy)
fprintf('EIzz:\t%12f\tkip-in^2\n',EIzz)
fprintf('EIyz:\t%12f\tkip-in^2\n',EIyz)
fprintf('GJ:\t%12f\tkip-in^2\n',GJ)

%% Section 2: Stress Resultant Diagrams
% input
L;
c;
w;
W = n*w; % W is distributed weight [lb/in]

%% try
% % Distributed load transformed to beam axis
syms x
px = 0;
py = @(x) pyo + py2*(x/L).^2; % (lb/in), Distributed drag load
pz = @(x) pzo + pz2*(x/L).^2 + pz4*(x/L).^4 - W; % Lift distribution act at c/4
% mx = @(x) mxo - W*(-yi(1)+c/2-yc); % moment due to point load only
mx = @(x) mxo - (pzo + pz2*(x/L).^2 + pz4*(x/L).^4)*(-yi(1)+yc-c/4) - W*(yi(1)+c/2-yc);
my = 0;
mz = 0;

%% Ro and Mo
Rxo = 0;
Ryo = integral(py, 0, L);
Rzo = integral(pz, 0, L);

Mxo = integral(mx, 0, L);
pz_times_x = @(x) x.*(pzo + pz2*(x/L).^2 + pz4*(x/L).^4 - W);
Myo = -1*integral(pz_times_x, 0, L);
py_times_x = @(x) x.*(pyo + py2*(x/L).^2);
Mzo = integral(py_times_x, 0, L);

% % correct values
% Rzo = 3.0151e3;
% Mxo = -1.7668e4;
% Myo = -3.5189e5;

%% V(x) and M(x)
xvalue = 0:L;
figure(1)

Vx(x) = x*0; % correct
subplot(2,3,1)
plot(xvalue,Vx(xvalue))
xlabel('x[in]'); ylabel('Vx[lb]')
set(gca,'XLim',[0 L])
grid on

Vy(x) = Ryo - int(py(x)); % correct
subplot(2,3,2)
plot(xvalue,Vy(xvalue))
xlabel('x[in]'); ylabel('Vy[lb]')
title('Stress Resultants Vs. Length')
set(gca,'XLim',[0 L])
grid on

Vz(x) = Rzo - int(pz(x));
subplot(2,3,3)
plot(xvalue,Vz(xvalue))
xlabel('x[in]'); ylabel('Vz[lb]')
set(gca,'XLim',[0 L])
grid on

Mx(x) = Mxo - int(mx(x));
subplot(2,3,4)
plot(xvalue,Mx(xvalue))
xlabel('x[in]'); ylabel('Mx[lb-in]')
set(gca,'XLim',[0 L])
grid on

My(x) = Myo + int(Vz(x));
subplot(2,3,5)
plot(xvalue,My(xvalue))
xlabel('x[in]'); ylabel('My[lb-in]')
set(gca,'XLim',[0 L])
grid on

Mz(x) = Mzo - int(Vy(x)); % correct
subplot(2,3,6)
plot(xvalue,Mz(xvalue))
xlabel('x[in]'); ylabel('Mz[lb-in]')
set(gca,'XLim',[0 L])
grid on

% syms x
% m_x(x) = mxo - (yc-c/4)*(pzo + pz2*(x/L).^2 + pz4*(x/L).^4) - W*(c/2-yc)
% M_x(x) = Mxo - int(m_x(x))
% 
% V_z = Rzo - int(mx);


%% Section 3: Stress and Margin of Safety in Wing Stringers

% Plot of stress vs. wing position for each stringer
% 1) calculate Kyy, Kzz, Kyz
Kyy = EIyy/(EIyy*EIzz-EIyz^2);
Kzz = EIzz/(EIyy*EIzz-EIyz^2);
Kyz = EIyz/(EIyy*EIzz-EIyz^2);
% 2) find y and z centeoidal coordinate for each stringers
for i = 1:6
    yST(i) = yi(i) - yc; % ys = y for string in centroidal coordinate
    zST(i) = zi(i) - zc; % same as above
end
yST(i);
zST(i);

% 3) construct matrix and plot
figure(2)
m2 = [1/EA 0 0;
    0 Kyy -Kyz;
    0 -Kyz Kzz];
m3 = [Vx(x);Mz(x);-My(x)];
% for i = 1:6
%     m1 = [1 -yST(i) -zST(i)];
%     sigmaxx(x) = (E(i)*10^3)*m1*m2*m3;
%     plot(xvalue,sigmaxx(xvalue))
%     hold on
% end
% 
% % sigmaxx(x)
% 
% hold off
% grid on
% title('Stringer Axial Stress Vs. Length')
% legend('Stringer 1','Stringer 2','Stringer 3','Stringer 4' ...
%     ,'Stringer 5','Stringer 6')

% correct till here

% make a new one without using for loop
m1_1 = [1 -yST(1) -zST(1)];
sigmaxx1(x) = (E(1)*10^3)*m1_1*m2*m3;
plot(xvalue,sigmaxx1(xvalue))
hold on

m1_2 = [1 -yST(2) -zST(2)]; %
sigmaxx2(x) = (E(2)*10^3)*m1_2*m2*m3;
plot(xvalue,sigmaxx2(xvalue))


m1_3 = [1 -yST(3) -zST(3)];
sigmaxx3(x) = (E(3)*10^3)*m1_3*m2*m3;
plot(xvalue,sigmaxx3(xvalue))


m1_4 = [1 -yST(4) -zST(4)];
sigmaxx4(x) = (E(4)*10^3)*m1_4*m2*m3;
plot(xvalue,sigmaxx4(xvalue))


m1_5 = [1 -yST(5) -zST(5)];
sigmaxx5(x) = (E(5)*10^3)*m1_5*m2*m3;
plot(xvalue,sigmaxx5(xvalue))

m1_6 = [1 -yST(6) -zST(6)]; %
sigmaxx6(x) = (E(6)*10^3)*m1_6*m2*m3;
plot(xvalue,sigmaxx6(xvalue))
hold off

grid on
title('Stringer Axial Stress Vs. Length')
xlabel('x [in]'); ylabel('Stringer \sigma_x_x [psi]')
legend('Stringer 1','Stringer 2','Stringer 3','Stringer 4' ...
    ,'Stringer 5','Stringer 6')

% for i = 1:6
%     m1 = [1 -yST(i) -zST(i)];
%     for x=0:L
%         sigmaxx(x) = (E(i)*10^3)*m1*m2*m3;
%         plot(xvalue,sigmaxx(xvalue))
%         hold on
%     end
%     sigmaxx(i)
% end


% Plot MS Vs. Position (x) of all stringers (y-range to MS = +50~-50)
% recall sigma_star = Min(sigmaY/FSy, sigmaU/FSu)
% recall MS = (sigma_star/sigma_applied) - 1
figure(3)
subplot(2,1,1)
sigma_star1 = min(Sy(1)/FSyST(1), Su(1)/FSuST(1))*1e3; %allowable stress
MS_sigmaxx1 = abs(sigma_star1/sigmaxx1(x)) - 1;
MS_sx_plot1 = matlabFunction(MS_sigmaxx1);
plot(xvalue,MS_sx_plot1(xvalue))
hold on

sigma_star2 = min(Sy(2)/FSyST(2), Su(2)/FSuST(2))*1e3; %allowable stress
MS_sigmaxx2 = abs(sigma_star2/sigmaxx2(x)) - 1; %
MS_sx_plot2 = matlabFunction(MS_sigmaxx2);
plot(xvalue,MS_sx_plot2(xvalue))

sigma_star3 = min(Sy(3)/FSyST(3), Su(3)/FSuST(3))*1e3; %allowable stress
MS_sigmaxx3 = abs(sigma_star3/sigmaxx3(x)) - 1;
MS_sx_plot3 = matlabFunction(MS_sigmaxx3);
plot(xvalue,MS_sx_plot3(xvalue))

sigma_star4 = min(Sy(4)/FSyST(4), Su(4)/FSuST(4))*1e3; %allowable stress
MS_sigmaxx4 = abs(sigma_star4/sigmaxx4(x)) - 1;
MS_sx_plot4 = matlabFunction(MS_sigmaxx4);
plot(xvalue,MS_sx_plot4(xvalue))

sigma_star5 = min(Sy(5)/FSyST(5), Su(5)/FSuST(5))*1e3; %allowable stress
MS_sigmaxx5 = abs(sigma_star5/sigmaxx5(x)) - 1;
MS_sx_plot5 = matlabFunction(MS_sigmaxx5);
plot(xvalue,MS_sx_plot5(xvalue))

sigma_star6 = min(Sy(6)/FSyST(6), Su(6)/FSuST(6))*1e3; %allowable stress
MS_sigmaxx6 = abs(sigma_star6/sigmaxx6(x)) - 1; %
MS_sx_plot6 = matlabFunction(MS_sigmaxx6);
plot(xvalue,MS_sx_plot6(xvalue))
hold off

set(gca,'YLim',[-50 50])

hold off
grid on
xlabel('x [in]'); ylabel('Margin of Safety')
title('Stringer Margin of Safety Vs. Length')
legend('Stringer 1','Stringer 2','Stringer 3','Stringer 4' ...
    ,'Stringer 5','Stringer 6','location','southeast')

% % set y lim [-2 5]

subplot(2,1,2)

% plot(xvalue,MS_sx_plot1(xvalue))
% hold on
% plot(xvalue,MS_sx_plot2(xvalue))
% plot(xvalue,MS_sx_plot3(xvalue))
% plot(xvalue,MS_sx_plot4(xvalue))
% plot(xvalue,MS_sx_plot5(xvalue))
% plot(xvalue,MS_sx_plot6(xvalue))
sigma_star1 = min(Sy(1)/FSyST(1), Su(1)/FSuST(1))*1e3; %allowable stress
MS_sigmaxx1 = abs(sigma_star1/sigmaxx1(x)) - 1;
MS_sx_plot1 = matlabFunction(MS_sigmaxx1);
plot(xvalue,MS_sx_plot1(xvalue))
hold on

sigma_star2 = min(Sy(2)/FSyST(2), Su(2)/FSuST(2))*1e3; %allowable stress
MS_sigmaxx2 = abs(sigma_star2/sigmaxx2(x)) - 1; %
MS_sx_plot2 = matlabFunction(MS_sigmaxx2);
plot(xvalue,MS_sx_plot2(xvalue))

sigma_star3 = min(Sy(3)/FSyST(3), Su(3)/FSuST(3))*1e3; %allowable stress
MS_sigmaxx3 = abs(sigma_star3/sigmaxx3(x)) - 1;
MS_sx_plot3 = matlabFunction(MS_sigmaxx3);
plot(xvalue,MS_sx_plot3(xvalue))

sigma_star4 = min(Sy(4)/FSyST(4), Su(4)/FSuST(4))*1e3; %allowable stress
MS_sigmaxx4 = abs(sigma_star4/sigmaxx4(x)) - 1;
MS_sx_plot4 = matlabFunction(MS_sigmaxx4);
plot(xvalue,MS_sx_plot4(xvalue))

sigma_star5 = min(Sy(5)/FSyST(5), Su(5)/FSuST(5))*1e3; %allowable stress
MS_sigmaxx5 = abs(sigma_star5/sigmaxx5(x)) - 1;
MS_sx_plot5 = matlabFunction(MS_sigmaxx5);
plot(xvalue,MS_sx_plot5(xvalue))

sigma_star6 = min(Sy(6)/FSyST(6), Su(6)/FSuST(6))*1e3; %allowable stress
MS_sigmaxx6 = abs(sigma_star6/sigmaxx6(x)) - 1; %
MS_sx_plot6 = matlabFunction(MS_sigmaxx6);
plot(xvalue,MS_sx_plot6(xvalue))

hold off
set(gca,'YLim',[-2 5])
set(gca,'XLim',[0 150])
grid on
xlabel('x [in]'); ylabel('Margin of Safety')
legend('Stringer 1','Stringer 2','Stringer 3','Stringer 4' ...
    ,'Stringer 5','Stringer 6','location','southeast')


% for i = 1:6
% %     sigma_star = min(abs(Sy(i)/FSyST(i)),abs(Su(i)/FSuST(i))) * 1e3
% %     MS_ST(x) = sigma_star/sigmaxx - 1
% %     plot(xvalue,MS_ST(xvalue))
%     m1_2 = [1 -yST(i) -zST(i)];
% %     for i = 1:6
% %         MS_ST(x) = (min(abs(Sy(i)/FSyST(i)),abs(Su(i)/FSuST(i))) * 1e3)./sigmaxx2(x) - 1;
% %         plot(xvalue,MS_ST(xvalue))
% %         hold on
% %     end
%     sigma_star = min(abs(Sy(i)/FSyST(i)),abs(Su(i)/FSuST(i))) * 1e3 %[lb/in^2]
% %     sigma_star
%     sigmaxx2(x) = (E(i)*10^3)*m1_2*m2*m3;
%     sigmaxx2_inverse(x) = finverse(sigmaxx2,x);
%     MS_ST(x) = sigma_star*(sigmaxx2_inverse(x)) - 1;
%     MS_ST_plotvalue = MS_ST(xvalue);
%     plot(xvalue,MS_ST(xvalue))
%     hold on
% end

%                           START HERE NEXT TIME
% for Stringer 1
% m1_2 = [1 -yST(1) -zST(1)];
% sigmaxx2(x) = (E(1)*10^3)*m1_2*m2*m3; %design stress
% 
% sigma_star2 = min((Sy(1)/FSyST(1)),(Su(1)/FSuST(1))) * 1e3 %constant(allowable stress)
% MS_ST(x) = abs(sigma_star2./(sigmaxx2(x))) - 1; % M.S.
% 
% 
% plot(xvalue,MS_ST(xvalue))

% sigmaSTAR(i) = sigma_star

% try
% sigma_star(i) = min(abs(Sy(i)/FSyST(i)),abs(Su(i)/FSuST(i))) * 1e3 %[lb/in^2]
% for i = 1:6
%     m1 = [1 -yST(i) -zST(i)];
%     sigmaxx(x) = (E(i)*10^3)*m1*m2*m3;
%     for j = 1:6
%         MS_ST(x) = abs(sigma_star(i)/sigmaxx(x)) - 1
%         plot(xvalue,MS_ST(xvalue))
%         hold on
%     end
% end

%% Section 4: Stress and Maragin of Safety in Wing Skin Panels
Vy;
Vz;

% if set point A

% % As12 = L12*t(10);
% % As23 = L23*t(11);
% % As34 = L34*t(12);
% % As45 = L45*t(13);
% % As56 = L56*t(14);
% % As61 = L61*t(15);
% As12 = (pi*a*b)/4
% As23 = (yc-yi(2))*b*0.5
% As34 = (yi(4)-yi(3))*b*0.5+(yi(3)-yc)*b*0.5
% As45 = As34
% As56 = As23
% As61 = As12
As12 = (pi*a*b)/4;
As23 = (yi(3)-yi(2))*b*0.5;
As34 = (yi(4)-yi(3))*b*0.5+(yi(3)-yi(2))*b*0.5;
As45 = As34;
As56 = As23;
As61 = As12;

syms qo 
% q12 = qo - 1e3*A(2)./(EIzz/E(2)).*Vy*(yi(2)-yc) - 1e3*A(2)./(EIyy/E(2)).*Vz*(zi(2)-zc);
% q23 = q12 - 1e3*A(3)./(EIzz/E(3)).*Vy*(yi(3)-yc) - 1e3*A(3)./(EIyy/E(3)).*Vz*(zi(3)-zc);
% q34 = q23 - 1e3*A(4)./(EIzz/E(4)).*Vy*(yi(4)-yc) - 1e3*A(4)./(EIyy/E(4)).*Vz*(zi(4)-zc);
% q45 = q34 - 1e3*A(5)./(EIzz/E(5)).*Vy*(yi(5)-yc) - 1e3*A(5)./(EIyy/E(5)).*Vz*(zi(5)-zc);
% q56 = q45 - 1e3*A(6)./(EIzz/E(6)).*Vy*(yi(6)-yc) - 1e3*A(6)./(EIyy/E(6)).*Vz*(zi(6)-zc);
% q61 = qo;

q12 = qo - 1e3*A(1)./(EIzz/E(1)).*Vy*(yi(1)-yc) - 1e3*A(1)./(EIyy/E(1)).*Vz*(zi(1)-zc);
q23 = q12 - 1e3*A(2)./(EIzz/E(2)).*Vy*(yi(2)-yc) - 1e3*A(2)./(EIyy/E(2)).*Vz*(zi(2)-zc);
q34 = q23 - 1e3*A(3)./(EIzz/E(3)).*Vy*(yi(3)-yc) - 1e3*A(3)./(EIyy/E(3)).*Vz*(zi(3)-zc);
q45 = q34 - 1e3*A(4)./(EIzz/E(4)).*Vy*(yi(4)-yc) - 1e3*A(4)./(EIyy/E(4)).*Vz*(zi(4)-zc);
q56 = q45 - 1e3*A(5)./(EIzz/E(5)).*Vy*(yi(5)-yc) - 1e3*A(5)./(EIyy/E(5)).*Vz*(zi(5)-zc);
q61 = qo;

% q12 = qo - 1e3*A(2)./(EIzz/E(2)).*Vy*(yc-yi(2)) - 1e3*A(2)./(EIyy/E(2)).*Vz*(zc-zi(2));
% q23 = q12 - 1e3*A(3)./(EIzz/E(3)).*Vy*(yi(3)-yc) - 1e3*A(3)./(EIyy/E(3)).*Vz*(zc-zi(3));
% q34 = q23 - 1e3*A(4)./(EIzz/E(4)).*Vy*(yi(4)-yc) - 1e3*A(4)./(EIyy/E(4)).*Vz*(zi(4)-zc);
% q45 = q34 - 1e3*A(5)./(EIzz/E(5)).*Vy*(yi(5)-yc) - 1e3*A(5)./(EIyy/E(5)).*Vz*(zi(5)-zc);
% q56 = q45 - 1e3*A(6)./(EIzz/E(6)).*Vy*(yc-yi(6)) - 1e3*A(6)./(EIyy/E(6)).*Vz*(zi(6)-zc);
% q61 = qo;
%External Moment about point A (center of ellipse)
M_ext = Vz*(yc-yi(2)) + Mx;
% Internal Moment aboyt point A
% % A12 = (pi*a*b)/4;
% % A23 = b*(yi(3)-yi(2))*0.5;
% % A34 = b*(yi(4)-yi(3))*0.5 + A23;
% % A45 = A34;
% % A56 = A23;
% % A61 = A12;
M_int = 2*As12*q12 + 2*As23*q23 + 2*As34*q34 + 2*As45*q45 + ...
    2*As56*q56 + 2*As61*q61;
% % M_int(x) = 2*A12*q12 + 2*A23*q23 + 2*A34*q34 + 2*A45*q45 + ...
% %     2*A56*q56 + 2*A61*q61;

% Let M_int = M_ext to solve q0, note q0 = qo
q0 = solve(M_ext == M_int, qo);

% q12 = q0 - 1e3*A(2)./(EIzz/E(2)).*Vy*(yi(2)-yc) - 1e3*A(2)./(EIyy/E(2)).*Vz*(zi(2)-zc);
% q23 = q12 - 1e3*A(3)./(EIzz/E(3)).*Vy*(yi(3)-yc) - 1e3*A(3)./(EIyy/E(3)).*Vz*(zi(3)-zc);
% q34 = q23 - 1e3*A(4)./(EIzz/E(4)).*Vy*(yi(4)-yc) - 1e3*A(4)./(EIyy/E(4)).*Vz*(zi(4)-zc);
% q45 = q34 - 1e3*A(5)./(EIzz/E(5)).*Vy*(yi(5)-yc) - 1e3*A(5)./(EIyy/E(5)).*Vz*(zi(5)-zc);
% q56 = q45 - 1e3*A(6)./(EIzz/E(6)).*Vy*(yi(6)-yc) - 1e3*A(6)./(EIyy/E(6)).*Vz*(zi(6)-zc);
% q61 = q0;

q121 = q0 - 1e3*A(1)./(EIzz/E(1)).*Vy*(yi(1)-yc) - 1e3*A(1)./(EIyy/E(1)).*Vz*(zi(1)-zc);
q231 = q121 - 1e3*A(2)./(EIzz/E(2)).*Vy*(yi(2)-yc) - 1e3*A(2)./(EIyy/E(2)).*Vz*(zi(2)-zc);
q341 = q231 - 1e3*A(3)./(EIzz/E(3)).*Vy*(yi(3)-yc) - 1e3*A(3)./(EIyy/E(3)).*Vz*(zi(3)-zc);
q451 = q341 - 1e3*A(4)./(EIzz/E(4)).*Vy*(yi(4)-yc) - 1e3*A(4)./(EIyy/E(4)).*Vz*(zi(4)-zc);
q561 = q451 - 1e3*A(5)./(EIzz/E(5)).*Vy*(yi(5)-yc) - 1e3*A(5)./(EIyy/E(5)).*Vz*(zi(5)-zc);
q611 = q0;

% q12 = q0 - 1e3*A(2)./(EIzz/E(2)).*Vy*(yc-yi(2)) - 1e3*A(2)./(EIyy/E(2)).*Vz*(zc-zi(2));
% q23 = q12 - 1e3*A(3)./(EIzz/E(3)).*Vy*(yi(3)-yc) - 1e3*A(3)./(EIyy/E(3)).*Vz*(zc-zi(3));
% q34 = q23 - 1e3*A(4)./(EIzz/E(4)).*Vy*(yi(4)-yc) - 1e3*A(4)./(EIyy/E(4)).*Vz*(zi(4)-zc);
% q45 = q34 - 1e3*A(5)./(EIzz/E(5)).*Vy*(yi(5)-yc) - 1e3*A(5)./(EIyy/E(5)).*Vz*(zi(5)-zc);
% q56 = q45 - 1e3*A(6)./(EIzz/E(6)).*Vy*(yc-yi(6)) - 1e3*A(6)./(EIyy/E(6)).*Vz*(zi(6)-zc);
% q61 = q0;

%           ``````````if set point at center
% syms qo
% As12 = (pi*a*b)/4+0.5*b*(yc-yi(2));
% As23 = (yc-yi(2))*b*0.5+(yi(3)-yc)*b*0.5;
% As34 = (yi(4)-yi(3))*b*0.5+(yi(3)-yc)*b*0.5;
% As45 = As34;
% As56 = As23;
% As61 = As12
% 
% q12(x) = qo - 1e3*A(2)/EIzz*E(2)*Vy*(yi(2)-yc) - 1e3*A(2)/EIyy*E(2)*Vz*abs(zi(2));
% q23(x) = qo - 1e3*A(3)/EIzz*E(3)*Vy*(yi(3)-yc) - 1e3*A(3)/EIyy*E(3)*Vz*abs(zi(3));
% q34(x) = qo - 1e3*A(4)/EIzz*E(4)*Vy*(yi(4)-yc) - 1e3*A(4)/EIyy*E(4)*Vz*abs(zi(4));
% q45(x) = qo - 1e3*A(5)/EIzz*E(5)*Vy*(yi(5)-yc) - 1e3*A(5)/EIyy*E(5)*Vz*(zi(5));
% q56(x) = qo - 1e3*A(6)/EIzz*E(6)*Vy*(yi(6)-yc) - 1e3*A(6)/EIyy*E(6)*Vz*(zi(6));
% q61(x) = qo;
% 
%     % q12(x) = qo - Vy/EIzz*y_c(1) - Vy/EIyy*z_c(1)
%     % q23(x) = qo - Vy/EIzz*y_c(2) - Vy/EIyy*z_c(2)
%     % q34(x) = qo - Vy/EIzz*y_c(3) - Vy/EIyy*z_c(3)
%     % q45(x) = qo - Vy/EIzz*y_c(4) - Vy/EIyy*z_c(4)
%     % q56(x) = qo - Vy/EIzz*y_c(5) - Vy/EIyy*z_c(5)
%     % q61(x) = qo
% 
% % External Moment about centroid
% M_ext(x) =  Mx;
% 
% % Internal Moment about centroid
% 
% 
% M_int(x) = 2*As12*q12 + 2*As23*q23 + 2*As34*q34 + 2*As23*q45 + ...
%     2*As56*q56 + 2*As61*q61;
% 
% % Let M_int = M_ext
% q0(x) = solve(M_ext(x) == M_int(x), qo)
% % 
%     % q12(x) = q0 - Vy/EIzz*y_c(1) - Vy/EIyy*z_c(1)
%     % q23(x) = q0 - Vy/EIzz*y_c(2) - Vy/EIyy*z_c(2)
%     % q34(x) = q0 - Vy/EIzz*y_c(3) - Vy/EIyy*z_c(3)
%     % q45(x) = q0 - Vy/EIzz*y_c(4) - Vy/EIyy*z_c(4)
%     % q56(x) = q0 - Vy/EIzz*y_c(5) - Vy/EIyy*z_c(5)
%     % q61(x) = q0
% 
% q12(x) = q0 - 1e3*A(2)/EIzz*E(2)*Vy*abs(yi(2)-yc) - 1e3*A(2)/EIyy*E(2)*Vz*abs(zi(2));
% q23(x) = q0 - 1e3*A(3)/EIzz*E(3)*Vy*abs(yi(3)-yc) - 1e3*A(3)/EIyy*E(3)*Vz*abs(zi(3));
% q34(x) = q0 - 1e3*A(4)/EIzz*E(4)*Vy*abs(yi(4)-yc) - 1e3*A(4)/EIyy*E(4)*Vz*abs(zi(4));
% q45(x) = q0 - 1e3*A(5)/EIzz*E(5)*Vy*abs(yi(5)-yc) - 1e3*A(5)/EIyy*E(5)*Vz*(zi(5));
% q56(x) = q0 - 1e3*A(6)/EIzz*E(6)*Vy*abs(yi(6)-yc) - 1e3*A(6)/EIyy*E(6)*Vz*(zi(6));
% q61(x) = q0;


% shear stress tau
tau12(x) = q121/t(10);
tau23(x) = q231/t(11);
tau34(x) = q341/t(12);
tau45(x) = q451/t(13);
tau56(x) = q561/t(14);
tau61(x) = q611/t(15);

figure(4)
plot(xvalue,tau12(xvalue))
hold on
plot(xvalue,tau23(xvalue))

plot(xvalue,tau34(xvalue))

plot(xvalue,tau45(xvalue))

plot(xvalue,tau56(xvalue))

plot(xvalue,tau61(xvalue))
hold off
grid on
title('Skin Shear Stress Vs. Length')
xlabel('x [in]'); ylabel('Skin \tau_x_s [psi]')
legend('Skin 1-2','Skin 2-3','Skin 3-4','Skin 4-5','Skin 5-6','Skin 6-1')

% Plot MS
% recall tau_star = Min(tau_Y/FSySK, tau_U/FSuSK)
% recall MS = (tau_star/tau_applied) - 1
tau_star1 = min(tau_y(10)/FSySK(10), tau_u(10)/FSuSK(10))*1e3; %allowable stress
tau_star2 = min(tau_y(11)/FSySK(11), tau_u(11)/FSuSK(11))*1e3; %allowable stress
tau_star3 = min(tau_y(12)/FSySK(12), tau_u(12)/FSuSK(12))*1e3; %allowable stress
tau_star4 = min(tau_y(13)/FSySK(13), tau_u(13)/FSuSK(13))*1e3; %allowable stress
tau_star5 = min(tau_y(14)/FSySK(14), tau_u(14)/FSuSK(14))*1e3; %allowable stress
tau_star6 = min(tau_y(15)/FSySK(15), tau_u(15)/FSuSK(15))*1e3; %allowable stress

MS_tau12(x) = abs(tau_star1/tau12(x)) - 1;
MS_tau23(x) = abs(tau_star2/tau23(x)) - 1;
MS_tau34(x) = abs(tau_star3/tau34(x)) - 1;
MS_tau45(x) = abs(tau_star4/tau45(x)) - 1;
MS_tau56(x) = abs(tau_star5/tau56(x)) - 1;
MS_tau61(x) = abs(tau_star6/tau61(x)) - 1;
figure(5)
subplot(2,1,1)
plot(xvalue,MS_tau12(xvalue))
hold on
plot(xvalue,MS_tau23(xvalue))
plot(xvalue,MS_tau34(xvalue))
plot(xvalue,MS_tau45(xvalue))
plot(xvalue,MS_tau56(xvalue))
plot(xvalue,MS_tau61(xvalue))
hold off
grid on
legend('Skin 1-2','Skin 2-3','Skin 3-4','Skin 4-5','Skin 5-6','Skin 6-1','location','southeast')
title('Skin Shear Margin of Safety Vs. Length')
xlabel('x [in]'); ylabel('Margin of Safety')
set(gca,'YLim',[-50 50])

subplot(2,1,2)
plot(xvalue,MS_tau12(xvalue))
hold on
plot(xvalue,MS_tau23(xvalue))
plot(xvalue,MS_tau34(xvalue))
plot(xvalue,MS_tau45(xvalue))
plot(xvalue,MS_tau56(xvalue))
plot(xvalue,MS_tau61(xvalue))
hold off
grid on
xlabel('x [in]'); ylabel('Margin of Safety')
legend('Skin 1-2','Skin 2-3','Skin 3-4','Skin 4-5','Skin 5-6','Skin 6-1','location','southeast')
set(gca,'YLim',[-2 5])


%% Section 5: Wing Displacement and Twist
% Plot of wing's in-plane displacement (v) in inch
dbo_Int_V(x) = Mz(x); % double integral of V(x)
slg_Int_V(x) = int(dbo_Int_V(x)); % integral of V(x)
V_x(x) = int(slg_Int_V(x)) / (EIzz*1e3); % EIzz was calculated in kip-in^2
figure(6)
subplot(3,1,1)
plot(xvalue,V_x(xvalue))
title('Wing Displacement Distribution and Teist vs. Length')
xlabel('x [in]'); ylabel('V_x [in]');
grid on

% Plot of wing's transverse displacement (w) in inch
dbo_Int_W(x) = My(x); % double integral of W(x); here W(x) is not weight
slg_Int_W(x) = int(dbo_Int_W(x)); % integral of W(x)
W_x(x) = int(slg_Int_W(x)) / (-EIyy*1e3); % EIyy was calculated in kip-in^2
subplot(3,1,2)
plot(xvalue,W_x(xvalue))
xlabel('x [in]'); ylabel('W_x [in]');
grid on

% Plot wing twist (degree) vs. wing position (x)
A12 = (pi*a*b)/4;
A23 = (yi(3)-yi(2))*b*0.5;
A34 = (yi(4)-yi(3))*b*0.5+(yi(3)-yi(2))*b*0.5;
A45 = A34;
A56 = A23;
A61 = A12;

L12 = 1/4*(pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b))));
L23 = yi(3)-yi(2);
L34 = sqrt((yi(4)-yi(3))^2+b^2);
L45 = L34;
L56 = L23;
L61 = L12;
sum_A = A12+A23+A34+A45+A56+A61;
sum_qiSi_overGt(x) = L12*(q121)/(G(10)*1e6*t(10)) + L23*(q231)/(G(11)*1e6*t(11)) ...
    + L34*(q341)/(G(12)*1e6*t(12)) + L45*(q451)/(G(13)*1e6*t(13)) ... 
    + L56*(q561)/(G(14)*1e6*t(14)) + L61*(q611)/(G(15)*1e6*t(15));
Der_delta(x) = 1/(2*sum_A).*(sum_qiSi_overGt(x));
delta(x) = int(Der_delta(x)) * 180/pi;
subplot(3,1,3)
plot(xvalue,delta(xvalue))
xlabel('x [in]'); ylabel('theta_x [degree]');
grid on

%% Output 6
Result = [yc-(yi(1)) zc EA*1e6 EIyy*1e3 EIzz*1e3 EIyz*1e3 GJ*1e3]

% Use this command
%Project3_SE160A('SE160A_Project3_InputFile.xlsx')