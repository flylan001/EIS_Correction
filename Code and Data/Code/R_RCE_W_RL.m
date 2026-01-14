function [x,Zs1,Zs2] = R_RCE_W_RL(eisdata)
%%%%% 输入待辨识参数
%   此处显示详细说明
freq = eisdata(:,1);
Zre = eisdata(:,2);
Zim = eisdata(:,3);
f=@(x)Battery_Para_identification(x,freq,Zre,Zim);

options = optimoptions('fmincon','ConstraintTolerance',1e-25,'Algorithm','active-set');
x0=[8.53e-08 0.000419 0.000147 	39 0.93 0.0000657 0.8 0.000365 8.18e-9];
% x0(6) = - sin(pi/4)/( (2*pi*freq(end))^0.5 * Zim(end) );
[lb,ub] = boundry(x0,10);
x = fmincon(f,x0,[],[],[],[],lb,ub,[],options);
[Zs1,Zs2,MAE]=Battery_Para_identification_result(x,freq,Zre,Zim);
end

function [lb,ub] = boundry(x0,a)

lb(1)=x0(1)/a;lb(2)=x0(2)/a;lb(3)=x0(3)/a;lb(4)=x0(4)/a;lb(5)=x0(5)/a;lb(6)=x0(6)/a;lb(7)=0;lb(8)=x0(8)/a;lb(9)=x0(9)/a;
ub(1)=x0(1)*a;ub(2)=x0(2)*a;ub(3)=x0(3)*a;ub(4)=x0(4)*a;ub(5)=1;ub(6)=x0(6)*a;ub(7)=1;ub(8)=x0(8)*a;ub(9)=x0(9)*a;
end

function [Zs1,Zs2,MAE]=Battery_Para_identification_result(x,freq,Zre,Zim)
L = x(1);
R0 = x(2);
Rct = x(3);
CPE1 = x(4);
CPEn = x(5);
W = x(6);
Wn = x(7);
R1 = x(8);
L1 = x(9);
w = 2*pi.*freq;
Zs1 = 1j*w*L+ R0 + Rct./( 1+(1j*w).^CPEn.*Rct.*CPE1 ) + 1./( (1j*w).^Wn).*W + (1j.*w.*R1.*L1)./( R1 + 1j*w*L1 );
Zs2 =  R0 + Rct./( 1+(1j*w).^CPEn.*Rct.*CPE1 ) + 1./( (1j*w).^Wn).*W;
Z_error = abs(Zre-real(Zs1)) + abs(Zim-imag(Zs1));
MAE=mean(Z_error);
end

function f=Battery_Para_identification(x,freq,Zre,Zim)
L = x(1);
R0 = x(2);
Rct = x(3);
CPE1 = x(4);
CPEn = x(5);
W = x(6);
Wn = x(7);
R1 = x(8);
L1 = x(9);
w = 2*pi.*freq;
Zf = 1j*w*L+ R0 + Rct./( 1+(1j*w).^CPEn.*Rct.*CPE1 ) + 1./( (1j*w).^Wn).*W + (1j.*w.*R1.*L1)./( R1 + 1j*w*L1 );
Z_error = abs(Zre-real(Zf)) + abs(Zim-imag(Zf));
f=sum(Z_error)*1000000;
end