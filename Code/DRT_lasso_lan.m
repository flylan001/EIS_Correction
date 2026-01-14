function [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan(impdata,num,lambda0)
% This part of the code follows the approach proposed in Electrochimica Acta 184 (2015) 483–499
freq=impdata(:,1)';
Zre=impdata(:,2);
Zim=impdata(:,3);  %% 为正值
if Zim(end)>0
    Zim=-Zim;
end
% % 阻抗数据截断
% minR=find(Zre==min(Zre));
% freq = freq(minR(1):end);
% Zre = Zre(minR(1):end);
% Zim = Zim(minR(1):end);

% 使用Gaussian函数
coeff=1;
lambda=lambda0;
rbf_gaussian_4_FWHM = @(x) exp(-(x).^2);  %设置高斯函数
FWHM_coeff = 2*fzero(@(x) rbf_gaussian_4_FWHM(x)-1/2, 1);  %求解高斯函数位于半高处的宽度，即半高宽，1为初始值
delta = mean(diff(log(1./freq)));
epsilon  = coeff*FWHM_coeff/delta;

[Are,freq_M]=cal_A(freq,epsilon,0);

Aim=cal_A(freq,epsilon,1);
Are(:,2) = 1;
Aim(:,1) = 2*pi*freq(:);
%计算M矩阵
M_temp=eye(length(freq),length(freq));
M=zeros(length(freq)+2,length(freq)+2);
M(3:end, 3:end)=M_temp;

%% 权重设置
z = abs(Zre+1j*Zim);
weight0 = 1./z;
weight = weight0 ./(mean(weight0));
M_temp=eye(length(freq),length(freq));
for i =1:1:length(freq)
    M_temp(i,i)=M_temp(i,i)*weight(i);
end

[H_combined,f_combined] = quad_format_combined(Are, Aim, Zre, Zim, M, lambda,M_temp);
options = optimset('algorithm','interior-point-convex','Display','off','TolFun',1e-15,'TolX',1e-10,'MaxFunEvals', 1E5);
lb = 0*ones(numel(freq_M)+2,1);
ub = Inf*ones(numel(freq_M)+2,1);
x_0 = -ones(size(lb));

x_ridge = quadprog(H_combined, f_combined, [], [], [], [], lb, ub, x_0, options);

mu_Z_re = Are*x_ridge;
mu_Z_im = Aim*x_ridge;
x_ridge_modify=x_ridge;
x_ridge_modify(1)=0;
mu_Z_re2 = Are*x_ridge_modify;
mu_Z_im2 = Aim*x_ridge_modify;
solve_Z=[mu_Z_re,mu_Z_im,mu_Z_re2,mu_Z_im2];
%下舍上入
taumax = ceil(max(log10(1./freq_M)))+0.5;
taumin = floor(min(log10(1./freq_M)))-0.5;
freq_fine = logspace(-taumin, -taumax, 20*numel(freq_M));
[gamma_ridge_fine,freq_fine] = map_array_to_gamma(freq_fine, freq_M, x_ridge(3:end), epsilon);

gamma=[1./freq_fine', gamma_ridge_fine];
%%
% 绘图部分
if num==0

else
    if num==1
        figure,
        set(gcf, 'Position',[100 100 1200 900/2]);
    end
    subplot(1,2,2);
    plot(1./freq_fine, gamma_ridge_fine, '-k', 'LineWidth', 3);hold on
    y_min = 0;
    y_max = max(gamma_ridge_fine);
    set(gca,'xscale','log','xlim',[10^taumin, 10^taumax],'ylim',[y_min, 1.1*y_max],'Fontsize',20,'xtick',10.^[-10:2:10],'TickLabelInterpreter','latex')
    % hold off
    xlabel('$\tau/s$', 'Interpreter', 'Latex','Fontsize',24)
    ylabel('$\gamma(\tau)/\Omega$','Interpreter', 'Latex','Fontsize',24);
    subplot(1,2,1);
    plot(mu_Z_re,-mu_Z_im,'ok', 'LineWidth', 3);hold on
    plot(Zre,-Zim,'*r', 'LineWidth', 3);
    axis equal
    set(gca,'Fontsize',20,'TickLabelInterpreter','latex');
    ylabel("-Imag/\Omega");xlabel("Real/\Omega");
    legend('计算结果','实验数据','Location','north');
    % plot(Zre,Zim);hold on
    % plot(mu_Z_re,mu_Z_im)
end

end

%line_num 行数
function  [impdata,iniU,F,linenum]=dataread_imp(path)
Testdata=importdata(path,',',20);
[linenum,~]=size(Testdata.data);
impdata(1:linenum,1:4)=Testdata.data(1:linenum,2:5);
F(1:linenum)=Testdata.data(1:linenum,1);
for ii=6:1:9
    if contains(Testdata.textdata(ii),'Init')
        break;
    end
end
text_inu=Testdata.textdata(ii);
text_iniU=text_inu{1};
iniU=str2num(text_iniU(14:length(text_iniU)));
end

function [A,freq_M]=cal_A(freq,epsilon,label)
N=length(freq);
freq_log=log10(freq);
num_mag=round(1/((freq_log(1)-freq_log(end))/(N-1)));
add_num=0;
freq_add=freq(1:add_num*num_mag)*10;
freq_M=[freq_add freq];
freq_add=freq((N-add_num*num_mag+1):end)/10;
freq_M=[freq_M freq_add];
A_temp=zeros(N,N+add_num*num_mag*2);  %构建矩阵A
A= zeros(N, N+add_num*num_mag*2+2);
% compute using brute force
for iter_freq_n = 1: N
    for iter_freq_m = 1: N+add_num*num_mag*2
        freq_n = freq(iter_freq_n);
        freq_m = freq_M(iter_freq_m);
        % compute all RBF terms
        A_temp(iter_freq_n, iter_freq_m) = g_i(freq_n, freq_m, epsilon,label);
    end
A(:, 3:end) = A_temp;
end
end


function out_val = g_i(freq_n, freq_m, epsilon,label)
%这个函数会生成A_re的元素
alpha = 2*pi*freq_n/freq_m;
rbf = @(x) exp(-(epsilon*x).^2);
if label==0
    integrand_g_i = @(x) 1./(1+alpha^2*exp(2*x)).*rbf(x);
elseif label==1
    integrand_g_i =@(x) -alpha./(1./exp(x)+alpha^2*exp(x)).*rbf(x);
end
out_val = integral(integrand_g_i, -inf, inf,'RelTol',1E-9,'AbsTol',1E-9);
end

function [H,c] = quad_format_combined(A_re,A_im,b_re,b_im,M,lambda,M_temp) 
% 这个函数将DRT回归重新格式化为一个二次规划问题-它同时使用re和im

    H = 2*((A_re'*M_temp*A_re+A_im'*M_temp*A_im));
    H = (H'+H)/2;
    c = -2*(b_im'*M_temp*A_im+b_re'*M_temp*A_re)+lambda*ones(1,size(H,1));

end

function [gamma_ridge_fine,freq_fine] = map_array_to_gamma(freq_fine, freq, x_ridge, epsilon)
rbf = @(y, y0) exp(-(epsilon*(y-y0)).^2);
y0 = -log(freq');
gamma_ridge_fine = zeros(size(freq_fine))';
for iter_freq_map = 1: numel(freq_fine)
    freq_loc = freq_fine(iter_freq_map);
    y = -log(freq_loc);
    gamma_ridge_fine(iter_freq_map) = x_ridge'*rbf(y, y0);
end
end
