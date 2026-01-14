%% ========================================================================
%  EIS Experimental Data Visualization and Correction
%  Author: LAN
%  Date: 2025-11-12
%  ------------------------------------------------------------------------
%  This script:
%   (1)  Plots raw EIS experimental data (µΩ)
%   (2) Scans parameters (a, b) for DRT-based reconstruction
%   (3) Evaluates MSE and visualizes parameter sensitivity
%  ========================================================================

clc; clear; close all;

load exp1_eis.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%% 10kHz-0.1Hz
%% 参数ab遍历 
start_freq_index=1;
n=0.0:0.1:2;
for j=1:1:length(n)
    for j1=1:1:length(n)
        n0 = [n(j),n(j1)]
        j00 = (j-1) * length(n) + j1;
        for i0=1:1:size(exp1_eis,1)
            eisdata=exp1_eis{i0,1};
            [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0);
            EIS_all0{i0,1} = solve_Z;
            EIS_all0{i0,2} = eisdata;
        end
        EIS_all{j00}=EIS_all0;
    end
end
%% 误差分布图
n0_all=[];error_data=[];
for j=1:1:length(n)
    for j1=1:1:length(n)
        n0_all = [n0_all;n(j),n(j1)];
        j00 = (j-1) * length(n) + j1;
        EIS_all0 = EIS_all{j00};
        for i0=1:1:size(exp1_eis,1)
            solve_Z = EIS_all0{i0,1};
            eisdata = EIS_all0{i0,2};
            error_data(:,i0) = [solve_Z(:,3)-solve_Z(1,3);solve_Z(:,4)];
        end
        error = std(error_data,0,2);
        error_sum(j00,1) = mean(error)*1000000;
    end
end

figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all=[];
x=[min(n0_all(:,1)),max(n0_all(:,1))];
y=[min(n0_all(:,2)),max(n0_all(:,2))];
imagesc(x,y,1./reshape(error_sum,[length(n),length(n)]));
colorbar('Ticks',[0.1,0.2,0.3,0.4,0.5]);
title("      Reciprocal of Average MSE(\rm\mu\Omega)")
% 反转Y轴坐标
set(gca, 'YDir', 'normal');
xlabel("\ita");ylabel('\itb');xtickformat('%.1f'); ytickformat('%.1f'); 
set(gca,'LineWidth',3,'position',[0.16,0.18,0.7,0.7],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
set(gca,'xtick',0.0:0.5:1000,'ytick',0.0:0.5:1000);
gca=boxoff(gca);clear gca

%% 精细划分
n1=1:0.01:1.2;
n2=0.3:0.02:0.7;
for j=1:1:length(n1)
    for j1=1:1:length(n2)
        n0 = [n1(j),n2(j1)]
        j00 = (j-1) * length(n2) + j1;
        for i0=1:1:size(exp1_eis,1)
            eisdata=exp1_eis{i0,1};
            [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0);
            EIS_all0{i0,1} = solve_Z;
            EIS_all0{i0,2} = eisdata;
        end
        EIS_all{j00}=EIS_all0;
    end
    j;
end
%% 误差分布图
n1=1:0.01:1.2;
n2=0.3:0.02:0.7;
eis_label=[2 4 6 8 10 12];
n0_all=[];error_data=[];
for j=1:1:length(n1)
    for j1=1:1:length(n2)
        n0_all = [n0_all;n1(j),n2(j1)];
        j00 = (j-1) * length(n2) + j1
        EIS_all0 = EIS_all{j00};
        % figure,
        for i0=1:1:size(exp1_eis,1)
            solve_Z = EIS_all0{i0,1};
            eisdata = EIS_all0{i0,2};
            error_data(:,i0) = [solve_Z(:,3)-solve_Z(1,3);solve_Z(:,4)];
            % plot(solve_Z(:,3)-solve_Z(1,3),-solve_Z(:,4),'-*','Color',color(i0,:));hold on
        end
        error = std(error_data,0,2);
        error_sum(j00,1) = mean(error)*1000000;
    end
end

figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all=[];
x=[min(n0_all(:,1)),max(n0_all(:,1))];
y=[min(n0_all(:,2)),max(n0_all(:,2))];
imagesc(x,y,1./reshape(error_sum,[length(n1),length(n2)]));
colorbar('Ticks',[0.3,0.4,0.5,0.6]);
title("      Reciprocal of Average MSE(\rm\mu\Omega)")
% 反转Y轴坐标
set(gca, 'YDir', 'normal');
xlabel("\ita");ylabel('\itb');xtickformat('%.1f'); ytickformat('%.1f'); 
set(gca,'LineWidth',3,'position',[0.16,0.18,0.7,0.7],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
set(gca,'xtick',0.1:0.1:1000,'ytick',0.1:0.1:1000);
gca=boxoff(gca);clear gca

%% 考虑耦合项 整体数据拟合
start_freq_index=1;
labelnn = find(error_sum == min(error_sum))
n0 = n0_all(labelnn,:);

figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all =[];
EIS_gather=[];EIS_gather_origin=[]; 
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0 );
    plo=plot(eisdata(:,2)*1000000,eisdata(:,3)*1000000,'o','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    plot_all=[plot_all plo];
    plo=plot(solve_Z(:,1)*1000000,-solve_Z(:,2)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,1)*1000000;-solve_Z(:,2)*1000000]'];
    EIS_gather_origin=[EIS_gather_origin;[eisdata(:,2)*1000000;eisdata(:,3)*1000000]'];
end
label={'+','x','s',"diamond","pentagram"};
eis_label0=[21 31 41 51];
for k=1:1:length(eis_label0)
    real0=EIS_gather(:,eis_label0(k));
    imag0=EIS_gather(:,51+eis_label0(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    real0=EIS_gather_origin(:,eis_label0(k));
    imag0=EIS_gather_origin(:,51+eis_label0(k));
    plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end
% plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
% plo_all=[plo_all plo];
% plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
% plo_all=[plo_all plo];
leg_all={'EIS1','DRT_RL1','EIS2','DRT_RL2','EIS3','DRT_RL3','EIS4','DRT_RL4','EIS5','DRT_RL5','EIS6','DRT_RL6','100 Hz','10 Hz','1 Hz','0.1 Hz'};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis equal
axis([0.369 0.669 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all,'Interpreter','none');
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',3,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%% 考虑耦合项 不考虑L部分
start_freq_index=1;
figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all =[];
EIS_gather=[];EIS_gather_origin=[]; 
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0 );
    plo=plot(solve_Z(:,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000;-solve_Z(:,4)*1000000]'];
    EIS_gather_origin=[EIS_gather_origin;[eisdata(:,3)*1000000;eisdata(:,4)*1000000]'];
end
label={'+','x','s',"diamond","pentagram"};
eis_label0=[21 31 41 51];
for k=1:1:length(eis_label0)
    real0=EIS_gather(:,eis_label0(k));
    imag0=EIS_gather(:,51+eis_label0(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'DRT_RL1-RL','DRT_RL2-RL','DRT_RL3-RL','DRT_RL4-RL','DRT_RL5-RL','DRT_RL6-RL','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis equal
axis([0.34 0.669 -0.079 0.24]*1000)
leg=legend(plot_all,leg_all,'Interpreter','none');
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%% 考虑耦合项 不考虑L R部分
start_freq_index=1;
figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all =[];
EIS_gather=[];EIS_gather_origin=[]; 
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0 );
    plo=plot(solve_Z(:,3)*1000000 - solve_Z(1,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000 - solve_Z(1,3)*1000000;-solve_Z(:,4)*1000000]'];
    EIS_gather_origin=[EIS_gather_origin;[eisdata(:,2)*1000000;eisdata(:,3)*1000000]'];
end
label={'+','x','s',"diamond","pentagram"};
eis_label0=[21 31 41 51];
for k=1:1:length(eis_label0)
    real0=EIS_gather(:,eis_label0(k));
    imag0=EIS_gather(:,51+eis_label0(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    % real0=EIS_gather_origin(:,eis_label0(k));
    % imag0=EIS_gather_origin(:,51+eis_label0(k));
    % plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'DRT_RL1-RL-R0','DRT_RL2-RL-R0','DRT_RL3-RL-R0','DRT_RL4-RL-R0','DRT_RL5-RL-R0','DRT_RL6-RL-R0','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis equal
axis([-0.01 0.24 -0.02 0.23]*1000)
leg=legend(plot_all,leg_all,'Interpreter','none');
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca


%% 1k-0.1Hz
%% 参数ab遍历
start_freq_index=11;
figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');
EIS_gather=[];
n=0.0:0.1:2;
for j=1:1:length(n)
    for j1=1:1:length(n)
        n0 = [n(j),n(j1)]
        j00 = (j-1) * length(n) + j1;
        for i0=1:1:size(exp1_eis,1)
            eisdata=exp1_eis{i0,1};
            [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0);
            EIS_all0{i0,1} = solve_Z;
            EIS_all0{i0,2} = eisdata;
        end
        EIS_all{j00}=EIS_all0;
    end
    j;
end
%% 误差分布图
n=00:0.1:2;
n0_all=[];error_data=[];
for j=1:1:length(n)
    for j1=1:1:length(n)
        n0_all = [n0_all;n(j),n(j1)];
        j00 = (j-1) * length(n) + j1;
        EIS_all0 = EIS_all{j00};
        for i0=1:1:size(exp1_eis,1)
            solve_Z = EIS_all0{i0,1};
            eisdata = EIS_all0{i0,2};
            error_data(:,i0) = [solve_Z(:,3)-solve_Z(1,3);solve_Z(:,4)];
        end
        error = std(error_data,0,2);
        error_sum(j00,1) = mean(error)*1000000;
    end
end

figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all=[];
x=[min(n0_all(:,1)),max(n0_all(:,1))];
y=[min(n0_all(:,2)),max(n0_all(:,2))];
imagesc(x,y,1./reshape(error_sum,[length(n),length(n)]));
colorbar('Ticks',[0.1,0.2,0.3,0.4,0.5,0.6]);
title("      Reciprocal of Average MSE(\rm\mu\Omega)")
% 反转Y轴坐标
set(gca, 'YDir', 'normal');
xlabel("\ita");ylabel('\itb');xtickformat('%.1f'); ytickformat('%.1f'); 
set(gca,'LineWidth',3,'position',[0.16,0.18,0.7,0.7],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
set(gca,'xtick',0.0:0.5:1000,'ytick',0.0:0.5:1000);
gca=boxoff(gca);clear gca


%% 考虑耦合项 整体数据拟合
start_freq_index=11;
labelnn = find(error_sum == min(error_sum))
n0 = n0_all(labelnn,:);
figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all =[];
EIS_gather=[];EIS_gather_origin=[]; 
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0 );
        
    plo=plot(eisdata(start_freq_index:end,2)*1000000,eisdata(start_freq_index:end,3)*1000000,'o','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    plot_all=[plot_all plo];
    plo=plot(solve_Z(:,1)*1000000,-solve_Z(:,2)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,1)*1000000;-solve_Z(:,2)*1000000]'];
    EIS_gather_origin=[EIS_gather_origin;[eisdata(start_freq_index:end,2)*1000000;eisdata(start_freq_index:end,3)*1000000]'];
end
label={'+','x','s',"diamond","pentagram"};
eis_label0=[11 21 31 41];
for k=1:1:length(eis_label0)
    real0=EIS_gather(:,eis_label0(k));
    imag0=EIS_gather(:,41+eis_label0(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    real0=EIS_gather_origin(:,eis_label0(k));
    imag0=EIS_gather_origin(:,41+eis_label0(k));
    plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end
% plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
% plo_all=[plo_all plo];
% plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
% plo_all=[plo_all plo];
leg_all={'EIS1','DRT_RL1','EIS2','DRT_RL2','EIS3','DRT_RL3','EIS4','DRT_RL4','EIS5','DRT_RL5','EIS6','DRT_RL6','100 Hz','10 Hz','1 Hz','0.1 Hz'};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis equal
axis([0.369 0.669 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all,'Interpreter','none');
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',3,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%% 考虑耦合项 不考虑L部分
figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all =[];
EIS_gather=[];EIS_gather_origin=[]; 
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0 );
    plo=plot(solve_Z(:,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000;-solve_Z(:,4)*1000000]'];
    EIS_gather_origin=[EIS_gather_origin;[eisdata(:,3)*1000000;eisdata(:,4)*1000000]'];
end
label={'+','x','s',"diamond","pentagram"};
eis_label0=[11 21 31 41];
for k=1:1:length(eis_label0)
    real0=EIS_gather(:,eis_label0(k));
    imag0=EIS_gather(:,41+eis_label0(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'DRT_RL1-RL','DRT_RL2-RL','DRT_RL3-RL','DRT_RL4-RL','DRT_RL5-RL','DRT_RL6-RL','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis equal
axis([0.34 0.669 -0.079 0.24]*1000)
leg=legend(plot_all,leg_all,'Interpreter','none');
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%% 考虑耦合项 不考虑L R部分
figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all =[];
EIS_gather=[];EIS_gather_origin=[]; 
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0 );
    plo=plot(solve_Z(:,3)*1000000 - solve_Z(1,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000 - solve_Z(1,3)*1000000;-solve_Z(:,4)*1000000]'];
    EIS_gather_origin=[EIS_gather_origin;[eisdata(:,2)*1000000;eisdata(:,3)*1000000]'];
end
label={'+','x','s',"diamond","pentagram"};
eis_label0=[11 21 31 41];
for k=1:1:length(eis_label0)
    real0=EIS_gather(:,eis_label0(k));
    imag0=EIS_gather(:,41+eis_label0(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'DRT_RL1-RL-R0','DRT_RL2-RL-R0','DRT_RL3-RL-R0','DRT_RL4-RL-R0','DRT_RL5-RL-R0','DRT_RL6-RL-R0','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis equal
axis([-0.01 0.24 -0.02 0.23]*1000)
leg=legend(plot_all,leg_all,'Interpreter','none');
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca