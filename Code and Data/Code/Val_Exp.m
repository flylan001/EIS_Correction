%% ========================================================================
%  EIS Experimental Data Visualization and Correction
%  Author: LAN
%  Date: 2025-11-12
%  ------------------------------------------------------------------------
%  This script:
%   (1) Plots raw EIS experimental data (µΩ)
%   (2) Scans parameters (a, b) for DRT-based reconstruction
%   (3) Evaluates RMSE and visualizes parameter sensitivity
%  ========================================================================

clc; clear; close all;

load val_exp_eis.mat
%%
figure,
color=hsv(11);
set(gcf, 'Position',[0 0 800 800],'Color','White');
for i=1:1:size(val_exp_eis,1)
    eisdata=val_exp_eis{i,1};
    plo=plot(eisdata(:,2)*1000-min(eisdata(:,2))*1000*0,eisdata(:,3)*1000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    axis equal
    eisdata0 = eisdata(:,1:3);
    eisdata0(:,3) = -eisdata(:,3);
end
xlabel("\itZ' \rm(m\Omega)");ylabel('-\itZ" \rm(m\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([-2.8 3.5 -6 0.3])
gca=boxoff(gca);clear gca
%% 
figure,
color=hsv(11);
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];plot_all =[];
for i=1:1:size(val_exp_eis,1)
    eisdata=val_exp_eis{i,1};
    plo=plot(eisdata(:,2)*1000000-min(eisdata(:,2))*1000000*0,eisdata(:,3)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    plot_all=[plot_all plo];
    axis equal
    EIS_gather=[EIS_gather;[eisdata(:,2)*1000000;eisdata(:,3)*1000000]'];
    eisdata0 = eisdata(:,1:3);
    eisdata0(:,3) = -eisdata(:,3);
end
label={'+','x','s',"diamond","pentagram"};
eis_label=[21 31 41 51];
for k=1:1:length(eis_label)
    real0=EIS_gather(:,eis_label(k));
    imag0=EIS_gather(:,51+eis_label(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'EIS1','EIS2','EIS3','EIS4','EIS5','EIS6','EIS7','EIS8','EIS9','EIS10','100 Hz','10 Hz','1 Hz','0.1 Hz',' ',' '};
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([0.22 0.44 -0.085 0.135]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',3,'box','off');
leg.ItemTokenSize=[60,18];  %
gca=boxoff(gca);clear gca
%%%%%%%%%%%%%%%%%%%%%%%%%%% 10kHz-0.1Hz
%% 参数ab遍历
start_freq_index=1;
figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');
n=0:0.1:2;
for j=1:1:length(n)
    for j1=1:1:length(n)
        n0 = [n(j),n(j1)]
        j00 = (j-1) * length(n) + j1;
        for i0=1:1:size(val_exp_eis,1)
            eisdata=val_exp_eis{i0,1};
            [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0);
            EIS_all0{i0,1} = solve_Z;
            EIS_all0{i0,2} = eisdata;
        end
        EIS_all{j00}=EIS_all0;
    end
    j;
end
%% 参数a,b范围内平均MSE
n0_all=[];error_data=[];
for j=1:1:length(n)
    for j1=1:1:length(n)
        n0_all = [n0_all;n(j),n(j1)];
        j00 = (j-1) * length(n) + j1;
        EIS_all0 = EIS_all{j00};
        for i0=1:1:size(val_exp_eis,1)
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
colorbar('Ticks',[0.3,0.6,0.9,1.0,1.1]);
title("      Reciprocal of Average RMSE(\rm\mu\Omega)")
a = reshape(error_sum,[length(n),length(n)]);
% 获取 colorbar 对象
cb = colorbar;
cb.TickLabels = compose('%.1f', cb.Ticks);

% 反转Y轴坐标
set(gca, 'YDir', 'normal');
xlabel("\ita");ylabel('\itb');xtickformat('%.1f'); ytickformat('%.1f'); 
set (gca,'LineWidth',3,'position',[0.16,0.18,0.7,0.7],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
set(gca,'xtick',0.0:0.5:1000,'ytick',0.0:0.5:1000,'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
gca=boxoff(gca);clear gca

%% 考虑耦合项 不考虑L R部分
start_freq_index=1;
labelnn = find(error_sum == min(error_sum));
n0 = n0_all(labelnn,:);
eis_label=[1 2 3 4 5 6 7 8 9 10];
figure,
color=hsv(11);
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all =[];
EIS_gather=[];EIS_exp_gather=[]; 
for i=1:1:size(val_exp_eis,1)
    eisdata=val_exp_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0 );
    plo=plot(solve_Z(:,3)*1000000 - solve_Z(1,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000 - solve_Z(1,3)*1000000;-solve_Z(:,4)*1000000]'];
    EIS_exp_gather=[EIS_exp_gather;[eisdata(:,2)*1000000;eisdata(:,3)*1000000]'];
end
label={'+','x','s',"diamond","pentagram"};
eis_label0=[11 21 31 41] + 10 ;
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
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'DRT_RL1-RL-R0','DRT_RL2-RL-R0','DRT_RL3-RL-R0','DRT_RL4-RL-R0','DRT_RL5-RL-R0','DRT_RL6-RL-R0','DRT_RL7-RL-R0','DRT_RL8-RL-R0','DRT_RL9-RL-R0','DRT_RL10-RL-R0','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  ','  ','  ','  ','  '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set (gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1,'FontName','Times New Roman');
set(gca,'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis equal
axis([-0.03 0.19 -0.025 0.195]*1000)
leg=legend(plot_all,leg_all,'Interpreter','none');
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%%%%%%%%%%%%%%%%%%%%%%%%%%% 1kHz-0.1Hz
%% 参数ab遍历
eis_label=[1 2 3 4 5 6 7 8 9 10];
start_freq_index=11;
figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');
EIS_gather=[];
n=0:0.1:2;
for j=1:1:length(n)
    for j1=1:1:length(n)
        n0 = [n(j),n(j1)]
        j00 = (j-1) * length(n) + j1;
        for i0=1:1:size(val_exp_eis,1)
            eisdata=val_exp_eis{i0,1};
            [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:1:end,1),eisdata(start_freq_index:1:end,2),eisdata(start_freq_index:1:end,3)],0,n0);
            EIS_all0{i0,1} = solve_Z;
            EIS_all0{i0,2} = eisdata;
        end
        EIS_all{j00}=EIS_all0;
    end
    j;
end
%% 参数a,b范围内平均MSE
n0_all=[];error_data=[];
for j=1:1:length(n)
    for j1=1:1:length(n)
        n0_all = [n0_all;n(j),n(j1)];
        j00 = (j-1) * length(n) + j1;
        EIS_all0 = EIS_all{j00};
        for i0=1:1:size(val_exp_eis,1)
            solve_Z = EIS_all0{i0,1};
            eisdata = EIS_all0{i0,2};
            error_data(:,i0) = [solve_Z(:,3)-solve_Z(1,3);solve_Z(:,4)];
        end
        error = std(error_data,0,2);
        error_sum(j00,1) = mean(error)*1000000;
    end
end
figure
plot3(n0_all(:,1)+2,n0_all(:,2),error_sum,'*')

figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all=[];
x=[min(n0_all(:,1)),max(n0_all(:,1))];
y=[min(n0_all(:,2)),max(n0_all(:,2))];
imagesc(x,y,1./reshape(error_sum,[21,21]));
colorbar('Ticks',[0.3,0.6,0.9,1.2,1.5]);
title("1/RMSE(\rm\mu\Omega) of \ita \rmand \itb")
% 反转Y轴坐标
set(gca, 'YDir', 'normal');
xlabel("\ita");ylabel('\itb');xtickformat('%.1f'); ytickformat('%.1f'); 
set (gca,'LineWidth',3,'position',[0.16,0.18,0.7,0.7],'FontSize',32,'LabelFontSizeMultiplier',1,'FontName','Times New Roman');
set(gca,'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');

set(gca,'xtick',0.0:0.5:1000,'ytick',0.0:0.5:1000,'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
gca=boxoff(gca);clear gca

%% %% 考虑耦合项 原始数据拟合
eis_label=[1 2 3 4 5 6 7 8 9 10];
start_freq_index=11;
labelnn = find(error_sum == min(error_sum))
n0 = n0_all(labelnn,:);
figure,
color=hsv(11);
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];plot_all =[];
EIS_gather=[];EIS_exp_gather=[]; 
for i=1:1:size(val_exp_eis,1)
    eisdata=val_exp_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:1:end,1),eisdata(start_freq_index:1:end,2),eisdata(start_freq_index:1:end,3)],0,n0);
    plo=plot(eisdata(start_freq_index:end,2)*1000000,eisdata(start_freq_index:end,3)*1000000,'o','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    plot_all=[plot_all plo];
    plo=plot(solve_Z(:,1)*1000000,-solve_Z(:,2)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,1)*1000000;-solve_Z(:,2)*1000000]'];
    EIS_exp_gather=[EIS_exp_gather;[eisdata(start_freq_index:end,2)*1000000;eisdata(start_freq_index:end,3)*1000000]'];
end
label={'+','x','s',"diamond","pentagram"};
eis_label0=[11 21 31 41];
for k=1:1:length(eis_label0)
    real0=EIS_gather(:,eis_label0(k));
    imag0=EIS_gather(:,41+eis_label0(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    real0=EIS_exp_gather(:,eis_label0(k));
    imag0=EIS_exp_gather(:,41+eis_label0(k));
    plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end
% plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
% plo_all=[plo_all plo];
% plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
% plo_all=[plo_all plo];
leg_all={'EIS1','DRT_RL1','EIS2','DRT_RL2','EIS3','DRT_RL3','EIS4','DRT_RL4','EIS5','DRT_RL5','EIS6','DRT_RL6','EIS7','DRT_RL7','EIS8','DRT_RL8','EIS9','DRT_RL9','EIS10','DRT_RL10','100 Hz','10 Hz','1 Hz','0.1 Hz'};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set (gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1,'FontName','Times New Roman');
set(gca,'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis equal
axis([0.22 0.44 -0.071 0.149]*1000)
leg=legend(plot_all,leg_all,'Interpreter','none');
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',3,'box','off');
leg.ItemTokenSize=[30,18];  %
gca=boxoff(gca);clear gca
save_path=strcat('EIS paraRL000 1k');
print(gcf, save_path,'-djpeg','-r1200');
print(gcf, save_path,'-dsvg','-r1200');

%% 考虑耦合项 不考虑L部分

figure,
color=hsv(11);
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all =[];
EIS_gather=[];EIS_exp_gather=[]; 
for i=1:1:size(val_exp_eis,1)
    eisdata=val_exp_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0 );
    plo=plot(solve_Z(:,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000;-solve_Z(:,4)*1000000]'];
    EIS_exp_gather=[EIS_exp_gather;[eisdata(:,2)*1000000;eisdata(:,3)*1000000]'];
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
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'DRT_RL1-RL-R0','DRT_RL2-RL-R0','DRT_RL3-RL-R0','DRT_RL4-RL-R0','DRT_RL5-RL-R0','DRT_RL6-RL-R0','DRT_RL7-RL-R0','DRT_RL8-RL-R0','DRT_RL9-RL-R0','DRT_RL10-RL-R0','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  ','  ','  ','  ','  '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set (gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1,'FontName','Times New Roman');
set(gca,'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis equal
axis([0.22 0.44 -0.03 0.19]*1000)
leg=legend(plot_all,leg_all,'Interpreter','none');
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca


%% 考虑耦合项 不考虑L R部分
figure,
color=hsv(11);
set(gcf, 'Position',[0 0 800 800],'Color','White');plot_all =[];
EIS_gather=[];EIS_exp_gather=[]; 
for i=1:1:size(val_exp_eis,1)
    eisdata=val_exp_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan_RL_n_dif([eisdata(start_freq_index:end,1),eisdata(start_freq_index:end,2),eisdata(start_freq_index:end,3)],0,n0 );
    plo=plot(solve_Z(:,3)*1000000 - solve_Z(1,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000 - solve_Z(1,3)*1000000;-solve_Z(:,4)*1000000]'];
    EIS_exp_gather=[EIS_exp_gather;[eisdata(:,2)*1000000;eisdata(:,3)*1000000]'];
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
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'DRT_RL1-RL-R0','DRT_RL2-RL-R0','DRT_RL3-RL-R0','DRT_RL4-RL-R0','DRT_RL5-RL-R0','DRT_RL6-RL-R0','DRT_RL7-RL-R0','DRT_RL8-RL-R0','DRT_RL9-RL-R0','DRT_RL10-RL-R0','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  ','  ','  ','  ','  '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set (gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1,'FontName','Times New Roman');
set(gca,'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis equal
axis([-0.03 0.19 -0.025 0.195]*1000)
leg=legend(plot_all,leg_all,'Interpreter','none');
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca