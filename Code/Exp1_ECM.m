%% ========================================================================
%  EIS Experimental Data Visualization and Correction
%  Author: LAN
%  Date: 2025-11-12
%  ------------------------------------------------------------------------
%  This script:
%   (1) Plots raw EIS experimental data (µΩ)
%   (2) Scans parameters (a, b) for DRT-based reconstruction
%   (3) Evaluates MSE and visualizes parameter sensitivity
%  ========================================================================

clc; clear; close all;

load exp1_eis.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%% 10kHz-0.1Hz
%% 原始数据 ECM拟合
color=hsv(8);
start_freq_index=1;
figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];EIS_exp_gather=[];Zs2_all = [];x_all= [];plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    plo=plot(eisdata(:,2)*1000000,eisdata(:,3)*1000000,'o','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    plot_all=[plot_all plo];
    axis equal
    EIS_gather=[EIS_gather;[eisdata(start_freq_index:end,2)*1000000;eisdata(start_freq_index:end,3)*1000000]'];
    eisdata0 = eisdata(start_freq_index:end,1:3);
    eisdata0(:,3) = -eisdata(start_freq_index:end,3);
    [x,Zs1,Zs2] = R_RCE_W_RL(eisdata0);
    x_all= [x_all;x];
    plo=plot(real(Zs1)*1000000,-imag(Zs1)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[real(Zs1)*1000000;-imag(Zs1)*1000000]'];
    Zs2_all{i} = Zs2;
    EIS_exp_gather=[EIS_exp_gather;[eisdata(:,3)*1000000;eisdata(:,4)*1000000]'];

end

label={'+','x','s',"diamond","pentagram"};
eis_label=[21 31 41 51];
for k=1:1:length(eis_label)
    real0=EIS_gather(:,eis_label(k));
    imag0=EIS_gather(:,51+eis_label(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    real0=EIS_exp_gather(:,eis_label(k));
    imag0=EIS_exp_gather(:,51+eis_label(k));
    plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end

leg_all={'EIS1','ECM1','EIS2','ECM2','EIS3','ECM3','EIS4','ECM4','EIS5','ECM5','EIS6','ECM6','100 Hz','10 Hz','1 Hz','0.1 Hz'};
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([0.37 0.67 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',3,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

a = x_all';
format short e
disp(a)


%% 原始数据 ECM拟合 除去RL部分
figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];Zs2_all = [];plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    eisdata0 = eisdata(start_freq_index:end,1:3);
    eisdata0(:,3) = -eisdata(start_freq_index:end,3);
    [x,Zs1,Zs2] = R_RCE_W_RL(eisdata0);
    plo=plot(real(Zs2)*1000000,-imag(Zs2)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    axis equal
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[real(Zs2)*1000000;-imag(Zs2)*1000000]'];
    Zs2_all{i} = Zs2;
end
% label={'^','+','x','s',"diamond","pentagram"};
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
leg_all={'ECM1-RL','ECM2-RL','ECM3-RL','ECM4-RL','ECM5-RL','ECM6-RL','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([0.37 0.67 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca


%% 原始数据 ECM拟合 除去R0部分
eis_label=[2 4 6 8 10 12];
figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];Zs2_all = [];plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    eisdata0 = eisdata(start_freq_index:end,1:3);
    eisdata0(:,3) = -eisdata(start_freq_index:end,3);
    [x,Zs1,Zs2] = R_RCE_W_RL(eisdata0);
    plo=plot(real(Zs2)*1000000-x(2)*1000000,-imag(Zs2)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    axis equal
    EIS_gather=[EIS_gather;[real(Zs2)*1000000-x(2)*1000000;-imag(Zs2)*1000000]'];
    Zs2_all{i} = Zs2;
end
% label={'^','+','x','s',"diamond","pentagram"};
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
leg_all={'ECM1-RL-R0','ECM2-RL-R0','ECM3-RL-R0','ECM4-RL-R0','ECM5-RL-R0','ECM6-RL-R0','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([-0.01 0.24 -0.02 0.23]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

a = EIS_gather(1:1:end,:);
error = std(a',0,2);
error_sum = mean(error);

%%%%%%%%%%%%%%%%%%%%%%%%%%% 1kHz-0.1Hz
%% 原始数据 ECM拟合
figure,
color=hsv(8);
start_freq_index=11;
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];EIS_exp_gather=[];Zs2_all = [];x_all= [];plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    plo=plot(eisdata(:,2)*1000000,eisdata(:,3)*1000000,'o','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    plot_all=[plot_all plo];
    axis equal
    % EIS_gather=[EIS_gather;[eisdata(start_num:end,2)*1000000;eisdata(start_num:end,3)*1000000]'];
    eisdata0 = eisdata(start_freq_index:end,1:3);
    eisdata0(:,3) = -eisdata(start_freq_index:end,3);
    [x,Zs1,Zs2] = R_RCE_W_RL(eisdata0);
    x_all= [x_all;x];
    plo=plot(real(Zs1)*1000000,-imag(Zs1)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[real(Zs1)*1000000;-imag(Zs1)*1000000]'];
    Zs2_all{i} = Zs2;
    EIS_exp_gather=[EIS_exp_gather;[eisdata(start_freq_index:end,2)*1000000;eisdata(start_freq_index:end,3)*1000000]'];
    % RMSE10(i) = rmse(eisdata(:,2)*1000000,real(Zs1)*1000000);
    % RMSE20(i) = rmse(eisdata(:,3)*1000000,-imag(Zs1)*1000000);
end
% label={'^','+','x','s',"diamond","pentagram"};
label={'+','x','s',"diamond","pentagram"};
eis_label=[21 31 41 51]-10;
for k=1:1:length(eis_label)
    real0=EIS_gather(:,eis_label(k));
    imag0=EIS_gather(:,41+eis_label(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    real0=EIS_exp_gather(:,eis_label(k));
    imag0=EIS_exp_gather(:,41+eis_label(k));
    plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end

leg_all={'EIS1','ECM1','EIS2','ECM2','EIS3','ECM3','EIS4','ECM4','EIS5','ECM5','EIS6','ECM6','100 Hz','10 Hz','1 Hz','0.1 Hz'};
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([0.37 0.67 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',3,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

a = x_all';
format short e
disp(a)


%% 原始数据 ECM拟合 除去RL部分

figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];Zs2_all = [];    plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    eisdata0 = eisdata(start_freq_index:end,1:3);
    eisdata0(:,3) = -eisdata(start_freq_index:end,3);
    [x,Zs1,Zs2] = R_RCE_W_RL(eisdata0);
    plo=plot(real(Zs2)*1000000,-imag(Zs2)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    axis equal
    EIS_gather=[EIS_gather;[real(Zs2)*1000000;-imag(Zs2)*1000000]'];
    Zs2_all{i} = Zs2;
end
% label={'^','+','x','s',"diamond","pentagram"};
label={'+','x','s',"diamond","pentagram"};
eis_label=[21 31 41 51]-10;
for k=1:1:length(eis_label)
    real0=EIS_gather(:,eis_label(k));
    imag0=EIS_gather(:,41+eis_label(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'ECM1-RL','ECM2-RL','ECM3-RL','ECM4-RL','ECM5-RL','ECM6-RL','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([0.37 0.67 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%% 原始数据 ECM拟合 除去R0部分
figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];Zs2_all = [];plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    eisdata0 = eisdata(start_freq_index:end,1:3);
    eisdata0(:,3) = -eisdata(start_freq_index:end,3);
    [x,Zs1,Zs2] = R_RCE_W_RL(eisdata0);
    plo=plot(real(Zs2)*1000000-x(2)*1000000,-imag(Zs2)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    axis equal
    EIS_gather=[EIS_gather;[real(Zs2)*1000000-x(2)*1000000;-imag(Zs2)*1000000]'];
    Zs2_all{i} = Zs2;
end
% label={'^','+','x','s',"diamond","pentagram"};
label={'+','x','s',"diamond","pentagram"};
eis_label=[21 31 41 51]-10;
for k=1:1:length(eis_label)
    real0=EIS_gather(:,eis_label(k));
    imag0=EIS_gather(:,41+eis_label(k));
    plo=plot(real0,imag0,label{k},'MarkerSize',16,'MarkerEdgeColor','k');
    plot_all=[plot_all plo];
end
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'ECM1-RL-R0','ECM2-RL-R0','ECM3-RL-R0','ECM4-RL-R0','ECM5-RL-R0','ECM6-RL-R0','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([-0.01 0.24 -0.02 0.23]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

a = EIS_gather(1:1:end,:);
error = std(a',0,2);
error_sum = mean(error);