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

load exp1_eis.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%% 10kHz-0.1Hz
%% DRT
color=hsv(8);
start_num = 1;
figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];EIS_exp_gather=[]; plot_all =[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan([eisdata(start_num:end,1),eisdata(start_num:end,2),eisdata(start_num:end,3)],0,1e-8);
    plo=plot(eisdata(:,2)*1000000,eisdata(:,3)*1000000,'o','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    plot_all=[plot_all plo];
    plo=plot(solve_Z(:,1)*1000000,-solve_Z(:,2)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,1)*1000000;-solve_Z(:,2)*1000000]'];
    EIS_exp_gather=[EIS_exp_gather;[eisdata(:,2)*1000000;eisdata(:,3)*1000000]'];
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
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'EIS1','DRT1','EIS2','DRT2','EIS3','DRT3','EIS4','DRT4','EIS5','DRT5','EIS6','DRT6','100 Hz','10 Hz','1 Hz','0.1 Hz',' ',' '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([0.369 0.669 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',3,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%% DRT-L
figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan([eisdata(start_num:end,1),eisdata(start_num:end,2),eisdata(start_num:end,3)],0,0);
    plo=plot(solve_Z(:,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000;-solve_Z(:,4)*1000000]'];
end
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

leg_all={'DRT1-RL','DRT2-RL','DRT3-RL','DRT4-RL','DRT5-RL','DRT6-RL','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([0.369 0.669 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%% DRT-L-R

figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan([eisdata(start_num:end,1),eisdata(start_num:end,2),eisdata(start_num:end,3)],0,0);
    plo=plot(solve_Z(:,3)*1000000-solve_Z(1,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000-solve_Z(1,3)*1000000;-solve_Z(:,4)*1000000]'];
end
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

leg_all={'DRT1-RL-R0','DRT2-RL-R0','DRT3-RL-R0','DRT4-RL-R0','DRT5-RL-R0','DRT6-RL-R0','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([-0.01 0.24 -0.02 0.23]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%%%%%%%%%%%%%%%%%%%%%%%%%%% 1000Hz-0.1Hz
%% DRT
start_num=11;
figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];EIS_exp_gather=[]; plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan([eisdata(start_num:end,1),eisdata(start_num:end,2),eisdata(start_num:end,3)],0,1e-8);
    plo=plot(eisdata(:,2)*1000000,eisdata(:,3)*1000000,'o','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    plot_all=[plot_all plo];
    plo=plot(solve_Z(:,1)*1000000,-solve_Z(:,2)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,1)*1000000;-solve_Z(:,2)*1000000]'];
    EIS_exp_gather=[EIS_exp_gather;[eisdata(start_num:end,2)*1000000;eisdata(start_num:end,3)*1000000]'];
end
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
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
plo=plot(-1,1,label{k},'MarkerSize',16,'MarkerEdgeColor','w');
plot_all=[plot_all plo];
leg_all={'EIS1','DRT1','EIS2','DRT2','EIS3','DRT3','EIS4','DRT4','EIS5','DRT5','EIS6','DRT6','100 Hz','10 Hz','1 Hz','0.1 Hz',' ',' '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([0.369 0.669 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',3,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%% DRT-L

figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan([eisdata(start_num:end,1),eisdata(start_num:end,2),eisdata(start_num:end,3)],0,0);
    plo=plot(solve_Z(:,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000;-solve_Z(:,4)*1000000]'];
end
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

leg_all={'DRT1-RL','DRT2-RL','DRT3-RL','DRT4-RL','DRT5-RL','DRT6-RL','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([0.369 0.669 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca

%% DRT-L-R
figure,
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    [gamma,x_ridge,solve_Z,epsilon]=DRT_lasso_lan([eisdata(start_num:end,1),eisdata(start_num:end,2),eisdata(start_num:end,3)],0,0);
    plo=plot(solve_Z(:,3)*1000000-solve_Z(1,3)*1000000,-solve_Z(:,4)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on
    plot_all=[plot_all plo];
    EIS_gather=[EIS_gather;[solve_Z(:,3)*1000000-solve_Z(1,3)*1000000;-solve_Z(:,4)*1000000]'];
end
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

leg_all={'DRT1-RL-R0','DRT2-RL-R0','DRT3-RL-R0','DRT4-RL-R0','DRT5-RL-R0','DRT6-RL-R0','100 Hz','10 Hz','1 Hz','0.1 Hz','  ','  '};xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set(gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([-0.01 0.24 -0.02 0.23]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[45,18];  %
gca=boxoff(gca);clear gca