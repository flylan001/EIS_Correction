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

%%
figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    plo=plot(eisdata(:,2)*1000-min(eisdata(:,2))*1000*0,eisdata(:,3)*1000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    axis equal
    eisdata0 = eisdata(:,1:3);
    eisdata0(:,3) = -eisdata(:,3);
    % [x,Zs1,Zs2] = R_RCE_W_RL(eisdata0);
    % plo=plot(real(Zs1)*1000000,-imag(Zs1)*1000000,'o','Color',color(i,:),'LineWidth',3,'MarkerSize',8);
    % Zs2{i} = Zs2;
end
% label={'^','+','x','s',"diamond","pentagram"};

xlabel("\itZ' \rm(m\Omega)");ylabel('-\itZ" \rm(m\Omega)');
set (gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1,'FontName','Times New Roman');
set(gca,'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([-2.8 3.5 -6 0.3])
gca=boxoff(gca);clear gca

%%
figure,
color=hsv(8);
set(gcf, 'Position',[0 0 800 800],'Color','White');EIS_gather=[];plot_all=[];
for i=1:1:size(exp1_eis,1)
    eisdata=exp1_eis{i,1};
    plo=plot(eisdata(:,2)*1000000,eisdata(:,3)*1000000,'-*','Color',color(i,:),'LineWidth',3,'MarkerSize',8);hold on;
    plot_all=[plot_all plo];
    axis equal
    EIS_gather=[EIS_gather;[eisdata(:,2)*1000000;eisdata(:,3)*1000000]'];
    eisdata0 = eisdata(:,1:3);
    eisdata0(:,3) = -eisdata(:,3);
    % [x,Zs1,Zs2] = R_RCE_W_RL(eisdata0);
    % plo=plot(real(Zs1)*1000000,-imag(Zs1)*1000000,'o','Color',color(i,:),'LineWidth',3,'MarkerSize',8);
    % Zs2{i} = Zs2;
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
leg_all={'EIS1','EIS2','EIS3','EIS4','EIS5','EIS6','100 Hz','10 Hz','1 Hz','0.1 Hz',' ',' '};
xlabel("\itZ' \rm(\mu\Omega)");ylabel('-\itZ" \rm(\mu\Omega)');
set (gca,'LineWidth',3,'position',[0.18,0.18,0.8,0.8],'FontSize',32,'LabelFontSizeMultiplier',1,'FontName','Times New Roman');
set(gca,'FontSize',32,'LabelFontSizeMultiplier',1.2,'FontName','Times New Roman');
axis([0.37 0.67 -0.06 0.24]*1000)
leg=legend(plot_all,leg_all);
set(leg,'Location','North','FontSize',24,'FontName','Times New Roman','Orientation','v','NumColumns',2,'box','off');
leg.ItemTokenSize=[60,18];  %
gca=boxoff(gca);clear gca
