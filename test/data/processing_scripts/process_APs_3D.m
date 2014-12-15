% A script to process some of the results from pilot sims
%
% In these files the data columns are:
%
% 1. time 
% 2. centre_node (scar)
% 3. top_of_centre (scar)
% 4. tissue_beside_scar (tissue)
% 5. top_of_tissue_beside_scar (tissue)
% 6. left_of_scar  (tissue)
% 7. top_of_left_of_scar (tissue)
% 8. right_of_scar (tissue)
% 9. top_of_right_of_scar (tissue)
% 10. left_in_scar (scar)
% 11. right in scar (scar)
% 12. side of scar (scar)
% 13. base of tissue (remotest tissue)

close all
clear all

%thicknesses = 0.05:-0.01:0.01;
thicknesses = 0.01;
capacitances = [0.1 0.01 0.001];

index_scar = 0; % 0 in the C++ file - centre of scar
index_tissue = 11; % 7 in the C++ file - middle of bottom edge in 2D

%% First let's have a look at the normal Myocyte activity
base_filename = 'pilot_results_3D/Myocyte_scar_cond_1_cap_1_boundary_1_period_100_end_80_thickness_';

data = load([base_filename num2str(thicknesses(1)) '_med_res.dat']);
% Synchronise the plotted times with T=0 as stimulus time.
times = data(:,1) - 10.0;

a = figure;
plot(times, data(:,index_scar+2),'k-','LineWidth',2.0) % +2 goes from C++ to matlab (plus a time column) indexing.
hold all
legend_entries{1} = ['Control: myocytes'];

% This is like the 7 day one.
early = importdata('pilot_results_3D/Neutral_scar_cond_0.001_cap_1_boundary_1_period_100_end_80_thickness_0.01_med_res_APs.dat');

for i=1:3
    figure(a)
    base_filename = ['pilot_results_3D/Neutral_scar_cond_0.1_cap_' num2str(capacitances(i)) '_boundary_1_period_100_end_80_thickness_0.01_med_res'];
    data = load([base_filename '.dat']);
    plot(times, data(:,index_scar+2))
    legend_entries{i+1} = ['Neutral capacitance ' num2str(100*capacitances(i)) '%'];
    
    % This is like the 30 day one.
    recovering = importdata([base_filename '_APs.dat']);

    % Plus ones to go from C to matlab indexing
    full_data = [early.data([index_tissue+1 index_scar+1],:); ...
                 recovering.data([index_tissue+1 index_scar+1],:)]; 

    % Action potential amplitudes
    amplitudes = full_data(:,7);
    resting = full_data(:,8);

    figure
    subplot(1,4,1)
    bar(resting,'grouped')
    ylabel('RMP (mV)')

    subplot(1,4,2)
    bar(amplitudes,'grouped')
    ylabel('APA (mV)')

    % dV/dt max
    upstrokes = full_data(:,5);
    subplot(1,4,3)
    bar(upstrokes,'grouped')
    ylabel('dV/dt (V/s)')
    title(['Scar capacitance of ' num2str(100*neutral_capacitances(i)) '% conductivity ' num2str(100*conductivities(c)) '%'])

    apds = full_data(:,[4 3 2]);
    subplot(1,4,4)
    bar(apds','grouped')

    ylabel('APD (ms)')
    set(gca,'XTickLabel',{'APD50','APD70','APD90'})   
end

figure(a)
xlabel('Time (ms)','FontSize',16)
xlim([0 70])
ylim([-90 40])
ylabel('Voltage (mV)','FontSize',16)
title('Action potential in the centre of the 3D region, conductivity 10%','FontSize',18)
set(gca,'FontSize',14)
legend(legend_entries,'Location','NorthEast')

% %% Myocytes Right hand node
% 
% data = load([base_filename num2str(thicknesses(1)) '_med_res.dat']);
% % Synchronise the plotted times with T=0 as stimulus time.
% times = data(:,1) - 10.0;
% 
% figure
% plot(times, data(:,8),'k-','LineWidth',2.0)
% hold all
% legend_entries{1} = ['Scar thickness ' num2str(thicknesses(1))];
% 
% for i=2:5
% data = load([base_filename num2str(thicknesses(i)) '_med_res.dat']);
% plot(times, data(:,8))
% legend_entries{i} = ['Scar thickness ' num2str(thicknesses(i))];
% end
% 
% xlabel('Time (ms)','FontSize',16)
% xlim([0 80])
% ylim([-90 40])
% 
% ylabel('Voltage (mV)','FontSize',16)
% title('Action potential after the central region','FontSize',18)
% set(gca,'FontSize',14)
% legend(legend_entries,'Location','East')
% 
% %% Neutral Centre
% base_filename = 'pilot_results_3D/Neutral_scar_1_boundary_1_period_100_end_100_thickness_';
% 
% data = load([base_filename num2str(thicknesses(1)) '_med_res.dat']);
% % Synchronise the plotted times with T=0 as stimulus time.
% times = data(:,1) - 10.0;
% 
% figure
% plot(times, data(:,2),'k-','LineWidth',2.0)
% hold all
% legend_entries{1} = ['Scar thickness ' num2str(thicknesses(1))];
% 
% for i=2:5
% data = load([base_filename num2str(thicknesses(i)) '_med_res.dat']);
% plot(times, data(:,2))
% legend_entries{i} = ['Scar thickness ' num2str(thicknesses(i))];
% end
% 
% xlabel('Time (ms)','FontSize',16)
% xlim([0 80])
% ylim([-90 40])
% ylabel('Voltage (mV)','FontSize',16)
% title('Action potential in the centre of neutral scar','FontSize',18)
% set(gca,'FontSize',14)
% legend(legend_entries,'Location','East')
% 
% %% Neutral After
% data = load([base_filename num2str(thicknesses(1)) '_med_res.dat']);
% % Synchronise the plotted times with T=0 as stimulus time.
% times = data(:,1) - 10.0;
% 
% figure
% plot(times, data(:,8),'k-','LineWidth',2.0)
% hold all
% legend_entries{1} = ['Scar thickness ' num2str(thicknesses(1))];
% 
% for i=2:5
% data = load([base_filename num2str(thicknesses(i)) '_med_res.dat']);
% plot(times, data(:,8))
% legend_entries{i} = ['Scar thickness ' num2str(thicknesses(i))];
% end
% 
% xlabel('Time (ms)','FontSize',16)
% xlim([0 80])
% ylim([-90 40])
% ylabel('Voltage (mV)','FontSize',16)
% title('Action potential after neutral scar','FontSize',18)
% set(gca,'FontSize',14)
% legend(legend_entries,'Location','East')
% 
% %% Neutral Centre 0.4 conductivity
% base_filename = 'pilot_results_3D/Neutral_scar_0.4_boundary_1_period_100_end_100_thickness_';
% 
% data = load([base_filename num2str(thicknesses(1)) '_med_res.dat']);
% % Synchronise the plotted times with T=0 as stimulus time.
% times = data(:,1) - 10.0;
% 
% figure
% plot(times, data(:,2),'k-','LineWidth',2.0)
% hold all
% legend_entries{1} = ['Scar thickness ' num2str(thicknesses(1))];
% 
% for i=2:5
% data = load([base_filename num2str(thicknesses(i)) '_med_res.dat']);
% plot(times, data(:,2))
% legend_entries{i} = ['Scar thickness ' num2str(thicknesses(i))];
% end
% 
% xlabel('Time (ms)','FontSize',16)
% xlim([0 80])
% ylim([-90 40])
% ylabel('Voltage (mV)','FontSize',16)
% title('Action potential in the centre of neutral scar 0.4','FontSize',18)
% set(gca,'FontSize',14)
% legend(legend_entries,'Location','East')
% 
% %% Neutral After 0.4 conductivity
% data = load([base_filename num2str(thicknesses(1)) '_med_res.dat']);
% % Synchronise the plotted times with T=0 as stimulus time.
% times = data(:,1) - 10.0;
% 
% figure
% plot(times, data(:,8),'k-','LineWidth',2.0)
% hold all
% legend_entries{1} = ['Scar thickness ' num2str(thicknesses(1))];
% 
% for i=2:5
% data = load([base_filename num2str(thicknesses(i)) '_med_res.dat']);
% plot(times, data(:,8))
% legend_entries{i} = ['Scar thickness ' num2str(thicknesses(i))];
% end
% 
% xlabel('Time (ms)','FontSize',16)
% xlim([0 80])
% ylim([-90 40])
% ylabel('Voltage (mV)','FontSize',16)
% title('Action potential after neutral scar 0.4','FontSize',18)
% set(gca,'FontSize',14)
% legend(legend_entries,'Location','East')
% 
