% A script to process some of the results from pilot sims
%
% In these files the data columns are:
%
% 1. time 
% 2. centre_node 
% 3. bottom_of_centre (outside scar)
% 4. left tissue (outside scar)
% 5. right tissue (outside scar)
% 6. left scar (inside scar)
% 7. right scar (inside scar)
% 8. side scar (inside scar)
% 9. base of tissue (far from scar, same activation time as centre).

close all
clear all

%%
% A script to process some of the results from pilot sims

%% Now load from central node in neutral cells with altered capacitance

conductivities = [0.1];
%neutral_capacitances = [1; 0.1; 0.01; 0.001; 0.0001];
neutral_capacitances = [1; 0.1; 0.01; 0.001];

index_scar = 0; % 0 in the C++ file - centre of scar
index_tissue = 7; % 7 in the C++ file - middle of bottom edge in 2D

cell_types = {'Neutral','Fibroblast'};

for t = 1:length(cell_types)
    
    for c = 1:length(conductivities)
        
        % This is from a simulation with myocytes everywhere and no change to
        % conductivity etc.
        data = load('results_capacitance/Myocyte_scar_cond_1_cap_1_boundary_cond_1_period_150_end_150.dat');
        times = data(:,1) - 5;
        
        a = figure;
        plot(times, data(:,index_scar+2),'k-','LineWidth',2.0)
        hold all
        legend_entries = {'Control: myocytes'};
        
        b=figure;
        
        
        for i=1:length(neutral_capacitances)
            figure(a)
            data = load(['results_capacitance/' cell_types{t}  '_scar_cond_' num2str(conductivities(c))...
                '_cap_' num2str(neutral_capacitances(i)) '_boundary_cond_1_period_150_end_150.dat']);
            plot(times, data(:,index_scar+2))
            legend_entries = {legend_entries{:} [cell_types{t} ' p =  ' num2str(conductivities(c)/neutral_capacitances(i))] };
            
            figure(b)
            subplot(1,length(neutral_capacitances),i)
            lw = {2,1,1,2};
            ls = {'-','--','-.',':'};
            temp = 1;
            for j=[0 4 5 6]
                plot(times, data(:,j+2),ls{temp},'LineWidth',lw{temp})
                hold all
                temp = temp+1;
            end
            xlim([0 75])
            
            ylim([-90 10])
            legend({'center','left','right','bottom'},'Location','NorthEast')
            xlabel('Time (ms)', 'FontSize',16)
            ylabel('Voltage (mV)', 'FontSize',16)
            title(['p = ' num2str(conductivities(c)/neutral_capacitances(i))], 'FontSize',16)
            set(gca,'FontSize',14)
        end
        
        figure(a)
        xlim([0 75])
        xlabel('Time (ms)')
        ylabel('Voltage (mV)')
        
        legend(legend_entries,'Location','NorthEast')
        title('Action potentials recorded in the center of the lesion region')
        
        
        % Now load up summary stats
        
        % For the centre of the scar for now and 'side' of normal tissue
        % Nodes 0 and 1 in the summary stats
        
        % The initial scar is same as zero conductivty from:
        
        for i=1:length(neutral_capacitances)
            recovering = importdata(['results_capacitance/' cell_types{t} '_scar_cond_' num2str(conductivities(c))...
                '_cap_' num2str(neutral_capacitances(i)) '_boundary_cond_1_period_150_end_150_APs.dat']);
            
            % The scar region in early recovery - let's say conductance is
            % 0.001
            early = importdata(['results_capacitance/' cell_types{t} '_scar_cond_0.001_cap_1_boundary_cond_1_period_150_end_150_APs.dat']);
            
            % Plus ones to go from C to matlab indexing
            full_data = [early.data([index_tissue+1 index_scar+1],:); ...
                recovering.data([index_tissue+1 index_scar+1],:)];
            
            % Action potential amplitudes
            amplitudes = full_data(:,7);
            
            % We don't show resting potential for the destroyed cells
            % at early time (it would stay at whatever value we set it
            % in the simulation as diffusion is switched off in this
            % region).
            resting = [full_data(1,8); 0; full_data([3 4],8)];
            
            figure
            subplot(1,4,1)
            bar(-resting,'grouped')
            ylabel('-RMP (mV)','FontSize',16)
            set(gca,'XTickLabel',{'R7','L7','R30','L30'},'FontSize',14)
               
            
            subplot(1,4,2)
            bar(amplitudes,'grouped')
            ylabel('APA (mV)','FontSize',16)
            set(gca,'XTickLabel',{'R7','L7','R30','L30'},'FontSize',14)
               
            % dV/dt max
            upstrokes = full_data(:,5);
            subplot(1,4,3)
            bar(upstrokes,'grouped')
            ylabel('dV/dt (V/s)','FontSize',16)
            title([cell_types{t} 's: p = ' num2str(conductivities(c)/neutral_capacitances(i))])
            set(gca,'XTickLabel',{'R7','L7','R30','L30'},'FontSize',14)
                     
            apds = full_data(:,[4 3 2]);
            subplot(1,4,4)
            bar(apds','grouped')
            
            ylabel('APD (ms)','FontSize',16)
            set(gca,'XTickLabel',{'APD50','APD70','APD90'},'FontSize',14)
            
            if (strcmp(cell_types{t},'Fibroblast') && (conductivities(c)/neutral_capacitances(i))==100)
                resting
                amplitudes
                upstrokes
                apds
            end
            
        end
    end
    
end