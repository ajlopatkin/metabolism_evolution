clear
close all
clc

global u ks k d eff A1 A2 Nm

%% Base model params
% base parameters
u_cnt = 0.5;                        % max growth rate
ks_cnt = 0.004;                     % substrate halfmax (ug/L)
k_cnt = 3.0;                        % ic50
d_cnt = 0.025;                      % death rate
eff_cnt = 5e-7;                     % conversion efficiency
Nm_cnt = 1e7;                       % carrying capacity
S = 1;                              % initial substrate

% evolving parameters
u = u_cnt;
ks = ks_cnt;
k = k_cnt;
d = d_cnt;
eff = eff_cnt;
Nm = Nm_cnt;
clrs=[36 129 237;42 245 110;222 216 49;247 143 40;208 0 0; 255 255 255]./255;

%% Evolutionary params
% NOISE LEVELS
pct_std_u = 0.05;                    % standard deviation of parameter distributions (assuming all are normally distro except ks)
pct_std_d = 0.2;                     % standard deviation of parameter distribution (assuming all are normally distro except ks)
pct_std_k = 0.2;                     % standard deviation of parameter distributions (assuming all are normally distro except ks)
pct_std_KS = -0.07;                  % the approx <lognormal std dev> of observed lognormal KS values at D0

% script parameters
num_seasons = 2+21;                 % number of seasons to evolve the cells
time_per_season = 24;               % length of time between seasons (in hours)
tspan = 0:0.1:time_per_season;      % ODE time range
iters = 50;                         % number of iterations

% flags
high_res_flag = 0;                  
static_flag = 1;                    % change to 1 to run with static drug
cidal_flag = 1;                     % change to 1 to run with cidal drug

if high_res_flag
    ause = [0,logspace(log10(.1),log10(4),30)];
else
    ause = [0 1 2 3 4];
end

savename = "evol_MAIN_A1_"+ static_flag + "_A2_" + cidal_flag + "_" + date;


% Population size and mutations
mutation_rate = 0.001;              % mutations per genome per cell generation
pop_after_dilution = 1e4;           % number of cells to keep after dilution
starting_n = 1;                     % number of populations to start each round
starting_cells = 1e4;               % number of cells per population, season 1
cap_cutoff = 10^.5;                 % the fold below season1 density at which an iteration must be replaced

% Intialize collection matrices
weighted_u = zeros(iters,num_seasons,length(ause));
weighted_ks = zeros(iters,num_seasons,length(ause));
weighted_k = zeros(iters,num_seasons,length(ause));
weighted_d = zeros(iters,num_seasons,length(ause));
weighted_eff = zeros(iters,num_seasons,length(ause));
weighted_dens = zeros(iters,num_seasons,length(ause));
num_mutants = zeros(iters,num_seasons,length(ause));

% loop parameters
collect_all_dens = cell(iters,num_seasons,length(ause));
collect_all_u = cell(iters,num_seasons,length(ause));
collect_all_k = cell(iters,num_seasons,length(ause));
collect_all_ks = cell(iters,num_seasons,length(ause));
collect_all_d = cell(iters,num_seasons,length(ause));
collect_all_eff = cell(iters,num_seasons,length(ause));

plot_flag = 0;                      % change to 1 if you want to see all cell density
save_flag = 1;


% loop through antibiotics
for aa = 1:length(ause)
    aa
    % loop through iterations
    ii = 1;
    while ii <= iters
        
        % this flag will break us out of an iter when all dens==0
        ALLDEAD_FLAG = 0;
        all_cell_dens = []; % re-initialize all density collection matrix
        all_cell_dens2 = [];
        tic;
        
        % initialize parameters
        dens = starting_cells.*ones(1,starting_n);
        u = repmat(u_cnt,1,starting_n);
        ks = repmat(ks_cnt,1,starting_n);
        k = repmat(k_cnt,1,starting_n);
        d = repmat(d_cnt,1,starting_n);
        eff = repmat(eff_cnt,1,starting_n);
        
        % loop through seaosns
        target_dens = [];
        plateaus = zeros(1,num_seasons);
        end_dens = zeros(1,num_seasons);
        
        for q = 1:num_seasons
            
            % set A conditions
            if q <= 2
                A = 0;
            else
                A = ause(aa);
            end
            
            
            % set static and cidal drugs
            if static_flag
                A1 = A; 
            else
                A1 = 0;
            end
            
            if cidal_flag
                A2 = A;
            else
                A2 = 0;
            end
            
            y0 = [dens S];
            [time, y] = ode45(@ddt_EVO,tspan,y0);
            dens = ceil(y(end,1:end-1));
            dens(dens<1) = 0;
            
            % get the final density of season 1
            if(q==1),target_dens=dens(end);end
            
            % collect all populations that grow
            if size(y,2) == 2
                sumY = y(:,1:end-1);
            else
                sumY = sum(y(:,1:end-1),2);
            end
            
            % did this season reach a plateau?
            tmp = diff(sumY);
            if(tmp(end)<1),plateaus(q)=1;end
            
            % save the end density of this season
            end_dens(q) = sumY(end);
            collect_all_dens{ii,q,aa} = y(1,1:(end-1));
            collect_all_u{ii,q,aa} = u; collect_all_k{ii,q,aa} = k; collect_all_ks{ii,q,aa} = ks;
            collect_all_d{ii,q,aa} = d; collect_all_eff{ii,q,aa} = eff;
            all_cell_dens2 = [all_cell_dens2;sumY]; 
            
            % generate new populations from mutations
            for qq = 1:length(dens)
                
                % generate binomial random based on pop size and mutation prob
                num_gen = log2(dens(qq)./pop_after_dilution);           % calculate number of generations for that population
                mutation_prob = mutation_rate .* num_gen;               % calculate mutation probability for that number of generations
                mutation_count = my_binornd(dens(qq), mutation_prob);   % get expected number of mutations in population
                
                % assume one new mutant population is generated per
                % sub-population
                if mutation_count > 0
                    dens = [dens, mutation_count];
                    
                    % resistance parameters
                    d = [d abs(normrnd(d(qq),pct_std_d*d_cnt))];
                    k = [k,abs(normrnd(k(qq),pct_std_k*k_cnt))];
                    
                    
                    % metabolism parameters                                  
                    u = [u,abs(normrnd(u(qq),pct_std_u*u_cnt))];
                    ks = [ks,lognrnd(log(ks(qq)),abs(pct_std_KS*log(ks_cnt)))];
                                        
                    eff = [eff, (ks(end)/ks_cnt)^2 * eff_cnt];
                                        
                end
            end
            
            % collect weighted cell densities
            ind = dens>=1;
            weighted_dens(ii,q,aa) = sum(dens(ind));
            weighted_u(ii,q,aa) = sum(u(ind).*dens(ind))./sum(dens(ind));
            weighted_ks(ii,q,aa) = sum(ks(ind).*dens(ind))./sum(dens(ind));
            weighted_k(ii,q,aa) = sum(k(ind).*dens(ind))./sum(dens(ind));
            weighted_d(ii,q,aa) = sum(d(ind).*dens(ind))./sum(dens(ind));
            
            weighted_eff(ii,q,aa) = sum(eff(ind).*dens(ind))./sum(dens(ind));
            
            
            
            %%%
            xx = find(cumsum(sort(dens./sum(dens).*100,'descend'))>80,1,'first');

            
            num_mutants(ii,q,aa) = find(cumsum(sort(dens./sum(dens).*100,'descend'))>80,1,'first');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Randomly choose cells from each population to dilute:
            
            % Determine the target size based on if the cells grew or if
            % they are dying:
            target_dilution_size = min(pop_after_dilution, sum(dens));
            
            % get percentage for each population
            pop_pct = dens/sum(dens);
            
            % create CDF and prob function for population
            pop_cdf = cumsum(pop_pct);
            pop_prob = @(r) find(r<pop_cdf,1,'first');
            
            % generate a unit random vector of with `init` elements
            R = rand(target_dilution_size, 1);
            
            % run the population prob function on R
            popR = arrayfun(pop_prob,R);
            
            for i = 1:length(dens)
                dens(i) = sum(popR == i);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end
               
        % if none of the plateaus reached fell below the season1 cutoff,
        % keep this iteration, else repeat
        if ~any((end_dens<target_dens/cap_cutoff) & plateaus)
            
            disp("Time taken for iteration "+ii+": "+toc)
            ii = ii + 1;
            if plot_flag
                figure(1), hold on
                plot(all_cell_dens2,'color',clrs(aa,:),'linewidth',4), hold on
                set(gca,'yscale','log')
                title('Total cell density')
                set(gca,'xtick',-1:size(tspan,2)*2:size(all_cell_dens2,1),'xticklabel',[-1:2:23])
                set(gca,'fontsize',40,'linewidth',4.0)
                xlim([241,size(all_cell_dens2,1)])
                title(''), ylim([10^4 10^8])
            end
        end
        
    end
end

if save_flag
    save(savename)
end



% Comparing the final average parameter values:
figure(1); hold on
for q = 1:length(ause)
    subplot(1,4,1), hold on
    bar(q,mean(squeeze(weighted_u(:,end,q))),'facecolor',clrs(q,:),'linewidth',2.0)
    line([0 10],[u_cnt u_cnt],'color','r','linewidth',3.0)
    errorbar(q,mean(squeeze(weighted_u(:,end,q))),std(squeeze(weighted_u(:,end,q))),'.','linewidth',3.0,'color','k')
    ylabel('u'), xlabel('[A_t_r_e_a_t]'), xlim([0,length(ause)+1])
    set(gca,'xtick',[1:length(ause)],'xticklabel',cellstr(string(ause)),'fontsize',14)
    axis square, box on, set(gca,'linewidth',2.0)

    subplot(1,4,2), hold on
    bar(q,mean(squeeze(weighted_k(:,end,q))),'facecolor',clrs(q,:),'linewidth',2.0)
    errorbar(q,mean(squeeze(weighted_k(:,end,q))),std(squeeze(weighted_k(:,end,q))),'.','linewidth',3.0,'color','k')
    ylabel('k'), xlabel('[A_t_r_e_a_t]'), xlim([0,length(ause)+1])
    line([0 10],[k_cnt k_cnt],'color','r','linewidth',3.0)
    set(gca,'xtick',[1:length(ause)],'xticklabel',cellstr(string(ause)),'fontsize',14)
    axis square, box on, set(gca,'linewidth',2.0)

    subplot(1,4,3), hold on
    bar(q,mean(squeeze(weighted_ks(:,end,q))),'facecolor',clrs(q,:),'linewidth',2.0)
    errorbar(q,mean(squeeze(weighted_ks(:,end,q))),std(squeeze(weighted_ks(:,end,q))),'.','linewidth',3.0,'color','k')
    ylabel('K_s'), xlabel('[A_t_r_e_a_t]'), xlim([0,length(ause)+1])
    line([0 10],[ks_cnt ks_cnt],'color','r','linewidth',3.0)
    set(gca,'xtick',[1:length(ause)],'xticklabel',cellstr(string(ause)),'fontsize',14)
    axis square, box on, set(gca,'linewidth',2.0)

    subplot(1,4,4), hold on
    bar(q,mean(squeeze(weighted_d(:,end,q))),'facecolor',clrs(q,:),'linewidth',2.0)
    errorbar(q,mean(squeeze(weighted_d(:,end,q))),std(squeeze(weighted_d(:,end,q))),'.','linewidth',3.0,'color','k')
    ylabel('d'), xlabel('[A_t_r_e_a_t]'), xlim([0,length(ause)+1])
    line([0 10],[d_cnt d_cnt],'color','r','linewidth',3.0)
    set(gca,'xtick',[1:length(ause)],'xticklabel',cellstr(string(ause)),'fontsize',14)
    axis square, box on, set(gca,'linewidth',2.0)

  


end

