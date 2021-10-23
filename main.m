    % Simulation of a learning tempotron

   

    % Computation and cognition undergrad - ex3
    % See PDF document for instructions

    clear; close all; clc;

%% Declare simulation parameters
gL      = 1e-4;     % conductance (S/cm^2), 1/R
C       = 1e-6;     % capacitance (F/cm^2)
tau_m   = C/gL;     % sec - time constant, equal to RC
tau_s   = tau_m/4;  % sec - time constant, equal to RC
TH      = 30*1E-3;  % firing threshold, in mV

dt      = 0.0001;       % time step for numerical integration
t_final = 0.5;          % sec, duration of numerical simulation
t       = 0:dt:t_final; % time vector

%% Define the (normalized) kernel function
alpha   = tau_m/tau_s;
kappa   = (alpha^(-1/(alpha - 1)) - alpha^(-alpha/(alpha - 1)))^(-1);
K       = @(x) (x > 0).*kappa.*(exp(-x/tau_m) - exp(-x/tau_s));

%% Declare the learning parameters
% TODO 1: Change the learning rate, as instructed in the PDF document. Only
%         try this after the learning already works. 
eta      	= 1e-3;	% Learning rate
max_steps   = 1000;	% Maximal number of learning steps

%% Load inputs
% TODO 2: After you implement all of the missing code, load the inputs (one
%         by one) and answer the questions in the PDF file. 
%load X_2SDIW;   % X
%load X_2SDIF;   % X
load X_2PGN;    % X
%load X_2PVTG;   % X
%load X_PRND;    % X

%% Deduce the number of input neurons from the input
N   = length(X{1}.x);          % number of input neurons

%% Initialize the tempotron's weights
W   = TH.*0.5.*(rand(N, 1));    % tempotron's weights
W0  = W;                        % initial weights (for future reference)

%% Learning loop
for learning_step = 1:max_steps
    
    % Choose a random sample
    sample = X{randi(length(X))};
    
    % Integrate & Fire neuron's simulation
    V	= zeros(size(t));	% the tempotron's voltage
    spk	= 0;                % spike occurance flag
    
    % Spikes' data per neuron
    % TODO 3: Describe the data structure created by the following command,
    %         and the data it contains (if you do not understand the
    %         fields' meanings that is OK, a technical description is 
    %         enough).
    neurs = cellfun(@(a, b) struct('inp', a, ...              % a copy of the input
                                   'idx', 1, ...              % next spike's index
                                   'N_spk', length(a), ...    % total number of spikes
                                   'id', b), ...              % neuron's id
                    sample.x, num2cell(1:N));
    
    % Ignore input neurons which never spike
    k = 1;
    while k <= length(neurs)
        if isempty(neurs(k).inp)
            neurs(k) = [];
        else
            k = k + 1;
        end
    end
    
    % While there are any non-processed spikes
    while ~isempty(neurs)
        
        % Get the next spike for each input neuron
        next_spks = arrayfun(@(a) a.inp(a.idx), neurs);
        
        % Select the next spike
        [t_i_j, k]	= min(next_spks);	% next spike's time & neuron
        i           = neurs(k).id;
        
        % Handle the input left for processing
        neurs(k).idx = neurs(k).idx + 1;	% update the relevant spike index
        if neurs(k).idx > neurs(k).N_spk
            neurs(k) = [];
        end
        
        % TODO 4: Describe the order in which the input spikes are reaching
        %         this point (using the debugger), and explain why this 
        %         particular order is important. 
        %Done
        
        % Update the voltage
        V = V + W(i).*K(t - t_i_j);
        
        % Check if the voltage threshold has been crossed
        % TODO 5: Numerically find the maximum voltage and the time at 
        %         which said maximum is achieved, and store the results in 
        %         the variables `V_max` and `t_max` (respectively).
        %Done
        V_max = max(V);
        t_max = t(V == max(V));
        
        if V_max > TH
            spk = 1;	% set the spike flag
            break;      % shunting
        end
        
    end
    
    % If there is no need to learn anything, skip the learning
    % TODO 6: Check if there is no need to update the synaptic weights
    %         (according to the tempotron's learning rules), and SKIP the 
    %         current sample if necessary. 
    
    if spk ~= sample.y0
    
        % Get the gradient of the weights
        W_grad = zeros(size(W));
        for i = 1:N                         % loop over neurons
            for j = 1:length(sample.x{i}) % loop over spikes
                t_i_j = sample.x{i}(j);% next spike time
                if t_i_j < t_max    % get only inputs prior to t_max
                    % TODO 7: Use the derivation from class to update the
                    %         gradient iteratively. You may use the
                    %         precalculated `t_max`, `t_i_j` and the kernel
                    %         function `K`. 
                    
                    W_grad(i) = W_grad(i) + K(t_max-t_i_j); % update the gradient
                end
            end
        end

        % Update the weights
        % TODO 8: Use the learning rate `eta`, teacher's signal `sample.y0` and 
        %         the gradient `W_grad` to update the synaptic weights `W`. 
        if sample.y0 == 1
        W   = W + eta*W_grad;
        else, W   = W - eta*W_grad;
        end
        

    end
   
end
%% Plots

% Set the subplots grid
n_plots	= 3*length(X);
n_rows 	= 3*ceil(length(X)/4);
n_cols 	= ceil(n_plots/n_rows);

% Create the figure
figure('Name', 'Tempotron''s results', ...
    'Units', 'normalized', ...
    'Position', [0, 0, 1, 1]);

% For each sample, draw 3 different plots
for n_sample = 1:length(X)
    
    % Define the sample's plot color according to the teacher's signal
    if X{n_sample}.y0
        sample_color = 'r';
    else
        sample_color = 'b';
    end
    
    % Input's raster plot
    subplot(n_rows, n_cols, 3*n_cols*floor((n_sample-1)/n_cols)+rem((n_sample-1),n_cols)+1);
    bonus2(t, X{n_sample}.x, sample_color);
    title(['Sample #' num2str(n_sample) ': Input']);
    
    % Tempotron's response - before learning
    subplot(n_rows, n_cols, 3*n_cols*floor((n_sample-1)/n_cols)+rem((n_sample-1),n_cols)+1+n_cols);
    plot_postsynaptic_voltage(t, X{n_sample}.x, K, W0, TH, tau_m, sample_color);
    title(['Sample #' num2str(n_sample) ': Voltage before learning']);
    
    % Tempotron's response - after learning
    subplot(n_rows, n_cols, 3*n_cols*floor((n_sample-1)/n_cols)+rem((n_sample-1),n_cols)+1+2*n_cols);
    plot_postsynaptic_voltage(t, X{n_sample}.x, K, W, TH, tau_m, sample_color);
    title(['Sample #' num2str(n_sample) ': Voltage after learning']);
    
end
