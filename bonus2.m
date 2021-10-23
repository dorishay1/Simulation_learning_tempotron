    % TODO: Student1, Name and ID
    % TODO: Student2, Name and ID

    % Computation and cognition undergrad - ex3 - BONUS
    % See pdf document for instructions

function bonus2(t, sample, sample_color)
%PLOT_INPUT_SPIKES Generate a raster plot of the input spikes for a given
%sample
%   t               - a times vector
%   sample          - the input neurons' spike times
%   sample_color    - the plot's color

N       = length(sample);   % number of input neurons
t_final = t(end);

% Get the data for the raster plots
inp = cellfun(@(a, b) [a; b.*ones(size(a))], ...
    sample, num2cell(1:N), ...
    'UniformOutput', false);
inp = cell2mat(inp);
inp_x = inp(1,:);   % spike times
inp_y = inp(2,:);   % neuron's id


% TODO: BONUS
%       Replace the following scatter command with a SINGLE plot command,
%       so that the resulting plot will be displayed in the same way as the
%       example in the PDF file. 
    

plot(inp_x*1000,inp_y,".",...
    'LineWidth',2,...
    'MarkerSize',5,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
    

xlim([0, t_final].*1000);
ylim([1 - 0.8, N + 0.8]);
xlabel('time [ms]');
ylabel('Neuron No.');
end

