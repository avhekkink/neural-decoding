clear;
% QUESTION 1

% Case with 5 stimuli
% Set up parameter values
T = 1000;                                                % 1000 ms = 1 sec
rateStim = [0.01, 0.02, 0.03, 0.04, 0.05];           % Rates in spikes/ ms
numStimPres = [100, 100, 100, 100, 100]; % All stimuli presented 100 times
binSize = 1;                                              % 1 ms time bins
seed = 1;

[trains] = generatetrains(T, binSize, numStimPres, rateStim, seed);


% QUESTION 2

% Set up an array to store no. of spikes in each spike train
numberOfStimuli = length(numStimPres);
r = numberOfStimuli;        % no. of rows = no. of stimuli
c = max(numStimPres);       % no. cols = max no. of stimuli presentations
spikeSums = zeros(r,c);

% Calculate spike count for each spike train
for i = 1:r
    for j = 1:c
        spikeSums(i,j) = sum(trains(i,j,:));
    end
end

% Determine max no. spikes occuring in any train,so find max of whole array
maxN = max(spikeSums, [], 'all');



% QUESTION 3

% Set up array pns to store the probabilities of observing a particular no.
% of spikes, n, given a particular stimulus, s
pns = zeros(numberOfStimuli, (maxN +1));

for i = 1:numberOfStimuli
    for n = 1:(maxN +1)
        % Look at row corresponding to stimulus i
        x = spikeSums(i,:);
        % Find all the values which equal n-1
        % (Note use of n-1 instead of n to account for possibility of zero)
        y = find( x == (n-1));
        % Probability of n = no. times n occured / total no. stim pres 
        p = length(y)/numStimPres(i);
        
        % Record this probability in pns
        pns(i, n) = p;
    end
end

figure(1)
plot(0:maxN, pns, 'LineWidth', 1)
title('Probability distribution for observing a particular number of spikes given a particular stimulus')
lgd = legend({'1','2','3','4','5'});
title(lgd, 'Stimulus')
xlabel('Number of spikes observed')
ylabel('Probability')



% QUESTION 4

% Set up vector for ML Estimate
mlEstimate = zeros(1,(maxN +1));

% ML estimate S_ml is the stimulus that maximises p(n|S_ml)
% So for each column of pns, we find the max value and return the
% corresponding stimulus.

% Note that the (m+1)th element of mlEstimate contains the MLE of the 
% stimulus,  given that m spikes occured.

for n = 1:maxN+1
    % Find max in column n and return index i
    % Index i corresponds to the ML estimate, S_ml
    % (Note that if two stimuli have the same probability, matlab 
    % automatically picks the smaller stimulus, as requested in assignment)
    [M, i] = max(pns(:,n));
    if M == 0
        mlEstimate(n) = NaN;
    else
    mlEstimate(n) = i;
    end
end

% Plot of the MLE for the stimulus, given the observed spike count

% Note: I am choosing to plot the data as points rather than a line because
% I think it displays the information more effectively and makes sure no 
% data is lost due to the fact some inputs are NaN

figure(2)
subplot(1,1,1)
plot(0:maxN, mlEstimate, 'Marker', '.', 'LineStyle', 'none')
title('ML estimate of the stimulus, given the observed spike count')
ylabel('ML estimate of the stimulus')
xlabel('Number of spikes observed')
xlim([0, maxN])
ylim([0.5, 5.5])
yticks([1, 2, 3, 4, 5])
    


% QUESTION 5

% Set up vector for MAP Estimate
mapEstimate = zeros(1,(maxN +1));

% To calculate the MAP Estimate, we need the probs of each stimulus, P(s)
ps = numStimPres./sum(numStimPres);

% Set up array for storing P(s_map|n) that is P(n|s_map)*P(s_map)
% i.e. psn = pns * ps
psn = zeros(numberOfStimuli, (maxN +1));
for i = 1:numberOfStimuli
    for n = 1:(maxN+1)
        psn(i,n) = pns(i,n) * ps(i);
    end
end

% MAP estimate S_map is the stimulus that maximises p(S_map|n)
% So for each column of psn, we find the max value and return the
% corresponding stimulus.

for n = 1:maxN+1
    % Find max in column n and return index i
    % Index i corresponds to the MAP estimate, S_map
    [M, i] = max(psn(:,n));
    if M == 0
        mapEstimate(n) = NaN;
    else
        mapEstimate(n) = i;
    end
end

% Plot of the MAP estimate for the stimulus, given the observed spike count
figure(3)
% Again, plotting just the points
subplot(1,1,1)
plot(0:maxN, mapEstimate, 'Marker', '.', 'LineStyle', 'none')
title('MAP estimate of the stimulus, given the observed spike count')
ylabel('MAP estimate of the stimulus')
xlabel('Number of spikes observed')
xlim([0, maxN])
ylim([0.5, 5.5])
yticks([1, 2, 3, 4, 5])

% Notice that this is the same plot as the one for the ML estimate, since 
% in this case all the stimuli are equally likely



%QUESTION 6

% See separate file for code

% Include plots from this code, and add in an explanation of the
% differences between these plots and those when all stimuli were presented
% the same number of times




% QUESTION 7

% Entropy in spike count response , H , is equal to:
% H = - sum( p(n) * log2( p(n) ), with n = 0:maxN

% To find p(n), we use p(n) = sum( p(n|s) * p(s) ), with s = 1:5

% Set up vector to store p(n)
pn = zeros(1,(maxN + 1));
for n = 1: (maxN+1)
    A = zeros(1, numberOfStimuli);
    for i = 1:5
        A(i) = ps(i) * pns(i, n);
    end
    pn(n) = sum(A);
end

% Set up vector for the values which need to be summed
B = zeros(1,(maxN+1));

for n = 1: (maxN+1)
    % Preventing NaN from log2(0)
    if pn(n) == 0
        B(n) = 0;
    else
        B(n) = pn(n) * log2( pn(n) );
    end
end

% Compute Entropy, H
H = (-1) * sum(B)
% Output:   H = 5.7485 bits

%H_noise = sum( sum( p(s) * p(n|s) * log2(p(n|s)) )), with s=1:5, n=0:maxN
a = zeros(numberOfStimuli, 1);

for i = 1:5
    b = zeros((maxN+1),1);
    for n = 1:maxN
        if pns(i,n) == 0
            b(n) = 0;
        else 
            b(n) = ps(i) * pns(i,n) * log2(pns(i,n));
        end
    end
    a(i) = sum(b);
end

% Compute noise entropy, H_noise
H_noise = -1 * sum(a)
% Output: H_noise = 4.1530


% Calculate mutual information, I = H - H_noise
I = H - H_noise
% Output: I = 1.5955

