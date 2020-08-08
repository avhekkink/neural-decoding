function [trains] = generatetrains(T, binSize, numStimPres, rateStim, seed)
  % [trains] = generatetrains(T, binSize, numStimPres, rateStim, seed)
  %
  % Returns a 3 dimensional array of poisson spiketrains.  The first
  % dimension is of length length(numStimPres) and indicates the stimulus
  % applied to produce the spike train.  The second dimension is of length
  % max(numStimPres) and indicates which presentation of the given stimulus
  % produced the spike train.  The third dimension of this array is of
  % length floor(T/binSize) and represents elements of a single spike train.
  % Note that each stimulus may be presented a different number of times,
  % and thus some elements of the array may not contain data --- these
  % elements are marked 'NaN'.
  %
  % INPUT VARIABLES:
  %   - T           : total time of each stimulus presentation
  %   - binSize     : length of each time bin
  %   - numStimPres : a vector indicating the number of times each
  %                   stimulus is to be presented.
  %   - rateStim    : a vector indicating the average firing rate produced
  %                   by each stimulus.
  %   - seed        : the random number seed
    
  s = RandStream('mt19937ar', 'Seed', seed);
  
  numberOfStimuli = length(numStimPres);
  maxNumStimPres = max(numStimPres);
  numBins = floor(T/binSize);
    
  trains = NaN(numberOfStimuli, maxNumStimPres, numBins);
    
  for iStim = 1:numberOfStimuli
     trains(iStim,1:numStimPres(iStim),:) = rand(s,1,numStimPres(iStim),numBins) ...
       < rateStim(iStim)*binSize;
  end
  
end