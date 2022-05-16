

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% Runs a set of experiments
% 
% Author: Gerard serra, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_exp(total,maxGen,totalPop,cutFreq)

currentExp = 1;
endExperiment = false;

experimentBestProposals = {};
experimentBestFits = {};
experimentBestShapes = {}
experimentCandidates = {}
experimentFitness = {}

while not(endExperiment)
    [experimentBestProposals{end+1},experimentBestFits{end+1},experimentBestShapes{end+1},experimentCandidates{end+1},experimentFitness{end+1}]= ABH_test_exp(40,maxGen,totalPop,cutFreq);
    currentExp = currentExp + 1;
       
    if(currentExp>total)
        endExperiment = true;
    end 
   
end 
 save(datestr(now))