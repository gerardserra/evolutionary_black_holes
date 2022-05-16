

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% Runs a set of experiments
% 
%   Example: 1 exp ,10 generations, 20 members, freq < 500
%    - test_exp(1,20,10,500) 
% 
%   Author: Gerard serra, 2022
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