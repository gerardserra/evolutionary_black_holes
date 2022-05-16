
function [bestProposals,bestFits,s]=ABH_test(numRings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% Returns the best proposal of each generation and its fitness
% 
% Examples:
% >> call with 2 rings [proposals, fitness] = ABH_test(2);
% >> call with 40 rings [proposals, fitness] = ABH_test(40);
% proposal = (20:1.5:80) is solution linear;
% Author: Gerard serra, 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAX_GENERATIONS = 1;
MAX_VALUE = 2500;

MAX_VALUE = 2720; 

POPULATION = 1;
%POPULATION = int32(4+3*log(numRings))
disp("Population:"+POPULATION)
L=0.5;
hring=0.001/1000;
xl=1.0e-3;
dx = (L-xl)/numRings;

bestProposals = {};
bestFits = {};
endExperiment = false;
totalRings = int32(numRings);

moduleSolver = py.importlib.import_module('solver');
                        
s = moduleSolver.Solver(numRings,POPULATION);
candidates = s.ask();
gen = 0;
    while not(endExperiment)
    disp("Generation "+gen)
    rewardList = zeros(1,POPULATION);
    iProposal = 1;
    for c = candidates
        y = [moduleSolver.npArray2MatlabList(c{1})];
        proposal = zeros(size(y));
            i = 1;
        for p = y
            proposal(i) = (p{1});
            i= i+1;
        end
       [R,f] = ABH_Optimitzation(proposal);
      
       reward = 0;
       for iR = 1:length(R)
           %if( abs(R(iR)) > 0.2)
           if isnan(R(iR))
           else
            reward = reward + abs(R(iR))*(iR/MAX_VALUE);
           end
           %end
       end

       %reward = sum(abs(R(FREQ_CUT:end)).^2);
       %if(reward < 1200)
       %if(reward < MAX_VALUE*10*0.5)
       %     plot(f,abs(R));
       %     disp(proposal);
       %     figure;
       %end
       disp("c:"+iProposal+" r:"+reward);
       rewardList(iProposal) = reward;
       iProposal = iProposal+1;
    end
    
	proposals = moduleSolver.simResults2List(rewardList,candidates);
    disp(proposals(1))
     disp(proposals(2))
    s.tell(proposals);
    candidates = s.ask();
    bestProposal = proposals(1);
    for n = bestProposal
        maxFit = n{1}{'fit'};
        p = cell(moduleSolver.npArray2MatlabList(n{1}{'value'}));
        bestProposals{end+1} = p;
        bestFits{end+1} = maxFit;
    end
    fprintf("Fitness " + maxFit +"\n\n");
    if(gen > MAX_GENERATIONS)
        endExperiment = true;
    end
    gen = gen + 1;
end











