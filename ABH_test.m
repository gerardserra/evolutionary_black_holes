
function [bestProposals,bestFits,bestShapes,all,allFitness]=ABH_test(numRings,maxGen,totalPop,cutFreq)

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

MAX_GENERATIONS = maxGen;

POPULATION = totalPop;
%POPULATION = int32(4+3*log(numRings))
disp("Population:"+POPULATION)

bestProposals = {};
bestFits = {};
bestShapes = {};
allFitness = {};
all = {};

endExperiment = false;

moduleSolver = py.importlib.import_module('solver_grow');
                        
s = moduleSolver.Solver(numRings,POPULATION);
candidates = s.ask();
gen = 0;

while not(endExperiment)
    disp("Generation "+gen)
    disp("Freq :"+cutFreq)
    rewardList = zeros(1,POPULATION);
    iProposal = 1;
    bestScore = 100000;
    allCandidates = {};
    allFitness = {};
    for c = candidates
       y = [moduleSolver.getShapeFromArray(c{1})];
       proposal = zeros(size(y));
       i = 1;
       for p = y
           proposal(i) = (p{1});
           i= i+1;   
       end
       
       allCandidates{end+1} = proposal;
       [R,f] = ABH_Optimitzation(proposal,'vec');
       sizeR = length(R);

       weights = stepWeights05(sizeR,75000);
       maxScore = sum(weights);
       disCount = 0;
       for iR = 2:sizeR
           if(f(iR) < 2000 )
               currentCoef = abs(R(iR));
               %disCount = disCount + currentCoef;
               disCount = disCount + currentCoef*weights(iR);
           end
       end
        %disCount = disCount/sizeR*100;
       disCount = disCount/maxScore*100;
       
       if(disCount < 0)
           disCount = 0.01
       end
       if(disCount < bestScore)
            bestPenalty = disCount;
            bestR = R;
            bestShape = proposal;
       end
       allFitness{end+1} = disCount;
       rewardList(iProposal) = disCount;
       iProposal = iProposal+1;
    end
    
	proposals = moduleSolver.simResults2List(rewardList,candidates);
    s.tell(proposals);
    candidates = s.ask();
    bestProposal = proposals(1);
    for n = bestProposal
        maxFit = n{1}{'fit'};
        p = cell(moduleSolver.npArray2MatlabList(n{1}{'value'}));
        bestProposals{end+1} = p;
        bestShapes{end+1} = bestShape;
        bestFits{end+1} = maxFit;
    end
    
    all{end+1} = allCandidates;
    
    fprintf("Fitness " + maxFit +"\n\n");
   
    if(gen > MAX_GENERATIONS)
        endExperiment = true;
    end
    gen = gen + 1;
end











