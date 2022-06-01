
function [bestProposals,bestFits,bestShapes,all,allFitness]=ABH_test_exp(numRings,maxGen,totalPop,cutFreq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% Returns the best proposal of each generation and its fitness
% 
% Examples:
% >> call with 2 rings [proposals, fitness] = ABH_test_exp(2,10,100);
% >> call with 40 rings [proposals, fitness] = ABH_test_exp(40,10,100);
% proposal = (20:1.5:80) is solution linear;
% Author: Gerard serra, 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAX_GENERATIONS = maxGen;
MAX_VALUE = 25000-cutFreq;

POPULATION = totalPop;
%POPULATION = int32(4+3*log(numRings))
disp("Population:"+POPULATION)

bestProposals = {};
bestFits = {};
bestShapes = {};
allFitness = {};
all = {};

endExperiment = false;

moduleSolver = py.importlib.import_module('solver');
                        
s = moduleSolver.Solver(numRings,POPULATION);
candidates = s.ask();
gen = 1;
while not(endExperiment)
    disp("Generation "+gen)
    disp("Freq :"+cutFreq)
    rewardList = zeros(1,POPULATION);
    iProposal = 1;
    bestScore = 100000;
    allCandidates = {};
    allFitness = {};
    for c = candidates
        y = [moduleSolver.generateRingsFromnpArray(c{1})];
        proposal = zeros(size(y));
            i = 1;
        for p = y
            proposal(i) = (p{1});
            i= i+1;
        end
        allCandidates{end+1} = proposal;
       [R,f] = ABH_Optimitzation(proposal,'vec');
       disCount = 0;
 
       for iR = 2:length(R)
           if(f(iR) < cutFreq)
               if(abs(R(iR))>0.5)
                    disCount = disCount+1;
               end
           end
       end

       disCount = disCount/MAX_VALUE*100;
       
       if(disCount < 0)
           disCount = 0.1
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
    
    disp(bestProposal(1))
    fprintf("Fitness " + maxFit +"\n\n");
   
    if(gen > MAX_GENERATIONS)
        endExperiment = true;
    end
    gen = gen + 1;
  end











