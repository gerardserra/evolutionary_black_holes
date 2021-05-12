function [bestProposals,bestFits]=ABH_test(numRings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% call with 2 rings [proposals, fitness] = ABH_test(2);
% call with 40 rings [proposals, fitness] = ABH_test(40);
% Returns the best proposal of each generation and its fitness
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAX_GENERATIONS = 100;
bestProposals = {};
bestFits = {};
endExperiment = false;
totalRings = int32(numRings);

moduleSolver = py.importlib.import_module('solver');
targetSolution = moduleSolver.linearEquation(totalRings);

disp('Solution'); 
disp(targetSolution);
s = moduleSolver.Solver(numRings);
candidates = s.ask();

gen=0;

while not(endExperiment)
 
    proposals = moduleSolver.evaluateMinimal(candidates,targetSolution); 
    bestProposal = proposals(1);
    for n = bestProposal
        maxFit = n{1}{'fit'};
        p = cell(moduleSolver.npArray2MatlabList(n{1}{'value'}));
        bestProposals{end+1} = p;
        bestFits{end+1} = maxFit;
    end

    endExperiment = s.tell(proposals);
    candidates = s.ask();
    
    disp("Generation "+gen)
    fprintf("Fitness " + maxFit +"\n\n")
    if(gen > MAX_GENERATIONS)
        endExperiment = true;
    end
    gen = gen + 1;
 
end

fprintf("\n Final result \n");
disp(bestProposals(end))
disp(bestFits(end))
