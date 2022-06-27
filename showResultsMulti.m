%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Shows R from a given profile comparing it to profile of N=2
%
%   From an optimization process
%       load('test_exp(1,20,10,500).mat')
%       showResults(experimentBestShapes{1, 1}{1, 51},experimentBestProposals{1, 1}{1, 51})
% 
%   From given profiles
%       load('profiles.mat')
%       showResults(profileN1,{0,1}) 
%       showResults(profileN2,{0,2}) 
%
%   Author: Gerard serra, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function showResultsMulti(shapes,lastExperiment)

    RlastArray = {}
    for i = 1:length(shapes)
        shape = shapes{1,i}{1, lastExperiment};
        [RLast, f] = ABH_Optimitzation(shape,'vec');
        RlastArray{i} = RLast;
     
    end
    [Rlinear, f] = ABH_Optimitzation(20:1.5:80,'qua');
     sizeR = length(Rlinear);
     weights = stepWeights05(sizeR,75000);
       maxScore = sum(weights);
       disCount = 0;
       for iR = 2:sizeR
           if(f(iR) < 2000 )
               currentCoef = abs(Rlinear(iR));
               disCount = disCount + currentCoef;
               %disCount = disCount + currentCoef*weights(iR);
           end
       end

      % disCount = disCount/maxScore*100;
       disCount = disCount/sizeR*100;
       
    

    disp(disCount)

    figure
    hold all
    
    for i = 1:length(RlastArray)
        if(RlastArray{1,i}(70000) < 0.2)
            plot(f/1000,abs(RlastArray{i}));
        end
    end
    plot(f/1000,abs(Rlinear),'LineWidth',1.5);
    
end