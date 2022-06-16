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

function showResultsMulti(shapes)
    RlastArray = {}
    for i = 1:length(shapes)
        shape = shapes{1,i}{1, 51};
        [RLast, f] = ABH_Optimitzation(shape,'vec');
        RlastArray{i} = RLast;
     
    end
    [Rlinear, f] = ABH_Optimitzation(20:1.5:80,'qua');
    figure
    hold all
    
    for i = 1:length(RlastArray)
        plot(f/1000,abs(RlastArray{i}),'--b');
    end
    plot(f/1000,abs(Rlinear),'LineWidth',1.5);
    
end