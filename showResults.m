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

function showResults(shape,exp)
 figure(100);
    
   [Rlinear, f] = ABH_Optimitzation(20:1.5:80,'qua');
    figure(200);
    [RLast, f] = ABH_Optimitzation(shape,'vec');
    figure;
    hold on
    plot(f/1000,abs(RLast));
    plot(f/1000,abs(Rlinear));
    %legend(strcat('Agent proposal ', num2str(exp{2})),'Quadratic')
%figure(200);
%L = 0.5
%numRings = 40
%xl=1.0e-3;
%hring=0.001/1000;
%dx = (L-hring*numRings)/numRings;
%x = xl - (0:numRings-1)*dx;
% bar(x,shape,1)

end