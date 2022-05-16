%
% showResults(experimentBestShapes{1, 1}{1, 51},experimentBestProposals{1, 1}{1, 51})
%
function showResults(shape,exp)
 figure(100);
    
   [Rlinear, f] = ABH_Optimitzation(20:1.5:80,'qua');
    figure(200);
    [RLast, f] = ABH_Optimitzation(shape,'vec');
    figure;
    hold on
    plot(f/1000,abs(RLast));
    plot(f/1000,abs(Rlinear));
    legend(strcat('Agent proposal ', num2str(exp{2})),'Quadratic')
end