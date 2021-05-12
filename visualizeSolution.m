function visualizeSolution(vec)
L=0.5;
hring=0.001/1000;
xl=1.0e-3;
numRings = 40;
dx = (L-xl)/numRings;

proposal = zeros(numRings);
i = 1;
for v = vec
    proposal(i) = (v{1});
    i= i+1;
end

 
[R,f] = ABH_Optimitzation(proposal);
plot(f,R)
figure

plot(f,abs(R))

figure

y = proposal;
x = xl + (0:numRings-1)*dx;
bar(x,y,1)


        


