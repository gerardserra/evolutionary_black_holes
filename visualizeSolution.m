function visualizeSolution(vec)
L=0.5;
hring=0.001/1000;
xl=1.0e-3;
numRings = 40;
dx = (L-xl)/numRings;


 
[R,f] = ABH_Optimitzation(vec);
plot(f,R)
figure

plot(f,abs(R))

figure

y = vec;
x = xl + (0:numRings-1)*dx;
bar(x,y,1)


        


