
import math
import numpy as np
from operator import itemgetter
import cma
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

class Solver(object):
    def __init__(self,rings,population):
        self.population = population
        self.rings = rings
        self.initialMeans = np.random.uniform(low=1, high=9.5, size=(2,))
        self.initialSigma = 5
        #self.es = cma.CMAEvolutionStrategy(self.initialMeans, self.initialSigma, {'popsize': self.population,'CMA_on':0 ,'bounds': [1, 9.5]})
        self.es = cma.CMAEvolutionStrategy(self.initialMeans, self.initialSigma, {'popsize': self.population,'bounds': [1, 9.5]})
        self.candidates = self.es.ask()
        print('Created solver')

    def ask(self):
        return self.es.ask()
    
    def showEvolution(self,target):
        x = []
        y = []
        fig, ax = plt.subplots()
        fig.suptitle(title)
        for candidate in self.candidates:
            ax.plot(candidate, color='tab:blue', linestyle='dashed',linewidth=0.5)
        ax.plot(self.candidates[0], color='tab:green', linestyle='dashed', linewidth=3, marker='o')
        ax.plot(target, color='tab:red', linestyle='dashed', linewidth=2, marker='o')
        ax.set(xlabel='Ring position', ylabel='Radius (R)', title='Quadratic ABH')
        ax.set_xlim(-1, len(target))
        ax.set_ylim((-1, int(target[len(target)-1]+target[len(target)-1]*0.5)))
        fig.show()

    def tell(self,proposals):
        best = []
        fitness = []
        for i in range(0,len(proposals)):
            fitness.append(proposals[i]["fit"])
            best.append(proposals[i]["value"])
        self.es.tell(best,fitness)
        if np.mean(fitness) < 0.05:
            end = True
            plot = True
        else:
            end = False
        return end

def simResults2List(rewardList, candidates):
    solutions = []
    fitNormal = []
    for i in range(len(candidates)):
        element = {
            "value":candidates[i],
            "fit": rewardList[i]
        }
        solutions.append(element)
    newList = sorted(solutions, key=itemgetter('fit'), reverse = False)
    return newList

def evaluateMinimal(candidates, target):
    solutions = []
    fitNormal = []
    for i in range(len(candidates)):
        fit = np.square(np.subtract(target,candidates[i])).mean()  
        element = {
            "value":candidates[i],
            "fit": fit
        }
        fitNormal.append(fit);
        solutions.append(element)
    newList = sorted(solutions, key=itemgetter('fit'), reverse = False)
    return newList

def npArray2MatlabList(npArray):
    proposal = []
    for element in npArray:
        proposal.append(element)
    return proposal

def generateRingsFromnpArray(npArray):
    profile = generateABHProfile(npArray)
    profilePerCent = radiusPerCent(profile,0.5)
    proposal = []    
    for element in profilePerCent:
        proposal.append(element)
    return proposal

def list2MatlabVector(x):
    eng = matlab.engine.start_matlab()
    columnVector = eng.fitdist(matlab.double(list(x), (x.size, 1)), 'stable')
    return columnVector

def generateABHProfile(arr):
    numRings = 40
    arrayRadius = []
    R = 0.23
    L = 0.5
    xl=1.0e-3;
    hring=0.001/1000;
    dx = (L-xl-hring*numRings)/numRings;
    n = 2
    arrayRadiusPositions = []
    pos = -xl
    for j in range(numRings): 
        pos = pos - dx
        arrayRadiusPositions.append(abs(pos))
    
    for i in range(numRings):
        currentRadius = R*(arrayRadiusPositions[i]/L)**n
        arrayRadius.append(currentRadius)
    return arrayRadius

def radiusPerCent(arrayRadius,L):
    return [ i / L *100 for i in arrayRadius]

def linearEquation(parameters):
    solution = []
    L=0.5;
    hring=0.001/1000;
    xl=1.0e-3;
    for i in range(parameters):
        solution.append(2*i)
    return solution

def quadraticEquation(parameters):
    distTruncation = 3
    ringSeparation = 3
    solution = []
    for i in range(parameters):
        solution.append(i*i*ringSeparation+distTruncation)
    return solution