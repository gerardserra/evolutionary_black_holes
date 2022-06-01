
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
        self.initialMeans = np.random.uniform(low=1, high=1.25, size=(self.rings,))
        self.initialSigma = 0.05
        #self.es = cma.CMAEvolutionStrategy(self.initialMeans, self.initialSigma, {'popsize': self.population,'CMA_on':0 ,'bounds': [1, 9.5]})
        self.es = cma.CMAEvolutionStrategy(self.initialMeans, self.initialSigma, {'popsize': self.population,'bounds': [1, 1.25]})
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

def npArray2MatlabList(npArray):
    proposal = []
    for element in npArray:
        proposal.append(element)
    return proposal


def getShapeFromArray(npArray):
    proposal = []
    previousElement = 1
    for element in npArray:
        currentElement = element*previousElement
        if(currentElement > 99):
            currentElement = 99
        if(currentElement < 0.05):
            currentElement = 0.05
        proposal.append(currentElement)
        previousElement = currentElement
    return proposal
 
    
