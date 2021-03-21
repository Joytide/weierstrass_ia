from math import sqrt, cos, pi
from random import randint,uniform

class iv():
    def __init__(self,a,b,c):
        self.a=a
        self.b=b
        self.c=c

    def __str__(self):
        return "("+str(self.a)+","+str(self.b)+","+str(self.c)+")"

    def Weierstrass(self,i):
        sum=0
        for n in range(self.c+1):
            sum+=pow(self.a,n)*cos(pow(self.b,n)*pi*i)
        return sum

    def fitness(self,samples): 
        score=0
        for i,t in samples.items():
            estimation=self.Weierstrass(i)
            score+=abs(estimation-t)
        return score

    def fit(self):
        return False if (self.a<=0 or self.a>=1 or self.b<1 or self.b>20 or self.c<1 or self.c>20) else True




def cross_over(pop_parents):
    #We choose 2 random parents in the most fit elements of the pop
    mum,dad = pop_parents[randint(0, len(pop_parents)-1)], pop_parents[randint(0, len(pop_parents)-1)]
    child=iv((mum.a+dad.a)/2,round((mum.b+dad.b)/2),round((mum.c+dad.c)/2)) #For continuous values (since mutations will only be discretes)
    return child



def mutation(pop, gen, samples):
    while True:
        parent=pop[randint(0, len(pop)-1)]
        child=iv(uniform(0, 1),randint(1, 20),randint(1, 20)) #Change -1%/0%/1%
        if parent.fitness(samples)<0.736:
            child=iv(parent.a+((randint(0,2)-1)*parent.a/(gen**2)),parent.b,parent.c)
        if child.fit: #Check for correct a,b,c
            return child



def parents(pop,n_parents,samples):
    sorted_pop=sorted(pop, key=lambda iv: iv.fitness(samples))
    print("Best parent fitness=",sorted_pop[0].fitness(samples),sorted_pop[0])
    #print(pop_dict,pop_dict.keys())
    #We only take a certain amount of the best parents
    return [sorted_pop[i] for i in range(len(sorted_pop)) if i<n_parents]




def load_samples():
    import os
    #print(os.getcwd())
    with open("./PROJET/temperature_sample.csv","r+") as sample:
        data=[line.strip().split(";") for line in sample.readlines()][1:]
        samples={}
        for line in data:
            samples[float(line[0])]=float(line[1])
    return samples

def main():
    samples=load_samples()
    pop_size=200
    parent_size=pop_size/10
    
    pop=[iv(uniform(0, 1),randint(1, 20),randint(1, 20)) for k in range(pop_size)]
    best_fitness=1
    for gen in range(1,10000):
        print("GEN",gen)

        pop_parents=parents(pop, parent_size, samples)
        childs=pop_parents
        while len(childs)<pop_size:
            childs.append(mutation(pop,gen,samples)) if not randint(0,1) else childs.append(cross_over(pop_parents))
        pop=childs


main()