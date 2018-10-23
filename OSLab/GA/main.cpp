#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "GA.h"

using namespace std;

typedef vector<GENE>    GENEList;
typedef vector<double>  DOUBLE1;
typedef vector<DOUBLE1> DOUBLE2;

int main()
{
    srand((unsigned)time(NULL));

    GA gene("dataset.txt", 3);
    GENE BEST;

    gene.initial();
    gene.fitness();
    BEST = gene.getBest();
    for(int i=0; i<10000; ++i) {
        
        gene.select();
        gene.crossover(0.95);
        gene.mutation(0.20);
        gene.fitness();

        if(gene.getBest().fitness < BEST.fitness){
            BEST = gene.getBest();
            cout<< BEST.fitness <<endl;
        }
    }

    return 0;
}
