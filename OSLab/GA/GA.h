#ifndef GA_H_INCLUDED
#define GA_H_INCLUDED

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>

using namespace std;

struct GENE {
    vector<int> glist;
    double      fitness;
};

typedef vector<GENE>    GENEList;
typedef vector<double>  DOUBLE1;
typedef vector<DOUBLE1> DOUBLE2;

class GA {
  private:
    DOUBLE2  DATA;
    GENEList GENES, POOL;
    int TYPE, GENES_CNT;
  public:
    GA(string fileName, int type);
    void initial();
    void fitness();
    void select();
    void crossover(double CRATE);
    void mutation(double MRATE);
    void caculateFitness(GENE& gene);
    double caculateSSE(GENE& gene, DOUBLE2& G);
    void readFile(string filename);
    GENE getBest();

};

    GA::GA(string fileName, int type) {
        TYPE = type;
        GENES_CNT = 20;
        readFile(fileName);
    }
    void GA::initial() {
        for(int i=0; i<GENES_CNT; ++i) {
            GENE gene, pool;
            int cnt = DATA.size();
            for(int j=0; j<cnt; ++j) {
                gene.glist.push_back( rand() % TYPE );
                pool.glist.push_back(0);
            }
            GENES.push_back(gene);
            POOL.push_back(pool);
        }
    }

    void GA::fitness() {
        for(int i=0; i<GENES.size(); ++i) {
            caculateFitness(GENES[i]);
        }

        bool change = true;
        while(change){
            change = false;
            for(int i=1; i<GENES.size(); ++i) {
                if(GENES[i].fitness > GENES[i-1].fitness) {
                    swap(GENES[i], GENES[i-1]);
                    change = true;
                }
            }
        }
    }

    void GA::select() {
        int sum = GENES.size()*GENES.size()/2;

        for(int i=0; i<POOL.size(); ++i) {
            // 輪盤挑選
            int number = rand()%sum, index=0;
            for( ; index<GENES.size()-1; ++index){
                number -= index;
                if(number <= 0) break;
            }
            // 放進交配池中
            POOL[i] = GENES[index];
        }
    }

    void GA::crossover(double CRATE) {
        int cnt = GENES.size(),
            len = GENES[0].glist.size();

        for(int i=0; i<cnt; ) {
            // 隨機決定兩條染色體
            int x = rand()%cnt,
                y = rand()%cnt;
            while(x == y) y = rand()%cnt;// 避免重複

            // 決定是否交配
            double rate = rand()%1000;
            if ( rate< CRATE*1000 ) {

                int a = rand()%len,
                    b = rand()%len;
                while(a == b) b = rand()%len; // 避免重複

                if(a > b) swap(a, b); // a < b

                // 交配並放入母體中
                for(int j=0; j<len; ++j) {
                    if( j >= a && j<=b )
                        GENES[i].glist[j] = POOL[y].glist[j];
                    else
                        GENES[i].glist[j] = POOL[x].glist[j];
                }++i;
                if(i>=cnt) break;
                for(int j=0; j<len; ++j) {
                    if( j >= a && j<=b )
                        GENES[i].glist[j] = POOL[x].glist[j];
                    else
                        GENES[i].glist[j] = POOL[y].glist[j];
                }++i;
            } else {
                // 不交配直接放入母體中
                GENES[i] = POOL[x]; ++i;
                if(i>=cnt) break;
                GENES[i] = POOL[y]; ++i;
            }
        }
    }

    void GA::mutation(double MRATE) {
        double rate = rand()%1000;
        int cnt = GENES.size(),
            len = GENES[0].glist.size();

        for(int i=0; i<cnt; ++i){
            // 決定是否發生突變
            
            if ( rate< MRATE*1000 ) {
                int a   = rand()%len;
                GENES[i].glist[a] = (GENES[i].glist[a]+1)%TYPE;
            }
        }
    }

    void GA::caculateFitness(GENE& gene) {
        DOUBLE2 G;
        DOUBLE1 cnt;
        // 初始化重心
        for(int i=0; i<TYPE; ++i) {
            DOUBLE1 point;
            for(int j=0; j<DATA[0].size(); ++j) {
                point.push_back(0);
            }
            G.push_back(point);
            cnt.push_back(0);
        }
        // 計算重心
        for(int i=0; i<gene.glist.size(); ++i){
            int index = gene.glist[i];
            for(int j=0; j<DATA[0].size(); ++j){
                G[index][j] += DATA[i][j];
            }
            ++cnt[index];
        }
        for(int i=0; i<TYPE; ++i){
            for(int j=0; j<DATA[0].size(); ++j){
                G[i][j] /= cnt[i];
            }
        }
        // 計算SSE
        gene.fitness = caculateSSE(gene, G);
    }

    double GA::caculateSSE(GENE& gene, DOUBLE2& G) {
        double sse = 0;
        for(int i=0; i<gene.glist.size(); ++i) {
            int index = gene.glist[i];
            for(int j=0; j<DATA[0].size(); ++j) {
                sse += ( (DATA[i][j]-G[index][j]) * (DATA[i][j]-G[index][j]) );
            }
        }
        return sse;
    }

    void GA::readFile(string filename) {
        string inputStr, token;
        ifstream readFile(filename);
        while(getline(readFile, inputStr)) {
            DOUBLE1 data;
            istringstream splitStr(inputStr);
            bool index = true;
            while(getline(splitStr, token, ',')) {
                data.push_back(stod(token));
            }
            DATA.push_back(data);
        }
        readFile.close();
    }

    GENE GA::getBest(){
        return GENES[GENES_CNT-1];
    }


#endif // GA_H_INCLUDED
