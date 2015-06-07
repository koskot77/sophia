#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <list>
using namespace std;

// g++ -g -Wall -o cluster cluster.cc

class UnionFind {
private:
    map<int,int>         node2cluster;
    map<int, list<int> > cluster2nodes;

public:
    int joinClusters(int cluster1, int cluster2){
        if( cluster1 == cluster2 ) return cluster1;
        list<int> &nodes1 = cluster2nodes[cluster1];
        list<int> &nodes2 = cluster2nodes[cluster2];
        if( nodes1.size()==0 || nodes2.size()==0 ) return 0;
        int newCluster = 0;
        if( nodes1.size() < nodes2.size() ){
            newCluster = cluster2;
            for(list<int>::const_iterator n = nodes1.begin(); n != nodes1.end(); n++)
                node2cluster[*n] = newCluster;
            nodes2.insert(nodes2.end(),nodes1.begin(),nodes1.end());
            cluster2nodes.erase(cluster1);
        } else {
            newCluster = cluster1;
            for(list<int>::const_iterator n = nodes2.begin(); n != nodes2.end(); n++)
                node2cluster[*n] = newCluster;
            nodes1.insert(nodes1.end(),nodes2.begin(),nodes2.end());
            cluster2nodes.erase(cluster2);
        }
        return 0;
    }

    int findCluster(int node) { return node2cluster[node]; }

    int nClusters(void) const { return cluster2nodes.size(); }

    const map<int, list<int> >& clusters(void) const { return cluster2nodes; }

    UnionFind(int maxNodes){
        for(int i=1; i<=maxNodes; i++){
             node2cluster [i] = i;
             cluster2nodes[i].push_back(i);
        }
    }
};

int main(int argc, char *argv[]){
    unsigned int nReads = 84205;

    int cutoffDistance = atoi(argv[2]);

    FILE *input;
    if( (input=fopen(argv[1],"rt")) == NULL ){
        printf("Cannot open %s\n",argv[1]);
        return 0;
    }

    multimap<int, pair<int,int> > edges;
    map<int, map<int,int> > edges2;

    for(int origin=0, dest=0, distance=0; !feof(input) && fscanf(input,"%d,%d,%d\n",&origin,&dest,&distance)==3;){
        if( distance > cutoffDistance ) continue;
        edges.insert( pair< int, pair<int,int> >(distance, pair<int,int>(origin,dest)) );
        edges2[origin][dest] = distance;
        edges2[dest][origin] = distance;
    }
    if( edges.size() == 0 ) return 0;

    UnionFind uf(nReads);

    size_t count = 0;
    for(multimap<int, pair<int,int> >::const_iterator edge=edges.begin(); edge!=edges.end(); edge++,count++){
        if( edge->first > cutoffDistance ) break;
        int node1 = edge->second.first;
        int node2 = edge->second.second;
        int cluster1 = uf.findCluster(node1);
        int cluster2 = uf.findCluster(node2);
        if( cluster1 != cluster2 ) uf.joinClusters(cluster1,cluster2);
    }

    cout<<"Found "<<uf.nClusters()<<" clusters"<<endl;

    ofstream output("clusters.csv");
    if( !output ){ cout<<"Cannot open clusters.csv"<<endl; return 0; }
    output<<"read1,read2,distance"<<endl;
    const map<int, list<int> > &clusters = uf.clusters();
    for(map<int, list<int> >::const_iterator iter = clusters.begin(); iter != clusters.end(); iter++){
        int seed = iter->first;
        for(list<int>::const_iterator node = iter->second.begin(); node != iter->second.end(); node++){
            output<<seed<<","<<*node<<","<<edges2[seed][*node]<<endl;
        }
    }
    output.close();

    return 0;
}
