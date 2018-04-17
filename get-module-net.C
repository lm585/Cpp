/*
1 4
2 3
3 2
3 4
4 1
4 3


8 9
8 11
9 8
10  11
11 8
11 10

given the above network (node in int values)
clustered into the 2 components

1 2 3 4 
8 9  10  11

3(maximal degree)  2(degree)   4(module size)  1,2,3,4
8(maximal degree)  2(degree)   4(module size)  8,9,10,11

class link {int, int}
map<int, bool> (node id, used)
vector<link>
map<int, set<int>>  m[1] = {4}
m[2] = {3}
m[3] = {2,4}
m[4] = {1, 3}

populate edge list, 
for each node, build it's neighbors, and node-list map (used = false), 

for each node in the map
  if node not used
    call function getModule
    find elem in the module with the largest degree
    print node degree module-size module-elem


function getModule  //seed node (int), node map, neighbor map, return a set<int> (module)
set<int> module_copy
    module clear
    add node to the module set
    do { module set old size
         module_copy = module
         for each elem in the module_copy
           add elem's neighbors into the module
           elem used = true in the node map
         module's new size
       }  
    while(old size != new size)
 */

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>
#include <sstream>
#include <map>

using namespace std;

class netLink
{
 public:
 int a, b;
};

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void maxDegreeNode(int & node, int & degree, map<int, set<int> > & neighbMap, const set<int> &  module);
void getModule(int seed,  map<int, bool> & nodeMap, map<int, set<int> > & neighbMap, set<int> &  module);

int main(int argc, char *argv[])
{
 if(argc != 3)
 {
  cerr << argv[0] << " input-file(network) output-file(module) " << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 int edge1, edge2;
 netLink  myLink;
 vector<netLink> linkVect;
 map<int, set<int> > neighbMap;
 map<int, bool> nodeMap;
 map<int, bool>::iterator nodeIter;
 set<int> module;

 readFile(argv[1], input);
 input >> edge1 >> edge2;
 while(!input.eof())
 {
  myLink.a = edge1;
  myLink.b = edge2;
  linkVect.push_back(myLink);
  input >> edge1 >> edge2;
 }
 input.close();
 writeFile(argv[2], output);
 cout << "Reading network file " << argv[1] << " ...done" << endl;
 for(int i = 0; i < linkVect.size(); i++)
 {
  int m = linkVect[i].a;
  int n = linkVect[i].b;
  nodeMap[m] = false;
  nodeMap[n] = false;
  neighbMap[m].insert(n);
 }
 cout << "# of nodes " << nodeMap.size() << endl;
 for(int i = 0; i < linkVect.size() && i < 6; i++)
 {
  int m = linkVect[i].a;
  int n = linkVect[i].b;
  cout << m << " has " << neighbMap[m].size() << " neighbors " << endl;
  cout << n << " has " << neighbMap[n].size() << " neighbors " << endl;
 }
 for(nodeIter = nodeMap.begin(); nodeIter != nodeMap.end(); nodeIter++)
 {
  if(nodeIter->second == false)
  {
   int node, degree;
   getModule(nodeIter->first, nodeMap, neighbMap, module);
   maxDegreeNode(node, degree, neighbMap, module);
   output << node << '\t' << degree << '\t' << module.size() << '\t';
   for(set<int>::const_iterator it = module.begin(); it != module.end(); it++)
   {
    output << *it << ',';
   }
   output << endl;
  }
 }
 
 output.close();
 return 0;
}

void readFile(const char * argv, ifstream & input)
{
 input.open(argv, ios::in);
 if(!input)
 {
  cerr << argv << " can not be opened for reading.\n";
  exit(1);
 }
}

void maxDegreeNode(int & node, int & degree, map<int, set<int> > & neighbMap, const set<int> &  module)
{
 int d, n;
 set<int>::const_iterator it = module.begin();
 node = *it;
 degree = neighbMap[node].size();
 for(it = module.begin(); it != module.end(); it++)
 {
  n = *it;
  d = neighbMap[n].size(); 
  if(d > degree)
  {
   degree = d;
   node = n;
  }
 }
}

/*
set<int> module_copy
    module clear
    add node to the module set
    do { module set old size
         module_copy = module
         for each elem in the module_copy
           add elem's neighbors into the module
           elem used = true in the node map
         module's new size
       }
    while(old size != new size)
*/ 
void getModule(int seed,  map<int, bool> & nodeMap, map<int, set<int> > & neighbMap, set<int> &  module)
{
 set<int> copy;
 int oldSize, newSize;
 module.clear();
 module.insert(seed);
 do
 {
  oldSize = module.size();
  copy = module;
  for(set<int>::iterator it = copy.begin(); it != copy.end(); it++)
  {
   int node = *it;
   for(set<int>::const_iterator ptr = neighbMap[node].begin(); ptr != neighbMap[node].end(); ptr++)
   {
    module.insert(*ptr);
   }
   nodeMap[node] = true;
  }
  newSize = module.size();
 }
 while(oldSize != newSize);
}

void writeFile(const char * argv, ofstream & output)
{
 output.open(argv, ios::out); 
 if(!output)
 {
  cerr << argv << " can not be opened for writing.\n";
  exit(1);
 }
}


