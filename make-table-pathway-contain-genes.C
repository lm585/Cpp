/*
input 1: a list of files, each file a pathway

input 2: a list of recurrent genes

op:
genes	path1	path2	path3	path4
g1	1	0	1	0
g2	0	0	0	1
g3	1	1	1	0


0: gene is NOT contained in the pathway LEGs
1: gene is contained in the pathway LEGs

algorithm:
vector<string> pathways, genes;
struct pathGenes {
string path;
map<string, int> LEGs;
vector<int> isGeneInPath;
};
vector<pathGenes> pathList;

read a pathway LEG file,
if LEGs contain the gene of interest
store LEGs in struct.LEGs (map)
if a gene in the "genes" in the map, =1, otherwise = 0

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
#include <cmath>

using namespace std;

struct pathGenes {
string path;
map<string, int> LEGs;
vector<int> isGeneInPath;
void clear();
};

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);

void pathGenes::clear()
{
 path = "";
 LEGs.clear();
 isGeneInPath.clear();
}

int main(int argc, char *argv[])
{
 if(argc != 4)
 {
  cerr << argv[0] << " arg1[ a list of files, each file a pathway]	arg2[a list of recurrent genes]  op-file" << endl;
  cerr << "see the header of " << argv[0] << ".C file for description" << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 vector<string> pathways, genes;
 vector<pathGenes> pathList;
 int i, j;

 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t', line);
  if(lineFields.size() > 0)
    pathways.push_back(lineFields[0]);

  getline(input, line);
 }
 input.close();

 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t', line);
  if(lineFields.size() > 0)
    genes.push_back(lineFields[0]);

  getline(input, line);
 }
 input.close();
 
 writeFile(argv[3], output);

 for(i = 0; i < pathways.size(); i++)
 {
  pathGenes pG;
  pG.clear();
  pG.path = pathways[i];
  readFile(pathways[i].c_str(), input);
  getline(input, line);
  while(!input.eof())
  {
   getFieldContent(lineFields, '\t', line);
   if(lineFields.size() > 0)
     pG.LEGs[lineFields[0] ] = 1;

   getline(input, line);
  }
  input.close();
  //after filling LEGs of a pathway
  for(j = 0; j < genes.size(); j++) //a        b       c       d
  {
   pG.isGeneInPath.push_back(0);
  }
  for(j = 0; j < genes.size(); j++)//a-0       b-0     c-1     d-1
  {
   if(pG.LEGs.count(genes[j]) > 0)
   {
    pG.isGeneInPath[j] = 1;
   }
  }
 pathList.push_back(pG);
 }
 
 output << "gene";
 for(i = 0; i < pathList.size(); i++)
 {
  output << '\t' << pathList[i].path;
 }
 output << endl;

 for(j = 0; j < genes.size(); j++)
 {
  output << genes[j];
  for(i = 0; i < pathList.size(); i++)
  {
   output << '\t' << pathList[i].isGeneInPath[j];
  }
  output << endl;
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

void writeFile(const char * argv, ofstream & output)
{
 output.open(argv, ios::out);
 if(!output)
 {
  cerr << argv << " can not be opened for writing.\n";
  exit(1);
 }
}

void getFieldContent(vector<string> & fields, char del, string line)
{
 vector<int> pos;
 string str;
 int len, size;

 fields.clear();
 for(int i = 0; i < line.length(); i++)
 {
  if(del == line[i])
    pos.push_back(i);
 }
 if(pos.size() > 0)
 {
  len = pos[0];
  str = line.substr(0, len);
  if(len > 0)
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   str = line.substr(pos[i] + 1, len);
   fields.push_back(str);
  }
  size = pos.size();
  if(pos[size-1] < line.length() - 1) //not at the end of line
  {
   str = line.substr(pos[size-1] + 1);
   fields.push_back(str);
  }
 }
 else
 {
  fields.push_back(line);
 }
}

void getFieldContent2(vector<string> & fields, string del, string line)
{
 vector<int> pos;
 string str;
 int len, size;

 fields.clear();
 for(int i = 0; i < line.length(); i++)
 {
  if(del.find(line[i]) != string::npos)
    pos.push_back(i);
 }
 if(pos.size() > 0)
 {
  len = pos[0];
  str = line.substr(0, len);
  if(len > 0)
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   if(len > 0)
   {
    str = line.substr(pos[i] + 1, len);
    fields.push_back(str);
   }
  }
  size = pos.size();
  if(pos[size-1] < line.length() - 1) //not at the end of line
  {
   str = line.substr(pos[size-1] + 1);
   fields.push_back(str);
  }
 }
 else
 {
  fields.push_back(line);
 }
}


