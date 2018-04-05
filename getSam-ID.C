#include <cstdlib>
#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

void readFile(const char * argv, ifstream & input);
void buildMap(string line, map<string, string> & samMap);

int main(int argc, char *argv[])
{
 if(argc != 3)
 {
  cerr << argv[0] << "read-ID-file    sam-file" << endl;
  return 1;
 }

 string line, id;
 map<string, string> samMap;
 ifstream input;

 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  buildMap(line, samMap);
  getline(input, line);
 }
 input.close();

 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  if(samMap.count(line) > 0)
    cout << samMap[line] << endl;
  getline(input, line);
 }
 input.close();
 
 return 0;
}

void buildMap(string line, map<string, string> & samMap)
{
 string id;
 int pos;

 pos = line.find('\t');
 if(pos == string::npos)
   return;
 // pos = 9; 0,1,2,...9; line[9] = '\t'; len = 9 (0,1,..., 8)
 id = line.substr(0, pos );
 if(samMap.count(id) == 0)
   samMap[id] = line;
 else // id already exist; mate read
 {
  string secLine("\n");
  secLine += line;
  samMap[id] += secLine;
 }
 return;
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


