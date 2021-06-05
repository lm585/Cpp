/*
BCL3	1.214417084	-13.9281(log10-P)
GZMB	1.863512174	-12.3335
SOCS3	1.396471393	-12.0778


vector<string> accessionVect;      //input file containing 126 accession names
construct map<string, bool> acc;
output header
for each line in conserv.txt file  //representative conserv.txt file
  parse line into fields
  call getAccPresent()
  output selected fields
  for each elem in accessionVector //output file
    if(acc[elem] == true)
      print "\t" 1
    else
      print "\t" 0


getAccPresent( map<string, bool> & acc, string col1, string field)
{
 set each member of acc to be false
 if(acc.count(col1) == 0)  //a new accession not in the input file
   cerr << new accession 
   exit(1)
 acc[col1] = 'true';
 parse field into a list
 for each member m in the list
   if(acc.count(col1) == 0)
     cerr << new accession
     exit(1)
   acc[m] = 'true'
}
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

struct twoInt {
 double sign;
 double p;
 };


void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);

int main(int argc, char *argv[])
{
 if(argc != 3)
 {
  cerr << argv[0] << "input   output" << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 map<string, twoInt> rnkMap;

 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t', line);
  string gene = lineFields[0];
  if(rnkMap.count(gene) == 0)
  {
   rnkMap[gene].sign = atof(lineFields[1].c_str());
   rnkMap[gene].p = atof(lineFields[2].c_str());
  }
  else
  {
   double d1, d2;
   d1 = atof(lineFields[1].c_str());
   d2 = atof(lineFields[2].c_str());
   if(d2 < rnkMap[gene].p)
   {
    rnkMap[gene].p = d2;
    rnkMap[gene].sign = d1;
   }
  }
  getline(input, line);
 }
 input.close();
 
 writeFile(argv[2], output);
 for(map<string, twoInt>::iterator it=rnkMap.begin(); it!=rnkMap.end(); ++it)
 {
  double t;
  if((it->second).sign < 0)
    t = 1;
  else
    t = -1;
  output << it->first << "\t" << (it->second).p * t << '\n';
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


