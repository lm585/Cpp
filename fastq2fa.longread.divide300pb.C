/*

100     Chr1    21089763        21089921        Chr1    21090056        21090218        3       F,      R,      INSERTION       380.333333333333        UNBAL   478    13       10      101,204,3,34,41,48,89,96,97,G13,

output 
Chr1 21089921 21090056  INSERTION     0	1	0 0 0 0 0 1 1 1 0 1 0 ....

126 accessions, if presented 1; otherwise 0;


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

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);

int main(int argc, char *argv[])
{
 if(argc != 4)
 {
  cerr << argv[0] << "	fastq(file)  outputFile   300(long read -> 300bp each)" << endl;
  
  return 1;
 }

 ifstream input;
 ofstream output;
 string line, readName;
 vector<string> lineFields;
 int len;

 len=atoi(argv[3]);
 writeFile(argv[2], output);

 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  string seq, str;
  int i, k,m;

  getFieldContent2(lineFields," \t\n", line);
  readName=lineFields[0];
  getline(input, line);
  seq=line;
  k=seq.size() / len; //read 100bp, 0block; 1000bp, 3 block; 900bp, 3 block
  if(k >= 1)
  {
   for(i = 1; i <= k; i++)
   {
    str= readName + "-" + std::to_string(i);
    output << ">" << str << endl;
    m=(i-1) * len; //0-299, 300-599, 600-899, 900
    output << seq.substr(m,len) << endl;
   }
  } 
  getline(input, line); //3rd line
  getline(input, line); //4th line

  getline(input, line); //5th line
 }
 input.close();
 
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
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   str = line.substr(pos[i] + 1, len);
   fields.push_back(str);
  }
  size = pos.size();
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


