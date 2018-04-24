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
void getAccPresent( map<string, bool> & acc, string col1, string columnLast);

int main(int argc, char *argv[])
{
 if(argc != 4)
 {
  cerr << argv[0] << " input1(file containing 126 accession names) input2(representative conserv.txt file) output " << endl;
  cerr << "this version: keep left and right coord for both anchoring windows. " << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields, accessVect;
 map<string, bool> accessMap;

 readFile(argv[1], input);
 input >> line;
 while(!input.eof())
 {
  accessVect.push_back(line);
  input >> line;
 }
 input.close();
 cout << accessVect.size() << " accessions: " << endl;
 for(int i = 0; i < accessVect.size(); i++)
 {
  cout << accessVect[i] << ',';
  accessMap[accessVect[i]] = false;
 }
 cout  << endl;

 char t = '\t';
 writeFile(argv[3], output);
 readFile(argv[2], input);
 getline(input, line);
 output << "left window chr" << t << "start" << t << "end" << t;
 output << "right window chr" << t << "start" << t << "end" << t;
 output << "# abn pairs" << t;
 output << "left window strand" << t;
 output << "right window strand" << t;
 output << "sv type" << t;
 output << "avg abn pair insert size" << t;
 output << "balance" << t;
 output << "avg norm pair insert size" << t;
 output << "norm pair insert size std" ;
 for(int i = 0; i < accessVect.size(); i++)
   output << t << accessVect[i];
 output << endl;
 while(!input.eof())
 {
  int left = 0, right = 0;
  getFieldContent(lineFields, '\t', line);
  getAccPresent(accessMap, lineFields[0], lineFields.back());
  if(lineFields[10] ==  "DELETION" || lineFields[10] ==  "INSERTION")
  {
   left = atoi(lineFields[3].c_str());
   right = atoi(lineFields[5].c_str());
  }
  else if(lineFields[10] ==  "INVERSION" || lineFields[10] ==  "LARGE_DUPLI")
  {
   left = atoi(lineFields[2].c_str());
   right = atoi(lineFields[6].c_str());
  }
  for(int i = 1; i <= 13; i++)
    output << lineFields[i] << t;
  output <<  lineFields[14]; 
  for(int i = 0; i < accessVect.size(); i++)
  {
   if(accessMap[accessVect[i]] == false)
     output << t << 0;
   else
     output << t << 1;
  }
  output << endl;
  getline(input, line);
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

void getAccPresent( map<string, bool> & acc, string col1, string columnLast)
{
 for(map<string, bool>::iterator it = acc.begin(); it != acc.end(); it++)
 {
  it->second = false;
 }
 if(acc.count(col1) == 0)
 {
  cerr << col1 << " not in the accession file. exit... " << endl;
  exit(1);
 }
 acc[col1] = true;
 
 vector<string> accList;
 getFieldContent(accList, ',', columnLast);
 for(int i = 0; i < accList.size(); i++)
 {
  if(acc.count(accList[i]) == 0)
  {
   cerr << accList[i]  << " not in the accession file. exit... " << endl;
   exit(1);
  }
  acc[ accList[i] ] = true;
 }
 return;
}
