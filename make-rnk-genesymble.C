/*
with one header line
geneID	logFC	logCPM	F	PValue	FDR	geneID	name
ENSG00000130844	1.71560015558198	6.80485947158675	23.1667928311584	7.14598479636303e-05	0.260209359869216	ENSG00000130844	ZNF331
ENSG00000047634

field 2 logfc, sign
field 5, p vale
filed 8, gene symbol
get signed log10 p-value
if logfc > 0, sign +
if logfc < 0, sign -
if logfc == 0, signed-log10-p = 0
if p == 0, signed-log10-p=notanumber
if one gene symbol has multiple gene IDs, signed-log10-p will be avg of all signed-log10-p

map<string, vector<double> > symbolValues
vector<string> symbolInorder 
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

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);

int main(int argc, char *argv[])
{
 if(argc != 3)
 {
  cerr << argv[0] << " diff-DEG-op-by-edgeR	op-file" << endl;
  cerr << "see the header of " << argv[0] << ".C file for description" << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 map<string, vector<double> > symbolValues;
 vector<string> symbolInorder;

 readFile(argv[1], input);
 getline(input, line); //skip 1st header line
 getline(input, line);
 while(!input.eof())
 {
  double signedlog10p;
  getFieldContent(lineFields, '\t', line);
  double sign = atof(lineFields[1].c_str());
  double p = atof(lineFields[4].c_str());
  if(sign > 0)
  {
   signedlog10p = -log10(p);
   if(signedlog10p == 0)
     signedlog10p = 0;
   //cout << p << endl;
   //cout << signedlog10p << endl;
  }
  else if(sign < 0)
  {
   signedlog10p = log10(p);
  }
  else
  {
   signedlog10p = 0;
  }
  symbolValues[lineFields[7]].push_back(signedlog10p);
  if(symbolValues[lineFields[7]].size() == 1)
  {
   symbolInorder.push_back(lineFields[7]);
  }

  getline(input, line);
 }
 input.close();
 
 writeFile(argv[2], output);
 {
  double sum;
  for(int i = 0; i < symbolInorder.size(); i++)
  {
   string str = symbolInorder[i];
   output << str << "\t";
   if(symbolValues[str].size() == 1)
   {
    output << symbolValues[str][0] << endl;
   }
   else
   {
    sum = symbolValues[str][0] ;
    for(int j = 1; j < symbolValues[str].size(); j++)
    {
     sum += symbolValues[str][j];
    }
    output << sum/(double)symbolValues[str].size() << endl;
   }
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


