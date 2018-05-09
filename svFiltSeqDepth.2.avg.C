/*
input:
pileup-file
chr-length (chr1 12345678; chr2 34567890; chr3 678987123)

SV-input-file

param:
copy numver ratio cutoff (>1.5)
ratio between unseq region and unseq chr ( < 1.0)

Calc over the chrs provided in length file, average, medium of seq depth, % of chr not sequenced
For large-dup sv, based on pileup file, calc average, medium of seq depth between [l2, r1], % of the region not sequenced
output discarded dup
output kept dup in a file
output kept dup copy number in 2nd file

map<string, vector<int> > g_pile;


g_pile["chr1"][0] = 0
g_pile["chr1"][1] = 1
g_pile["chr1"][2] = 2
g_pile["chr1"][3] = 3
g_pile["chr1"][4] = 3
g_pile["chr1"][5] = 3
...
g_pile["chr1"][12345678] = 1
g_pile["chr567"][34567890] = 0

 */

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <set>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>

using namespace std;

struct param
{
 double copyRatio;
 double unseqRatio;
};

map<string, vector<int> > g_pile;
param g_param;

void isHighDepth(const vector<string> & lineFields, const string & line, const double * chrStat, ofstream & output);
void calcRegionStat(const string & chr, int begin, int end, double stat[]);
void calcSeqDepStat(double chrStat[]);
void assign_g_pile(const vector<string> & lineFields);
void allocChrMemory(const vector<string> & lineFields);

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);

int main(int argc, char *argv[])
{
 if(argc != 6)
 {
  cerr << argv[0] << " pileup-file  chrLen-file  sv-file copyNumber-ratio(> 1.5) unSeq-ratio-vs-chr(< 1.0)";
  cerr << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line, str;
 vector<string> lineFields;
 double chrStat[10]; //[0] medium, [1] avg, [2] std, [3] unseq fraction, 0.42

 g_param.copyRatio = atof(argv[4]);
 g_param.unseqRatio = atof(argv[5]);
 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t', line);
  allocChrMemory(lineFields);
  getline(input, line);
 }
 input.close();

 cerr << "memory allocated for all chr" << endl;

 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t', line);
  assign_g_pile(lineFields);
  getline(input, line);
 }
 input.close();

 calcSeqDepStat(chrStat); 
 
 str = argv[3];
 str = str + ".seqDepth.filt";
 writeFile(str.c_str(), output);

 readFile(argv[3], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t', line);
  isHighDepth(lineFields, line, chrStat, output);
  getline(input, line);
 }
 input.close();
 
 output.close();
 return 0;
}

void isHighDepth(const vector<string> & lineFields, const string & line, const double * chrStat, ofstream & output)
{
 if(lineFields.size() < 20)
 {
  output << line << endl;
  return;
 }
 if(lineFields[16] != "LARGE_DUPLI")
 {
  output << line << endl;
  return;
 }
 // is large dup
 string chr = lineFields[0];
 double stat[10];
 int begin = atoi(lineFields[2].c_str());
 int end = atoi(lineFields[4].c_str());
 calcRegionStat(chr, begin, end, stat);
 if(stat[1] > g_param.copyRatio * chrStat[1]) 
 {
  //if((chrStat[3] == 0 && stat[3] == 0) || (stat[3] < g_param.unseqRatio * chrStat[3]))
  {
   output << line << endl;
   return;
  }
 }
 //large del, but filtered out
 for(int i = 0; i <= 6; i++)
 {
  cout << lineFields[i] << '\t';
 }
 cout << lineFields[8] << '\t' << lineFields[9] << '\t' << lineFields[16] << '\t';
 cout << lineFields[17] << '\t' << lineFields[19] << '\t';
 cout << chrStat[0] << '\t'  <<  chrStat[1] << '\t' << chrStat[2] << '\t' << chrStat[3]  << '\t';
 cout << stat[0] << '\t'  <<  stat[1] << '\t' << stat[3] <<  endl;
}

void calcRegionStat(const string & chr, int begin, int end, double stat[])
{
 map<string, vector<int> >::const_iterator iter;
 int totalLen = 0, mapLen = 0, unSeqLen;
 long long int totalReads = 0 ;
 int  j;
 double s = 0;
 vector<int> intVect;

 if(begin > end) 
 {
  j = end;
  end = begin;
  begin = j;
 }
 totalLen = end - begin;
 if(g_pile.count(chr) == 0) //SV not in the chr-seq-length file, no memory allocated
 {
  stat[1] = 0; //avg
  stat[3] = 1; //unmapped fraction
  stat[0] = 0; //median
  return;
 }

 for(int i = begin; i < end; i++)
 {
  j= g_pile[chr][i-1];
  if(j > 0)
  {
   mapLen++;
   totalReads += j;
  }
  intVect.push_back(j);
 }
 if(mapLen == 0)
 {
  stat[1] = 0; //avg
  stat[3] = 1; //unmapped fraction
  stat[0] = 0; //median
  return;
 }
 stat[1] = (double) totalReads / (double) totalLen;
 stat[3] = (double) (totalLen - mapLen) / (double) totalLen;
 j = (intVect.size() ) / 2;
 nth_element(intVect.begin(), intVect.begin() + j , intVect.end());
 stat[0] = intVect[j];
 return;
}
void calcSeqDepStat(double chrStat[])
{
 map<string, vector<int> >::const_iterator iter;
 int totalLen = 0, mapLen = 0, unSeqLen;
 long long int totalReads = 0 ; 
 int j;
 double s = 0;
 vector<int> intVect;

 for(iter = g_pile.begin(); iter != g_pile.end(); iter++)
 {
  totalLen += iter->second.size();
 } 
 cout << "total chr length = " << totalLen << endl;
 for(iter = g_pile.begin(); iter != g_pile.end(); iter++)
 {
  for(int i = 0; i < iter->second.size(); i++)
  {
   if(iter->second[i] > 0)
   {
    mapLen++;
    totalReads += iter->second[i];
   }
  }
 }
 chrStat[1] = (double) totalReads / (double) totalLen;
 chrStat[3] = (double) (totalLen - mapLen) / (double) totalLen;
 for(iter = g_pile.begin(); iter != g_pile.end(); iter++)
 {
  for(int i = 0; i < iter->second.size(); i++)
  {
   {
    s += ((double)iter->second[i] - chrStat[1]) * ((double)iter->second[i] - chrStat[1]);
   }
  }
 }
 chrStat[2] =  sqrt(s / (double) (totalLen - 1));
 cout << "mapped seq length = " << mapLen << endl;
 cout <<"un-mapped chr fraction = " << chrStat[3] << endl;
 cout << "whole-chr-seq-depth-avg +/ std " << chrStat[1] << "+/" << chrStat[2] << endl;
 
 for(iter = g_pile.begin(); iter != g_pile.end(); iter++)
 {
  for(int i = 0; i < iter->second.size(); i++)
  {
   {
    intVect.push_back(iter->second[i]);
   }
  }
 }
 j = (intVect.size() ) / 2;
 nth_element(intVect.begin(), intVect.begin() + j , intVect.end());
 chrStat[0] = intVect[j];
 cout << "whole-chr-seq-depth-median " << chrStat[0] << endl; 
}

void assign_g_pile(const vector<string> & lineFields)
{
//seq1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
 if(lineFields.size() < 6) 
   return;

 string chr = lineFields[0];
 int pos = atoi(lineFields[1].c_str());
 int reads = atoi(lineFields[3].c_str());
 
//only anal chr in seqLen file
 if(g_pile.count(chr) > 0)
 {
  g_pile[chr][pos - 1] = reads; //272 24 => seq1[271] = 24
 }
}

void allocChrMemory(const vector<string> & lineFields)
{
 if(lineFields.size() < 2)
   return;

 string chr = lineFields[0];
 int len;
 bool isChr;
 
 if(isdigit(lineFields[1][0]))
 {
  len = atoi(lineFields[1].c_str());
 }
 else
  return;
 //chr1   123456789
 if(g_pile.count(chr) == 0)
 {
  g_pile[chr].resize(len, 0);
 }
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

