/*
Chr6    22517399        22517399        Chr6    22517568        22517568        3       F,      R,      DELETION        567     UNBAL   478     13      N   NN       N       N       N       N       N       N       N       N       N       N       N       N       N       N       N       N       N       N       N   NN       N       1       N       N       N       N       N       N       N       N       N       N       N       N       N       N       1       1       1   1N       1       N       N       N       N       N       N       1       1       N       N       N       N       N       N       N       N       N       N   NN       N       N       N       N       N       N       N       N       N       N       N       N       N       N       N       1       N       N       N   NN       N       N       N       N       N       N       N       N       N       N       N       0       N       0       N       N       1       N       N   NN       N       1       N       N       N       0       1       N       0       N       N       0       N       1       1       1       1       1       0   NN       N       1       N       N       1       1       N       N

input 1: svdetect del table
input 2: pindel del table

foreach svd in svdetect
1) if >1 pindel del overlap with svd, the overlapping pindel del discarded
2) no overlap, do nothing
3) if == 1 pindel del overlap, if not sig overlapping, discard the pindel del
4) if == 1 pindel del overlap, if sig  overlapping, discard svdetect one, 
	keep pindel one
	if(svd != pindel)
         if pindel == 'N' pindel=svd (0/1)
         if pindel == '0' && svd == '1'  pindel=svd='1' (0->1)
         if pindel == '1' doing nothing
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

struct SV
{
 bool include;
 vector<string> fields;
 int start;
 int end;
};

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void writeSv(const vector<SV> & sv, ofstream & output);
void mergeSv();
void overlap(vector<int> & ovIds, bool & sigOv, int q);

vector<SV> svd_sv_g, pindel_sv_g;

int main(int argc, char *argv[])
{
 if(argc != 4)
 {
  cerr << argv[0] << " file1(svdetect table)  file2(pindel table) output-file" << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 SV svt;

 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t', line);
  svt.include = true;
  svt.fields = lineFields;
  svt.start = atoi(lineFields[2].c_str());
  svt.end = atoi(lineFields[4].c_str());
  svd_sv_g.push_back(svt);
  getline(input, line);
 }
 input.close();

 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t', line);
  svt.include = true;
  svt.fields = lineFields; //vector<string> copy
  svt.start = atoi(lineFields[2].c_str()); //Chr6    22517399        22517399        Chr6    22517568        22517568   
  svt.end = atoi(lineFields[4].c_str());
  pindel_sv_g.push_back(svt);
  getline(input, line);
 }
 input.close();

 writeFile(argv[3], output);
 mergeSv();
 writeSv(svd_sv_g , output);
 writeSv(pindel_sv_g,  output);
 output.close();
 return 0;
}

void writeSv(const vector<SV> & sv, ofstream & output)
{
 for(int i = 0; i < sv.size(); i++)
 {
  if(sv[i].include == true)
  {
   for(int j = 0; j < sv[i].fields.size(); j++)
   { 
    output << sv[i].fields[j];
    if( j == sv[i].fields.size() - 1)
    {
     output << endl;
    }
    else
    {
     output << "\t";
    }
   }
  }
 }
}

void mergeSv()
{
 vector<int> ovIds;
 bool sigOv;
 string strS, strP;
 for(int i = 0; i < svd_sv_g.size(); i++)
 {
  overlap( ovIds, sigOv, i);
  if(ovIds.size() >= 1)
  {
   if(ovIds.size() > 1)
   {
    for(int j = 0; j < ovIds.size(); j++)
    {
     pindel_sv_g[ovIds[j] ].include = false;
    }
   }
   else //ovIds.size() == 1
   {
    if(sigOv == false)
    {
     pindel_sv_g[ovIds[0] ].include = false;
    }
    else //single 1, sig overlap
    {
     svd_sv_g[i].include = false;
     int p = ovIds[0];
     pindel_sv_g[p].include = true;
     for(int j = 14; j <  svd_sv_g[i].fields.size(); j++)
     {
      strS = svd_sv_g[i].fields[j];
      strP = pindel_sv_g[p].fields[j];
      if(strS != strP)
      {
       if(strP == "N")
       {
        pindel_sv_g[p].fields[j] = strS; //should be 0 or 1
       }
       else if (strP == "0" )
       {
        if(strS == "1")  
        {
         pindel_sv_g[p].fields[j] = strS; //0 -> 1
        }
       }
       else;
      }
     }
    }
   }
  }
 }
 return;
}

void overlap(vector<int> & ovIds, bool & sigOv, int q)
{
 int l1, l2, l;
 ovIds.clear();
 sigOv = false;
 for(int i = 0; i <  pindel_sv_g.size(); i++)
 {
  if( svd_sv_g[q].fields[0] == pindel_sv_g[i].fields[0] && svd_sv_g[q].fields[9] == pindel_sv_g[i].fields[9]  //same chr, same "deletion"
&& !(svd_sv_g[q].start > pindel_sv_g[i].end || pindel_sv_g[i].start > svd_sv_g[q].end)) //coord overlapping
  {
   ovIds.push_back(i);
   l1 = max (pindel_sv_g[i].start, svd_sv_g[q].start);
   l2 =  min(pindel_sv_g[i].end,   svd_sv_g[q].end);
   l = max(pindel_sv_g[i].end - pindel_sv_g[i].start  , svd_sv_g[q].end - svd_sv_g[q].start );
   if(pindel_sv_g[i].start - svd_sv_g[q].start >= -20 && pindel_sv_g[i].start - svd_sv_g[q].start <= 40
      && pindel_sv_g[i].end - svd_sv_g[q].end >= -40 && pindel_sv_g[i].end - svd_sv_g[q].end <= 20
      && (double) (l2 - l1) / (double) (l) >= 0.8 )
   {
    sigOv = true;
   }
  }
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

