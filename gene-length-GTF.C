/*
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

class chrPos 
{
 protected:
 string chr;
 int pos;

 public:
 chrPos()
 {return;}
 chrPos(string c, int l):chr(c), pos(l)
 {return;}
 string getChr() const
 {
  return chr;
 }
 int getPos() const
 {
  return pos;
 }
 bool operator<(const chrPos& rh) const
 {
  if(chr < rh.chr)
    return true;
  if(chr > rh.chr)
    return false;
  //same chromsome
  return (pos < rh.pos);
 }
};

class geneLenth
{
 //1 gene, 1 map, usually one chromosome with many positions
 map<chrPos, int> gL;

 public:
 geneLenth()
 {return;}
 void add1pos(const chrPos & cp)
 {
  //only keep unique chr-position
  gL[cp] = 1;
 }
 int getGeneLength()
 {
  int ll=0;
  
  ll = gL.size();
  return ll;
 }
};

map<string, geneLenth> glob_GL;

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);
string getID(string str);

int main(int argc, char *argv[])
{
 if(argc != 3)
 {
  cerr << argv[0] << " GTF_file(1)	op_file(2) " << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;

 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields,'\t', line);
  if(lineFields.size() > 6) //usually 9 fields in GTF
  {
   if(lineFields[2] == "exon" || lineFields[2] == "three_prime_utr" || lineFields[2] == "five_prime_utr")
   {
    string geneID;
    geneID = getID(lineFields[8]);
    for(int i = atoi(lineFields[3].c_str()); i <= atoi(lineFields[4].c_str()); i++)
    {
     //'geneA', exon, 200-250;
     //'geneA', exon, 210-250; from a diff transcript
     //'geneB'
     glob_GL[geneID].add1pos(chrPos(lineFields[0], i));
    }
   }
  }
  getline(input, line);
 }
 input.close();
 
 writeFile(argv[2], output);
 output << "geneID\tlength\n";
 for(map<string, geneLenth>::iterator iter=glob_GL.begin(); iter != glob_GL.end(); iter++)
 {
  output << iter->first << '\t';
  output << (iter->second).getGeneLength() << endl;
 }
 output.close();
 return 0;
}

string getID(string str)
{
 vector<int> pos;

 for(int i = 0; i < str.length(); i++)
 {
  if(str[i] == '"')
    pos.push_back(i);
 }
 if(pos.size() >= 2)
 {
  return str.substr(pos[0]+1,pos[1]-pos[0]-1); 
 }
 else
  return "notAppicable";
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


