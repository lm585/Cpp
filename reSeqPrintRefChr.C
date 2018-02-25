#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>
#include <sstream>

using namespace std;

struct chrSeq
{
 string name;
 vector<char> seq;
};

void setChrSeq(ifstream & fasta);

vector<chrSeq> glob_chr;

int main(int argc, char * argv[])
{
 if(argc != 3)
 {
  cerr << argv[0] << "    ref-genome-fasta   output-file-name \n";
  return 1;
 }

 ifstream input;
 ofstream output;

 input.open(argv[1], ios::in);
 if(!input)
   {
    cerr << argv[1] << " cannot be opened for reading!" << endl;
    return 1;
   }
 setChrSeq(input);
 input.close();

 for(int i = 0; i < glob_chr.size(); i++)
 {
  cout << glob_chr[i].name << '\t' << glob_chr[i].seq.size() << endl;
 }
 
 output.open(argv[2], ios::out);
 if(!output)
 {
  cerr << argv[2] << " cannot be written to.\n";
  return 1;
 }
 output <<  "chr" << '\t' << "pos"  <<  '\t' << "base" << endl;
 for(int i = 0; i < glob_chr.size(); i++)
 {
  for(int j = 0; j < glob_chr[i].seq.size(); j++)
    output <<  glob_chr[i].name << '\t' << (j+1) <<  '\t' << glob_chr[i].seq[j] << endl;
 }

 output.close();
 
 return 0;
}

void setChrSeq(ifstream & fasta)
{
 string line, str;
 char * unitStr, * lineCstr;
 vector<string> lineStrList;
 chrSeq myChrSeq;

 getline(fasta, line);
 while(!fasta.eof())
 {
  lineStrList.clear();
  lineCstr = new char [line.length() + 10];
  strcpy(lineCstr, line.c_str());
  unitStr = strtok(lineCstr, " \t\n");
  while(unitStr != NULL)
  {
   str = unitStr;
   lineStrList.push_back(str);
   unitStr = strtok(NULL, " \t\n");
  }

  if(lineStrList.size() > 0)
  {
   if(lineStrList[0][0] == '>')
   {
    if(myChrSeq.seq.size() > 0)
    {
      glob_chr.push_back(myChrSeq);
    }
    myChrSeq.name = lineStrList[0].substr(1);
    myChrSeq.seq.clear();
   }
   else
   {
    for(int i = 0; i < lineStrList.size(); i++)
      for(int j = 0; j < lineStrList[i].length(); j++)
        myChrSeq.seq.push_back(lineStrList[i][j]); 
   }
  }
  getline(fasta, line);
 }
 if(myChrSeq.seq.size() > 0)
 {
  glob_chr.push_back(myChrSeq);
 }
 return;
}


