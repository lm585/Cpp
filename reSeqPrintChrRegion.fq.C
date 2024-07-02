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

struct chrSeq
{
 string name;
 vector<char> seq;
};

void setChrSeq(ifstream & fasta);

vector<chrSeq> glob_chr;

int main(int argc, char * argv[])
{
 if(argc != 6)
 {
  cerr << argv[0] << "    ref-genome-fasta   output-file-name chr-name start end  \n";
  return 1;
 }

 ifstream input;
 ofstream output;
 map<string, int> chr2i;
 struct chrSeg 
 {
  string chr;
  int begin;
  int end;
 };

 chrSeg segment;

 segment.chr = argv[3];
 segment.begin = atoi(argv[4]);
 segment.end = atoi(argv[5]);
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
  //cout << glob_chr[i].name << '\t' << glob_chr[i].seq.size() << endl;
  chr2i[ glob_chr[i].name ] = i;
 }
 
 {
  int index;
  index = chr2i[ segment.chr];
  cout << ">" << glob_chr[index].name << endl;
  for(int i = segment.begin -1; i < segment.end; i++) //100-200; 0-based index 99-199
    cout << glob_chr[index].seq[i];
  cout << endl;
 }

 output.open(argv[2], ios::out);
 if(!output)
 {
  cerr << argv[2] << " cannot be written to.\n";
  return 1;
 }
/*
 output <<  "chr" << '\t' << "pos"  <<  '\t' << "base" << endl;
 for(int i = 0; i < glob_chr.size(); i++)
 {
  for(int j = 0; j < glob_chr[i].seq.size(); j++)
    output <<  glob_chr[i].name << '\t' << (j+1) <<  '\t' << glob_chr[i].seq[j] << endl;
 }
 */
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


