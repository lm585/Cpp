/*
input 1:
chr     pos     base    DNA-2228
chr1    1       A       N
chr1    2       C       N
chr1    3       C       N
chr1    4       C       N
chr1    5       T       N
chr1    6       A       N
chr1    7       A       N
chr1    8       A       N
chr1    9       C       N


input 2:
readID  chr     pos     DNA-1380        RNA-P1  DNA-2228
>HWI-ST1085:86:D0U1MACXX:2:1101:10520:24699     chr15   18324777        T       C
>HWI-ST1085:86:D0U1MACXX:2:1101:10650:30062     chr13   10379787        G       A
>HWI-ST1085:86:D0U1MACXX:2:1101:10714:57654     chr7    2571824 G       A
>HWI-ST1085:86:D0U1MACXX:2:1101:10974:44114     chr11   2866781 G       C
>HWI-ST1085:86:D0U1MACXX:2:1101:11561:20451     chr12   20596522        A       G
>HWI-ST1085:86:D0U1MACXX:2:1101:11561:20451     chr12   20596546        A       G

output
readID  chr     pos     DNA-1380        RNA-P1  DNA-2228
>HWI-ST1085:86:D0U1MACXX:2:1101:10520:24699     chr15   18324777        T       C	C
>HWI-ST1085:86:D0U1MACXX:2:1101:10650:30062     chr13   10379787        G       A	A
>HWI-ST1085:86:D0U1MACXX:2:1101:10714:57654     chr7    2571824 G       A	AG
>HWI-ST1085:86:D0U1MACXX:2:1101:10974:44114     chr11   2866781 G       C	C
>HWI-ST1085:86:D0U1MACXX:2:1101:11561:20451     chr12   20596522        A       G	G+AA
>HWI-ST1085:86:D0U1MACXX:2:1101:11561:20451     chr12   20596546        A       G	*

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
void getSeqs(string chr, int i, int j,  vector<string> & seqs);
string getName(string chr, string i, string j, int k);
map< string, vector< string > > g_chrGenotype;

int main(int argc, char *argv[])
{
 if(argc != 4)
 {
  cerr << argv[0] << "  ref3colum-sampGenotype-file(must have a header)  5-columns(readID, chr, pos, scion-DNA, scion-RNA. Must have a header)   outputFile(add one column of rootstock-DNA " << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;

 readFile(argv[1], input);
 getline(input, line);//skip the header
 getline(input, line);
 while(!input.eof())
 {
  string chr;
  getFieldContent(lineFields,'\t', line); //chr1    9       C       N
  chr = lineFields[0];
  g_chrGenotype[chr].push_back(lineFields[3]);
  getline(input, line);
 }
 input.close();
 
 writeFile(argv[3], output);
 readFile(argv[2], input);
 getline(input, line);//skip the header
 output << line <<  endl;
 getline(input, line);
 while(!input.eof())
 {
  string chr, gen;
  int pos;
  getFieldContent(lineFields,'\t', line); //>HWI-ST1085:86:D0U1MACXX:2:1101:10520:24699     chr15   18324777        T       C
  chr = lineFields[1];
  pos = atoi(lineFields[2].c_str());//pos 100, index [100-1]
  gen = g_chrGenotype[chr][pos-1];
  output << line << '\t' << gen << endl;
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

string getName(string chr, string i, string j, int k)
{
 string str = chr + '_' + i  + '_' + j  + '_';
 string Result; 
 ostringstream convert;
 convert << k;
 Result = convert.str();
 str = str + Result;
 return str;
}

void recGetSeqs(string chr, int i, int j, int k, vector<string> & seqs, string & fa)
{
 if(k > j)
 {
  seqs.push_back(fa); 
  return;
 }
 //chr1 1 1000 size==1000; [0] refers to position 1, [999] refers to position 1000
 if(g_chrGenotype.count(chr) == 0)
   return;
 if(i < 1 || j > g_chrGenotype[chr].size())
   return;

 string gen = g_chrGenotype[chr][k-1];
 int pos,pos2, t = g_chrGenotype[chr][k-1].length();
 if(gen.find('*') != string::npos) //*   A*
 {
  recGetSeqs(chr, i, j, k+1,  seqs, fa);
  return;
 }
 if(gen.find('+') != string::npos) //A+AAC
 {
  fa += gen[0]; //'A'
  pos = gen.find_last_of('+');
  fa += gen.substr(pos + 1); //"AAC"
  recGetSeqs(chr, i, j, k+1,  seqs, fa);
  return;
 }
 if( t == 1) //A  N
 {
  fa += gen[0];
  recGetSeqs(chr, i, j, k+1,  seqs, fa);
  return;
 }
 // AC ACG
 fa += '~'; //place holder
 pos2 = fa.length() - 1;
 for(int m = 0; m < t; m++)
 {
  fa[pos2] = gen[m];
  recGetSeqs(chr, i, j, k+1,  seqs, fa);
  fa = fa.substr(0, pos2 + 1);
 }
}

void getSeqs(string chr, int i, int j,  vector<string> & seqs)
{
 string fa = "";
 recGetSeqs(chr, i, j, i, seqs, fa);
}
