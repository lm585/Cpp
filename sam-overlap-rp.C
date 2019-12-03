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

struct intChar 
{
 int num;
 char mid; //'M' 'I' 'D'
};
 
void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);
bool isSAMheader(const char *f);
bool isCigarMID_Read(string cigar);
bool parseCigar(string cigar, vector<intChar> & vC);
int parseReadMap(string line, string & res);
string int2str(int n);

int main(int argc, char *argv[])
{
 if(argc != 3)
 {
  cerr << argv[0] << "	SAM-file-ip(1)	output(2)	" << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line, res;
 vector<string> lineFields;
 int parseRes;

 writeFile(argv[2], output);
 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  parseRes = parseReadMap(line, res);
  if(parseRes == 1)
    output << line << endl;
  else if(parseRes == -1)
    output << res << endl;
  else;
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

bool isSAMheader(const char *f)
{
 bool res = false;
 if(strlen(f) == 3)
   if(f[0] == '@')
     if(isalpha(f[1]) && isalpha(f[2]))
       res = true;

 return res;
}

bool isCigarMID_Read(string cigar)
{
 string chars = "0123456789MID";
 for(int i = 0; i < cigar.length(); i++)
 {
  if(chars.find(cigar[i]) == string::npos) //cigar contains non \dMID characters
  {
    return false;
  }
 }
 return true;
}

bool parseCigar(string cigar, vector<intChar> & vC)
{
 string numOfOp = "";
 intChar cigarPair;
 vC.clear();

 if(isCigarMID_Read(cigar))
 {
  for(int i = 0; i < cigar.length(); i++)
  {
   if(isdigit(cigar[i]))
   {
    numOfOp += cigar[i];
   }
   else //MID
   {
    cigarPair.num = atoi(numOfOp.c_str());
    cigarPair.mid = cigar[i]; //M I D
    vC.push_back(cigarPair);
    numOfOp = "";
   }
  }
  return true;
 }
 else
   return false;
}

//SAM header line, return 1;
//read line, cannot recognize cigar string, read kept, return 1;
//read; cigar recognized; no overlapping, return 1;
//read; cigar recognized; overlappiong; return -1 and a new line string
//if entire read is trimmed, return 0
int parseReadMap(string line, string & res)
{
 res = "";
 vector<string>  fields;
 bool bo;
 vector<intChar> vC;
 int a, prev, ic, start, mateStart, keptBases;
 string newCigar = "";

 getFieldContent( fields, '\t', line);
 if(fields.size() == 0)
   return 1;
 if(isSAMheader(fields[0].c_str()))
   return 1;
//read line
 if(fields.size() < 11)
   return 1;
 if(atoi(fields[8].c_str()) <= 0) //read in downstream 
   return 1; 
 bo = parseCigar(fields[5], vC);
 if(!bo)//bo == false,   cannot recognize cigar string
   return 1;
 //cigar recognized; read upstream
 a = 0;
 for(int i = 0; i < vC.size(); i++)
 {
  if(vC[i].mid == 'M' || vC[i].mid == 'D')
    a += vC[i].num;
 }
 start = atoi(fields[3].c_str());
 mateStart = atoi(fields[7].c_str());
 if(!(start - 1 + a >= mateStart))
   return 1; //two reads no overlapping
 a = 0; prev = 0;
 keptBases = 0;
 newCigar = "";
 for(int i = 0; i < vC.size(); i++)
 {
  prev = a;
  if(vC[i].mid == 'M')
  {
    a += vC[i].num; //101M; 30M+2D; curr extension L from start
    if(start + a >= mateStart) //reached mate start position
    {
     ic = mateStart - start - prev;
     newCigar = newCigar + int2str(ic) + 'M';
     keptBases += ic;
     break;
    }
    else
    {
     newCigar = newCigar + int2str(vC[i].num) + 'M';
     keptBases += vC[i].num;
    }
  }
  if(vC[i].mid == 'D')
  {
    a += vC[i].num; //101M; 30M+2D; curr extension L from start
    if(start + a >= mateStart) //reached mate start position
    {
     ic = mateStart - start - prev;
     break;
    }
    else
    {
     newCigar = newCigar + int2str(vC[i].num) + 'D';
    }
  }
  if(vC[i].mid == 'I')
  {
    if(start + a >= mateStart) //reached mate start position
    {
     ic = mateStart - start - prev;
     break;
    }
    else
    {
     newCigar = newCigar + int2str(vC[i].num) + 'I';
     keptBases += vC[i].num;
    }
  }
 } //end for(int i = 0; i < vC.size(); i++)
 if(keptBases == 0)
   return 0;
 res = "";
 for(int i = 0; i < fields.size(); i++)
 {
  string temps;
  if(i == 5)
  {
   temps = newCigar; 
  }  
  else if (i == 9 || i == 10)
  {
   temps = fields[i].substr(0, keptBases);
  }
  else
    temps = fields[i];

  res += temps;
  if(i != fields.size() -1 )
    res += "\t";
 }
 cout << fields[0] << '\t' << fields[3] << '\t' << fields[5] << '\t' << fields[7] << '\t';
 cout << newCigar << '\t' << keptBases << endl;
 return -1;
}

string int2str(int n)
{
 stringstream ss;
 ss << n;
 return ss.str();
}
