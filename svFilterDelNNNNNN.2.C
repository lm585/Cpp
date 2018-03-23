/*
SL2.40ch01	1333706	1333861	SL2.40ch07	52060515	52060662
SL2.40ch01	1881316	1881482	SL2.40ch10	4700284	4700541


ref. genome seq
>SL2.40ch00
AATAATAATAATAATAATAATAATAAATAAATAAATAAATAATAATAATAATAATAATAATAAATAAATAAATAAATAAA
TAAATAAATAAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATA
>SL2.40ch01
GGCTTGTCCGACCCATTTGGAAAGGTCAAACGAGCCTCAAAGCGAGCATACCCCTCATTTCGACGATTTTCGTGTGCTAT


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

using namespace std;

struct chrSeq
{
 string name;
 vector<char> seq;
};

struct param
{
 int leftWinDist;
 int rightWinDist;
 int readLen;
 int numOfN;
 string eValue;
 double identity;
 string pairEnd;
};

class SV
{
 protected:
 string chr[2];
 int window[4];
 string strand[2];
 bool isBal;
 string starts[2];
 
 public:
 void setValues(string ch1, int a, int b, string chr2, int c, int d, string str1, string str2, bool isB, string readStarts1, string readStarts2);
 void getMinMaxSubWin(int minSubWin[], int maxSubWin[]) const;
 bool filtSv(const param & par) ;
 bool peBalFilt(const param & par) ;
 bool peFiltFR(const param & par) const;
 bool peFiltRF(const param & par) const;
 bool peFiltFF(const param & par) const;
 bool peFiltRR(const param & par) const;
 bool filtUnit(int newWindow[], const param & par, string strand2blast) const;
 friend ostream& operator<<(ostream & o, const SV & a); 
};


void setChrSeq(ifstream & fasta);
string coord2seq(string chrName, int a, int b);
void readFile(char * argv, ifstream & input);
void write2file(string s, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);
bool blastMap(string seq1, string seq2Window, const param & par, string strand2blast, bool isBal, string & blastRes, double ACGTratio);

vector<chrSeq> glob_chr;
int glob_lineNum;

ostream& operator<<(ostream & ost, const SV & a)
{
 ost << a.chr[0] << '\t' << a.window[0] << '\t' << a.window[1] << '\t';
 ost << a.chr[1] << '\t' << a.window[2] << '\t' << a.window[3] << '\t'; 
 ost << a.strand[0] << '\t' << a.strand[1] << '\t';
 if(a.isBal) 
   ost << "BAL"  ;
 else
   ost << "UNBAL"  ;
 return ost;
}

void SV::setValues(string chr1, int a, int b, string chr2, int c, int d, string str1, string str2, bool isB, string start1, string start2)
{
 chr[0] = chr1; 
 chr[1] = chr2;
 window[0] = a;
 window[1] = b;
 if(a > b)
 {
  window[0] = b;
  window[1] = a;
 }
 window[2] = c;
 window[3] = d;
 if(c > d)
 {
  window[2] = d;
  window[3] = c;
 }
 strand[0] = str1;
 strand[1] = str2;
 isBal = isB;
 starts[0] =  start1;
 starts[1] =  start2;
}

//(a,b)(c,d)	(a',b')(c',d')
//starts[0]	starts[1]
//a, subwind [0] min, b max
//c  subwind [1] min, d max
//a' subwind [2] min, b' max
//c' subwind [3] min, d' max
void SV::getMinMaxSubWin(int minSubWin[], int maxSubWin[]) const
{
 if(isBal)
 {
  string del = " \t\n(,)";
  vector<string> strVect;
  string line = starts[0] + starts[1];
  getFieldContent2(strVect, del, line);
  for(int i = 0; i <= 6; i = i+2)
  {
   minSubWin[i/2] = atoi(strVect[i].c_str());
   maxSubWin[i/2] = atoi(strVect[i+1].c_str());
  }
 } 
 return;
}

bool SV::filtUnit(int newWindow[4], const param & par, string strand2blast) const
{
 bool res = false;
 string seq1,  seq2Window;

 seq2Window = coord2seq(chr[1], window[1] + 1, window[2] - 1);
 string blastRes1;
 
 double ratio, cutoffStart;

 bool b1 = blastMap(seq1, seq2Window, par, strand2blast, isBal, blastRes1, ratio);

 if( b1 ) //
 {
  cout << glob_lineNum << '\t';
  cout << "NNNinDel\t" ;
  cout << *this;
  cout <<  '\t' << blastRes1 << endl;
  return true;
 }
 return res;
}

bool SV::peFiltFR(const param & par) const
{
 int newWindow[4];
 newWindow[0] = window[0] + par.leftWinDist;
 newWindow[1] = window[1] + par.rightWinDist;
 newWindow[2] = window[2] - par.rightWinDist;
 newWindow[3] = window[3] - par.leftWinDist;
 return filtUnit(newWindow,  par, " -S 1 ");
}

bool SV::peFiltRF(const param & par) const
{
 int newWindow[4];
 newWindow[0] = window[0] - par.rightWinDist;
 newWindow[1] = window[1] - par.leftWinDist;
 newWindow[2] = window[2] + par.leftWinDist;
 newWindow[3] = window[3] + par.rightWinDist;
 return filtUnit(newWindow,  par, " -S 1 ");
}

bool SV::peFiltFF(const param & par) const
{
 int newWindow[4];
 newWindow[0] = window[0] + par.leftWinDist;
 newWindow[1] = window[1] + par.rightWinDist;
 newWindow[2] = window[2] + par.leftWinDist;
 newWindow[3] = window[3] + par.rightWinDist;
 return filtUnit(newWindow,  par, " -S 2 ");
}

bool SV::peFiltRR(const param & par) const
{
 int newWindow[4];
 newWindow[0] = window[0] - par.rightWinDist;
 newWindow[1] = window[1] - par.leftWinDist;
 newWindow[2] = window[2] - par.rightWinDist;
 newWindow[3] = window[3] - par.leftWinDist;
 return filtUnit(newWindow,  par, " -S 2 ");
}

bool SV::peBalFilt(const param & par) 
{
 bool res = false, b1, b2;
 int minSubWin[4], maxSubWin[4];
 getMinMaxSubWin(minSubWin, maxSubWin);
 //F,R,      R,F,
 //R,F,      F,R,
 //F,R,      F,R,
 //R,F,      R,F,
 if(strand[0].find("F,R")  != string::npos && strand[1].find("R,F") != string::npos)
 {
  //(a,b,c),(d,e)	(a',b',c'),(d',e')
  window[0] = minSubWin[0];   //FR	(a,b,c)(a',b',c')  [0][2]
  window[1] = maxSubWin[0] + par.readLen - 1;
  window[2] = minSubWin[2] - par.readLen + 1;
  window[3] = maxSubWin[2];
  b1 = peFiltFR(par); //two sub windows, one F, one R. If true, one connection gone
  //RF (d,e)(d',e')	sub-window[1][3]
  window[0] = minSubWin[1] - par.readLen + 1;
  window[1] = maxSubWin[1];
  window[2] = minSubWin[3];
  window[3] = maxSubWin[3]  + par.readLen - 1; 
  b2 = peFiltRF(par); //two sub windows, one R, one F. If true, the 2nd connection gone
 }
 else if (strand[1].find("F,R")  != string::npos && strand[0].find("R,F") != string::npos)
 {
  //RF (a,b,c)(a',b',c')  [0][2]
  window[0] = minSubWin[0] - par.readLen + 1;
  window[1] = maxSubWin[0];
  window[2] = minSubWin[2];
  window[3] = maxSubWin[2]  + par.readLen - 1;
  b1 = peFiltRF(par); //two sub windows, one R, one F. If true, the 1st connection gone
  //FR (d,e)(d',e')     sub-window[1][3]
  window[0] = minSubWin[1];
  window[1] = maxSubWin[1] + par.readLen - 1;
  window[2] = minSubWin[3] - par.readLen + 1;
  window[3] = maxSubWin[3];
  b2 = peFiltFR(par); //two sub windows, one R, one F. If true, the 2nd connection gone
 }
 else if(strand[0].find("F,R")  != string::npos && strand[1].find("F,R") != string::npos)
 {
  //FF (a,b,c)(a',b',c')  [0][2]
  window[0] = minSubWin[0];
  window[1] = maxSubWin[0] + par.readLen - 1;
  window[2] = minSubWin[2];
  window[3] = maxSubWin[2] + par.readLen - 1;
  b1 = peFiltFF(par); //two whole windows, one F, one F. If true, one connection gone
  //RR (d,e)(d',e')  [1][3]
  window[0] = minSubWin[1] - par.readLen + 1;
  window[1] = maxSubWin[1];
  window[2] = minSubWin[3] - par.readLen + 1;
  window[3] = maxSubWin[3];
  b2 = peFiltRR(par); //two whole windows, one R, one R. If true, the 2nd connection gone
 }
 else if(strand[0].find("R,F")  != string::npos && strand[1].find("R,F") != string::npos)
 {
  //RR (a,b,c)(a',b',c')  [0][2]
  window[0] = minSubWin[0] - par.readLen + 1;
  window[1] = maxSubWin[0];
  window[2] = minSubWin[2] - par.readLen + 1;
  window[3] = maxSubWin[2];
  b1 = peFiltRR(par); 
  //FF  (d,e)(d',e')   [1][3]
  window[0] = minSubWin[1];
  window[1] = maxSubWin[1] + par.readLen - 1;
  window[2] = minSubWin[3];
  window[3] = maxSubWin[3] + par.readLen - 1;
  b2 = peFiltFF(par);
 }
 else
 {
  b1 = false;
  b2 = false;
 }
 if(b1 && b2)
    return true; //filter SV
  else
    return false;

 return res;
}

bool SV::filtSv(const param & par) 
{
 bool res = false;

 if( par.pairEnd == "FR")//paired-end
 {
  if(true)
  {
   if(!isBal) //unBal same chr
   {
    //FR, FF, RF, RR
    if(strand[0].find('F') != string::npos && strand[1].find('R')  != string::npos)
    {
     return peFiltFR(par);
    }
    if(strand[0].find('R') != string::npos && strand[1].find('F')  != string::npos)
    {
     return peFiltRF(par);
    }
    if(strand[0].find('F') != string::npos && strand[1].find('F')  != string::npos)
    {
     return peFiltFF(par);
    }
    if(strand[0].find('R') != string::npos && strand[1].find('R')  != string::npos)
    {
     return peFiltRR(par);
    }
   }
   else // bal same chr
   {
    return peBalFilt(par);
   }
  }
 }
 else //mate-pair
 {
 }

 return res;
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

string coord2seq(string chrName, int a, int b)
{
 string res = "";
 int start, end, len;

  start = a;
  end = b;
 
 for(int i = 0; i < glob_chr.size(); i++)
 {
  if(glob_chr[i].name == chrName)
  {
   for(int j = start; j <= end; j++)
   {
    if(j > 0 && j <= glob_chr[i].seq.size())
      res += glob_chr[i].seq[j-1];
   }
   break;
  }
 }
 return res;
}

void readFile(char * argv, ifstream & input)
{
 input.open(argv, ios::in);
 if(!input)
 {
  cerr << argv << " can not be opened for reading.\n";
  exit(1);
 }
}

void write2file(string s, ofstream & output)
{
 string str = s + ".DELnnn.filt";
 output.open(str.c_str(), ios::out);
 if(!output)
 {
  cerr << str << " cannot be written to.\n";
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

//ratio: perc of ACGT needed in a window
bool blastMap(string seq1, string seq2Window, const param & par, string strand2blast, bool isBal, string & blastRes, double ACGTratio)
{
 bool res = false, firstHit = false;
 int  numOfN = 0, requireN = 1;
 double ratio;
 string  blastCmd = "", line;
 stringstream ss;

 blastRes = "";
 //cerr << glob_lineNum << '\t' << "seq1" << '\t' << seq1.length() << '\t' << "seq2Window" << '\t' << seq2Window.length() << endl;
 for(int i = 0; i < seq2Window.length(); i++)
 {
  if(toupper(seq2Window[i]) == 'N')
    numOfN++;  
 }
 if(isBal)
 {
  ratio = numOfN / (1.0 * seq2Window.length());
 }
 else
 {
  ratio = numOfN / (double) seq2Window.length(); 
 }
 ss <<  numOfN << '\t' << seq2Window.length()  << '\t' << ratio;
 blastRes += ss.str();
 
 if(numOfN >= par.numOfN && ratio >= par.identity)
   res = true;
 return res;
}

int main(int argc, char * argv[])
{
 if(argc != 9)
 {
  cerr << argv[0] << "  sv-file-simple(filtered)   ref-genome-fasta    leftWinDist(mu-read_length)   rightWinDist(mu-read_length) numOfN_cutoff(depending on window size, e.g. 6)    fraction-of-window-with-N(0.10) FR  read-length(80)\n";
  cerr << "Format of sv-file-simple format: \n";
  cerr << "SL2.40ch11      35988688        35989527        SL2.40ch12      22323933        22324032        7       (HWUSI-EAS517:6:1:182:683#0,HWUSI-EAS517:6:59:1655:429#0,HWUSI-EAS517:6:60:361:1239#0,HWUSI-EAS517:6:68:1433:1185#0),(HWUSI-EAS517:6:43:206:1374#0,HWUSI-EAS517:6:110:392:1655#0,HWUSI-EAS517:6:44:93:407#0)     F,R,	R,F, 	(2,1,1,1),(1,1,1)        (1,2,2,2),(2,2,2)       (1,2,3,4),(5,6,7)       (2,4,1,3),(6,7,5)       (a,b)(c,d)	(a',b')(c',d')        TRANSLOC                7/7     BAL     7/7             (35988885,35989426)     (22324031,22323933)     1       7" << endl;
  cerr << "(a,b)(c,d)     (a',b')(c',d')" << "   a: subwind min; b: subwind max. min & max refer to read's start positions. F read, leftmost; R read, rightmost" << endl;
  return 1;
 }
 
 param myParam;
 myParam.leftWinDist = atoi(argv[3]);
 myParam.rightWinDist = atoi(argv[4]);
 myParam.numOfN = atoi(argv[5]);
 myParam.identity = atof(argv[6]);
 myParam.pairEnd = argv[7];
 myParam.readLen = atoi(argv[8]);
 ifstream input;
 input.open(argv[2], ios::in);
 if(!input)
   {
    cerr << argv[2] << " cannot be opened for reading!" << endl;
    return 1;
   }
 setChrSeq(input);
 input.close();
 for(int i = 0; i < glob_chr.size(); i++)
 {
  cout << glob_chr[i].name << '\t' << glob_chr[i].seq.size() << endl;
 }

 ofstream output;
 write2file(string(argv[1]), output);

 readFile(argv[1], input);
 SV mySv;
 string line;
 vector<string> lineFields;
 glob_lineNum = 1;
 bool bal, filt;

 getline(input, line);
 while(!input.eof())
 {
  if(glob_lineNum % 100 == 0)
    cout << "line # " << glob_lineNum << endl;
  getFieldContent(lineFields, '\t', line);
  if(lineFields[19] == "BAL")
    bal = true;
  else
    bal = false;
  mySv.setValues(lineFields[0], atoi(lineFields[1].c_str()), atoi(lineFields[2].c_str()), lineFields[3],  atoi(lineFields[4].c_str()),  atoi(lineFields[5].c_str()),
     lineFields[8], lineFields[9], bal, lineFields[14], lineFields[15]);

  if(true) //lineFields[16] == "DELETION")
  {
   filt = mySv.filtSv(myParam);
   if(!filt)
     output << line << endl;
  }
  else
  {
   output << line << endl; 
  }
  getline(input, line);
  glob_lineNum++;
 }
 input.close();
 output.close();
 return 0;
}
