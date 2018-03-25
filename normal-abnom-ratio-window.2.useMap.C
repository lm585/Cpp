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

struct param
{
 double normAbnRatio;
 int readLen;
 string pairEnd;
};

class mappedRead
{
 protected:
 string chr;
 int leftPos;
 int rightPos;
 bool forward; //true forward; false reverse

 public:
 bool parseOneRead(string line); //if false, head line; It sets all member values in the class
 mappedRead()
 {return;}
 mappedRead(string c, int l, int r, bool f):chr(c), leftPos(l), rightPos(r), forward(f)
 {return;}
 string getChr() const { return chr;}
 int getLeftPos() const { return leftPos;}
 int getRightPos() const {return rightPos;}
 bool getStrand() const { return forward;}
 bool operator<(const mappedRead& rh) const
 {
  if(chr < rh.chr)
    return true;
  if(chr > rh.chr)
    return false;
  if(forward && !rh.forward) //true && false)
    return true;
  if(!forward && rh.forward) //false && true)
    return false;
  //same chr, same strand
  if(leftPos < rh.leftPos)
    return true;
  if(leftPos > rh.leftPos)
    return false;
  //equivalent mapped read
  return false;
 }
};

class compare
{
 public: 
 bool operator() (const mappedRead& lh, const mappedRead& rh) const
 {
  return lh < rh;
 }
};

//read list per chromosome-strand. Thus millions of reads will be grouped into chromosomes.
//It will increase speed . It will also save memory.

class SV
{
 protected:
 string chr[2];
 int window[4];
 int numOfPairs;
 string strand[2];
 bool isBal;
 string starts[2];
 string rank[2];
 
 public:
 void setValues(string ch1, int a, int b, string chr2, int c, int d, string str1, string str2, bool isB, string readStarts1, string readStarts2, int n, string rank1, string rank2);
 void getMinMaxSubWin(int minSubWin[], int maxSubWin[]) const;
 bool filtSv(const param & par, const map<mappedRead, int, compare> & allReads);
 bool peBalFilt(const param & par, const map<mappedRead, int, compare> & allReads) ;
 bool peFiltFR(const param & par, const map<mappedRead, int, compare> & allReads) const;
 bool peFiltRF(const param & par, const map<mappedRead, int, compare> & allReads) const;
 bool peFiltFF(const param & par, const map<mappedRead, int, compare>& allReads) const;
 bool peFiltRR(const param & par, const map<mappedRead, int, compare>& allReads) const;
 bool filtUnit( const param & par, bool dir[], const map<mappedRead, int, compare>& allReads) const;
 friend ostream& operator<<(ostream & o, const SV & a); 
};


void parseSAMflag(int flag, int & strand, bool & isProperPair,  bool & is2ndRead, bool & isMapped);
int coordReverseRead(int pos, string cigar);
void insertReads(const mappedRead & read, map<mappedRead, int, compare> & allReads);
int getNormPairs(string chr, int w0,  int w1, bool dir, const map<mappedRead, int, compare>& allReads);
void getBalNumReads(string rank, int num[]);
void print(const map<mappedRead, int, compare> & allReads);
void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);

param glob_param;
int glob_lineNum;

bool mappedRead::parseOneRead(string line)
{
 if(line.length() == 0 || line[0] == '@') 
   return false;  

 char * linePtr = new char [line.length()+ 10];
 char *f;
 int field = 1, flag, strand, pos;
 bool isProperPair, is2ndRead, isMapped;
 string cigar;

 strcpy(linePtr, line.c_str());
 f = strtok(linePtr, "\t");
 while(f != NULL)
 {
  if(field == 2)
    flag = atoi(f);
  if(field == 3)
    this->chr =  f;
  if(field == 4)
    pos = atoi(f);
  if(field == 6)
    cigar = f;
  f = strtok(NULL, "\t");
  field++;
 }
 parseSAMflag(flag, strand, isProperPair,  is2ndRead, isMapped);
 if(strand == 1) //reverse strand
 {
  this->forward = false;
 }
 else //forward strand
 {
  this->forward = true;
 }
 this->leftPos = pos;
 this->rightPos = this->leftPos + glob_param.readLen - 1; //coordReverseRead(pos, cigar);

 delete [] linePtr;

 if(field > 10) //at least 11 mandatory sam fields 
   return true;
 else
   return false;
}

void SV::setValues(string chr1, int a, int b, string chr2, int c, int d, string str1, string str2, bool isB, string start1, string start2, int n, string rank1, string rank2 )
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

 numOfPairs = n;
 rank[0] = rank1;
 rank[1] = rank2;
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

//Window A (w[0], w[1], chrA, F), window B (w[2], w[3], chrB, R), #ofpairs
//bool [2] for two windows' strands, respectively
//all normal reads list
bool SV::filtUnit( const param & par, bool dir[], const map<mappedRead, int, compare>& allReads) const
{
 int norm_1 = getNormPairs(chr[0], window[0],  window[1], dir[0], allReads);
 int norm_2 = getNormPairs(chr[1], window[2],  window[3], dir[1],  allReads);
 double r1 = norm_1 / (double) numOfPairs;
 double r2 = norm_2 / (double) numOfPairs;
 if(r1 > par.normAbnRatio || r2 > par.normAbnRatio) //too many normal pairs
 {
  cout << glob_lineNum << '\t';
  cout << *this;
  cout << numOfPairs << '\t' << norm_1 << '\t' << r1 << '\t' << norm_2 << '\t' << r2 << endl;
  return true; //SV link removed
 }
 return false;

}

bool SV::peFiltFR(const param & par, const map<mappedRead, int, compare>& allReads) const
{
 bool dir[2];
 dir[0] = true;
 dir[1] = false;
 return filtUnit(par, dir, allReads);
}

bool SV::peFiltRF(const param & par, const map<mappedRead, int, compare>& allReads) const
{
 bool dir[2];
 dir[0] = false;
 dir[1] = true;
 return filtUnit(par, dir, allReads);
}

bool SV::peFiltFF(const param & par, const map<mappedRead, int, compare>& allReads) const
{
 bool dir[2];
 dir[0] = true;
 dir[1] = true;
 return filtUnit(par, dir, allReads); 
}

bool SV::peFiltRR(const param & par, const map<mappedRead, int, compare>& allReads) const
{
 bool dir[2];
 dir[0] = false;
 dir[1] = false;
 return filtUnit(par, dir, allReads);
}

bool SV::peBalFilt(const param & par, const map<mappedRead, int, compare>& allReads)
{
 bool res = false, b1, b2;
 int minSubWin[4], maxSubWin[4];
 int num[2];
 getMinMaxSubWin(minSubWin, maxSubWin);
 getBalNumReads(rank[0], num);
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
  numOfPairs = num[0];
  b1 = peFiltFR(par, allReads); //two sub windows, one F, one R. If true, one connection gone
  //RF (d,e)(d',e')	sub-window[1][3]
  window[0] = minSubWin[1] - par.readLen + 1;
  window[1] = maxSubWin[1];
  window[2] = minSubWin[3];
  window[3] = maxSubWin[3]  + par.readLen - 1; 
  numOfPairs = num[1];
  b2 = peFiltRF(par, allReads); //two sub windows, one R, one F. If true, the 2nd connection gone
 }
 else if (strand[1].find("F,R")  != string::npos && strand[0].find("R,F") != string::npos)
 {
  //RF (a,b,c)(a',b',c')  [0][2]
  window[0] = minSubWin[0] - par.readLen + 1;
  window[1] = maxSubWin[0];
  window[2] = minSubWin[2];
  window[3] = maxSubWin[2]  + par.readLen - 1;
  numOfPairs = num[0];
  b1 = peFiltRF(par, allReads); //two sub windows, one R, one F. If true, the 1st connection gone
  //FR (d,e)(d',e')     sub-window[1][3]
  window[0] = minSubWin[1];
  window[1] = maxSubWin[1] + par.readLen - 1;
  window[2] = minSubWin[3] - par.readLen + 1;
  window[3] = maxSubWin[3];
  numOfPairs = num[1];
  b2 = peFiltFR(par, allReads); //two sub windows, one R, one F. If true, the 2nd connection gone
 }
 else if(strand[0].find("F,R")  != string::npos && strand[1].find("F,R") != string::npos)
 {
  //FF (a,b,c)(a',b',c')  [0][2]
  window[0] = minSubWin[0];
  window[1] = maxSubWin[0] + par.readLen - 1;
  window[2] = minSubWin[2];
  window[3] = maxSubWin[2] + par.readLen - 1;
  numOfPairs = num[0];
  b1 = peFiltFF(par, allReads); //two whole windows, one F, one F. If true, one connection gone
  //RR (d,e)(d',e')  [1][3]
  window[0] = minSubWin[1] - par.readLen + 1;
  window[1] = maxSubWin[1];
  window[2] = minSubWin[3] - par.readLen + 1;
  window[3] = maxSubWin[3];
  numOfPairs = num[1];
  b2 = peFiltRR(par, allReads); //two whole windows, one R, one R. If true, the 2nd connection gone
 }
 else if(strand[0].find("R,F")  != string::npos && strand[1].find("R,F") != string::npos)
 {
  //RR (a,b,c)(a',b',c')  [0][2]
  window[0] = minSubWin[0] - par.readLen + 1;
  window[1] = maxSubWin[0];
  window[2] = minSubWin[2] - par.readLen + 1;
  window[3] = maxSubWin[2];
  numOfPairs = num[0];
  b1 = peFiltRR(par, allReads); 
  //FF  (d,e)(d',e')   [1][3]
  window[0] = minSubWin[1];
  window[1] = maxSubWin[1] + par.readLen - 1;
  window[2] = minSubWin[3];
  window[3] = maxSubWin[3] + par.readLen - 1;
  numOfPairs = num[1];
  b2 = peFiltFF(par, allReads);
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

bool SV::filtSv(const param & par, const map<mappedRead, int, compare> & allReads) 
{
 bool res = false;
 if( true ) //par.pairEnd == "FR")//paired-end
 {
  if(true)
  {
   if(!isBal) //unBal same chr
   {
    //FR, FF, RF, RR
    if(strand[0].find('F') != string::npos && strand[1].find('R')  != string::npos)
    {
     return peFiltFR(par, allReads);
    }
    else if(strand[0].find('R') != string::npos && strand[1].find('F')  != string::npos)
    {
     return peFiltRF(par, allReads);
    }
    else if(strand[0].find('F') != string::npos && strand[1].find('F')  != string::npos)
    {
     return peFiltFF(par, allReads);
    }
    else if(strand[0].find('R') != string::npos && strand[1].find('R')  != string::npos)
    {
     return peFiltRR(par, allReads);
    }
    else
      return false;
   }
   else // bal same chr
   {
    return peBalFilt(par, allReads);
   }
  }
 }
 else //mate-pair
 {
 }

 return res;
}


ostream& operator<<(ostream & ost, const SV & a)
{
 ost << a.chr[0] << '\t' << a.window[0] << '\t' << a.window[1] << '\t';
 ost << a.chr[1] << '\t' << a.window[2] << '\t' << a.window[3] << '\t'; 
 ost << a.strand[0] << '\t' << a.strand[1] << '\t';
 if(a.isBal) 
   ost << "BAL"  << '\t';
 else
   ost << "UNBAL"  << '\t';
 return ost;
}


void parseSAMflag(int flag, int & strand, bool & isProperPair,  bool & is2ndRead, bool & isMapped)
{
 int q, r;
 int time = 0 ;
 int pp, secondR, map;

 q = flag / 2;
 r = flag % 2;
 while(true)
 {
  time++;
  if(time == 2)
    pp = r;
  if(time == 3)
    map = r;
  if(time == 5)
    strand = r;
  if(time == 8)
    secondR = r;
  if(q == 0)
    break;
  r = q % 2;
  q /=  2;
 }

 if(pp == 0)
   isProperPair = false;
 else
   isProperPair = true;

 if(secondR == 0)
   is2ndRead = false;
 else
   is2ndRead = true;

 if(map == 1)
   isMapped = false;
 else
   isMapped = true;
}

int coordReverseRead(int pos, string cigar)
{
 int end = pos;
 int length = 0;
 string numOfOp;

 for(int i = 0; i < cigar.length(); i++)
 {
  if(isdigit(cigar[i]))
  {
   numOfOp += cigar[i];
  }
  else //MID
  {
   if(cigar[i] == 'M' || cigar[i] == 'D')
     length += atoi(numOfOp.c_str());
   numOfOp = "";
  }
 }
 int coord = end - 1 + length;
 //cout << cigar << "\t" << length << "\t" << end << "\t" << ss.str() << endl;
 return coord;
}

void getBalNumReads(string rank, int num[])
{
 string del1 = "()", del2 = ",";
 vector<string> strVect1, strVect2;
 int count = 0;

 getFieldContent2(strVect1, del1 , rank);
 for(int i = 0; i < strVect1.size(); i++) //strVect1[1] == "," for BAL
 {
  getFieldContent2(strVect2, del2, strVect1[i]);
  count = 0;
  for(int j = 0; j < strVect2.size(); j++) //2556496 2556528  2556529  2556533   2557713$
  {
   if(strVect2[j].find('$') == string::npos && strVect2[j].find('*') == string::npos && strVect2[j].find('@') == string::npos)
   {
    count++;
   }
  }
  if(i == 0)
    num[0] = count;
  if(i == 2)
    num[1] = count;
 }
}


void insertReads(const mappedRead & read, map<mappedRead, int, compare> & allReads)
{
 if(allReads.count(read) == 0) //a new read
   allReads[read] = 1;
 else //already exist
   allReads[read]++;
 return;
}

int getNormPairs(string chr, int w0,  int w1, bool dir, const map<mappedRead, int, compare>& allReads)
{
 mappedRead read(chr, w0, 0, dir);
 map<mappedRead, int, compare>::const_iterator iter;
 int res = 0;

 iter = allReads.lower_bound(read);
 //normal reads >= the w0-chr-dir read; however they must be same chr and dir (then leftPos >= W0)
 while(iter != allReads.end())
 {
  string c = (*iter).first.getChr();
  int r = (*iter).first.getRightPos();
  bool b = (*iter).first.getStrand();
  if(c == chr && b == dir && r <= w1)
  {
   res += (*iter).second;
   iter++;
  }
  else
    break;
 } 
 return res;  
}


void print(const map<mappedRead, int, compare> & allReads)
{
 char t = '\t';
 cout << "# of unique reads in the map" << t << allReads.size() << endl;
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


int main(int argc, char *argv[])
{
 if(argc != 5)
 {
  cerr << "executable SAM-normal-paired-end-input-file  sv-file-simple   normAbnRatio-cutoff(0.25) read-length(80)" << endl;
  cerr << "Format of sv-file-simple format: \n";
  cerr << "SL2.40ch11      35988688        35989527        SL2.40ch12      22323933        22324032        7       (HWUSI-EAS517:6:1:182:683#0,HWUSI-EAS517:6:59:1655:429#0,HWUSI-EAS517:6:60:361:1239#0,HWUSI-EAS517:6:68:1433:1185#0),(HWUSI-EAS517:6:43:206:1374#0,HWUSI-EAS517:6:110:392:1655#0,HWUSI-EAS517:6:44:93:407#0)     F,R,	R,F, 	(2,1,1,1),(1,1,1)        (1,2,2,2),(2,2,2)       (1,2,3,4),(5,6,7)       (2,4,1,3),(6,7,5)       (a,b)(c,d)	(a',b')(c',d')        TRANSLOC                7/7     BAL     7/7             (35988885,35989426)     (22324031,22323933)     1       7" << endl;
  cerr << "(a,b)(c,d)     (a',b')(c',d')" << "   a: subwind min; b: subwind max. min & max refer to read's start positions. F read, leftmost; R read, rightmost" << endl;

  return 1;
 }

 ifstream input;
 ofstream output;
 string line1, line2, line, chr;
 mappedRead read;
 map<mappedRead, int, compare> allReads;
 map<mappedRead, int, compare>::iterator mapIter;
 vector<string> lineFields;
 bool h, forward;
 int left, right, oneMil = 1000000;
 SV mySv;
 glob_lineNum = 1;
 bool bal, filt;


 glob_param.normAbnRatio = atof(argv[3]);
 glob_param.readLen = atoi(argv[4]);
 readFile(argv[1], input);
 //processing normal paired sam file
 getline(input, line1);
 while(!input.eof())
 {
  h = read.parseOneRead(line1); 
  if(h) //not head, has at least 11 fields
  {
   insertReads(read, allReads);
   if(allReads.size() % oneMil == 0)
     cout << allReads.size() << " unique normal reads added to the map." << endl;                                
  }
  else
  {
   cout << "head\tline1" << endl;
  }
  getline(input, line1);
 }
 input.close(); 

 print(allReads);
 
 readFile(argv[2], input);
 line2 = argv[2];
 line2 += ".normPair.filt";
 writeFile(line2.c_str(), output);
 //processing SV link file
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
     lineFields[8], lineFields[9], bal, lineFields[14], lineFields[15],
     atoi(lineFields[6].c_str()), lineFields[10], lineFields[11]);
  filt = mySv.filtSv(glob_param, allReads);
  if(!filt) //not true, no filter, keep
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
