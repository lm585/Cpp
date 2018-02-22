/*
HWUSI-EAS517_0047:1:1:1445:950#0        83      SL2.40ch09      116875  9       101M    =       116425  -551    TCTTGATGGACGTCCACGAAAAAATTTGGCGTTTTTGAAGTCGGAATCCGGATCACCCAAAAAATCATGTGCTATAGCACACGAAAATCGTCGAAATGAGN   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB_QQ__VVTPVVVTVTNHLKNT[VY[______Y___NNMNNJLOKD   XT:A:R  NM:i:4  SM:i:0  AM:i:0  X0:i:3  X1:i:0  XM:i:4  XO:i:0  XG:i:0  MD:Z:30A7C50A10G0       XA:Z:SL2.40ch09,-116875,101M,4;SL2.40ch09,-118319,101M,4;
HWUSI-EAS517_0047:1:1:1445:950#0        163     SL2.40ch09      116425  9       101M    =       116875  551     AAATGAGGGGTATGCTCGCTTCGGGGCTCGTTTGACCTTCCAAATGGGTCGACCAAGCCGTGATGGCCAACCGTATGCATAGACAAGGTCTTGCGGGACGC   gagggggggggagggcfcdfdffff_cfffcffff]eeWdbb[d]addQ^SFFFPQUXKPW`YY```^^\edaaV[QWK]V[Y^BBBBBBBBBBBBBBBBB   XT:A:U  NM:i:5  SM:i:0  AM:i:0  X0:i:1  X1:i:22 XM:i:5  XO:i:0  XG:i:0  MD:Z:51G0A40A0C5T0

Given SAM paired-end output, the C++ program will prepare the input to the BAM_preprocessingPairs.pl program in the SVDetect package.
The perl program generates average and std_dev for insertion size. It also generates anomalous
paired-ends and normal paired-ends.

The C++ program has the following functions:
1. Remove the header lines of SAM file. So that multiple SAM files can be concatenated later.
2. Properly paired reads appear before other not properly paired reads in the output file. 
3. both reads in a pair need to be uniquely mapped reads
4. Remove duplicated read pairs.

Read_pair_1, read_pair_2

pair_1 and pair_2 are duplicated if:
1) 1st read in pair_1 and 1st read in pair_2 are duplicated
2) 2nd read in pair_1 and 2nd read in pair_2 are also duplicated.
3) two reads are duplicated if they mapped to the same chromosome, same 
coordinates (either starting or ending base), same orientations

1st or 2nd read judged from SAM bitwise flag.
How to determine the two  first reads are duplicated, they have same chromosomes,
same coordinates, same orientation (from bitwise flag).

this version,
some SAM flag (e.g. 167), the read is unique, NM < 4, proper paired, but unmapped, 
- strange flag
It will remove a pair in which a read has unmapped bit set.


*/

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>
#include <map>

using namespace std;


class readPairKey
{
 public:
  void setMemberValues(string str1, string str2, int p1, int p2, int s1, int s2);
  void setTwoReadsPos(int p1, int p2);
  bool operator< (const  readPairKey & x) const;
  friend ostream& operator<< (ostream& ostr, const readPairKey & x); 
 private:
  string chr1, chr2;
  int pos1, pos2;
  int strand1, strand2;
};

void readPairKey::setMemberValues(string str1, string str2, int p1, int p2, int s1, int s2)
{
 chr1 = str1;
 chr2 = str2;
 pos1 = p1;
 pos2 = p2;
 strand1 = s1;
 strand2 = s2;
}

void readPairKey::setTwoReadsPos(int p1, int p2)
{
 pos1 = p1;
 pos2 = p2;
}

bool readPairKey::operator< (const  readPairKey & rhs  ) const
{
 bool res;

 if(chr1 < rhs.chr1)
   return true;
 if(chr1 > rhs.chr1)
   return false;
 //chr1 == rhs.chr1
  if(pos1 < rhs.pos1)
    return true;
  if(pos1 > rhs.pos1)
    return false;
  //pos1 == rhs.pos1
  if(strand1 < rhs.strand1)
    return true;
  if(strand1 > rhs.strand1)
    return false;
  //1st read duplicated
  if(chr2 < rhs.chr2)
    return true;
  if(chr2 > rhs.chr2)
    return false;
  //chr2 == rhs.chr2
  if(pos2 < rhs.pos2)
    return true;
  if(pos2 > rhs.pos2)
    return false;
  //pos2 == rhs.pos2
  if(strand2 < rhs.strand2)
    return true;
  if(strand2 > rhs.strand2)
    return false;
  //duplicated pair
  return false;
}

ostream& operator<< (ostream& ostr, const readPairKey & x)
{
 ostr << "Chr1:" << x.chr1 << "\t";
 ostr << "pos1:" << x.pos1 << "\t";
 ostr << "strand:" << x.strand1 << "\t";
 ostr << "Chr2:" << x.chr2 << "\t";
 ostr << "pos2:" << x.pos2 << "\t";
 ostr << "strand:" << x.strand2 << "\t";
 return ostr;
}

struct classcomp
{
  bool operator() (const readPairKey & lhs, const readPairKey & rhs) const
  {return lhs<rhs;}
};


class readPairProperties
{
 public:
  string read1, read2;
  int read1Pos[2], read2Pos[2];
  bool isProperPair;
  int sumNM;
};

void readHeader(string line1, vector<string> &  headLines,  bool & firstRead);
void populateRP(string line1, string line2,  vector<readPairKey> & rpKeys,  vector<readPairProperties> & rpProperties, int nm_cutoff);
void write2file(ofstream & output, const  map<readPairKey, int, classcomp>  & readlocLinenumber,  const vector<readPairProperties> & rpProperties);
void removeRedun( map<readPairKey, int, classcomp> & readlocLinenumber,const vector<readPairKey> &  rpKeys,vector<readPairProperties> & rpProperties);
void removeRedunMultiScans(map<readPairKey, int, classcomp> & readlocLinenumber, vector<readPairKey> &  rpKeys,vector<readPairProperties> & rpProperties);

int main(int argc, char *argv[])
{
 if(argc != 4)
 {
  cout << "executable SAM-paired-end-mapping-input-file  output-file-name NM-cutoff" << endl;
  cout << "Reads with NM > cut0ff will be discarded" << endl;
  return 1;
 }

 vector<string> headLines;
 vector<readPairKey> rpKeys;
 vector<readPairProperties> rpProperties;
 map<readPairKey, int, classcomp> readlocLinenumber;
 map<readPairKey, int, classcomp>::iterator mapIt;
 string line1, line2;
 int readBeginLineNum = 0, nmCutoff = 100;
 bool firstRead = false;
 ifstream input;
 ofstream output;


 nmCutoff = atoi(argv[3]);
 input.open(argv[1], ios::in);
 if(!input)
   {
    cerr << argv[1] << "SAM file cannot be opened for reading!" << endl;
    return 1;
   }
 getline(input, line1);
 while(!input.eof())
 {
  readHeader(line1,  headLines,  firstRead);
  if(firstRead)
  {
   break;
  }
  else
  {
   readBeginLineNum++;
   getline(input, line1);
  }
 }
 input.close();

 cout << readBeginLineNum << endl;
 cout << "SVDetectInput.3.C.  MultiScans. Either ends of two reads same as duplicated\n";
 cout << "Remove reads with NM > 4 \n";
 cout << "remove read with contradictory bitwise flag \n";
 input.open(argv[1], ios::in);
 for(int i = 0; i < readBeginLineNum; i++)
 {
  getline(input, line1);
 }
 while(!input.eof())
 {
  getline(input, line1);
  getline(input, line2);
  populateRP(line1, line2,  rpKeys, rpProperties, nmCutoff);
 }
 input.close();

 removeRedunMultiScans(readlocLinenumber,  rpKeys, rpProperties);
 output.open(argv[2], ios::out);
 if(!output)
   {
    cerr << argv[2] << " cannot be opened for writing" << endl;
    return 1;
   }
 write2file(output, readlocLinenumber, rpProperties);
 output.close();

 return 0;
}

bool isSAMheader(char *f)
{
 bool res = false;
 if(strlen(f) == 3)
   if(f[0] == '@')
     if(isalpha(f[1]) && isalpha(f[2]))
       res = true;

 return res;
}


void readHeader(string line1, vector<string> &  headLines,  bool & firstRead)
{
 char * linePtr = new char [line1.length()+ 10];
 strcpy(linePtr, line1.c_str());
 char *f = strtok(linePtr, "\t");
 
 if(f != NULL)
 {
  if(isSAMheader(f) && firstRead == false)
  {
   headLines.push_back(line1);
  }
  else 
  {
   firstRead = true;
  }
 }
 delete [] linePtr;
}

void parseMisMatch(char * f, int & nm)
{
 string str = f;
 int t = 0;
 int pos[2];

 for(int i = 0; i < str.length(); i++)
 {
  if(str[i] == ':')
  {
   t++;
   if(t == 1)
     pos[0] = i;
  }
  if(t == 2)
    {
     pos[1] = i;
     break;
    }
 }
 if(t == 2)
 {
  if(str.substr(0, pos[0]) == "NM")
    {
      nm = atoi(str.substr(pos[1] + 1).c_str());
    }
 }
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

void parseOneRead(string line, int & pos, int & end, int & strand, string & chr,
bool & isProperPair, int & nm, bool & is2ndRead, bool & isMapped)
{
 char * linePtr = new char [line.length()+ 10];
 int field = 1, flag;
 string cigar;

 strcpy(linePtr, line.c_str());
 char *f = strtok(linePtr, "\t");
 while(f != NULL)
 {
  if(field == 2)
    flag = atoi(f);
  if(field == 3)
    chr =  f;
  if(field == 4)
    pos = atoi(f);
  if(field == 6)
    cigar = f;
  if(field > 11)
    {
     parseMisMatch(f, nm);
    }
  f = strtok(NULL, "\t");
  field++;
 }
 parseSAMflag(flag, strand, isProperPair,  is2ndRead, isMapped);
 if(strand == 1) //reverse strand
 {
  end = pos;
  pos = coordReverseRead(pos, cigar);
 }
 else //forward strand
 {
  end = coordReverseRead(pos, cigar);
 }

 delete [] linePtr;
}

bool isUniqAndCigarMID_Read(string str)
{
 bool res = false;
 vector<int> fieldPos;

 for(int i = 0; i < str.length(); i++)
 {
  if(str[i] == '\t')
    fieldPos.push_back(i);
 }
 if(fieldPos.size() <= 10)
   return false;
 string option = str.substr(fieldPos[10], str.length());
 if(option.find("\tXT:A:U\t") == string::npos)
   return false;
 //a unique read
 //CIGAR string must be in field #6 (1-based), left is fP[4], right is fP[5]
 int start = fieldPos[4] + 1;
 int  n =  fieldPos[5] - start;
 string cigar = str.substr(start, n);
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

void populateRP(string line1, string line2,  vector<readPairKey> & rpKeys,  vector<readPairProperties> & rpProperties,
int cutoff)
{
 int pos[2], strand[2], end[2];
 string chr[2];
 bool isProperPair, is2ndRead;
 bool map[2];
 int NM[2];
 readPairKey rpk;
 readPairProperties rpp;

 if(isUniqAndCigarMID_Read(line1) && isUniqAndCigarMID_Read(line2))
 {
  parseOneRead(line1, pos[0], end[0], strand[0], chr[0], isProperPair, NM[0], is2ndRead, map[0]);
  parseOneRead(line2, pos[1], end[1], strand[1], chr[1], isProperPair, NM[1], is2ndRead, map[1]);
  if(NM[0] > cutoff || NM[1] > cutoff)
    return;
  if(map[0] == false ||  map[1] == false)
    return;
  if(is2ndRead) //line2 is 2nd read, line1 is 1st read
  {
   rpk.setMemberValues(chr[0],chr[1],  pos[0], pos[1],  strand[0], strand[1]); 
   rpp.read1 = line1;
   rpp.read1Pos[0] =  pos[0];
   rpp.read1Pos[1] =  end[0];
   rpp.read2 = line2;
   rpp.read2Pos[0] =  pos[1];
   rpp.read2Pos[1] =  end[1];
  }
  else
  {
   rpk.setMemberValues(chr[1], chr[0],  pos[1],  pos[0], strand[1], strand[0]);
   rpp.read1 = line2;
   rpp.read1Pos[0] =  pos[1];
   rpp.read1Pos[1] =  end[1];
   rpp.read2 = line1;
   rpp.read2Pos[0] =  pos[0];
   rpp.read2Pos[1] =  end[0];
  }
  rpKeys.push_back(rpk);
  rpp.isProperPair = isProperPair;
  rpp.sumNM = NM[0] + NM[1];
  rpProperties.push_back(rpp);
 }
}

void removeRedun( map<readPairKey, int, classcomp> & readlocLinenumber,const vector<readPairKey> &  rpKeys,vector<readPairProperties> & rpProperties)
{
 map<readPairKey, int, classcomp>::iterator iter;
 pair<map<readPairKey, int, classcomp>::iterator, bool> ret;
 int index, preMM, curMM;
 bool  prePP, curPP;

 if( rpKeys.size() == 0)
   return;
 vector<bool> readsToBeUsed;
 //in the beginning, the map is empty
 //after 1st, 2nd, 3rd, ..., scan, at least 1 element in the map. 
 //in the next scan, only focus on the reads already in the map
 if(readlocLinenumber.size() == 0) //at the beginning
 {
  for(int i = 0; i < rpKeys.size(); i++)
    readsToBeUsed.push_back(true);
 }
 else //only focus on the reads already in the map
 {
  for(int i = 0; i < rpKeys.size(); i++)
    readsToBeUsed.push_back(false);
  map<readPairKey, int, classcomp>::const_iterator iter;
  int index;

  for( iter = readlocLinenumber.begin(); iter != readlocLinenumber.end(); iter++)
  {
   index = (*iter).second;
   readsToBeUsed[index] = true;
  }
 }
 readlocLinenumber.clear(); 
 for(int i = 0; i < rpKeys.size(); i++)
 {
  if(readsToBeUsed[i])
  {
    ret = readlocLinenumber.insert(pair<readPairKey, int>(rpKeys[i], i) );
    if(ret.second == false)
    {
     index = (*(ret.first)).second;
     preMM = rpProperties[index].sumNM;
     prePP = rpProperties[index].isProperPair;
     curMM =  rpProperties[i].sumNM;
     curPP = rpProperties[i].isProperPair;
     if(prePP == curPP)
     {
      if(preMM > curMM)
      {
        (*(ret.first)).second = i;
      }
     }
     else
     {
      if(curPP)
        (*(ret.first)).second = i;
     }
    }
  }
 }
}

void removeRedunMultiScans(map<readPairKey, int, classcomp> & readlocLinenumber,vector<readPairKey> &  rpKeys,vector<readPairProperties> & rpProperties)
{
 int p1, p2;

 //1st read start, 2nd read start as the key
 removeRedun( readlocLinenumber,  rpKeys, rpProperties);
 cout << "1st scan completed\n";
 //1st read start, 2nd read end as the key
 for(int i = 0; i < rpKeys.size(); i++)
 {
  p1 = rpProperties[i].read1Pos[0];
  p2 = rpProperties[i].read2Pos[1];
  rpKeys[i].setTwoReadsPos(p1, p2); 
 }
 removeRedun( readlocLinenumber,  rpKeys, rpProperties);
 cout << "2nd scan completed\n";
 //1st read end, 2nd read start as the key
 for(int i = 0; i < rpKeys.size(); i++)
 {
  p1 = rpProperties[i].read1Pos[1];
  p2 = rpProperties[i].read2Pos[0];
  rpKeys[i].setTwoReadsPos(p1, p2);
 }
 removeRedun( readlocLinenumber,  rpKeys, rpProperties);
 cout << "3rd scan completed\n";
  //1st read end, 2nd read end as the key
 for(int i = 0; i < rpKeys.size(); i++)
 {
  p1 = rpProperties[i].read1Pos[1];
  p2 = rpProperties[i].read2Pos[1];
  rpKeys[i].setTwoReadsPos(p1, p2);
 }
 removeRedun( readlocLinenumber,  rpKeys, rpProperties);
 cout << "4th scan completed\n";
}

void write2file(ofstream & output, const  map<readPairKey, int, classcomp>  & readlocLinenumber,  const vector<readPairProperties> & rpProperties)
{
 map<readPairKey, int, classcomp>::const_iterator iter;
 int index;
 int normPair = 0;

 for( iter = readlocLinenumber.begin(); iter != readlocLinenumber.end(); iter++)
 {
  index = (*iter).second;
  if(rpProperties[index].isProperPair)
  {
   output <<  rpProperties[index].read1 << endl;
   output <<  rpProperties[index].read2 << endl;
   normPair++;
  }
 }
 for( iter = readlocLinenumber.begin(); iter != readlocLinenumber.end(); iter++)
 {
  index = (*iter).second;
  if(!rpProperties[index].isProperPair)
  {
   output <<  rpProperties[index].read1 << endl;
   output <<  rpProperties[index].read2 << endl;
  }
 }
 
 cout << "normPair: " << normPair << endl;
}
