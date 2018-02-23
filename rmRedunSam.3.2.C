/*
the program selects uniquley mapped reads, removes redundent reads that are mapped to the same genome location.

Given SAM output by BWA single-end mapping

@SQ     SN:SL2.40ch11   LN:53386025
@SQ     SN:SL2.40ch12   LN:65486253
@PG     ID:bwa  PN:bwa  VN:0.5.9-r16
HWUSI-EAS517_0047:1:1:1445:950#0        16      SL2.40ch09      124733  0       101M    *       0       0       TCTTGATGGACGTCCACGAAAAAATTTGGCGTTTTTGAAGTCGGAATCCGGATCACCCAAAAAATCATGTGCTATAGCACACGAAAATCGTCGAAATGAGN   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB_QQ__VVTPVVVTVTNHLKNT[VY[______Y___NNMNNJLOKD   XT:A:R  NM:i:4  X0:i:3  X1:i:0  XM:i:4  XO:i:0  XG:i:0  MD:Z:30A7C1C59G0        XA:Z:SL2.40ch09,-116875,101M,4;SL2.40ch09,-118319,101M,4;
HWUSI-EAS517_0047:1:1:1564:935#0        16      SL2.40ch08      21741247        37      101M    *       0       0       CACCAATCAAGAATCCCCAAGGAATTCCAAGAATTCTTGGGATACATTTGTTCCAACCCTCAACCCTCTTTTCTAGATTCATGGATATTGAATCAATCGAN   BBBBBBBBBBBBBBBBB_________[[Y[[__V__YYYYYTVVVOQWVXW[[[[[WWXRXTRRRMQQQ[[________________W___MMNQMHIMIE   XT:A:U  NM:i:4  X0:i:1  X1:i:0  XM:i:4  XO:i:0  XG:i:0  MD:Z:40T33T14C10A0


Keep the headers.
The header line 1st field @SQ or one line after @SQ

For each read line, only the line with "\tXT:A:U\t" will be kept


If multiple reads mapped to the same genome location:
1. same chromosome
2. same starting coordinate
3. same strand
then, only one read with the mimimal mismatches will be kept.Read A and B are both uniquely mapped. If they are the same strand, then both A and B have identical SAM flag field value (0 or 16).  

create map<string, int> objects,
key is "chr-position-strand", the value is the index of uniqReadLines

If the key is not in the map, insert (key, value) pair
Otherwise, compare the mismatch of the read with the new one. If the mismatch decreases, update the map element.

Output header lines and non-redundant uniquely mapped reads 

v2.0
How to determine the genome coordinate for the starting base of a read given SAM output?
Method:
Select uniquely mapped reads & the reads with CIGAR string that contains only 'MID'. Summation of M and I operations should be 101 for each qualified read (has been checked for 296117 reads).
Thus, for a forward read, the coordinate reported by SAM corresponds to the starting base of a read;  
For a reverse strand, the genome coordinate for the starting base is given by
SAM_coord - 1 + numOf_M_operations + numOf_D_operations 

v3.0
Read with NM < 5 (hard coded) will be kept
ONLY reads with identical end positions will be deemed as duplicates.

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
#include <sstream>

using namespace std;

bool isSAMheader(char *f);
void populateReadLines(string line, vector<string> & headLines, vector<string> & uniqReadLines, vector<string> & chrPosStrand,  vector<int> & uniqReadMismatch, bool & firstRead);

void removeRedun( map<string, int> & readlocLinenumber, const vector<string> & chrPosStrand,  const vector<int> & uniqReadMismatch);

void write2file(ofstream & output, const vector<string> & headLines, const  map<string, int> & readlocLinenumber, const  vector<string> & uniqReadLines);

vector<int> gl_readsGaps;
int glo_maxNM;

int main(int argc, char *argv[])
{
 if(argc != 4)
 {
  cout << "executable SAM-input-file  output-file-name max_NM" << endl;
  return 1;
 }

 vector<string> headLines, uniqReadLines, chrPosStrand;
 vector<int> uniqReadMismatch;
 map<string, int> readlocLinenumber;
 map<string, int>::iterator mapIt;
 string line;
 bool firstRead = false;
 ifstream input;
 ofstream output;


 input.open(argv[1], ios::in);
 if(!input)
   {
    cerr << argv[1] << "SAM file cannot be opened for reading!" << endl;
    return 1;
   }
 glo_maxNM = atoi(argv[3]);
 glo_maxNM++;
 getline(input, line);
 while(!input.eof())
 {
  populateReadLines(line,  headLines, uniqReadLines, chrPosStrand,  uniqReadMismatch, firstRead);
  getline(input, line);
 }
 input.close();
 removeRedun( readlocLinenumber,  chrPosStrand,   uniqReadMismatch);

 output.open(argv[2], ios::out);
 if(!output)
   {
    cerr << argv[2] << " cannot be opened for writing" << endl;
    return 1;
   }
 write2file(output, headLines, readlocLinenumber, uniqReadLines);
 output.close();

 cout << "rmRedunSam.3.C\n";
 cout << "Only retain reads with NM < " << glo_maxNM  << "\n";
 cout << "Only reads with identical end pos are deemed duplicated \n";
 return 0;
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

bool isSAMheader(char *f)
{
 bool res = false;
 if(strlen(f) == 3)
   if(f[0] == '@')
     if(isalpha(f[1]) && isalpha(f[2]))
       res = true;

 return res;
}

void extractFields(string str, vector<string> & list, char del)
{
 vector<int> fieldPos;
 int start, n;

 for(int i = 0; i < str.length(); i++)
 {
  if(str[i] == del)
    fieldPos.push_back(i);
 }
 for(int i = 0; i <= fieldPos.size(); i++)
 {
  if(i == 0)
  {
   start = 0;
   n = fieldPos[0];
  }
  else if(i == fieldPos.size())
  {
   start = fieldPos[i-1] + 1;
   n = str.length()- start;
  }
  else
  {
   start = fieldPos[i-1] + 1;
   n =  fieldPos[i] - start;
  }
  list.push_back(str.substr(start, n));
 }
}

bool isUniqRead(string str)
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
 if(option.find("\tXT:A:U\t") != string::npos)
   return true;
 
 return res;
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

int calcNumGaps(string cigar)
{
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
   if(cigar[i] == 'I' || cigar[i] == 'D')
     length += atoi(numOfOp.c_str());
   numOfOp = "";
  }
 }
 return length;
}

string coordReverseRead(string pos, string cigar)
{
 int end = atoi(pos.c_str());
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
 stringstream ss;
 ss << coord;
 //cout << cigar << "\t" << length << "\t" << end << "\t" << ss.str() << endl;
 return ss.str();
}


void populateReadLines(string line, vector<string> & headLines, vector<string> & uniqReadLines, vector<string> & chrPosStrand,  vector<int> & uniqReadMismatch, bool & firstRead)
{
 char * linePtr = new char [line.length()+ 10];
 strcpy(linePtr, line.c_str());
 char *f = strtok(linePtr, "\t");
 
 if(f != NULL)
 {
  if(isSAMheader(f) && firstRead == false)
  {
   headLines.push_back(line);
  }
  else 
  {
   firstRead = true;
  }

  if(firstRead)
  {
   if(isUniqAndCigarMID_Read(line))
    {
     int field = 1, nm;
     string temp, pos, chr, strand, cigar;
     while(f != NULL)
     {
      if(field == 2)
        strand = f;
      if(field == 3)
        chr =  f;
      if(field == 4)
        pos = f;
      if(field == 6)
        cigar = f;
      if(field > 11)
        {
         parseMisMatch(f, nm);
        }
      f = strtok(NULL, "\t");
      field++;
     }
     //forward read, calc. end position; reverse read, pos is already end pos 
     if(strand == "0")
     {
       pos = coordReverseRead(pos, cigar);
     }
     temp = chr  + "-" + pos + "-" + strand;
     if(nm < glo_maxNM)
     {
      uniqReadLines.push_back(line);
      chrPosStrand.push_back(temp);
      gl_readsGaps.push_back(calcNumGaps(cigar));
      uniqReadMismatch.push_back(nm);
     }
    }
  }
 }
 delete [] linePtr;
 return;
}

void removeRedun( map<string, int> & readlocLinenumber, const vector<string> & chrPosStrand,  const vector<int> & uniqReadMismatch)
{
 map<string, int>::iterator iter;
 pair<map<string, int>::iterator, bool> ret;
 int index, preMM, curMM, prevGap, currGap;

 for(int i = 0; i < chrPosStrand.size(); i++)
 {
  ret = readlocLinenumber.insert(pair<string, int>(chrPosStrand[i], i) );
  if(ret.second == false)
  {
   index = (*(ret.first)).second;
   preMM = uniqReadMismatch[index];
   curMM = uniqReadMismatch[i];
   prevGap = gl_readsGaps[index];
   currGap = gl_readsGaps[i];
   if(preMM > curMM)
   {
     (*(ret.first)).second = i;
   }
   else
   {
    if(preMM == curMM)
    {
     if(prevGap > currGap) //NM same, currGap smaller
       (*(ret.first)).second = i;
    }
   }
  }
 }
}


void write2file(ofstream & output, const vector<string> & headLines, const  map<string, int> & readlocLinenumber, const  vector<string> & uniqReadLines)
{
 for(int i = 0; i < headLines.size(); i++)
 {
  output <<  headLines[i] << endl;
 }

 int index;
 map<string, int>::const_iterator iter;

 for( iter = readlocLinenumber.begin(); iter != readlocLinenumber.end(); iter++)
 {
  index = (*iter).second;
  output <<  uniqReadLines[index] << endl;
 }
}
