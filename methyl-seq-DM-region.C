/*

input file:
-----------------
chr	pos	base	SCMRA	SCMRARO
chr10	73522	cg	0.565217	0.565217
chr10	73536	cg	-0.761905	0.428571
chr10	73555	cg	0.85	0.333333
chr10	74103	cg	NA	0.794872
chr10	74109	cg	0.185185	-0.025641
chr10	74132	cg	0.130435	NA
chr10	78342	CG	0.608696	0.448718
chr10	89821	CG	0.780059	0.632997
chr10	89867	CG	0.891892	0.9
chr10	92623	CG	0.916667	0.776515
chr10	93373	CG	0.961538	0.810811
chr10	93393	CG	0.677419	0.659574
chr10	93396	CG	0.882353	0.90566


> info <-c("chr1", 12345, 7, "less")  ## consecu 7 sites
> c1 <- c(0.1, 0.1, 0.2, 0.2, NA, 0.3, 0.3)
> c2 <- c(0.2, 0.2, 0.21, 0.2, 1, 0.4, 0.4)
> res <- t.test(c1, c2, paired = T, alternative = "less") #### “greater” 
> res

	Paired t-test

data:  c1 and c2
t = -3.4049, df = 5, p-value = 0.009575
alternative hypothesis: true difference in means is less than 0
95 percent confidence interval:
        -Inf -0.02789267
sample estimates:
mean of the differences 
            -0.06833333 

> res$p.value
[1] 0.009574844


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

struct cgMethyl
{
 string chr;
 int pos;
 string nts;
 bool hasAvalue;
 bool hasBeenUsed;
 double methA, methB, methylDiff; //0.3 0.5 0.3-0.5=-0.2
};

struct param
{
 int minNumSites ; // 6
 double mostNeg ; // -0.1, if meth < -0.1, then treated as NA
 ofstream opLess, opGreater;
};

map<string, vector<cgMethyl> > globMapChrCpGsites; // "chr1" -> a vector of CpG sites in "chr1"
//vector<string> globVectChrs; //
set<string> globSetChrs;
param globParam;

void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);
bool findDMR();
bool write2DMR(string chr, int begin, int end, string comp);
bool setglobMapChrCpGsites(vector<string>  & lineFields);

int main(int argc, char *argv[])
{
 if(argc != 6)
 {
  cerr << argv[0] << "	methyl-28.7million-sites-input(1)	output-greater-file(sampA > sampB, 2)";
  cerr << "	output-less-file(3)	minNumOfConsecuSites(e.g. 6, 4)	mostNegMeth(e.g. -0.1, if meth < -0.1, then treated as NA, 5)";
  cerr << "	\n";
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 cgMethyl myCG;

 globParam.minNumSites = atoi(argv[4]);
 globParam.mostNeg  = atof(argv[5]);
 writeFile(argv[2], globParam.opGreater);
 writeFile(argv[3], globParam.opLess);

 readFile(argv[1], input);
 getline(input, line);
 if(line.find("chr	pos	base") == string::npos) //tab header-line
 {
  cerr << argv[1] << " has no header line\n";
  cerr << "chr	pos	base	TN	TEM	SCMRA	SCMRARO" << endl;
  return 1;
 } 
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields,'\t', line);
  setglobMapChrCpGsites(lineFields);
  getline(input, line);
 }
 input.close();
 findDMR();
 globParam.opGreater.close();
 globParam.opLess.close();
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

bool findDMR()
{
 string chr;
 int numCpGsitesChr;
 for(set<string>::iterator it = globSetChrs.begin(); it != globSetChrs.end(); ++it)
 {
  chr = *it ;
  numCpGsitesChr = globMapChrCpGsites[chr].size();
  for(int j = 0; j < numCpGsitesChr;  )
  {
   int begin, end;
   bool bo;
 
   if(globMapChrCpGsites[chr][j].hasAvalue == false)
     j++;
   else if(globMapChrCpGsites[chr][j].methylDiff == 0)
     j++;
   else //#451
   {
    if(globMapChrCpGsites[chr][j].methylDiff > 0)
    {
     begin = j;
     bo = true;
     while(bo)
     {
      j++;
      if(j == numCpGsitesChr || globMapChrCpGsites[chr][j].hasAvalue == false
         || globMapChrCpGsites[chr][j].methylDiff <= 0)
      {
       end = j - 1;
       bo = false;
      }
     }
     if(end - begin + 1 >= globParam.minNumSites)
     {
      write2DMR(chr, begin, end, "greater");
     }
    }
    else if(globMapChrCpGsites[chr][j].methylDiff < 0)
    {
     begin = j;
     bo = true;
     while(bo)
     {
      j++;
      if(j == numCpGsitesChr || globMapChrCpGsites[chr][j].hasAvalue == false
         || globMapChrCpGsites[chr][j].methylDiff >= 0)
      {
       end = j - 1;
       bo = false;
      }
     }
     if(end - begin + 1 >= globParam.minNumSites)
     {
      write2DMR(chr, begin, end, "less");
     }
    }
    else
    {
     cerr <<"sth wrong(510): " << chr << " " << globMapChrCpGsites[chr][j].pos << endl;
     exit(1);
    }
   //update j
   } //END - else #451  globMapChrCpGsites[chr][j].methylDiff != 0
  } //END - for(int j = 0; j < numCpGsitesChr;  )
 } //END - for(set<string>::iterator it = globSetChrs.begin(); it != globSetChrs.end()
 return true;
}

bool write2DMR(string chr, int begin, int end, string comp)
{
 if(comp == "greater")
 {
  globParam.opGreater << "c(\"" << chr << "\", " << globMapChrCpGsites[chr][begin].pos << ", " << globMapChrCpGsites[chr][end].pos << ", " << end - begin + 1 << ", " << "\"greater\")" << endl; 
  globParam.opGreater << "c(";
  for(int j = begin; j <= end - 1; j++)
  {
   globParam.opGreater << globMapChrCpGsites[chr][j].methA << ", ";
  } 
  globParam.opGreater << globMapChrCpGsites[chr][end].methA << ")" << endl;
  globParam.opGreater << "c(";
  for(int j = begin; j <= end - 1; j++)
  {
   globParam.opGreater << globMapChrCpGsites[chr][j].methB << ", ";
  }
  globParam.opGreater << globMapChrCpGsites[chr][end].methB << ")" << endl;
 }
 else if(comp == "less")
 {
  globParam.opLess << "c(\"" << chr << "\", " << globMapChrCpGsites[chr][begin].pos << ", " << globMapChrCpGsites[chr][end].pos << ", " << end - begin + 1 << ", " << "\"less\")" << endl;
  globParam.opLess << "c(";
  for(int j = begin; j <= end - 1; j++)
  {
   globParam.opLess << globMapChrCpGsites[chr][j].methA << ", ";
  }
  globParam.opLess << globMapChrCpGsites[chr][end].methA << ")" << endl;
  globParam.opLess << "c(";
  for(int j = begin; j <= end - 1; j++)
  {
   globParam.opLess << globMapChrCpGsites[chr][j].methB << ", ";
  }
  globParam.opLess << globMapChrCpGsites[chr][end].methB << ")" << endl;
 }
 else;
 return true;
}

bool setglobMapChrCpGsites(vector<string>  & lineFields)
{
 cgMethyl cg;
 string chr;
 //chr10   73536   cg      -0.761905       0.428571
 if(lineFields.size() >= 5)
 {
  chr = lineFields[0];
  cg.chr = chr;
  cg.pos = atoi(lineFields[1].c_str());
  cg.nts = lineFields[2];
  if(lineFields[3].find('N') != string::npos || lineFields[4].find('N') != string::npos ||
     lineFields[3].find('n') != string::npos || lineFields[4].find('n') != string::npos )
  {
   cg.hasAvalue = false;
  }
  else if(atof(lineFields[3].c_str()) < globParam.mostNeg || atof(lineFields[4].c_str()) < globParam.mostNeg )
  { 
   cg.hasAvalue = false;
  }
  else
  {
   cg.hasAvalue = true;
   cg.methA = atof(lineFields[3].c_str());
   cg.methB = atof(lineFields[4].c_str());
   cg.methylDiff = cg.methA - cg.methB;
  }
  cg.hasBeenUsed = false;
  globMapChrCpGsites[chr].push_back(cg);
  globSetChrs.insert(chr);
 }
 return true;
}
