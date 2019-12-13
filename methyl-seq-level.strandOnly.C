/*
$ head hg19/hg19.CpG.sites.list
chr	pos	base
chr10	60025	CG
chr10	60089	CG
chr10	60109	CG

pileup:
chr1	9999	N	1	^~c	?
chr1	10000	N	9	c^~C^~C^~A^~A^~A^~A^~A^~A	8C>>&!>5!
chr1	10001	t	11	,.......C^~.^~.	1C?>&!<5!%/
chr1	10002	a	16	,...GT..G..^~T^~T^~.^~G^~.	?CCC&!?5!%/!!!$B
chr1	10003	a	17	,........TCTTG..^~.	9FFF0%D6+%/!!!0C!
chr1	10004	c	19	,......G.....T..A^~,^~,	DFFF=CD37=899!@C!A?
 
output - a list of methylation levels ordered by CpG sites in the hg19.CpG.sites.list file
T-cell-EM
NA (not available)
0.95
0.95
.9
0.8
0
1.1

Algorithm:
If sufficient read coverage
  how many T for C site position, 
  how many A for the next G site position
  unM = T/covC + A/covG (70% + 30%, unMethylated, C>T in forward strand, G>A in reverse strand,  could be > 1)
     30 read cover, 20 reads in forw, 10 reads in rev
     60% unmethyl (C>T, G>A), 12/20 T, 6/10 A, 12/30 + 6/30 = 18/30 = 60%
     a/b = c/d = x
     (a+c)/(b+d) = (bx + dx) / (b+d) = x(b+d) / (b+d) = x 
   met = 1 - unM (could be < 0)

read next line in pileup file
if the position is CpG site (cgSiteMap)
  calc T/covC
  read next line (G position)
  calc A/covG  
  calc met_level
  updateMethyList( chr, posAs-a-string, cgSiteMap, met_level)

struct cgMethyl
{
 string chr;
 int pos;
 string nts;
 bool hasAvalue;
 double methyl;
};

vector<cgMethyl> glob_cgMethyl_list;

map<string, int> cgSiteMap; //basically determine if the position is a CpG site based on the hg19.CpG.sites.list    <int> refers to the index in the glob_cgMethyl_list
struct param
{
 int coverage;
 int indBaseCover;
 double indBasePerc;
};

void updateMethyList( string chr, string pos, map<string, int> & indexMap, double m)
{
 string str = chr + "\t\t\t" + pos;
 int index = indexMap[str];
 glob_cgMethyl_list[index].hasAvalue = TRUE;
 glob_cgMethyl_list[i].methyl = m;
 return;
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

struct pileInfo
{
 int a, c, g, t;
 int ins, del;
};

struct cgMethyl
{
 string chr;
 int pos;
 string nts; //CG, cg
 bool hasAvalue;
 double methyl;
 int numOfCorG, total;
 pileInfo pI[2]; //C site, pI[0], A/C/G/T/ins/del read number; G site - pI[1] ......
};

struct pileup
{
 string chr;
 int pos;
 char base;
 int cover;
 string readStr;
 string quality;
 void clear();
};

struct param
{
 int coverage;
 int indBaseCover;
 double indBasePerc;
 bool forward, reverse;
};


vector<cgMethyl> glob_cgMethyl_list;

void updateMethyList( string key, map<string, int> & indexMap, double m, const pileInfo & piC, const pileInfo & piG)
{
 if(indexMap.count(key) > 0)
 {
  int index = indexMap[key];
  glob_cgMethyl_list[index].hasAvalue = true;
  glob_cgMethyl_list[index].pI[0] = piC; //piC assignment, six integers
  glob_cgMethyl_list[index].pI[1] = piG;
  return;
 }
}

bool line2pile(string line, pileup & pile,  string & str);
int parseInDelStr(string mapStr, int i, string & ins);  //+20TTTTTTTTTTTT....TG; -3AAC
void pile2ACGTcount(const pileup & pile, pileInfo & pI, int & cov);
bool calcMethyl(const pileup & pile_C, const pileup & pile_G, const param & par, double & meth, pileInfo & pIc, pileInfo & pIg);

int main(int argc, char * argv[])
{
 if(argc != 7)
 {
  cerr << argv[0] << "    CpG_site_list_file_1   output-file-name_2 pileup-file_3  sample-name_4(T-cell-EM)   coverage_cutoff_5(10)  forward/reverse()\n";
  return 1;
 }

 map<string, int> cgSiteMap; //basically determine if a chromosome position is a CpG site based on the hg19.CpG.sites.list    <int> refers to the index in the glob_cgMethyl_list for a quick access
 long long int lineNum = 0;
 ifstream input;
 ofstream output;
 param par;
 string line, baseRep;
 pileup pile, pileG;

 par.coverage = atoi(argv[5]);
 if(string(argv[6]) == "forward" || string(argv[6]) == "Forward")
 {
  par.forward = true;
  par.reverse = false;
 }  
 else
 {
  par.forward = false;
  par.reverse = true;
 }

 input.open(argv[1], ios::in);
 if(!input)
   {
    cerr << argv[1] << " cannot be opened for reading!" << endl;
    return 1;
   }
 int numTab = 0;
 vector<int> tabPos;
 int len, coverage;
 char refBase;
 string mapStr, str, str2;
 cgMethyl metStru;

 getline(input, line);
 if(line == "chr\tpos\tbase")
   getline(input, line);
 else
   {
    cerr << "Cannot recognize the header line in the " << argv[1] << endl;
    cerr << "chr\tpos\tbase" << endl;
    return 1;
   }
 while(!input.eof())
 {
  numTab = 0;
  tabPos.clear();
  for(int i = 0; i < line.length(); i++)
  {
   if(line[i] == '\t')
   {
    tabPos.push_back(i);
   }
  }
  if(tabPos.size() == 2)
  {
   len = tabPos[0] ;
   metStru.chr = line.substr(0, len);
   //chr position 2nd field, between 1st and 2nd tab
   len = tabPos[1] - tabPos[0] - 1;
   str = line.substr(tabPos[0] + 1, len);
   str2 = metStru.chr + "\t\t\t" + str;
   metStru.pos = atoi(str.c_str());
   metStru.nts = line.substr(tabPos[1] + 1, 2);
   metStru.hasAvalue = false;
   metStru.methyl = -999999;
   glob_cgMethyl_list.push_back(metStru);
   cgSiteMap[str2] = glob_cgMethyl_list.size() - 1; //size 10, index 9, glob_cgMethyl_list[9] 
  }
  getline(input, line);
 }
 input.close();

 bool b, t;
 double meth = -999999;
 pileInfo pIc, pIg;
 string chrPosStr;
 input.open(argv[3], ios::in);
 if(!input)
   {
    cerr << argv[3] << " pile-up-file cannot be opened for reading!" << endl;
    return 1;
   }
 getline(input, line);
 while(!input.eof())
 {
  if(lineNum % 10000000 == 0)
    cout <<  argv[3] << " pile-up-file line: " << lineNum << endl;
  lineNum++;
  b = line2pile(line, pile, str);
  if(b) // a genuine pileup line
  {
   chrPosStr = pile.chr + "\t\t\t" + str;
   if(cgSiteMap.count(chrPosStr) > 0) //a CpG site
   {
    getline(input, line);
    line2pile(line, pileG, str2);
    if(pileG.chr == pile.chr && pileG.pos - 1 == pile.pos) //pileup file has both C & G line;
    {
     t = calcMethyl(pile,pileG,par, meth, pIc, pIg);
     if(t) // sufficient coverage to calc methy level
     {
      updateMethyList(chrPosStr, cgSiteMap, meth, pIc, pIg);
     }
    }
   } //a CpG site
  } // a genuine pileup line
  getline(input, line);
 }
 input.close();

 output.open(argv[2], ios::out);
 if(!output)
 {
  cerr << argv[2] << " cannot be written to.\n";
  return 1;
 }
 output <<   argv[4] << "\t" << argv[6] << endl;
 output << "C-site A/C/G/T/ins/del\tG-site A/C/G/T/ins/del" << endl;
 output << "A\tC\tG\tT\tIns\tDel\tA\tC\tG\tT\tIns\tDel\n";
 {
  for(int j = 0; j < glob_cgMethyl_list.size(); j++)
  {
   if(glob_cgMethyl_list[j].hasAvalue == true)
   {
    for(int i = 0; i < 2; i++)
    {
     output << glob_cgMethyl_list[j].pI[i].a << "\t";
     output << glob_cgMethyl_list[j].pI[i].c << "\t";
     output << glob_cgMethyl_list[j].pI[i].g << "\t";
     output << glob_cgMethyl_list[j].pI[i].t << "\t";
     output << glob_cgMethyl_list[j].pI[i].ins << "\t";
     output << glob_cgMethyl_list[j].pI[i].del;
     if(i == 0)
       output << "\t";
     else
       output << "\n";
    }
   }
   else
   {
    for(int i = 1; i <= 11; i++)
       output << "0\t";
    output << "0\n";
   }
  }
 }

 output.close();

 return 0; //end of main()
}

bool line2pile(string line, pileup & pile, string & strPos)
{
 int numTab = 0;
 vector<int> tabPos;

 for(int i = 0; i < line.length(); i++)
 {
  if(line[i] == '\t')
  {
   tabPos.push_back(i);
  }
 }
 if(tabPos.size() < 5)
   return false;
 int len, coverage;
 char refBase;
 string mapStr, str;

 //chrosome 1st field
 len = tabPos[0] ;
 pile.chr = line.substr(0, len);
 //chr position 2nd field, between 1st and 2nd tab
 len = tabPos[1] - tabPos[0] - 1;
 str = line.substr(tabPos[0] + 1, len);
 strPos = str;
 pile.pos = atoi(str.c_str());
 //ref base is following the 2nd tab
 refBase = line[tabPos[1] + 1 ];
 pile.base = refBase;
 //coverage between 3rd tab and 4th tab
 len = tabPos[3] - tabPos[2] - 1;
 str = line.substr(tabPos[2] + 1, len);
 coverage = atoi(str.c_str());
 pile.cover = coverage;
 //length should be 4th tab and 5th tab
 len = tabPos[4] - tabPos[3] - 1;
 mapStr = line.substr(tabPos[3] + 1, len);
 pile.readStr = mapStr;
 return true;
}

bool calcMethyl(const pileup & pile_C, const pileup & pile_G, const param & par, double & meth, pileInfo & pIc, pileInfo & pIg)
{
 bool hasMethValue = false;
 int a,c,g,t,cov, a2,c2,g2,t2,cov2;

 pile2ACGTcount(pile_C, pIc,cov);
 pile2ACGTcount(pile_G, pIg,cov2);
 /*if(par.forward == true)
 {
  numOfC = c;
  cAndT = c + t;
 }
 else
 {
  numOfC = g2;
  cAndT = g2 + a2;
 }
 */
/* cov = a+c+g+t;
 cov2 = a2+c2+g2+t2;
 if(cov >= par.coverage && cov2 >= par.coverage)
 {
  double unm = (double)t / cov + (double) a2 / cov2;
  meth = 1.0 - unm;
  hasMethValue = true;
 } */
 hasMethValue = true;
 return hasMethValue; 
}

void pile2ACGTcount(const pileup & pile, pileInfo & pI, int & cov)
{
 char letter, refBase;
 string mapStr = "", insStr = "", acg = "ACGT";
 vector<string> insList;
 int numIn = 0, numDel = 0;
 int numMatches = 0, mutCover = pile.cover, insCover = pile.cover, delCover = pile.cover;
 map<char, int> baseCount;
 map<string, int> insDormMap;
 bool res = false;

 for(int i = 0; i < acg.length(); i++)
   baseCount[acg[i]] = 0;
 
 refBase = toupper(pile.base);
 
 for(int i = 0; i < pile.readStr.length(); i++)
   mapStr += toupper(pile.readStr[i]);

 for(int i = 0; i < mapStr.length(); i++)
 {
  letter = mapStr[i];
  if(letter == '^')
  {
   i++; //skip the char following ^
   delCover--;
  }
  else if (letter == '$')
  {
   insCover--;
   delCover--;
  }
  else if (letter == '+' ) //+20TTTTTTTTTTTT....TG
  {
   i = parseInDelStr(mapStr, i, insStr); //i points to the last 'agctn'
   insList.push_back(insStr);
   numIn++;
  }
  else if (letter == '-' )
  {
   i = parseInDelStr(mapStr, i, insStr); //i points to the last 'agctn'
  }
  else if (letter == '*' )
  {
   numDel++;
  }
  else if (letter == ',' || letter == '.')
  {
   numMatches++;
  }
  else if (letter == 'N' || letter == 'n')
  {
   mutCover--;
  }
  else if(acg.find(letter) != string::npos) //ACGT acgt
  {
   baseCount[letter]++;
  }
  else; //other symbol
 }
 baseCount[refBase] = numMatches;
 pI.a = baseCount['A'];
 pI.c = baseCount['C'];
 pI.g = baseCount['G'];
 pI.t = baseCount['T'];
 pI.ins = numIn;
 pI.del = numDel;
}

int parseInDelStr(string mapStr, int i, string & ins)  //+20TTTTTTTTTTTT....TG; -3AAC
{
   int begin ;
   i++; //the first digit
   string strIn(1, mapStr[i]);
   i++;
   while(isdigit(mapStr[i]))
   {
    strIn += mapStr[i];
    i++;
   }
   //exit while loop, i points to the first 'acgtn'
   begin = i;
   int numIn = atoi(strIn.c_str());
   i = i - 1 + numIn; //i points to the last 'agctn'
   int length = i + 1 - begin;
   ins = mapStr.substr(begin, length);
   return i;
}

