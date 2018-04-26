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

/*
Chr1    TAIR9   chromosome      1       30427671        .       .       .       ID=Chr1;Name=Chr1
Chr1    TAIR9   gene    3631    5899    .       +       .       ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
 */

struct gff
{
 string chr;
 string type;
 string desc;
 char strand;
 int left;
 int right;

 bool setValues(string line);
 bool operator< (const gff &) const;
};

struct snp
{
 string chr;
 int pos;
 vector<string> fields;
 bool setValues(string line);
};

class compare
{
 public: 
 bool operator() (const gff & lh, const gff & rh) const
 {
  return lh < rh;
 }
};

void getFieldContent(vector<string> & fields, char del, string line);
void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
string parseSnpGff(const snp & mySnp, const map<gff, vector<gff>, compare> & gffMap, int dist, string line);
string snpType(const vector<gff> & gffVect, const string & prom);
string int2str(int n);

int glob_param_prom;
vector<gff> glob_gff_vect;

int main (int argc, char *argv[])
{
 if(argc != 6)
 {
  cerr << argv[0] << " snp-file  gff-file output-file  distance(1000000)  upstream(1000)" << endl;
  cerr << "distance: if the left of a gff element to the snp position is greater than distance, then snp is not within the element." << endl;
  cerr << argv[5] << " bp upstream of a gene start as promoter" << endl;
  return 1;
 }
 
 int dist = atoi(argv[4]);
 glob_param_prom = atoi(argv[5]);
 ifstream input;
 ofstream output;
 string line, res;
 gff myGff;
 snp mySnp;
 map<gff, vector<gff>, compare > gffMap;

 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  if(myGff.setValues(line) == false) return 1;
  gffMap[myGff].push_back(myGff);
  glob_gff_vect.push_back(myGff);
  getline(input, line);
 }
 input.close(); 

 writeFile(argv[3], output);
 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  if(mySnp.setValues(line)== false) return false;
  res = parseSnpGff(mySnp, gffMap, dist, line);
  output << res;
  getline(input, line);
 }
 input.close();
 output.close();

 return 0;
}

bool gff::setValues(string line)
{
 bool res = false;
 vector<string> fields;

 getFieldContent(fields, '\t', line);
 if(fields.size() != 9)
 {
  cerr << "**********************************************" << endl;
  cerr << "the following line in gff file has fewer or larger than 9 field" << endl;
  cerr << line << endl;
  return false;
 }
//scaffold_1      phytozomev10    gene    3909585 3920476    .       +       .       ID=Si016125m.g.version2.1;Name=Si016125m.g
 chr = fields[0];
 type = fields[2];
 desc = fields[8];
 left = atoi(fields[3].c_str());
 right = atoi(fields[4].c_str());
 strand = fields[6][0];
 return true;
}

bool snp::setValues(string line)
{
 bool res = false;
 vector<string> fields;

 getFieldContent(fields, '\t', line);
 if(fields.size() != 9)
 {
  cerr << "**********************************************" << endl;
  cerr << "the following line in snp file has fewer or larger than 10 fields" << endl;
  cerr << line << endl;
  return false;
 }
//M       T->G;   Chr1    892     T       30      ,$..,,..,c,.,.,,,,,,,,,.,,..,., 10      gggggggggg
 chr = fields[2];
 pos = atoi(fields[3].c_str());
 this->fields.clear();
 for(int i = 0; i < fields.size(); i++)
 {
  this->fields.push_back(fields[i]);
 }
 //if(this->fields[0] == "MD") this->fields[0] = "D"; 
 return true;
}

bool gff::operator< (const gff & r) const
{
 if(chr < r.chr)
   return true;
 if(chr > r.chr)
   return false;
 if(left < r.left)
   return true;
 if(left > r.left)
   return false;
 //same chr, same left
 return false;
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

void readFile(const char * argv, ifstream & input)
{
 input.open(argv, ios::in);
 if(!input)
 {
  cerr << argv << " can not be opened for reading.\n";
  exit(1);
 }
}

string isInPromoter(const snp & mySnp)
{
 string type = "", res = "";
 stringstream ss;

 for(int i = 0; i < glob_gff_vect.size(); i++)
 {
  if( glob_gff_vect[i].chr == mySnp.chr)
  {
   type = "";
   for(int j = 0; j < glob_gff_vect[i].type.length(); j++)
     type = type + (char) toupper(glob_gff_vect[i].type[j]);
   if(type.find("GENE") !=string::npos)
   {
    if(glob_gff_vect[i].strand == '+' && glob_gff_vect[i].left > mySnp.pos && glob_gff_vect[i].left - mySnp.pos <= glob_param_prom )
    {
     ss << mySnp.pos - glob_gff_vect[i].left;
     res = string("promoter of the gene; ") + ss.str() + " bp;" +  glob_gff_vect[i].desc;
     return res;
    }  
    else if (glob_gff_vect[i].strand == '-' && glob_gff_vect[i].right < mySnp.pos && mySnp.pos - glob_gff_vect[i].right <= glob_param_prom )
    {
     ss << glob_gff_vect[i].right -  mySnp.pos;
     res = string("promoter of the gene; ") + ss.str() + " bp;" +  glob_gff_vect[i].desc;
     return res;
    }
    else; 
   }
  }
 } // next gff element
 return res;
}

string parseSnpGff(const snp & mySnp, const map<gff, vector<gff>, compare> & gffMap, int dist, string line)
{
 string res, annot = "", prom = "";
 vector<gff> gffVect;
 map<gff, vector<gff>, compare>::const_iterator iter;
 gff g;
 
 g.chr = mySnp.chr;
 g.left = mySnp.pos;
 iter = gffMap.upper_bound(g);
 //iter: the first gff elem. greater than g.
 //iter--: last elem <= g. 
 // 1) same chr, left pos of elem <=  mySnp.pos
 // 2) diff chr which is < mySnp.chr 
 if(iter == gffMap.begin())
 {// no gff elem can contain snp
 }
 else
 {
  while(true)
  {
   iter--;
   string c = iter->first.chr;
   int left = iter->first.left;
   if(c == mySnp.chr && mySnp.pos - left <= dist)
   {
    for(int i = 0; i < iter->second.size(); i++)
    {
     if(iter->second[i].right >= mySnp.pos) //same chr, left <= mySnp.pos <= right
       gffVect.push_back(iter->second[i]);
    }
   }
   else //diff chr OR same chr but too far away from snp.Pos
     break; //no need to iter-- again
   
   if(iter == gffMap.begin())
     break;
  } //while(true)
 } //if(iter == gffMap.begin())   else

 prom = isInPromoter(mySnp);
 res =  snpType(gffVect, prom);
 res =  mySnp.fields[0] + '\t' + mySnp.fields[1] + '\t' + mySnp.fields[2] + '\t'
        + mySnp.fields[3] + '\t'
        + mySnp.fields[4]  + '\t' + mySnp.fields[5]+ '\t' + mySnp.fields[6]
       + '\t' + mySnp.fields[7] + '\t' + mySnp.fields[8] + '\t' + res + '\n';
 return res;
} //end of parseSnpGff()

void writeFile(const char * argv, ofstream & output)
{
 output.open(argv, ios::out); 
 if(!output)
 {
  cerr << argv << " can not be opened for writing.\n";
  exit(1);
 }
}

string int2str(int n)
{
 stringstream ss;
 ss << n;
 return ss.str();
}

//snp not in gene, res = "."
//snp in gene, in cds, res = "CDS"
//snp in gene, not in cds res = "gene; non-coding region; "
string snpType(const vector<gff> & gffVect, const string & prom)
{
 string res = ".", type = "";
 
 for(int i = 0; i < gffVect.size(); i++)
 {
  type = "";
  for(int j = 0; j < gffVect[i].type.length(); j++)
    type = type + (char) toupper(gffVect[i].type[j]);
  if(type.find("CD") !=string::npos)
  {
   res = string("CDS; ") +  gffVect[i].desc;
   return res;
  }
  if(type.find("GENE") !=string::npos)
  {
   res = string("gene; non-coding region; ") +  gffVect[i].desc;
  }
 }
 if(res == "." && prom.find("promoter") !=string::npos)
 {
  res = prom;
 }
 return res;
}
