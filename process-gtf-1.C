/*
linyong@fry /lab/solexa_page/linyong$ cat  gencode/gencode.v42.primary_assembly.noPar.annotation.gtf | grep -iw ddx3x | awk '$3 == "transcript"' | grep -o "transcript_type.\{30\}" | sort| uniq -c
     14 transcript_type "nonsense_mediated_decay"; tr
      3 transcript_type "protein_coding_CDS_not_defin
     23 transcript_type "protein_coding"; transcript_
     29 transcript_type "retained_intron"; transcript

cat  gencode/gencode.v42.primary_assembly.noPar.annotation.gtf | grep -iw ddx3x | awk '$3 == "transcript"' | grep "transcript_type \"protein_coding\";"
chrX    HAVANA  gene    41333348        41364472        .       +       .       gene_id "ENSG00000215301.11"; gene_type "protein_coding"; gene_name "DDX3X"; level 2; hgnc_id "HGNC:2745"; havana_gene "OTTHUMG00000021369.22";

chrX    HAVANA  transcript      41333348        41350287        .       +       .       gene_id "ENSG00000215301.11"; transcript_id "ENST00000629496.3"; gene_type "protein_coding"; gene_name "DDX3X"; transcript_type "protein_coding"; transcript_name "DDX3X-217"; level 2; protein_id "ENSP00000487224.1"; transcript_support_level "5"; hgnc_id "HGNC:2745"; tag "alternative_5_UTR"; tag "non_submitted_evidence"; tag "basic"; tag "appris_principal_4"; tag "CCDS"; ccdsid "CCDS43931.1"; havana_gene "OTTHUMG00000021369.22"; havana_transcript "OTTHUMT00000481647.3";

chrX    HAVANA  exon    41333348        41333699        .       +       .       gene_id "ENSG00000215301.11"; transcript_id "ENST00000629496.3"; gene_type "protein_coding"; gene_name "DDX3X"; transcript_type "protein_coding"; transcript_name "DDX3X-217"; exon_number 1; exon_id "ENSE00003771363.1"; level 2; protein_id "ENSP00000487224.1"; transcript_support_level "5"; hgnc_id "HGNC:2745"; tag "alternative_5_UTR"; tag "non_submitted_evidence"; tag "basic"; tag "appris_principal_4"; tag "CCDS"; ccdsid "CCDS43931.1"; havana_gene "OTTHUMG00000021369.22"; havana_transcript "OTTHUMT00000481647.3"; 

chrX    HAVANA  UTR     41347717        41350287        .       +       .       gene_id "ENSG00000215301.11"; transcript_id "ENST00000629496.3"; gene_type "protein_coding"; gene_name "DDX3X"; transcript_type "protein_coding"; transcript_name "DDX3X-217"; exon_number 18; exon_id "ENSE00003830364.1"; level 2; protein_id "ENSP00000487224.1"; transcript_support_level "5"; hgnc_id "HGNC:2745"; tag "alternative_5_UTR"; tag "non_submitted_evidence"; tag "basic"; tag "appris_principal_4"; tag "CCDS"; ccdsid "CCDS43931.1"; havana_gene "OTTHUMG00000021369.22"; havana_transcript "OTTHUMT00000481647.3";

chrX    HAVANA  CDS     41341484        41341616        .       +       2       gene_id "ENSG00000215301.11"; transcript_id "ENST00000629496.3"; gene_type "protein_coding"; gene_name "DDX3X"; transcript_type "protein_coding"; transcript_name "DDX3X-217"; exon_number 5; exon_id "ENSE00003572608.1"; level 2; protein_id "ENSP00000487224.1"; transcript_support_level "5"; hgnc_id "HGNC:2745"; tag "alternative_5_UTR"; tag "non_submitted_evidence"; tag "basic"; tag "appris_principal_4"; tag "CCDS"; ccdsid "CCDS43931.1"; havana_gene "OTTHUMG00000021369.22"; havana_transcript "OTTHUMT00000481647.3";

linyong@fry /lab/solexa_page/linyong$ cat  gencode/gencode.v42.primary_assembly.noPar.annotation.gtf | grep  "gene_id.*gene_id" | wc
      0       0       0 
cat  gencode/gencode.v42.primary_assembly.noPar.annotation.gtf | grep  "transcript_id.*transcript_id" | wc
      0       0       0
linyong@fry /lab/solexa_page/linyong$ cat gencode.v42.primary_assembly.noPar.protein.xist.gtf | cut -f 9 | sed 's/peak_id=/     gene_id /' > temp-gene-ids-131pm

linyong@fry /lab/solexa_page/linyong$ cat  gencode/gencode.v42.primary_assembly.noPar.annotation.gtf | awk '$3 == "transcript"' | grep "transcript_type \"protein_coding\";"  | cut -f 9 | awk 'BEGIN {FS = ";"} {print $2}'  > temp-transcriptid-140pm

geneID-input-file 'gene_id "ENSG00000215301.11"'      protein coding genes, xist
transcript-ip-file 'transcript_id "ENST00000629496.3"'       transcript_type "protein_coding"; 

if field #3 is a gene, then gene_id xxxxxx contained in the geneID-input-file
else
  gene_id in the file && transcript_id in second file


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
void getFieldContent2(vector<string> & fields, string del, string line);

int main(int argc, char *argv[])
{
 if(argc != 5)
 {
  cerr << argv[0] << "	gtf-file	gene-ID-file	transcript-ID-file	output-gtf-file-name" << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 map<string, int> reads, transc;

 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  if(line.length() > 0)
    reads[line] = 1234; 
  getline(input, line);
 }
 input.close();

 readFile(argv[3], input);
 getline(input, line);
 while(!input.eof())
 {
  if(line.length() > 0)
    transc[line] = 4321;
  getline(input, line);
 }
 input.close();

 
 writeFile(argv[4], output);

 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent(lineFields, '\t',  line);
  if(lineFields.size() >= 9) //at least 9 fields in gtf file
  {
   string s, g, t;
   vector<string> lf1;
   s=lineFields[8];
   if(lineFields[2] == "gene")
   {
    getFieldContent(lf1, ';', s);
    g=string("	") + lf1[0];
    if(reads.count(g) > 0)
    {
     output << line << endl;
    }
   }
   else
   {
    getFieldContent(lf1, ';', s);
    if(lf1.size() > 2)
    {
     g=string("	") + lf1[0];
     t=lf1[1];
     if(reads.count(g) > 0 && transc.count(t) > 0)
     {
      output << line << endl;
     }
    }    
   }
  }
  getline(input, line);
 }
 output.close();
 input.close();
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


