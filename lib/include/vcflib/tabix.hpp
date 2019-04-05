#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include "htslib/bgzf.h"
#include "htslib/tbx.h"
#include "htslib/kseq.h"
#include <iostream>
#include <cstring>
#include <vector>


using namespace std;

class Tabix {

    htsFile* fn;
    tbx_t* tbx;
    kstring_t str;
    hts_itr_t* iter;
    const tbx_conf_t *idxconf;
    int tid, beg, end;
    string firstline;
    bool has_jumped;
    vector<string>::iterator current_chrom;

public:
    string filename;
    vector<string> chroms;

    Tabix(void);
    Tabix(string& file);
    ~Tabix(void);

    const kstring_t * getKstringPtr();
    void getHeader(string& header);
    bool setRegion(string& region);
    bool getNextLine(string& line);
    bool getNextLineKS();

};
