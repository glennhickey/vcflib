/*
    vcflib C++ library for parsing and manipulating VCF files

    This software is published under the MIT License. See the LICENSE file.
*/

// vcfleftshift: trim and left-align variants against the reference, like bcftools norm -f, but
// never move a variant onto or across a base that another variant on one of its haplotypes
// already describes.
//
// bcftools norm -f places each record by the reference alone.  In a phased VCF, that can slide an
// insertion onto the base of a SNP the same haplotype carries, or a deletion through a tandem
// copy and across that haplotype's own SNPs: the haplotype then has two records describing one
// base, or a SNP inside its own deletion.  Here, each haplotype keeps the end of the last variant
// it carries, and an indel only moves left while it stays clear of that end for every haplotype
// that carries it.  The bases it moves across are reference on all of those haplotypes, so the
// move changes none of their sequences.

#include <htslib/faidx.h>

#include <getopt.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using namespace std;

static void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " -f ref.fa [options] < in.vcf > out.vcf" << endl
         << endl
         << "Trim and left-align variants against the reference, like bcftools norm -f, but never" << endl
         << "move a variant onto or across a base that another variant on one of its haplotypes" << endl
         << "describes.  Each haplotype (one allele of a sample's GT) remembers where the last" << endl
         << "variant it carries ends, and an indel only moves left while it stays to the right of" << endl
         << "that for every haplotype that carries it.  The bases it moves across are then reference" << endl
         << "on all of them, so no haplotype's sequence changes." << endl
         << endl
         << "Records at the same position that share a haplotype are trimmed but not moved." << endl
         << "Multiallelic and symbolic records are passed through (split with bcftools norm -m -any" << endl
         << "first).  The input must be sorted by position.  Moved records can come out of order:" << endl
         << "sort the output afterwards.  Without genotypes nothing limits a move, as in bcftools norm." << endl
         << endl
         << "options:" << endl
         << "    -f, --fasta-reference FILE  reference (plain or bgzipped FASTA, with .fai)" << endl
         << "    -s, --stats                 print counts of what was trimmed and moved to stderr" << endl
         << "    -h, --help                  this message" << endl
         << endl
         << "Type: transformation" << endl;
}

// reference access, cached in windows since moves read the bases just left of each record
class RefCache {
public:
    RefCache(const string& path) {
        fai = fai_load(path.c_str());
        if (!fai) {
            cerr << "error: could not load reference " << path << endl;
            exit(1);
        }
    }
    ~RefCache() { fai_destroy(fai); }
    // 1-based inclusive range, upper case; empty if out of range or unknown sequence
    string get(const string& chrom, int64_t start, int64_t end) {
        if (start < 1 || end < start) return "";
        if (chrom != cacheChrom && !faidx_has_seq(fai, chrom.c_str())) {
            cerr << "error: sequence " << chrom << " is not in the reference" << endl;
            exit(1);
        }
        if (chrom != cacheChrom || start < cacheStart || end >= cacheStart + (int64_t)cache.size()) {
            int64_t fetchStart = max<int64_t>(1, start - window);
            int64_t fetchEnd = max(end, start + window);
            hts_pos_t len = 0;
            char* seq = faidx_fetch_seq64(fai, chrom.c_str(), fetchStart - 1, fetchEnd - 1, &len);
            if (!seq || len <= 0) {
                if (seq) free(seq);
                cacheChrom.clear();
                cache.clear();
                return "";
            }
            cache.assign(seq, len);
            free(seq);
            for (auto& c : cache) c = toupper(c);
            cacheChrom = chrom;
            cacheStart = fetchStart;
            if (end >= cacheStart + (int64_t)cache.size()) return "";
        }
        return cache.substr(start - cacheStart, end - start + 1);
    }
private:
    faidx_t* fai = nullptr;
    string cacheChrom;
    int64_t cacheStart = 0;
    string cache;
    static const int64_t window = 1 << 20;
};

struct Record {
    string line;                  // kept as is apart from POS, REF and ALT
    vector<size_t> tabs;          // positions of the first five tabs in line
    string chrom;
    int64_t pos = 0;
    string ref;
    vector<string> alts;
    vector<pair<int, int>> calls; // (haplotype, allele > 0) for each called non-reference allele
    bool editable = false;        // biallelic, plain bases
    bool changed = false;
};

static vector<string> split(const string& s, char d) {
    vector<string> out;
    size_t b = 0;
    for (size_t e = s.find(d); e != string::npos; b = e + 1, e = s.find(d, b)) out.push_back(s.substr(b, e - b));
    out.push_back(s.substr(b));
    return out;
}

static bool plainBases(const string& s) {
    if (s.empty()) return false;
    for (char c : s) {
        char u = toupper(c);
        if (u != 'A' && u != 'C' && u != 'G' && u != 'T' && u != 'N') return false;
    }
    return true;
}

// parse what's needed out of a data line
static bool parseRecord(const string& line, Record& rec) {
    rec.line = line;
    rec.tabs.clear();
    size_t p = 0;
    for (int i = 0; i < 5; ++i) {
        p = line.find('\t', p);
        if (p == string::npos) return false;
        rec.tabs.push_back(p++);
    }
    rec.chrom = line.substr(0, rec.tabs[0]);
    rec.pos = atoll(line.c_str() + rec.tabs[0] + 1);
    rec.ref = line.substr(rec.tabs[2] + 1, rec.tabs[3] - rec.tabs[2] - 1);
    rec.alts = split(line.substr(rec.tabs[3] + 1, rec.tabs[4] - rec.tabs[3] - 1), ',');
    rec.editable = rec.alts.size() == 1 && plainBases(rec.ref) && plainBases(rec.alts[0]);
    for (auto& c : rec.ref) c = toupper(c);
    if (rec.editable) for (auto& c : rec.alts[0]) c = toupper(c);

    // genotypes: FORMAT is column 9, samples follow
    rec.calls.clear();
    vector<string> rest = split(line.substr(rec.tabs[4] + 1), '\t');
    // rest: QUAL FILTER INFO FORMAT samples...
    if (rest.size() < 5) return true;
    vector<string> format = split(rest[3], ':');
    auto gtIt = find(format.begin(), format.end(), "GT");
    if (gtIt == format.end()) return true;
    size_t gtIndex = gtIt - format.begin();
    for (size_t s = 4; s < rest.size(); ++s) {
        const string& sample = rest[s];
        // the GT subfield
        size_t b = 0;
        for (size_t k = 0; k < gtIndex && b != string::npos; ++k) {
            b = sample.find(':', b);
            if (b != string::npos) ++b;
        }
        if (b == string::npos) continue;
        size_t e = sample.find(':', b);
        string gt = sample.substr(b, e == string::npos ? string::npos : e - b);
        int h = 0;
        size_t i = 0;
        while (i <= gt.size()) {
            size_t j = gt.find_first_of("|/", i);
            string a = gt.substr(i, j == string::npos ? string::npos : j - i);
            if (!a.empty() && a != "." && a != "0") {
                int allele = atoi(a.c_str());
                if (allele > 0) rec.calls.push_back(make_pair((int)((s - 4) * 2 + h), allele));
            }
            ++h;
            if (j == string::npos || h > 1) break;
            i = j + 1;
        }
    }
    return true;
}

// last base an allele describes: its REF minus the suffix it shares with REF
static int64_t alleleEnd(const Record& rec, int allele) {
    const string& r = rec.ref;
    if (allele < 1 || allele > (int)rec.alts.size()) return rec.pos + (int64_t)r.size() - 1;
    const string& a = rec.alts[allele - 1];
    size_t i = 0;
    while (i + 1 < r.size() && i + 1 < a.size() && r[r.size() - 1 - i] == a[a.size() - 1 - i]) ++i;
    return rec.pos + (int64_t)(r.size() - i) - 1;
}

struct Stats {
    uint64_t records = 0, trimmed = 0, indels = 0, moved = 0, limited = 0, samePos = 0, refMismatch = 0, skipped = 0;
};

// trim and left-align one biallelic record, not moving it to a position <= limit
static void normalize(Record& rec, RefCache& refs, int64_t limit, bool mayMove, Stats& stats) {
    string r = rec.ref;
    string a = rec.alts[0];
    int64_t pos = rec.pos;
    if (r == a) return;
    if (refs.get(rec.chrom, pos, pos + (int64_t)r.size() - 1) != r) {
        ++stats.refMismatch;
        return;
    }
    if (r.size() != a.size()) ++stats.indels;
    bool wanted = false; // would have moved further without the limit
    while (!r.empty() && !a.empty() && r.back() == a.back()) {
        if (r.size() > 1 && a.size() > 1) {
            // drop a shared last base: no move
            r.pop_back();
            a.pop_back();
            continue;
        }
        // one allele would be empty: move left one base, taking that base as the new anchor
        if (pos - 1 < 1) break;
        if (!mayMove || pos - 1 <= limit) {
            wanted = true;
            break;
        }
        string b = refs.get(rec.chrom, pos - 1, pos - 1);
        if (b.empty()) break;
        r = b + r.substr(0, r.size() - 1);
        a = b + a.substr(0, a.size() - 1);
        --pos;
    }
    // drop shared first bases, keeping one
    while (r.size() > 1 && a.size() > 1 && r[0] == a[0]) {
        r.erase(0, 1);
        a.erase(0, 1);
        ++pos;
    }
    if (wanted) ++stats.limited;
    if (pos < rec.pos) ++stats.moved;
    if (r != rec.ref || a != rec.alts[0] || pos != rec.pos) {
        if (pos >= rec.pos) ++stats.trimmed;
        rec.ref = r;
        rec.alts[0] = a;
        rec.pos = pos;
        rec.changed = true;
    }
}

static void writeRecord(const Record& rec) {
    if (!rec.changed) {
        cout << rec.line << '\n';
        return;
    }
    const string& l = rec.line;
    // CHROM, POS, the ID with the tabs around it, REF, ALT, then the rest from the fifth tab on
    cout << rec.chrom << '\t' << rec.pos << l.substr(rec.tabs[1], rec.tabs[2] - rec.tabs[1] + 1)
         << rec.ref << '\t' << rec.alts[0] << l.substr(rec.tabs[4]) << '\n';
}

// all records at one position: limits come from earlier positions only, so their order
// within the group does not matter
static void processGroup(vector<Record>& group, vector<int64_t>& lastEnd, RefCache& refs, Stats& stats) {
    vector<int> count;
    for (const auto& rec : group) {
        for (const auto& call : rec.calls) {
            if (call.first >= (int)lastEnd.size()) lastEnd.resize(call.first + 1, 0);
            if (call.first >= (int)count.size()) count.resize(call.first + 1, 0);
            ++count[call.first];
        }
    }
    for (auto& rec : group) {
        ++stats.records;
        if (!rec.editable) {
            ++stats.skipped;
            continue;
        }
        int64_t limit = 0;
        bool mayMove = true;
        for (const auto& call : rec.calls) {
            limit = max(limit, lastEnd[call.first]);
            if (count[call.first] > 1) mayMove = false;
        }
        if (!mayMove) ++stats.samePos;
        normalize(rec, refs, limit, mayMove, stats);
    }
    for (const auto& rec : group) {
        for (const auto& call : rec.calls) {
            lastEnd[call.first] = max(lastEnd[call.first], alleleEnd(rec, call.second));
        }
        writeRecord(rec);
    }
    group.clear();
}

int main(int argc, char** argv) {
    string fasta;
    bool printStats = false;
    static struct option longOptions[] = {
        {"fasta-reference", required_argument, 0, 'f'},
        {"stats", no_argument, 0, 's'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "f:sh", longOptions, nullptr)) != -1) {
        switch (c) {
        case 'f': fasta = optarg; break;
        case 's': printStats = true; break;
        case 'h': printSummary(argv); return 0;
        default: printSummary(argv); return 1;
        }
    }
    if (fasta.empty()) {
        printSummary(argv);
        return 1;
    }
    ios_base::sync_with_stdio(false);
    RefCache refs(fasta);
    Stats stats;
    vector<int64_t> lastEnd;
    vector<Record> group;
    set<string> doneChroms;
    string line;
    string chrom;
    int64_t lastPos = 0;
    while (getline(cin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            cout << line << '\n';
            continue;
        }
        Record rec;
        if (!parseRecord(line, rec)) {
            cerr << "error: could not parse VCF line: " << line.substr(0, 200) << endl;
            return 1;
        }
        if (rec.chrom != chrom) {
            if (!group.empty()) processGroup(group, lastEnd, refs, stats);
            if (!chrom.empty()) doneChroms.insert(chrom);
            if (doneChroms.count(rec.chrom)) {
                cerr << "error: input is not sorted: " << rec.chrom << " appears again after " << chrom << endl;
                return 1;
            }
            chrom = rec.chrom;
            lastPos = 0;
            fill(lastEnd.begin(), lastEnd.end(), 0);
        } else if (rec.pos < lastPos) {
            cerr << "error: input is not sorted: " << chrom << ":" << rec.pos << " after " << chrom << ":" << lastPos << endl;
            return 1;
        } else if (rec.pos != lastPos && !group.empty()) {
            processGroup(group, lastEnd, refs, stats);
        }
        lastPos = rec.pos;
        group.push_back(move(rec));
    }
    if (!group.empty()) processGroup(group, lastEnd, refs, stats);
    cout.flush();
    if (printStats) {
        cerr << "records\t" << stats.records << endl
             << "not biallelic or not plain bases (passed through)\t" << stats.skipped << endl
             << "REF does not match the reference (passed through)\t" << stats.refMismatch << endl
             << "indels\t" << stats.indels << endl
             << "moved left\t" << stats.moved << endl
             << "trimmed in place\t" << stats.trimmed << endl
             << "stopped short of a variant on one of their haplotypes\t" << stats.limited << endl
             << "shares its position with another record on one of its haplotypes (not moved)\t" << stats.samePos << endl;
    }
    return 0;
}
