#include <iostream>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <string>
#include <map>
#include <utility>
extern "C" {
    #include <htslib/sam.h>
}

struct Variant {
    std::string chrom;
    int pos;
};

bool operator<(const Variant &a, const Variant &b) {
    if (a.chrom == b.chrom) {
        return a.pos < b.pos;
    }
    return a.chrom < b.chrom;
}

using VariantPair = std::pair<Variant, Variant>;

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: program <variant_file> <cram_file>" << std::endl;
        return 1;
    }

    // Read variants
    std::ifstream file(argv[1]);
    std::vector<Variant> variants;
    std::string line;
    while (std::getline(file, line)) {
        // Assuming each line is: chrom pos
        Variant v;
        sscanf(line.c_str(), "%s %d", &v.chrom, &v.pos);
        variants.push_back(v);
    }

    // Generate pairs
    std::vector<VariantPair> variant_pairs;
    for (size_t i = 0; i < variants.size(); ++i) {
        for (size_t j = i + 1; j < variants.size(); ++j) {
            if (variants[i].chrom == variants[j].chrom && std::abs(variants[i].pos - variants[j].pos) < 500) {
                variant_pairs.push_back({variants[i], variants[j]});
            }
        }
    }

    // Open CRAM
    samFile *in = sam_open(argv[2], "r");
    bam_hdr_t *header = sam_hdr_read(in);
    hts_idx_t *idx = sam_index_load(in, argv[2]);
    if (!idx) {
        std::cerr << "Failed to open the index for the CRAM file." << std::endl;
        return 1;
    }

    // Check pairs in CRAM
    std::map<VariantPair, std::tuple<int, int, int>> results; // same-read, diff-read, total-read counts
    for (const auto& pair : variant_pairs) {
        int tid = bam_name2id(header, pair.first.chrom.c_str());
        hts_itr_t *iter = sam_itr_queryi(idx, tid, pair.first.pos, pair.second.pos + 1);
        bam1_t *read = bam_init1();

        std::unordered_set<std::string> read_names_first, read_names_second;
        while (sam_itr_next(in, iter, read) >= 0) {
            if (read->core.pos <= pair.first.pos && bam_endpos(read) > pair.first.pos) {
                read_names_first.insert(bam_get_qname(read));
            }
            if (read->core.pos <= pair.second.pos && bam_endpos(read) > pair.second.pos) {
                read_names_second.insert(bam_get_qname(read));
            }
        }
        bam_destroy1(read);
        hts_itr_destroy(iter);

        int same_read = 0, diff_read = 0;
        for (const auto& name : read_names_first) {
            if (read_names_second.count(name)) {
                same_read++;
            } else {
                diff_read++;
            }
        }
        results[pair] = std::make_tuple(same_read, diff_read, same_read + diff_read);
    }

    // Output results
    for (const auto& pair : results) {
        std::cout << pair.first.first.chrom << "\t" << pair.first.first.pos << "\t" 
                  << pair.first.second.chrom << "\t" << pair.first.second.pos << "\t"
                  << std::get<0>(pair.second) << "\t" << std::get<1>(pair.second) << "\t" << std::get<2>(pair.second) << std::endl;
    }

    // Cleanup
    bam_hdr_destroy(header);
    sam_close(in);
    hts_idx_destroy(idx);

    return 0;
}


