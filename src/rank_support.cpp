#include "rank_support.h"

static inline uint64_t popcount64(uint64_t x) {
    return (uint64_t)__builtin_popcountll(x);
}

void RankSupport64::build(const std::vector<uint64_t>& bits, uint64_t bit_len, uint64_t factor_words) {
    bits_ = &bits;
    bit_len_ = bit_len;

    factor_words_ = (factor_words == 0 ? 8 : factor_words);
    sample_bits_ = WORD_BITS * factor_words_;

    Rs_.clear();

    if (bit_len_ == 0) {
        return;
    }

    const uint64_t num_words = (bit_len_ + 63ULL) >> 6;
    const uint64_t num_sample_blocks = (num_words + factor_words_ - 1) / factor_words_;

    Rs_.resize(num_sample_blocks + 1, 0ULL);

    uint64_t total_ones = 0;

    for (uint64_t block = 0; block < num_sample_blocks; ++block) {
        Rs_[block] = total_ones;

        const uint64_t word_start = block * factor_words_;
        const uint64_t word_end = std::min(word_start + factor_words_, num_words);

        for (uint64_t w = word_start; w < word_end; ++w) {
            uint64_t word = (*bits_)[w];

            // mask garbage bits in the final word
            if (w == num_words - 1) {
                const uint64_t valid_bits = bit_len_ - (w << 6);
                if (valid_bits < 64) {
                    const uint64_t mask = (valid_bits == 0)
                        ? 0ULL
                        : ((1ULL << valid_bits) - 1ULL);
                    word &= mask;
                }
            }

            total_ones += popcount64(word);
        }
    }

    Rs_[num_sample_blocks] = total_ones;
}

uint64_t RankSupport64::rank1_before(uint64_t pos) const {
    if (!bits_ || bit_len_ == 0 || pos == 0) {
        return 0ULL;
    }

    if (pos > bit_len_) {
        pos = bit_len_;
    }

    const uint64_t block = pos / sample_bits_;
    uint64_t ans = Rs_[block];

    const uint64_t word_start = block * factor_words_;
    const uint64_t word_end = pos >> 6;   // full words before pos
    const uint64_t offset = pos & 63ULL;  // remaining bits in current word

    for (uint64_t w = word_start; w < word_end; ++w) {
        ans += popcount64((*bits_)[w]);
    }

    if (offset > 0 && word_end < ((bit_len_ + 63ULL) >> 6)) {
        const uint64_t mask = (1ULL << offset) - 1ULL;
        ans += popcount64((*bits_)[word_end] & mask);
    }

    return ans;
}

uint64_t RankSupport64::rank0_before(uint64_t pos) const {
    if (pos > bit_len_) {
        pos = bit_len_;
    }
    return pos - rank1_before(pos);
}

uint64_t RankSupport64::space_in_bits() const {
    return (uint64_t)Rs_.size() * 64ULL;
}