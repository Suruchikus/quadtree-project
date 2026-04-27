#pragma once

#include <cstdint>
#include <vector>

class RankSupport64 {
public:
    RankSupport64() = default;

    // factor_words = how many 64-bit words per checkpoint block
    // Example: factor_words = 8 means one checkpoint every 512 bits.
    void build(const std::vector<uint64_t>& bits, uint64_t bit_len, uint64_t factor_words = 8);

    // Number of 1s in [0, pos)
    uint64_t rank1_before(uint64_t pos) const;

    // Number of 0s in [0, pos)
    uint64_t rank0_before(uint64_t pos) const;

    uint64_t bit_length() const { return bit_len_; }

    uint64_t space_in_bits() const;

private:
    const std::vector<uint64_t>* bits_ = nullptr;
    uint64_t bit_len_ = 0;

    static constexpr uint64_t WORD_BITS = 64;

    // Number of 64-bit words per sampled block
    uint64_t factor_words_ = 8;

    // Bits per sampled block = 64 * factor_words_
    uint64_t sample_bits_ = 512;

    // Rs_[j] = total number of 1s before sampled block j
    std::vector<uint64_t> Rs_;

};