#pragma once

#include "shared.hpp"

#include <vector>
#include <array>
#include <compare>

namespace rez {
    std::vector<std::uint8_t> decompress_deflate(std::span<const std::uint8_t> deflate_data);
}

namespace rez::impl::deflate {
    // base lengths for length symbols (257~285)
    constexpr std::array<int, 29> length_bases {
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
        35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258
    };

    // extra bits for length symbols
    constexpr std::array<int, 29> length_extra_bits {
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0
    };

    // base distances for distance symbols (0~29)
    constexpr std::array<int, 30> distance_bases {
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
        257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145,
        8193, 12289, 16385, 24577
    };

    // extra bits for distance symbols
    constexpr std::array<int, 30> distance_extra_bits {
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13
    };

    struct Huffman_code {
        Huffman_code() noexcept {}
        // just for emplace_back
        Huffman_code(int a_code, int a_bit_length, int a_symbol) noexcept : code {a_code}, bit_length {a_bit_length}, symbol {a_symbol} {}

        int code {0};
        int bit_length {0};
        int symbol {0};
    };

    using Deflate_bitstream = Bitstream<Bitstream_format::gif>;

    void decompress_uncompressed(std::vector<std::uint8_t>& inflated_data, Deflate_bitstream& bitstream);
    void decompress_fixed(std::vector<std::uint8_t>& inflated_data, Deflate_bitstream& bitstream);
    void decompress_dynamic(std::vector<std::uint8_t>& inflated_data, Deflate_bitstream& bitstream);

    std::vector<Huffman_code> make_fixed_huffman_table(); // for literals and lengths
    // used in decompress_dynamic
    std::vector<Huffman_code> make_huffman_codes_from_bit_lengths(std::span<const int> bit_lengths);

    // about fetch_symbol_in_fixed_block: only for literals and lengths
    int fetch_symbol_in_fixed_block(const std::vector<Huffman_code>& huffman_codes, Deflate_bitstream& bitstream);
    int fetch_symbol_in_dynamic_block(const std::vector<Huffman_code>& huffman_codes, Deflate_bitstream& bitstream);

    void lz77_copy(std::vector<std::uint8_t>& inflated_data, const int length, const std::int32_t distance);
}