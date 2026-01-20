#include "deflate.hpp"

#include <algorithm>

std::vector<std::uint8_t> rez::decompress_deflate(std::span<const std::uint8_t> deflate_data)
{
    impl::deflate::Deflate_bitstream bitstream {deflate_data};
    std::vector<std::uint8_t> inflated_data;
    inflated_data.reserve(5000); // 5KB
    int bfinal {0};
    int btype {0};
    do {
        bfinal = bitstream.read_bits(1);
        btype = bitstream.read_bits(2);
        switch(btype) {
            case 0: // no compression
                bitstream.skip_until_next_byte_boundary();
                impl::deflate::decompress_uncompressed(inflated_data, bitstream);
                break;
            case 1: // fixed Huffman codes
                impl::deflate::decompress_fixed(inflated_data, bitstream);
                break;
            case 2: // dynamic Huffman codes
                impl::deflate::decompress_dynamic(inflated_data, bitstream);
                break;
            default:
                throw Exception {Error::bad_formed_data};
        }
    } while(not bfinal);

    return inflated_data;
}

void rez::impl::deflate::decompress_uncompressed(std::vector<std::uint8_t>& inflated_data, Deflate_bitstream& bitstream)
{
    int len {static_cast<int>(bitstream.read_bits(16))};
    int nlen {static_cast<int>(bitstream.read_bits(16))};
    if(~len != nlen) throw Exception {Error::bad_formed_data};
    if(len == 0) return; // zero length is allowed

    std::span<const std::uint8_t> uncompressed_data {bitstream.read_bytes(len)};
    inflated_data.insert(inflated_data.end(), uncompressed_data.begin(), uncompressed_data.end());
}

void rez::impl::deflate::decompress_fixed(std::vector<std::uint8_t>& inflated_data, Deflate_bitstream& bitstream)
{
    static const std::vector<Huffman_code> huffman_codes(make_fixed_huffman_table());
    int symbol {fetch_symbol_in_fixed_block(huffman_codes, bitstream)};
    while(symbol != 256) {
        if(symbol < 256) { inflated_data.push_back(static_cast<std::uint8_t>(symbol)); }
        else {
            if(symbol > 285) throw Exception {Error::bad_formed_data};
            symbol -= 257;
            const int length {length_bases[symbol] + static_cast<int>(bitstream.read_bits(length_extra_bits[symbol]))};

            // for fixed huffman blocks, a distance huffman code and its symbol have the same values
            symbol = bitstream.read_bits(5);
            symbol = reverse_bits_from_lsbit(symbol, 5);
            if(symbol > 29) throw Exception {Error::bad_formed_data};
            const std::int32_t distance {distance_bases[symbol] + bitstream.read_bits(distance_extra_bits[symbol])};
            if(distance > inflated_data.size() or distance > 32768) throw Exception {Error::bad_formed_data};

            lz77_copy(inflated_data, length, distance);
        }

        symbol = fetch_symbol_in_fixed_block(huffman_codes, bitstream);
    }
}

void rez::impl::deflate::decompress_dynamic(std::vector<std::uint8_t>& inflated_data, Deflate_bitstream& bitstream)
{
    const int hlit {static_cast<int>(bitstream.read_bits(5) + 257)};
    const int hdist {static_cast<int>(bitstream.read_bits(5) + 1)};
    const int hclen {static_cast<int>(bitstream.read_bits(4) + 4)};
    if(hlit > 286 or hdist > 30) throw Exception {Error::bad_formed_data};

    // the order of slots with which to place the bit-lengths of the codes of the code bit-length alphabet
    static constexpr std::array<int, 19> ordered_indexes {
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    };

    std::array<int, 19> code_length_alphabet_bit_lengths;
    for(int i = 0; i < hclen; ++i) {
        code_length_alphabet_bit_lengths[ordered_indexes[i]] = bitstream.read_bits(3);
    }
    for(std::uint32_t i = hclen; i < 19; ++i) {
        code_length_alphabet_bit_lengths[ordered_indexes[i]] = 0;
    }
    std::vector<Huffman_code> code_length_alphabet(make_huffman_codes_from_bit_lengths(code_length_alphabet_bit_lengths));

    // bit-lengths of both the literal+length alphabet and the distance alphabet
    std::vector<int> alphabets_bit_lengths;
    const int hlit_hdist {hlit + hdist};
    alphabets_bit_lengths.reserve(hlit_hdist);
    // cannot be like this: for(int i = 0; i < hlit_hdist; ++i)
    while(alphabets_bit_lengths.size() < hlit_hdist) {
        const int symbol {fetch_symbol_in_dynamic_block(code_length_alphabet, bitstream)};
        if(symbol < 16) { alphabets_bit_lengths.push_back(symbol); }
        else if(symbol == 16) {
            if(alphabets_bit_lengths.empty()) throw Exception {Error::bad_formed_data};
            const int value_to_copy {alphabets_bit_lengths.back()};
            const int times_to_copy {static_cast<int>(bitstream.read_bits(2)) + 3};
            alphabets_bit_lengths.insert(alphabets_bit_lengths.end(), times_to_copy, value_to_copy);
        }
        else if(symbol == 17) {
            const int times_to_copy {static_cast<int>(bitstream.read_bits(3)) + 3};
            alphabets_bit_lengths.insert(alphabets_bit_lengths.end(), times_to_copy, 0);
        }
        else if(symbol == 18) {
            const int times_to_copy {static_cast<int>(bitstream.read_bits(7)) + 11};
            alphabets_bit_lengths.insert(alphabets_bit_lengths.end(), times_to_copy, 0);
        }
        else { throw Exception {Error::bad_formed_data}; }
    }

    std::span<const int> literal_length_alphabet_bit_lengths(alphabets_bit_lengths.begin(), hlit);
    std::vector<Huffman_code> literal_length_alphabet(make_huffman_codes_from_bit_lengths(literal_length_alphabet_bit_lengths));
    /* this is so silly: the case in where the amount of bit-lengths for the distance alphabet is 1
    * and that lonely bit-length happens to be zero is valid, it means that the data to decompress
    * is all literals and there aren't length or distance codes. It's silly because the "no compression"
    * block (BTYPE == 0) already exists, this was completely unnecessary.
    */
    std::span<const int> distance_alphabet_bit_lengths(alphabets_bit_lengths.begin() + hlit, hdist);
    std::vector<Huffman_code> distance_alphabet;
    if(distance_alphabet_bit_lengths.size() == 1 and distance_alphabet_bit_lengths[0] == 0) { /* do nothing */ }
    else { distance_alphabet = make_huffman_codes_from_bit_lengths(distance_alphabet_bit_lengths); }

    /* loop copy-pasted from decompress_fixed, the differences are small and I could put this loop
    * in a single function to avoid code duplication, but in this case, I am going to allow the
    * duplication */
    int symbol {fetch_symbol_in_dynamic_block(literal_length_alphabet, bitstream)};
    while(symbol != 256) {
        if(symbol < 256) { inflated_data.push_back(static_cast<std::uint8_t>(symbol)); }
        else {
            if(distance_alphabet.empty()) throw Exception {Error::bad_formed_data};
            if(symbol > 285) throw Exception {Error::bad_formed_data};
            symbol -= 257;
            const int length {length_bases[symbol] + static_cast<int>(bitstream.read_bits(length_extra_bits[symbol]))};

            symbol = fetch_symbol_in_dynamic_block(distance_alphabet, bitstream);
            if(symbol > 29) throw Exception {Error::bad_formed_data};
            const std::int32_t distance {distance_bases[symbol] + bitstream.read_bits(distance_extra_bits[symbol])};
            if(distance > inflated_data.size() or distance > 32768) throw Exception {Error::bad_formed_data};

            lz77_copy(inflated_data, length, distance);
        }

        symbol = fetch_symbol_in_dynamic_block(literal_length_alphabet, bitstream);
    }
}

std::vector<rez::impl::deflate::Huffman_code> rez::impl::deflate::make_fixed_huffman_table()
{
    // bl_count[7 (for example)] == number of codes that have 7 bits
    //constexpr std::array<int, 10> bl_count { 0, 0, 0, 0, 0, 0, 0, 24, 152, 112 };

    std::vector<int> bit_lengths;
    bit_lengths.reserve(288);
    for(int i = 0; i < 144; ++i) { bit_lengths.push_back(8); }
    for(int i = 144; i < 256; ++i) { bit_lengths.push_back(9); }
    for(int i = 256; i < 280; ++i) { bit_lengths.push_back(7); }
    for(int i = 280; i < 288; ++i) { bit_lengths.push_back(8); }

    // generate the codes for each bit-length
    std::vector<Huffman_code> huffman_codes;
    huffman_codes.reserve(288);
    int code {0}; // the smallest valid code
    for(int i = 7; i < 10; ++i) {
        for(int j = 0; j < bit_lengths.size(); ++j) {
            if(bit_lengths[j] != i) continue;

            huffman_codes.emplace_back(code, i, j);
            // within a bit-length, the codes are assigned consecutive values
            ++code;
        }

        // an extra bit must be added just before going to the next bit-length
        code <<= 1;
    }

    return huffman_codes;
}

std::vector<rez::impl::deflate::Huffman_code> rez::impl::deflate::make_huffman_codes_from_bit_lengths(std::span<const int> bit_lengths)
{
    auto iter {std::max_element(bit_lengths.begin(), bit_lengths.end())};
    // bl_count[7 (for example)] == number of codes that have 7 bits
    std::vector<int> bl_count;
    bl_count.resize(*iter + 1); // + 1 because you must account for bit-length 0

    // count the number of codes for each bit-length
    for(int i = 0; i < bit_lengths.size(); ++i) {
        bl_count[bit_lengths[i]] += 1;
    }

    // generate the codes for each bit-length
    std::vector<Huffman_code> huffman_codes;
    int code {0}; // the smallest valid code
    for(int i = 1; i <= *iter; ++i) {
        // an extra bit must be added just before going to the next bit-length
        code <<= 1;
        /* the above "code <<= 1;" must be before the below
        * if(bl_count[i] == 0) continue;" because there can be bl_counts that
        * are equal to zero between bl_counts that aren't, and the left-shift
        * must still be done
        */
        if(bl_count[i] == 0) continue;

        for(int j = 0; j < bit_lengths.size(); ++j) {
            if(bit_lengths[j] != i) continue;

            huffman_codes.emplace_back(code, i, j);
            // within a bit-length, the codes are assigned consecutive values
            ++code;
        }
    }

    return huffman_codes;
}

int rez::impl::deflate::fetch_symbol_in_fixed_block(const std::vector<Huffman_code>& huffman_codes, Deflate_bitstream& bitstream)
{
    for(int i = 7; i < 10; ++i) {
        int code {static_cast<int>(bitstream.peek_bits(i))};
        code = reverse_bits_from_lsbit(code, i);

        for(const Huffman_code hc : huffman_codes) {
            if(hc.bit_length == i and hc.code == code) {
                bitstream.skip_bits(i);
                return hc.symbol;
            }
        }
    }

    throw Exception {Error::bad_formed_data};
}

int rez::impl::deflate::fetch_symbol_in_dynamic_block(const std::vector<Huffman_code>& huffman_codes, Deflate_bitstream& bitstream)
{
    for(int i = 1; i < 16; ++i) {
        int code {static_cast<int>(bitstream.peek_bits(i))};
        code = reverse_bits_from_lsbit(code, i);

        for(const Huffman_code hc : huffman_codes) {
            if(hc.bit_length == i and hc.code == code) {
                bitstream.skip_bits(i);
                return hc.symbol;
            }
        }
    }

    throw Exception {Error::bad_formed_data};
}

void rez::impl::deflate::lz77_copy(std::vector<std::uint8_t>& inflated_data, int length, const std::int32_t distance)
{
    const std::size_t beginning_of_copy {inflated_data.size() - distance};
    std::size_t copy_from {beginning_of_copy};

    inflated_data.reserve(inflated_data.size() + length);
    for(int i = 0; i < length; ++i) {
        inflated_data.push_back(inflated_data[copy_from]);
        ++copy_from;
        if(copy_from == inflated_data.size()) {
            copy_from = beginning_of_copy;
        }
    }
}
