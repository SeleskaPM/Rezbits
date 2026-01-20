#include "deflate.hpp"

#include <algorithm>

std::vector<std::uint8_t> rez::decompress_deflate(std::span<const std::uint8_t> deflate_data)
{
    impl::deflate::Deflate_bitstream bitstream {deflate_data};
    std::vector<std::uint8_t> inflated_data;
    inflated_data.reserve(10000); // 10KB
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
    const std::int32_t len {bitstream.read_bits(16)};
    const std::int32_t nlen {bitstream.read_bits(16)};

    if((~len & 0xFFFF) != nlen) throw Exception {Error::bad_formed_data};
    if(len == 0) return; // zero length is allowed

    std::span<const std::uint8_t> uncompressed_data {bitstream.read_bytes(len)};
    inflated_data.insert(inflated_data.end(), uncompressed_data.begin(), uncompressed_data.end());
}

void rez::impl::deflate::decompress_fixed(std::vector<std::uint8_t>& inflated_data, Deflate_bitstream& bitstream)
{
    // bl_count[7 (for example)] == number of codes that have 7 bits
    //constexpr std::array<int, 10> bl_count { 0, 0, 0, 0, 0, 0, 0, 24, 152, 112 };

    std::vector<int> code_lengths;
    code_lengths.reserve(288);
    for(int i = 0; i < 144; ++i) { code_lengths.push_back(8); }
    for(int i = 144; i < 256; ++i) { code_lengths.push_back(9); }
    for(int i = 256; i < 280; ++i) { code_lengths.push_back(7); }
    for(int i = 280; i < 288; ++i) { code_lengths.push_back(8); }
    Decoding_table literal_length_alphabet {make_decoding_table_from_code_lengths(code_lengths, 9, 512)};

    std::vector<int> code_lengths_2;
    code_lengths_2.resize(30, 5);
    Decoding_table distance_alphabet {make_decoding_table_from_code_lengths(code_lengths_2, 5, 32)};

    process_symbols(inflated_data, literal_length_alphabet, distance_alphabet, bitstream);
}

void rez::impl::deflate::decompress_dynamic(std::vector<std::uint8_t>& inflated_data, Deflate_bitstream& bitstream)
{
    const int hlit {static_cast<int>(bitstream.read_bits(5) + 257)};
    const int hdist {static_cast<int>(bitstream.read_bits(5) + 1)};
    const int hclen {static_cast<int>(bitstream.read_bits(4) + 4)};
    if(hlit > 286 or hdist > 30) throw Exception {Error::bad_formed_data};

    /* the order of slots with which to place the bit-lengths of the
    * codes of the code bit-length alphabet */
    static constexpr std::array<int, 19> ordered_indexes {
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    };

    std::array<int, 19> code_length_alphabet_code_lengths {};
    for(int i = 0; i < hclen; ++i) {
        code_length_alphabet_code_lengths[ordered_indexes[i]] = bitstream.read_bits(3);
    }
    Decoding_table code_length_alphabet {make_decoding_table_from_code_lengths(code_length_alphabet_code_lengths, 7, 128)};

    // code-lengths of both the literal+length alphabet and the distance alphabet
    std::vector<int> alphabets_code_lengths;
    const int hlit_hdist {hlit + hdist};
    alphabets_code_lengths.reserve(hlit_hdist);
    // cannot be like this: for(int i = 0; i < hlit_hdist; ++i)
    while(alphabets_code_lengths.size() < hlit_hdist) {
        const int symbol {fetch_symbol(code_length_alphabet, bitstream)};
        if(symbol < 16) { alphabets_code_lengths.push_back(symbol); }
        else if(symbol == 16) {
            if(alphabets_code_lengths.empty()) throw Exception {Error::bad_formed_data};
            const int value_to_copy {alphabets_code_lengths.back()};
            const int times_to_copy {static_cast<int>(bitstream.read_bits(2) + 3)};
            alphabets_code_lengths.insert(alphabets_code_lengths.end(), times_to_copy, value_to_copy);
        }
        else if(symbol == 17) {
            const int times_to_copy {static_cast<int>(bitstream.read_bits(3) + 3)};
            alphabets_code_lengths.insert(alphabets_code_lengths.end(), times_to_copy, 0);
        }
        else if(symbol == 18) {
            const int times_to_copy {static_cast<int>(bitstream.read_bits(7) + 11)};
            alphabets_code_lengths.insert(alphabets_code_lengths.end(), times_to_copy, 0);
        }
        else { throw Exception {Error::bad_formed_data}; }
    }
    if(alphabets_code_lengths.size() > hlit_hdist) throw Exception {Error::bad_formed_data};

    std::span<const int> literal_length_alphabet_code_lengths(alphabets_code_lengths.begin(), hlit);
    /* the bit-width of the first area is 9 because that is the value
    * used by zlib for the literal+length alphabet. I could use the
    * same value that zlib-ng uses, but I decided to be more
    * conservative with memory. */
    Decoding_table literal_length_alphabet(make_decoding_table_from_code_lengths(literal_length_alphabet_code_lengths, 9, 852));
    /* this is so silly: the case in where the amount of code-lengths
    * for the distance alphabet is 1 and that lonely bit-length happens
    * to be zero is valid, it means that the data to decompress is all
    * literals and there aren't length or distance codes. It's silly
    * because the "no compression" block (BTYPE == 0) already exists,
    * this was completely unnecessary */
    std::span<const int> distance_alphabet_code_lengths(alphabets_code_lengths.begin() + hlit, hdist);
    Decoding_table distance_alphabet;
    if(distance_alphabet_code_lengths.size() == 1 and distance_alphabet_code_lengths[0] == 0) { /* do nothing */ }
    /* the bit-width of the first area is 6 because that is the value
    * used by zlib for the distance alphabet. I could use the same
    * value that zlib-ng uses, but I decided to be more conservative
    * with memory */
    else { distance_alphabet = make_decoding_table_from_code_lengths(distance_alphabet_code_lengths, 6, 592); }

    if(distance_alphabet.entries.empty()) {
        int symbol {fetch_symbol(literal_length_alphabet, bitstream)};
        while(symbol != 256) {
            if(symbol > 256) throw Exception {Error::bad_formed_data};
            inflated_data.push_back(static_cast<std::uint8_t>(symbol));
            symbol = fetch_symbol(literal_length_alphabet, bitstream);
        }    
    }
    else { process_symbols(inflated_data, literal_length_alphabet, distance_alphabet, bitstream); }
}

rez::impl::deflate::Decoding_table rez::impl::deflate::make_decoding_table_from_code_lengths(std::span<const int> code_lengths, int first_area_bitwidth, int allocation_size)
{
    /* The purpose of this function is to construct a decoding table
    * that can be indexed with raw bits from the DEFLATE stream and
    * return an entry that contains either a symbol or information
    * used to find another entry that will contain a symbol */

    // bl_count[7 (for example)] == number of codes that have 7 bits
    std::array<int, 16> bl_count {}; // [0,15] bits

    /* for this use, symbol_count can also be named "code_count" but
    * the name symbol_count is more general as the code uses the
    * variable multiple times */
    const int symbol_count {static_cast<int>(code_lengths.size())};
    for(int symbol = 0; symbol < symbol_count; ++symbol) {
        bl_count[code_lengths[symbol]] += 1;
    }

    // get the minimum and maximum bit-lengths for clamping
    int min {1};
    while(min < 16) {
        if(bl_count[min] != 0) break;
        ++min;
    }
    int max {15};
    while(max > 0) {
        if(bl_count[max] != 0) break;
        --max;
    }
    first_area_bitwidth = std::clamp(first_area_bitwidth, min, max);

    /* compute an offset for the first symbol of each code-length.
    * the offsets will be used to index into a symbol table for the
    * purpose of initialising it in a sorted manner. At the end,
    * a symbol that pertains to a code that has a bit-length of 1,
    * will have a smaller offset than a symbol that pertains to a
    * code that has a bit-length of 2 and so on */
    std::array<int, 16> offsets {};
    for(int code_length = min; code_length < max; ++code_length) {
        offsets[code_length + 1] = offsets[code_length] + bl_count[code_length];
    }

    /* initialise the symbols table. Thanks to the 'offsets' array, at
    * the end of the initialision, the 'symbols' table is ordered like
    * this:
    * [symbols that pertain to codes that have bit-length 1] and then
    * [symbols that pertain to codes that have bit-length 2] and then
    * [symbols that pertain to codes that have bit-length 3] and so on,
    * and within each set of symbols, the symbols are ordered in
    * sequential order, like this: 0, 1, 2, 3 and so on */
    std::vector<int> symbols;
    symbols.resize(symbol_count - bl_count[0]);
    for(int symbol = 0; symbol < symbol_count; ++symbol) {
        const int code_length {code_lengths[symbol]};
        if(code_length == 0) continue;
        symbols[offsets[code_length]] = symbol;
        offsets[code_length] += 1;
    }

    // now comes the loop that fills the decoding table with entries
    int current_area_bitwidth {first_area_bitwidth};
    int code {0}; // the smallest valid Huffman code (bit-reversed)
    unsigned int icode {0}; // will help in computing the next bit-reversed Huffman code. i stands for "intermediate"
    int current_code_bit_length = min;
    int area_offset {0};
    int bits_consumed_by_previous_area {0};
    /* impossible initial value on purpose to go to the next area of
    * the decoding table when the bit-length of the current code
    * becomes greater than the bit-width of the first area */
    int current_prefix {0x7FFF}; // impossible initial value on purpose
    const int prefix_mask {(1 << first_area_bitwidth) - 1};
    Decoding_table decoding_table;
    decoding_table.first_area_bitwidth = first_area_bitwidth;
    decoding_table.entries.resize(allocation_size);

    const int actual_symbol_count {static_cast<int>(symbols.size())};
    for(int i = 0; i < actual_symbol_count; ++i) {
        // create table entry
        Huffman_entry entry;
        entry.category = Huffman_entry::Category::symbol;
        entry.value = symbols[i];
        entry.bits_to_consume = current_code_bit_length - bits_consumed_by_previous_area;

        /* the following code computes all the decoding table indices
        * that can correspond to the current Huffman code (reversed
        * because that is how it comes in the DEFLATE stream). The
        * core of the idea is bit-patterns: imagine you need to read
        * a Huffman code that has a bit-length of 2, but that you
        * index the decoding table with 6 bits, so you grab 6 bits
        * from the DEFLATE stream; in that scenario, you have this
        * layout: 0000'00. The right part are the bits that pertain
        * to the Huffman code and the left part are bits that do not.
        * The bit-pattern of the left part could be whatever and due
        * to that, the same symbol must be assigned to all of those
        * bit-patterns.
        *
        * the code traverses the bit-patterns in reverse order, like
        * this:
        * decoding_table[1111'00]
        * decoding_table[1110'00]
        * decoding_table[1101'00]
        * until
        * decoding_table[0000'00]
        *
        * the code leverages the following dynamic: if you want to
        * have a certain number of high bits set to 1, and a certain
        * number of low bits set to 0, then you can do this:
        *
        * int a = 1 << 8; // 0b1'0000'0000
        * int b = 1 << 4; // 0b0'0001'0000
        * int c = a - b; //..0b0'1111'0000
        *
        * the previous numbers mean "I want 8 bits, but the first 4
        * low bits must be zeros". And if you keep subtracting b from
        * c, then you will traverse the following bit-patterns:
        *
        * c -= b; // 0b1110'0000
        * c -= b; // 0b1101'0000
        * c -= b; // 0b1100'0000
        * and so on */

        // the bits of the "right part" in the previous explanation
        const int code_bits {code >> bits_consumed_by_previous_area};
        int a {1 << current_area_bitwidth};
        const int b {1 << (current_code_bit_length - bits_consumed_by_previous_area)};
        do {
            a -= b;
            decoding_table.entries[area_offset + (a + code_bits)] = entry;
        } while(a != 0);

        /* this code goes to the next Huffman code (bit-reversed).
        * thanks to the way in which the 'symbols' vector is sorted,
        * "going to the next Huffman code" just means increasing the
        * Huffman code to the next valid value */
        icode += 1u << (16 - current_code_bit_length);
        code = static_cast<int>(reverse_bits(static_cast<std::uint16_t>(icode)));

        // do some book-keeping
        bl_count[current_code_bit_length] -= 1; // one code less to take care off
        if(bl_count[current_code_bit_length] == 0) {
            // go to the next bit-length
            if(current_code_bit_length == max) break;
            /* we cannot just write ++current_code_bit_length because
            * there can be bl_counts that are equal to zero between
            * bl_counts that aren't */
            current_code_bit_length = code_lengths[symbols[i + 1]];
        }

        /* go to the next area of the decoding table if needed. It is
        * needed every time the prefix changes and we cannot stay in
        * the first area */
        if(current_code_bit_length > first_area_bitwidth and (code & prefix_mask) != current_prefix) {
            // do some book-keeping
            current_prefix = code & prefix_mask;
            area_offset += 1 << current_area_bitwidth; // go to the next area
            // despite the variable name, it's always the bit-width of the first area
            bits_consumed_by_previous_area = first_area_bitwidth;

            // determine the bit-width of the next area
            current_area_bitwidth = current_code_bit_length - bits_consumed_by_previous_area;
            int free_nodes {1 << current_area_bitwidth}; // simulate a Huffman tree
            int bit_length {current_code_bit_length};
            while(bit_length < max) {
                free_nodes -= bl_count[bit_length];
                if(free_nodes <= 0) break;
                ++bit_length;
                ++current_area_bitwidth;
                free_nodes <<= 1;
            }

            // adjust the entry
            entry.category = Huffman_entry::Category::bridge;
            entry.value = area_offset;
            entry.bits_to_consume = current_area_bitwidth;
            decoding_table.entries[current_prefix] = entry;
        }
    }

    return decoding_table;
}

int rez::impl::deflate::fetch_symbol(const Decoding_table& decoding_table, Deflate_bitstream& bitstream)
{
    Huffman_entry entry {decoding_table.entries[bitstream.peek_bits(decoding_table.first_area_bitwidth)]};
    switch(entry.category) {
        case Huffman_entry::Category::symbol:
            bitstream.skip_bits(entry.bits_to_consume);
            return entry.value;
        case Huffman_entry::Category::bridge: {
            bitstream.skip_bits(decoding_table.first_area_bitwidth);
            Huffman_entry entry2 {decoding_table.entries[entry.value + bitstream.peek_bits(entry.bits_to_consume)]};
            bitstream.skip_bits(entry2.bits_to_consume);
            return entry2.value;
        }
        default:
            throw Exception {Error::bad_formed_data};
    }
}

void rez::impl::deflate::process_symbols(std::vector<std::uint8_t>& inflated_data, const Decoding_table& literal_length_alphabet, const Decoding_table& distance_alphabet, Deflate_bitstream& bitstream)
{
    int symbol {fetch_symbol(literal_length_alphabet, bitstream)};
    while(symbol != 256) {
        if(symbol < 256) { inflated_data.push_back(static_cast<std::uint8_t>(symbol)); }
        else {
            if(symbol > 285) throw Exception {Error::bad_formed_data};
            symbol -= 257;
            const int length {length_bases[symbol] + static_cast<int>(bitstream.read_bits(length_extra_bits[symbol]))};

            symbol = fetch_symbol(distance_alphabet, bitstream);
            if(symbol > 29) throw Exception {Error::bad_formed_data};
            const std::int32_t distance {distance_bases[symbol] + bitstream.read_bits(distance_extra_bits[symbol])};
            if(distance > inflated_data.size()) throw Exception {Error::bad_formed_data};

            lz77_copy(inflated_data, length, distance);
        }

        symbol = fetch_symbol(literal_length_alphabet, bitstream);
    }
}

void rez::impl::deflate::lz77_copy(std::vector<std::uint8_t>& inflated_data, const int length, const std::int32_t distance)
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
