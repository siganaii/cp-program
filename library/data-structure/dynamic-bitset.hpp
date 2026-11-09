#pragma once
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// 実行時に長さを変更できる bitset。添字 0 が最下位ビット。
// 二項演算は同じ長さ同士で行い、シフトでは長さを変えない。
// 1 ビットの操作は O(1)、全体の操作は O(ceil(size / 64))。
class DynamicBitset {
    using word_type = std::uint64_t;
    std::size_t n_ = 0;
    std::vector<word_type> words_;

    static std::size_t word_count(std::size_t n) { return n / 64 + (n % 64 != 0); }
    void trim() {
        if (n_ % 64 != 0) words_.back() &= (word_type(1) << (n_ % 64)) - 1;
    }
    template <class F>
    static void combine_words(word_type *__restrict dst,
                              const word_type *__restrict src,
                              std::size_t count, F operation) {
        const std::size_t blocks = count / 4;
        // 4 ワードずつ処理し、固定幅の SIMD 演算に変換しやすくする。
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
        for (std::size_t b = 0; b < blocks; ++b) {
            const std::size_t i = b * 4;
            dst[i] = operation(dst[i], src[i]);
            dst[i + 1] = operation(dst[i + 1], src[i + 1]);
            dst[i + 2] = operation(dst[i + 2], src[i + 2]);
            dst[i + 3] = operation(dst[i + 3], src[i + 3]);
        }
        for (std::size_t i = blocks * 4; i < count; ++i) dst[i] = operation(dst[i], src[i]);
    }
    template <class F>
    DynamicBitset &combine(const DynamicBitset &other, F operation) {
        assert(n_ == other.n_);
#if defined(__clang__)
        // Clang は通常のループを自動 SIMD 化できるため、その形を保つ。
        for (std::size_t i = 0; i < words_.size(); ++i)
            words_[i] = operation(words_[i], other.words_[i]);
#else
        // 小さい配列では SIMD 用のループ準備を省く。
        if (words_.size() <= 32) {
            for (std::size_t i = 0; i < words_.size(); ++i)
                words_[i] = operation(words_[i], other.words_[i]);
            return *this;
        }
        // 自己演算を分離し、SIMD 用ループでは領域が重ならないことを保証する。
        if (this == &other) {
            for (auto &word : words_) word = operation(word, word);
            return *this;
        }
        combine_words(words_.data(), other.words_.data(), words_.size(), operation);
#endif
        return *this;
    }

public:
    class reference {
        friend class DynamicBitset;
        word_type &word_;
        word_type mask_;
        reference(word_type &word, word_type mask) : word_(word), mask_(mask) {}
    public:
        operator bool() const { return (word_ & mask_) != 0; }
        reference &operator=(bool value) {
            if (value) word_ |= mask_;
            else word_ &= ~mask_;
            return *this;
        }
        reference &operator=(const reference &other) { return *this = bool(other); }
        bool operator~() const { return !bool(*this); }
        reference &flip() { word_ ^= mask_; return *this; }
    };

    // value は初期値。true の場合は全ビットを 1 にする。
    explicit DynamicBitset(std::size_t n = 0, bool value = false)
        : n_(n), words_(word_count(n), value ? ~word_type(0) : 0) { trim(); }
    DynamicBitset(const DynamicBitset &) = default;
    DynamicBitset &operator=(const DynamicBitset &) = default;
    DynamicBitset(DynamicBitset &&other) noexcept
        : n_(other.n_), words_(std::move(other.words_)) { other.clear(); }
    DynamicBitset &operator=(DynamicBitset &&other) noexcept {
        if (this != &other) {
            words_ = std::move(other.words_);
            n_ = other.n_;
            other.clear();
        }
        return *this;
    }

    std::size_t size() const { return n_; }
    bool empty() const { return n_ == 0; }

    // 下位のビットを保持し、増えた部分を value で埋める。
    // 容量の再確保がなければ、変更するワード数に比例する時間で処理する。
    // 再確保や縮小により、operator[] で得た参照が無効になる場合がある。
    void resize(std::size_t n, bool value = false) {
        std::size_t old = n_;
        words_.resize(word_count(n), value ? ~word_type(0) : 0);
        if (n > old && value && old % 64 != 0) {
            words_[old / 64] |= ~word_type(0) << (old % 64);
        }
        n_ = n;
        trim();
    }
    void reserve(std::size_t n) { words_.reserve(word_count(n)); }
    void clear() { words_.clear(); n_ = 0; }
    void push_back(bool value) {
        if (n_ % 64 == 0) words_.push_back(0);
        if (value) words_.back() |= word_type(1) << (n_ % 64);
        ++n_;
    }
    void pop_back() { assert(n_ != 0); resize(n_ - 1); }

    bool test(std::size_t pos) const {
        assert(pos < n_);
        return (words_[pos / 64] >> (pos % 64)) & 1;
    }
    bool operator[](std::size_t pos) const { return test(pos); }
    reference operator[](std::size_t pos) {
        assert(pos < n_);
        return reference(words_[pos / 64], word_type(1) << (pos % 64));
    }
    DynamicBitset &set(std::size_t pos, bool value = true) { (*this)[pos] = value; return *this; }
    DynamicBitset &reset(std::size_t pos) { return set(pos, false); }
    DynamicBitset &flip(std::size_t pos) {
        assert(pos < n_);
        words_[pos / 64] ^= word_type(1) << (pos % 64);
        return *this;
    }
    DynamicBitset &set() {
        std::fill(words_.begin(), words_.end(), ~word_type(0));
        trim();
        return *this;
    }
    DynamicBitset &reset() {
        std::fill(words_.begin(), words_.end(), 0);
        return *this;
    }
    DynamicBitset &flip() {
        for (auto &word : words_) word = ~word;
        trim();
        return *this;
    }
    std::size_t count() const {
        std::size_t result = 0;
        for (auto word : words_) result += __builtin_popcountll(word);
        return result;
    }
    bool any() const {
        for (auto word : words_) if (word) return true;
        return false;
    }
    bool none() const { return !any(); }
    bool all() const { return count() == n_; }

    DynamicBitset &operator&=(const DynamicBitset &other) {
        return combine(other, [](word_type a, word_type b) { return a & b; });
    }
    DynamicBitset &operator|=(const DynamicBitset &other) {
        return combine(other, [](word_type a, word_type b) { return a | b; });
    }
    DynamicBitset &operator^=(const DynamicBitset &other) {
        return combine(other, [](word_type a, word_type b) { return a ^ b; });
    }
    DynamicBitset &operator<<=(std::size_t shift) {
        if (shift >= n_) return reset();
        if (shift == 0) return *this;
        std::size_t whole = shift / 64, part = shift % 64;
        // 上位から更新し、読み出す元の値を保持する。
        for (std::size_t i = words_.size(); i-- > whole;) {
            word_type value = words_[i - whole] << part;
            if (part && i > whole) value |= words_[i - whole - 1] >> (64 - part);
            words_[i] = value;
        }
        std::fill(words_.begin(), words_.begin() + whole, 0);
        trim();
        return *this;
    }
    DynamicBitset &operator>>=(std::size_t shift) {
        if (shift >= n_) return reset();
        if (shift == 0) return *this;
        std::size_t whole = shift / 64, part = shift % 64;
        std::size_t remaining = words_.size() - whole;
        for (std::size_t i = 0; i < remaining; ++i) {
            word_type value = words_[i + whole] >> part;
            if (part && i + whole + 1 < words_.size()) value |= words_[i + whole + 1] << (64 - part);
            words_[i] = value;
        }
        std::fill(words_.begin() + remaining, words_.end(), 0);
        return *this;
    }
    friend DynamicBitset operator&(DynamicBitset a, const DynamicBitset &b) { return a &= b; }
    friend DynamicBitset operator|(DynamicBitset a, const DynamicBitset &b) { return a |= b; }
    friend DynamicBitset operator^(DynamicBitset a, const DynamicBitset &b) { return a ^= b; }
    friend DynamicBitset operator<<(DynamicBitset a, std::size_t shift) { return a <<= shift; }
    friend DynamicBitset operator>>(DynamicBitset a, std::size_t shift) { return a >>= shift; }
    DynamicBitset operator~() const { DynamicBitset result(*this); return result.flip(); }
    bool operator==(const DynamicBitset &other) const { return n_ == other.n_ && words_ == other.words_; }
    bool operator!=(const DynamicBitset &other) const { return !(*this == other); }

    // 見つからない場合は size()。find_next は pos より大きい添字を探す。
    std::size_t find_first() const {
        for (std::size_t i = 0; i < words_.size(); ++i) {
            if (words_[i]) return i * 64 + __builtin_ctzll(words_[i]);
        }
        return n_;
    }
    std::size_t find_next(std::size_t pos) const {
        if (n_ == 0 || pos >= n_ - 1) return n_;
        ++pos;
        std::size_t i = pos / 64;
        word_type word = words_[i] & (~word_type(0) << (pos % 64));
        if (word) return i * 64 + __builtin_ctzll(word);
        for (++i; i < words_.size(); ++i) {
            if (words_[i]) return i * 64 + __builtin_ctzll(words_[i]);
        }
        return n_;
    }
    // std::bitset と同様、最上位ビットから文字列にする。
    std::string to_string() const {
        std::string result(n_, '0');
        for (std::size_t i = 0; i < n_; ++i) if (test(i)) result[n_ - 1 - i] = '1';
        return result;
    }
};
