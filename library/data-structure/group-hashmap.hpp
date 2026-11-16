#pragma once
#if !defined(GROUP_HASHMAP_FORCE_SCALAR) && defined(__ARM_NEON) && (!defined(__BYTE_ORDER__) || __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
#include <arm_neon.h>
#define GROUP_HASHMAP_NEON
#elif !defined(GROUP_HASHMAP_FORCE_SCALAR) && defined(__SSE2__)
#include <emmintrin.h>
#define GROUP_HASHMAP_SSE2
#endif
// 整数キー用。16スロットのfingerprintをNEON/SSE2で同時比較する。
// 対応しない環境ではスカラー実装を使う。
// 最大負荷率7/8。再配置時は全iterator/referenceが無効になる。
// eraseは縮小を伴う場合がある。erase(key, false)は次iteratorの探索を省く。
// 全要素の走査が多い時速いらしい
template <typename Key, typename Val>
struct GroupHashMap {
    static_assert(is_integral<Key>::value, "GroupHashMap requires integral keys");
    using u64 = uint64_t;
    using Data = pair<Key, Val>;
    static inline u64 r = 0x9e3779b97f4a7c15ULL;
    size_t cap = 16;
private:
    static constexpr uint8_t EMPTY = 128, DELETED = 254;
    vector<Data> data;
    vector<uint8_t> ctrl;
    size_t s = 0, occupied = 0;
    static u64 hash(Key key) {
        u64 h = u64(key) + r;
        h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
        h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
        return h ^ (h >> 31);
    }
#if defined(GROUP_HASHMAP_NEON)
    using Group = uint8x16_t;
    Group group(size_t pos) const { return vld1q_u8(ctrl.data() + pos); }
    static u64 matches(Group word, uint8_t value) {
        auto eq = vceqq_u8(word, vdupq_n_u8(value));
        // 16個の比較結果を4bitずつに圧縮し、各要素につき1bitだけ残す。
        auto packed = vshrn_n_u16(vreinterpretq_u16_u8(eq), 4);
        return vget_lane_u64(vreinterpret_u64_u8(packed), 0) & 0x8888888888888888ULL;
    }
    static size_t first(u64 mask) { return __builtin_ctzll(mask) / 4; }
#elif defined(GROUP_HASHMAP_SSE2)
    using Group = __m128i;
    Group group(size_t pos) const { return _mm_loadu_si128(reinterpret_cast<const __m128i*>(ctrl.data() + pos)); }
    static u64 matches(Group word, uint8_t value) {
        return unsigned(_mm_movemask_epi8(_mm_cmpeq_epi8(word, _mm_set1_epi8(char(value)))));
    }
    static size_t first(u64 mask) { return __builtin_ctzll(mask); }
#else
    struct Group { uint8_t bytes[16]; };
    Group group(size_t pos) const { Group g; memcpy(g.bytes, ctrl.data() + pos, 16); return g; }
    static u64 matches(const Group& word, uint8_t value) {
        u64 mask = 0;
        for (int i = 0; i < 16; ++i) if (word.bytes[i] == value) mask |= u64(1) << i;
        return mask;
    }
    static size_t first(u64 mask) { return __builtin_ctzll(mask); }
#endif
    size_t lookup(Key key) const {
        u64 h = hash(key);
        size_t pos = (h >> 7) & (cap - 16);
        for (;;) {
            auto word = group(pos);
            u64 mask = matches(word, h & 127);
            while (mask) {
                size_t i = pos + first(mask);
                if (data[i].first == key) return i;
                mask &= mask - 1;
            }
            if (matches(word, EMPTY)) return cap;
            pos = (pos + 16) & (cap - 1);
        }
    }
    void rehash(size_t ncap) {
        auto old_data = move(data);
        auto old_ctrl = move(ctrl);
        cap = ncap;
        data = vector<Data>(cap);
        ctrl.assign(cap, EMPTY);
        for (size_t i = 0; i < old_ctrl.size(); ++i) if (old_ctrl[i] < EMPTY) {
            u64 h = hash(old_data[i].first);
            size_t pos = (h >> 7) & (cap - 16);
            while (!matches(group(pos), EMPTY)) pos = (pos + 16) & (cap - 1);
            size_t j = pos + first(matches(group(pos), EMPTY));
            data[j] = move(old_data[i]);
            ctrl[j] = h & 127;
        }
        occupied = s;
    }
    // 挿入位置と、新規挿入かどうかを返す。重複キーは上書きしない。
    pair<size_t, bool> prepare(Key key) {
        u64 h = hash(key);
        size_t pos = (h >> 7) & (cap - 16), dead = cap;
        for (;;) {
            auto word = group(pos);
            u64 mask = matches(word, h & 127);
            while (mask) {
                size_t i = pos + first(mask);
                if (data[i].first == key) return {i, false};
                mask &= mask - 1;
            }
            u64 deleted = matches(word, DELETED);
            if (dead == cap && deleted) dead = pos + first(deleted);
            u64 empty = matches(word, EMPTY);
            if (empty) {
                size_t i = dead == cap ? pos + first(empty) : dead;
                if (dead == cap && occupied + 1 > cap - cap / 8) {
                    rehash(s + 1 <= cap / 2 ? cap : cap * 2);
                    pos = (h >> 7) & (cap - 16);
                    dead = cap;
                    continue;
                }
                if (ctrl[i] == EMPTY) ++occupied;
                ctrl[i] = h & 127;
                ++s;
                return {i, true};
            }
            pos = (pos + 16) & (cap - 1);
        }
    }
public:
    template <bool Const> struct Iterator {
        using owner = conditional_t<Const, const GroupHashMap, GroupHashMap>;
        using value_type = Data;
        using difference_type = ptrdiff_t;
        using pointer = conditional_t<Const, const Data*, Data*>;
        using reference = conditional_t<Const, const Data&, Data&>;
        using iterator_category = bidirectional_iterator_tag;
        size_t i = 0;
        owner* p = nullptr;
        reference operator*() const { return p->data[i]; }
        pointer operator->() const { return &p->data[i]; }
        Iterator& operator++() {
            assert(i < p->cap);
            do { ++i; } while (i < p->cap && p->ctrl[i] >= EMPTY);
            return *this;
        }
        Iterator operator++(int) { auto old = *this; ++*this; return old; }
        Iterator& operator--() {
            do { assert(i > 0); --i; } while (p->ctrl[i] >= EMPTY);
            return *this;
        }
        Iterator operator--(int) { auto old = *this; --*this; return old; }
        bool operator==(const Iterator& other) const { return p == other.p && i == other.i; }
        bool operator!=(const Iterator& other) const { return !(*this == other); }
    };
    using iterator = Iterator<false>;
    using const_iterator = Iterator<true>;
    using itr = iterator;
    GroupHashMap() : data(cap), ctrl(cap, EMPTY) {}
    size_t size() const { return s; }
    bool empty() const { return s == 0; }
    iterator begin() { size_t i=0; while(i<cap && ctrl[i]>=EMPTY)++i; return {i,this}; }
    const_iterator begin() const { size_t i=0; while(i<cap && ctrl[i]>=EMPTY)++i; return {i,this}; }
    iterator end() { return {cap,this}; }
    const_iterator end() const { return {cap,this}; }
    iterator find(Key key) { return {lookup(key),this}; }
    const_iterator find(Key key) const { return {lookup(key),this}; }
    bool contain(Key key) const { return lookup(key) != cap; }
    int count(Key key) const { return contain(key); }
    iterator emplace(Key key, Val value) {
        auto [i, fresh] = prepare(key);
        if (fresh) data[i] = {key,move(value)};
        return {i,this};
    }
    iterator insert(const Data& value) { return emplace(value.first,value.second); }
    Val& operator[](Key key) {
        auto [i, fresh] = prepare(key);
        if (fresh) data[i] = {key,Val()};
        return data[i].second;
    }
    iterator erase(iterator it, bool get_next = true) {
        if (it == end()) return end();
        auto next = end();
        if (get_next) { next = it; ++next; }
        size_t pos = it.i & ~size_t(15);
        // 元から空きのあるgroupは探索の終端なので、空きに戻してよい。
        if (matches(group(pos), EMPTY)) { ctrl[it.i] = EMPTY; --occupied; }
        else ctrl[it.i] = DELETED;
        --s;
        if (cap > 16 && s <= cap / 10) {
            bool has_next = next != end();
            Key key = has_next ? next->first : Key();
            rehash(cap / 2);
            return has_next ? find(key) : end();
        }
        return next;
    }
    iterator erase(Key key, bool get_next = true) { return erase(find(key),get_next); }
    void clear() { fill(ctrl.begin(),ctrl.end(),EMPTY); s=occupied=0; }
    void reserve(int n) {
        if (n <= 0) return;
        size_t target = cap;
        while (size_t(n) > target - target / 8) target *= 2;
        if (target != cap) rehash(target);
    }
};

#undef GROUP_HASHMAP_NEON
#undef GROUP_HASHMAP_SSE2
