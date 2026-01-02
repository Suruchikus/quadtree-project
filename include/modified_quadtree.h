// #pragma once

// #include <cstdint>
// #include <vector>
// #include <cassert>
// #include <limits>

// #include "tree_base.h"
// #include "rect.h"
// #include "quadtree.h"


// // Rects are half-open: [xmin, xmax) x [ymin, ymax)
// // Node regions are NOT stored per node. Region is derived during traversal.


// ///////////////////////////////////////////////////////////////
// // Bit vector with rank1 support
// ///////////////////////////////////////////////////////////////
// class BitVectorRank {
// public:
//     BitVectorRank() = default;

//     void reset(size_t bit_len) {
//         bit_len_ = bit_len;
//         blocks_.assign((bit_len + 63) / 64, 0ULL);
//         super_.clear();
//         total_ones_ = 0;
//         built_ = false;
//     }

//     size_t size() const { return bit_len_; }

//     void set(size_t i, bool value) {
//         assert(i < bit_len_);
//         const size_t w = i >> 6;  // /64
//         const size_t b = i & 63;  // %64
//         const uint64_t mask = 1ULL << b;
//         if (value) blocks_[w] |= mask;
//         else       blocks_[w] &= ~mask;
//         built_ = false;
//     }

//     bool get(size_t i) const {
//         assert(i < bit_len_);
//         const size_t w = i >> 6;
//         const size_t b = i & 63;
//         return (blocks_[w] >> b) & 1ULL;
//     }

//     void build_rank() {
//         static constexpr size_t SUPERBLOCK_BITS = 512; // 8 words
//         static constexpr size_t WORDS_PER_SUPER = SUPERBLOCK_BITS / 64;

//         const size_t nwords = blocks_.size();
//         const size_t nsupers = (nwords + WORDS_PER_SUPER - 1) / WORDS_PER_SUPER;

//         super_.assign(nsupers + 1, 0);
//         uint32_t running = 0;

//         for (size_t s = 0; s < nsupers; ++s) {
//             super_[s] = running;
//             const size_t start = s * WORDS_PER_SUPER;
//             const size_t end = (start + WORDS_PER_SUPER < nwords) ? (start + WORDS_PER_SUPER) : nwords;
//             for (size_t w = start; w < end; ++w) {
//                 running += popcnt64(blocks_[w]);
//             }
//         }
//         super_[nsupers] = running;
//         total_ones_ = running;
//         built_ = true;
//     }

//     // rank1(i): ones in [0..i]
//     uint32_t rank1(size_t i) const {
//         assert(built_);
//         assert(i < bit_len_);

//         static constexpr size_t SUPERBLOCK_BITS = 512;
//         static constexpr size_t WORDS_PER_SUPER = SUPERBLOCK_BITS / 64;

//         const size_t word = i >> 6;
//         const size_t bit  = i & 63;

//         const size_t super_id = word / WORDS_PER_SUPER;
//         uint32_t ans = super_[super_id];

//         const size_t start_word = super_id * WORDS_PER_SUPER;
//         for (size_t w = start_word; w < word; ++w) ans += popcnt64(blocks_[w]);

//         const uint64_t mask = (bit == 63) ? ~0ULL : ((1ULL << (bit + 1)) - 1ULL);
//         ans += popcnt64(blocks_[word] & mask);

//         return ans;
//     }

//     uint32_t ones() const {
//         assert(built_);
//         return total_ones_;
//     }

//     size_t bytes_used() const {
//         // blocks_ is uint64_t, super_ is uint32_t
//         return blocks_.size() * sizeof(uint64_t) + super_.size() * sizeof(uint32_t);
//     }



// private:
//     static inline uint32_t popcnt64(uint64_t x) {
//         #if defined(__GNUG__) || defined(__clang__)
//                 return (uint32_t)__builtin_popcountll(x);
//         #else
//                 uint32_t c = 0;
//                 while (x) { x &= (x - 1); ++c; }
//                 return c;
//         #endif
//             }

//     size_t bit_len_ = 0;
//     std::vector<uint64_t> blocks_;
//     std::vector<uint32_t> super_;
//     uint32_t total_ones_ = 0;
//     bool built_ = false;
// };


// class ModQuadTree : public BaseTree {
// public:
//     enum Quadrant : uint8_t { NW=0, NE=1, SW=2, SE=3 };
//     static constexpr uint32_t NIL = std::numeric_limits<uint32_t>::max();

//     ModQuadTree() = default;

//     // ---- BaseTree API ----
//     bool member(const Point& p) const override;
//     void range_query(const Rect& r, std::vector<Point>& out) const override;
//     Stats stats() const override;

//     // ---- Build API ----
//     static ModQuadTree build_from_pointer_tree(QuadTree::Node* root,
//                                            const Rect& root_region,
//                                            int maxLevel);



//     // ---- Space usage reporting ----
//     size_t bytes_used_total() const;
//     size_t bits_used_total() const { return bytes_used_total() * 8ULL; }

//     // Includes tree structure arrays + stored midpoints + leaf offsets.
//     // Excludes the actual point coordinates arrays if you want structure-only.
//     size_t bits_used_structure_only() const;

//     double bits_per_node_total() const {
//         return (numNodes_ == 0) ? 0.0 : (double)bits_used_total() / (double)numNodes_;
//     }

//     double bits_per_node_structure_only() const {
//         return (numNodes_ == 0) ? 0.0 : (double)bits_used_structure_only() / (double)numNodes_;
//     }

//     double bits_per_point_total() const {
//     const size_t P = pts_x_.size();   // number of stored points
//     if (P == 0) return 0.0;
//     return (double)(bytes_used_total() * 8ULL) / (double)P;
//     }

//     double bits_per_point_payload_only() const {
//         const size_t P = pts_x_.size();
//         if (P == 0) return 0.0;
//         const size_t payload_bytes = (pts_x_.size() + pts_y_.size()) * sizeof(int);
//         return (double)(payload_bytes * 8ULL) / (double)P; // ~64 bits/point if int32
//     }

//     double bits_per_point_structure_only() const {
//         const size_t P = pts_x_.size();
//         if (P == 0) return 0.0;
//         return (double)bits_used_structure_only() / (double)P;
//     }




// private:
//     // Root info
//     Rect root_;
//     int maxLevel_ = 0;

//     // Node counts
//     uint32_t numNodes_ = 0;      // all nodes in BFS order
//     uint32_t numInternals_ = 0;  // internal nodes only

//     // nodeId -> internalIndex (or NIL if leaf)
//     std::vector<uint32_t> node_to_internal_;

//     // Topology bitmap: 4 bits per internal node
//     // B_[4*i + q] = 1 iff internal node i has child in quadrant q
//     BitVectorRank B_;

//     // Each 1-bit in B_ corresponds to an existing child slot.
//     // child_slot_to_node_[k] stores the child nodeId for the k-th 1-bit (k=rank1-1).
//     std::vector<uint32_t> child_slot_to_node_;

//     // Leaf marker over ALL nodes (nodeId space)
//     BitVectorRank isLeaf_;

//     // Mode bits over INTERNAL nodes (internalIndex space)
//     BitVectorRank mxStored_;
//     BitVectorRank myStored_;

//     // Stored midpoint values, packed by rank among stored bits
//     std::vector<int> mx_vals_;
//     std::vector<int> my_vals_;

//     // Leaf points storage: flattened points + offsets
//     std::vector<uint32_t> leaf_off_;
//     std::vector<uint32_t> leaf_len_;
//     std::vector<int> pts_x_;
//     std::vector<int> pts_y_;

// private:
//     // half-open geometric midpoint
//     static inline int geom_mid_x(const Rect& r) { return (r.xmin + r.xmax) / 2; }
//     static inline int geom_mid_y(const Rect& r) { return (r.ymin + r.ymax) / 2; }

//     // Determine child region from parent + (mx,my) + quadrant
//     static inline Rect child_region(const Rect& parent, int mx, int my, Quadrant q) {
//         Rect c = parent;
//         switch (q) {
//             case NW: c.xmax = mx; c.ymin = my; return c;
//             case NE: c.xmin = mx; c.ymin = my; return c;
//             case SW: c.xmax = mx; c.ymax = my; return c;
//             case SE: c.xmin = mx; c.ymax = my; return c;
//         }
//         return c;
//     }

//     // Navigation + payload access (we implement inline now)
//     inline bool is_leaf_node(uint32_t nodeId) const {
//         assert(nodeId < numNodes_);
//         return isLeaf_.get(nodeId);
//     }

//     inline uint32_t child_from_internal(uint32_t internalIndex, Quadrant q) const {
//         const size_t pos = 4ULL * internalIndex + (uint8_t)q;
//         if (!B_.get(pos)) return NIL;
//         const uint32_t k = B_.rank1(pos) - 1;
//         return child_slot_to_node_[k];
//     }

//     inline uint32_t child(uint32_t nodeId, Quadrant q) const {
//         const uint32_t internalIndex = node_to_internal_[nodeId];
//         if (internalIndex == NIL) return NIL;
//         return child_from_internal(internalIndex, q);
//     }

//     inline int get_mx(uint32_t nodeId, const Rect& region) const {
//         const uint32_t internalIndex = node_to_internal_[nodeId];
//         assert(internalIndex != NIL);
//         if (!mxStored_.get(internalIndex)) return geom_mid_x(region);
//         const uint32_t k = mxStored_.rank1(internalIndex) - 1;
//         return mx_vals_[k];
//     }

//     inline int get_my(uint32_t nodeId, const Rect& region) const {
//         const uint32_t internalIndex = node_to_internal_[nodeId];
//         assert(internalIndex != NIL);
//         if (!myStored_.get(internalIndex)) return geom_mid_y(region);
//         const uint32_t k = myStored_.rank1(internalIndex) - 1;
//         return my_vals_[k];
//     }

//     inline void get_leaf_points(uint32_t nodeId, std::vector<Point>& out) const {
//         out.clear();
//         const uint32_t leafId = isLeaf_.rank1(nodeId) - 1;
//         const uint32_t off = leaf_off_[leafId];
//         const uint32_t len = leaf_len_[leafId];
//         out.reserve(len);
//         for (uint32_t i = 0; i < len; ++i) {
//             out.push_back(Point{pts_x_[off + i], pts_y_[off + i]});
//         }
//     }
// };

