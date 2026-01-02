// #include "modified_quadtree.h"

// #include <iostream>
// #include <queue>
// #include <unordered_map>
// #include <utility>   

// ///////////////////////////////////////////////////////////////
// // 2) Build a succinct ModQuadTree from a pointer-based tree
// ///////////////////////////////////////////////////////////////
// ModQuadTree ModQuadTree::build_from_pointer_tree(QuadTree::Node* root,
//                                                  const Rect& root_region,
//                                                  int maxLevel)
// {
//     std::cout << "[succinct] build_from_pointer_tree: start\n";
//     std::cout.flush();

//     ModQuadTree out;
//     out.root_ = root_region;
//     out.maxLevel_ = maxLevel;

//     if (!root) {
//         out.numNodes_ = 0;
//         out.numInternals_ = 0;
//         return out;
//     }

//     // --- BFS enumerate nodes in BFS order and assign nodeId ---
//     std::vector<QuadTree::Node*> bfs_nodes;
//     bfs_nodes.reserve(1024);

//     std::queue<QuadTree::Node*> q;
//     q.push(root);

//     while (!q.empty()) {
//         QuadTree::Node* u = q.front(); q.pop();
//         u->tmp_id = (uint32_t)bfs_nodes.size();   // assign id here
//         bfs_nodes.push_back(u);

//         for (int k = 0; k < 4; ++k) {
//             if (u->child[k]) q.push(u->child[k]);
//         }
//     }
//     std::cout << "[succinct] BFS nodes = " << bfs_nodes.size() << "\n";
//     std::cout.flush();
//     out.numNodes_ = (uint32_t)bfs_nodes.size();


//     // Map pointer -> nodeId (needed to fill child_slot_to_node_)
//     std::unordered_map<const QuadTree::Node*, uint32_t> id_of;
//     id_of.reserve(bfs_nodes.size() * 2);

//     for (uint32_t i = 0; i < out.numNodes_; ++i) {
//         id_of[bfs_nodes[i]] = i;
//     }

//     // --- Determine internal nodes and create node_to_internal_ mapping ---
//     out.node_to_internal_.assign(out.numNodes_, NIL);

//     // Count internals first
//     for (uint32_t nodeId = 0; nodeId < out.numNodes_; ++nodeId) {
//         const QuadTree::Node* u = bfs_nodes[nodeId];
//         if (!u->is_leaf()) out.numInternals_++;
//     }

//     // Assign internal indices in BFS order
//     uint32_t internalIndex = 0;
//     for (uint32_t nodeId = 0; nodeId < out.numNodes_; ++nodeId) {
//         const QuadTree::Node* u = bfs_nodes[nodeId];
//         if (!u->is_leaf()) {
//             out.node_to_internal_[nodeId] = internalIndex++;
//         }
//     }
//     assert(internalIndex == out.numInternals_);

//     // --- Allocate bitvectors ---
//     out.isLeaf_.reset(out.numNodes_);
//     out.B_.reset(4ULL * out.numInternals_);
//     out.mxStored_.reset(out.numInternals_);
//     out.myStored_.reset(out.numInternals_);

//     out.child_slot_to_node_.clear();
//     out.child_slot_to_node_.reserve(out.numNodes_ * 2);

//     out.leaf_off_.clear();
//     out.leaf_len_.clear();
//     out.pts_x_.clear();
//     out.pts_y_.clear();

//     out.mx_vals_.clear();
//     out.my_vals_.clear();

//     uint32_t leafCount = 0;

//     for (uint32_t nodeId = 0; nodeId < out.numNodes_; ++nodeId) {
//         const QuadTree::Node* u = bfs_nodes[nodeId];

//         if (u->is_leaf()) {
//             out.isLeaf_.set(nodeId, true);
//             leafCount++;

//             const uint32_t off = static_cast<uint32_t>(out.pts_x_.size());
//             const uint32_t len = static_cast<uint32_t>(u->points.size());

//             out.leaf_off_.push_back(off);
//             out.leaf_len_.push_back(len);

//             out.pts_x_.reserve(out.pts_x_.size() + len);
//             out.pts_y_.reserve(out.pts_y_.size() + len);

//             for (const auto& p : u->points) {
//                 out.pts_x_.push_back(p.x);
//                 out.pts_y_.push_back(p.y);
//             }
//             continue;
//         }

//         // internal node
//         out.isLeaf_.set(nodeId, false);
//         const uint32_t iIdx = out.node_to_internal_[nodeId];
//         assert(iIdx != NIL);

//         // store mode bits (normal QuadTree has mode=0, but keep general)
//         const bool mxS = (u->mode & 1) != 0;
//         const bool myS = (u->mode & 2) != 0;

//         out.mxStored_.set(iIdx, mxS);
//         out.myStored_.set(iIdx, myS);

//         if (mxS) out.mx_vals_.push_back(u->mx_store);
//         if (myS) out.my_vals_.push_back(u->my_store);

//         // set topology bits and child_slot_to_node_
//         for (int k = 0; k < 4; ++k) {
//             if (u->child[k]) {
//                 const size_t pos = 4ULL * iIdx + (uint8_t)k;
//                 out.B_.set(pos, true);

//                 out.child_slot_to_node_.push_back(u->child[k]->tmp_id);

//             }
//         }
//     }

//     // --- Build rank support tables ---
//     out.B_.build_rank();
//     out.isLeaf_.build_rank();
//     out.mxStored_.build_rank();
//     out.myStored_.build_rank();

//     // sanity checks
//     assert(out.leaf_off_.size() == leafCount);
//     assert(out.leaf_len_.size() == leafCount);
//     assert(out.child_slot_to_node_.size() == out.B_.ones());
//     assert(out.mx_vals_.size() == out.mxStored_.ones());
//     assert(out.my_vals_.size() == out.myStored_.ones());

//     std::cout << "[succinct] build_from_pointer_tree: done\n";
//     std::cout.flush();

//     return out;
// }


// ///////////////////////////////////////////////////////////////
// // 3) member(p): descend using derived regions, then scan leaf points
// ///////////////////////////////////////////////////////////////
// bool ModQuadTree::member(const Point& p) const {
//     if (numNodes_ == 0) return false;
//     if (!root_.contains(p)) return false;

//     uint32_t nodeId = 0;
//     Rect region = root_;

//     while (!is_leaf_node(nodeId)) {
//         const int mx = get_mx(nodeId, region);
//         const int my = get_my(nodeId, region);

//         // Determine quadrant for half-open rectangles:
//         // x < mx => west; else east
//         // y < my => south; else north
//         const bool east  = (p.x >= mx);
//         const bool north = (p.y >= my);

//         Quadrant q;
//         if (!east && north) q = NW;
//         else if (east && north) q = NE;
//         else if (!east && !north) q = SW;
//         else q = SE;

//         const uint32_t nxt = child(nodeId, q);
//         if (nxt == NIL) return false;

//         region = child_region(region, mx, my, q);
//         nodeId = nxt;
//     }

//     // scan leaf points directly from flattened arrays
//     const uint32_t leafId = isLeaf_.rank1(nodeId) - 1;
//     const uint32_t off = leaf_off_[leafId];
//     const uint32_t len = leaf_len_[leafId];

//     for (uint32_t i = 0; i < len; ++i) {
//         if (pts_x_[off + i] == p.x && pts_y_[off + i] == p.y) return true;
//     }
//     return false;
// }

// ///////////////////////////////////////////////////////////////
// // 4) range_query(r): stack DFS, prune by intersection, scan leaves
// ///////////////////////////////////////////////////////////////
// void ModQuadTree::range_query(const Rect& r, std::vector<Point>& out) const {
//     out.clear();
//     if (numNodes_ == 0) return;
//     if (!root_.intersects(r)) return;

//     std::vector<std::pair<uint32_t, Rect>> st;
//     st.reserve(128);
//     st.push_back({0u, root_});

//     while (!st.empty()) {
//         auto [nodeId, region] = st.back();
//         st.pop_back();

//         if (!region.intersects(r)) continue;

//         if (is_leaf_node(nodeId)) {
//             const uint32_t leafId = isLeaf_.rank1(nodeId) - 1;
//             const uint32_t off = leaf_off_[leafId];
//             const uint32_t len = leaf_len_[leafId];

//             for (uint32_t i = 0; i < len; ++i) {
//                 Point p{pts_x_[off + i], pts_y_[off + i]};
//                 if (r.contains(p)) out.push_back(p);
//             }
//             continue;
//         }

//         const int mx = get_mx(nodeId, region);
//         const int my = get_my(nodeId, region);

//         // Push children that exist; you can push in reverse order if you want stable ordering
//         for (int k = 0; k < 4; ++k) {
//             const Quadrant qd = static_cast<Quadrant>(k);
//             const uint32_t cId = child(nodeId, qd);
//             if (cId == NIL) continue;

//             Rect cr = child_region(region, mx, my, qd);
//             if (cr.intersects(r)) st.push_back({cId, cr});
//         }
//     }
// }

// ///////////////////////////////////////////////////////////////
// // 5) stats(): depends on your Stats struct.
// // For now: return default object.
// // Paste stats.h and I’ll fill this properly.
// ///////////////////////////////////////////////////////////////
// Stats ModQuadTree::stats() const {
//     Stats s{};
//     // TODO: fill based on your Stats fields.
//     // Typical things people store:
//     // - s.numNodes = numNodes_;
//     // - s.numInternals = numInternals_;
//     // - s.numLeaves = isLeaf_.ones();
//     // - s.numPoints = pts_x_.size();
//     // - s.bytes = memory usage (sum of vectors/bitvectors)
//     return s;
// }



// size_t ModQuadTree::bytes_used_total() const {
//     size_t bytes = 0;

//     // bitvectors
//     bytes += B_.bytes_used();
//     bytes += isLeaf_.bytes_used();
//     bytes += mxStored_.bytes_used();
//     bytes += myStored_.bytes_used();

//     // structural arrays
//     bytes += node_to_internal_.size() * sizeof(uint32_t);
//     bytes += child_slot_to_node_.size() * sizeof(uint32_t);
//     bytes += leaf_off_.size() * sizeof(uint32_t);
//     bytes += leaf_len_.size() * sizeof(uint32_t);

//     // stored midpoint values
//     bytes += mx_vals_.size() * sizeof(int);
//     bytes += my_vals_.size() * sizeof(int);

//     // point arrays (data payload)
//     bytes += pts_x_.size() * sizeof(int);
//     bytes += pts_y_.size() * sizeof(int);

//     return bytes;
// }

// size_t ModQuadTree::bits_used_structure_only() const {
//     size_t bits = 0;

//     // raw bitmap bits (not counting rank tables) – useful for “pure succinct” view
//     bits += B_.size();
//     bits += isLeaf_.size();
//     bits += mxStored_.size();
//     bits += myStored_.size();

//     // packed-ish arrays (treat ints/uint32 as 32 bits)
//     bits += 32ULL * node_to_internal_.size();
//     bits += 32ULL * child_slot_to_node_.size();
//     bits += 32ULL * leaf_off_.size();
//     bits += 32ULL * leaf_len_.size();
//     bits += 32ULL * mx_vals_.size();
//     bits += 32ULL * my_vals_.size();

//     // exclude pts_x_/pts_y_ since that’s “stored data”, not structure
//     return bits;
// }

