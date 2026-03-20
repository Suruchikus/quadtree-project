#pragma once

#include<cstdint>
#include<vector>
#include<array>
#include<string>

#include "point.h"
#include "rect.h"

class MXQuadtreeBits{
    public:
        struct Params{     //All the universal construction parameter we will use, like D is the depth of the quadtree, which is a constant once we know the size of grid
            int D = 0;
        };

        struct Stats{      //I will be measuring the bits per point, so I will declare all the stat parameters here
            uint64_t points=0;
            uint64_t N=0;
            int D=0;
            uint64_t T_bits = 0;
            uint64_t UM_bits = 0;
            uint64_t ULL_bits = 0;   //length bits for unary to leaf
            uint64_t ULD_bits = 0;   //dir bits for unary to leaf
            uint64_t UML_bits = 0;   //length from unary to mixed
            uint64_t UMD_bits = 0;   //dir bits for unary to mixed

            uint64_t total_nodes = 0;
            uint64_t unary_to_leaf_nodes = 0;
            uint64_t unary_to_mixed_nodes = 0;
            uint64_t internal_nodes = 0;
            uint64_t fullblock_nodes = 0;
            uint64_t leaf_nodes = 0;
            uint64_t mixed_internal = 0;

            double bpp = 0.0;
        };

        MXQuadtreeBits() = default;

        void build(const Rect& region, const std::vector<Point>& pts, const Params& params);   //Entry point for construction of the data structure
        bool contains(const Point& q) const;   //Membership Queries

        const Stats& stats() const {
            return stats_;
        }

    private:
        //We first define a node and node analysis. This is to build teh logical structute of the tree which will help build the compressed version later. So, these two are more of helpers
        struct Node{
            Rect r;
            int depth = 0;
            std::vector<int> ids;
        };

        struct NodeAnalysis{
            bool is_unit_leaf = false;
            bool is_fullblock = false;
            bool expandable = false;

            uint8_t childMask = 0;
            uint8_t nonempty_children = 0;

            std::array<std::vector<int>, 4> child_ids;
            std::array<Rect, 4> child_rects;
        };

        //We will need a query state to know the node depth and region details as we navigate
        struct QueryState{
            Rect r;
            int depth = 0;
            uint64_t t_pos =0;  //starting bit position of this node in T
        };

        Node make_root() const;
        Node make_child(const NodeAnalysis& a, const Node& parent, int child_idx) const;

        //The next two functions give the midpoint and quadrant index for building the logical quadtree
        static inline int midpoint(int a, int b){
            return (a+b) >> 1;
        }

        static inline int quadrant_index(const Point& p, int xm, int ym){
            const bool east = (p.x>=xm);
            const bool north = (p.y>=ym);

            if(!east && !north) return 0; //SW
            if(east && !north) return 1; //SE
            if(!east && north) return 2; //NW
            return 3; //NE
        }

        NodeAnalysis analyze_node(const Rect& r, int depth, const std::vector<int>& ids) const;
        
        Params params_;
        Stats stats_;
        Rect region_{0,0,0,0};
        std::vector<Point> points_;
        bool root_is_fullblock_=false;
        bool root_is_unitleaf_=false;

        //Bit Packed storage
        std::vector<uint64_t> T_;
        std::vector<uint64_t> UM_;
        std::vector<uint64_t> ULL_;
        std::vector<uint64_t> ULD_;
        std::vector<uint64_t> UML_;
        std::vector<uint64_t> UMD_;

        uint64_t T_len_=0;
        uint64_t UM_len_=0;
        uint64_t ULL_len_=0;
        uint64_t ULD_len_=0;
        uint64_t UML_len_=0;
        uint64_t UMD_len_=0;

        //bit helpers
        static void push_bits(std::vector<uint64_t>& dst, uint64_t& bit_len, uint64_t value, int width);
        static uint64_t get_bit(const std::vector<uint64_t>& src, uint64_t bit_pos);
        static uint64_t read_bits(const std::vector<uint64_t>& src, uint64_t bit_pos, int width);
        Rect child_rect(const Rect& r,int child_idx) const;

        //unary skip support
        struct UnarySkipResult{
            uint8_t L=0;
            std::vector<uint8_t> dirs;
            Node endpoint;
            NodeAnalysis endpoint_analysis;
        };

        UnarySkipResult follow_unary_chain(const Node& start) const;

        //contracted tree builder
        void build_bfs_contracted();

        //Rank/Select helpers
        uint64_t rank1_T(uint64_t bit_pos) const;               //Number of 1 bit before bit_pos in T
        uint64_t state_at_S(uint64_t child_index) const;  //Read the 2 bit state
        uint64_t rank_unary_S(uint64_t child_index) const;  //Count earlier unary entires
        uint64_t rank_continue_S(uint64_t child_index) const;  //Count earlier children that continue in T
        uint64_t next_node_id_for_child(uint64_t child_bit_pos) const; //Find the point in T
        uint64_t unary_path_offset(uint64_t unary_index) const;    //Find where the directions given the unary entry in U begin in P


};