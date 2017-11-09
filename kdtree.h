#ifndef KDTREE_H_
#define KDTREE_H_

#include <memory>
#include <limits>
#include <queue>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>

namespace spatial_index {

namespace util {

// Some compile time recursion in order to get a dimension dynamically
template <typename Point, std::size_t Dimension, std::size_t Count>
struct dimension_extractor {
    static inline typename boost::geometry::default_distance_result<Point>::type subtract(const Point &p1, const Point &p2, std::size_t dim) {
        if (Dimension == dim) {
            return boost::geometry::get<Dimension>(p1) - boost::geometry::get<Dimension>(p2);
        }
        return dimension_extractor<Point, Dimension + 1, Count>::subtract(p1, p2, dim);
    }
};
// end recursion
template <typename Point, std::size_t Count>
struct dimension_extractor<Point, Count, Count> {
    static inline typename boost::geometry::default_distance_result<Point>::type subtract(const Point &p1, const Point &p2, std::size_t dim) {
    }
};

template <typename Point>
typename boost::geometry::default_distance_result<Point>::type subtract(const Point &p1, const Point &p2, std::size_t dim) {
    return dimension_extractor<Point, 0, boost::geometry::dimension<Point>::type::value>::subtract(p1, p2, dim);
}

} // namespace util

template <typename Data, typename Point = boost::geometry::model::d2::point_xy<double>, typename Range = boost::geometry::model::box<Point>>
class kdtree {
public:
    kdtree() {}
    virtual ~kdtree() {}

    void reserve(const std::size_t size){
        m_nodes.reserve(size);
    }

    void add(const Point point, const Data *data) {
        typename kdnode::ptr node = std::make_shared<kdnode>(point, data);
        m_nodes.push_back(node);
    }
    void build(const Range range){
        if (m_nodes.empty()) {
            return;
        }
        m_root = build(m_nodes, 0, range);
    }
    void clear() {
        m_root.reset();
        m_nodes.clear();
    }
    const Data *nearest_recursive(const Point &query) const {
        if (!m_root) {
            return NULL;
        }
        best_match best(m_root, std::numeric_limits<double>::max());
        nearest(query, m_root, best);
        return best.node->data;
    }
    void knearest(const Point &query, size_t k, std::vector<const Data*> &result) const {
        if (!m_root || k < 1) {
            return;
        }
        MaxPriorityQueue tmp;
        knearest(query, m_root, k, tmp);
        size_t size = tmp.size();
        result.resize(size);
        for (size_t i = 0; i < size; i++) {
            // Reverse order
            result[size - i - 1] = tmp.top().second->data;
            tmp.pop();
        }
    }
    const Data *nearest_iterative(const Point &query) const {
        if (!m_root) {
            return NULL;
        }
        MinPriorityQueue pq;
        best_match best(m_root, std::numeric_limits<double>::max());
        pq.push(DistanceTuple(0, m_root));
        while (!pq.empty()) {
            const auto current = pq.top();
            if (current.first >= best.distance) {
                return best.node->data;
            }
            pq.pop();
            auto currentNode = current.second;
            double d = boost::geometry::comparable_distance(query, currentNode->split); // no sqrt
            double dx = util::subtract(query, currentNode->split, currentNode->axis);
            if (d < best.distance) {
                best.node = currentNode;
                best.distance = d;
            }
            node_ptr near = dx <= 0 ? currentNode->left : currentNode->right;
            node_ptr far = dx <= 0 ? currentNode->right : currentNode->left;
            if (far) pq.push(DistanceTuple(dx * dx, far));
            if (near) pq.push(DistanceTuple(0, near));
        }
        return best.node->data;
    }

    const std::vector<const Data *> range_search(const Range &range) const{
        std::vector<const Data *> result;
        range_search(result, m_root, range);
        return result;
    }

private:
    struct kdnode {
        typedef std::shared_ptr<kdnode> ptr;
        ptr left;
        ptr right;
        int axis;
        const Point split;
        const Data *data;
        Range range;
        kdnode(const Point g, const Data *d) : axis(0), split(g), data(d) {}
    };
    typedef typename kdnode::ptr node_ptr; // get rid of annoying typename
    typedef std::vector<node_ptr> Nodes;
    typedef std::pair<double, node_ptr> DistanceTuple;
    struct SmallestOnTop {
        bool operator()(const DistanceTuple &a, const DistanceTuple &b) const {
            return a.first > b.first;
        }
    };
    struct LargestOnTop {
        bool operator()(const DistanceTuple &a, const DistanceTuple &b) const {
            return a.first < b.first;
        }
    };
    typedef std::priority_queue<DistanceTuple, std::vector<DistanceTuple>, SmallestOnTop> MinPriorityQueue;
    typedef std::priority_queue<DistanceTuple, std::vector<DistanceTuple>, LargestOnTop> MaxPriorityQueue;
    Nodes m_nodes;
    node_ptr m_root;

    template<typename NODE_TYPE>
    struct Sort : std::binary_function<NODE_TYPE, NODE_TYPE, bool> {
        Sort(std::size_t dim) : m_dimension(dim) {}
        bool operator()(const NODE_TYPE &lhs, const NODE_TYPE &rhs) const {
            return util::subtract(lhs->split, rhs->split, m_dimension) < 0;
        }
        std::size_t m_dimension;
    };
    struct best_match {
        node_ptr node;
        double distance;
        best_match(const node_ptr &n, double d) : node(n), distance(d) {}
    };

    node_ptr build(Nodes &nodes, int depth, Range range) {
        if (nodes.empty()) {
            return node_ptr();
        }
        int axis = depth % boost::geometry::dimension<Point>();
        size_t median = nodes.size() / 2;
        std::nth_element(nodes.begin(), nodes.begin() + median, nodes.end(), Sort<node_ptr>(axis));
        node_ptr node = nodes.at(median);
        node->axis = axis;
        node->range = range;

        Range leftRange, rightRange;
        if(axis == 0){
            leftRange = Range(range.min_corner(), Point(node->split.x(), range.max_corner().y()));
            rightRange = Range(Point(node->split.x(), range.min_corner().y()), range.max_corner());
        }else{
            leftRange = Range(range.min_corner(), Point(range.max_corner().x(), node->split.y()));
            rightRange = Range(Point(range.min_corner().x(), node->split.y()), range.max_corner());
        }

        Nodes left(nodes.begin(), nodes.begin() + median);
        Nodes right(nodes.begin() + median + 1, nodes.end());
        node->left = build(left, depth + 1, leftRange);
        node->right = build(right, depth + 1, rightRange);

        return node;
    }

    static void nearest(const Point &query, const node_ptr &currentNode, best_match &best) {
        if (!currentNode) {
            return;
        }
        double d = boost::geometry::comparable_distance(query, currentNode->split); // no sqrt
        double dx = util::subtract(query, currentNode->split, currentNode->axis);
        if (d < best.distance) {
            best.node = currentNode;
            best.distance = d;
        }
        node_ptr near = dx <= 0 ? currentNode->left : currentNode->right;
        node_ptr far = dx <= 0 ? currentNode->right : currentNode->left;
        nearest(query, near, best);
        if ((dx * dx) >= best.distance) {
            return;
        }
        nearest(query, far, best);
    }
    template <typename PriorityQueue>
    static void knearest(const Point &query, const node_ptr &currentNode, size_t k, PriorityQueue &result) {
        if (!currentNode) {
            return;
        }
        double d = boost::geometry::comparable_distance(query, currentNode->split); // no sqrt
        double dx = util::subtract(query, currentNode->split, currentNode->axis);
        if (result.size() < k or d <= result.top().first) {
            result.push(DistanceTuple(d, currentNode));
            if (result.size() > k) {
                result.pop();
            }
        }
        node_ptr near = dx <= 0 ? currentNode->left : currentNode->right;
        node_ptr far = dx <= 0 ? currentNode->right : currentNode->left;
        knearest(query, near, k, result);
        if ((dx * dx)  >= result.top().first) {
            return;
        }
        knearest(query, far, k, result);
    }

    void report_subtree(std::vector<const Data *> &result, const node_ptr &currentNode) const {
        if(currentNode->left.get() != NULL){
            result.push_back(currentNode->left->data);
            report_subtree(result, currentNode->left);
        }
        if(currentNode->right.get() != NULL){
            result.push_back(currentNode->right->data);
            report_subtree(result, currentNode->right);
        }
    }

    void range_search(std::vector<const Data *> &result, const node_ptr &currentNode, const Range &searchRange) const{
        if(boost::geometry::within(currentNode->split, searchRange))
            result.push_back(currentNode->data);
        if(boost::geometry::within(currentNode->range, searchRange))
            report_subtree(result, currentNode);
        else{
            auto overlap = boost::geometry::intersects(currentNode->range, searchRange);
            if(overlap){
                if(currentNode->left.get() != NULL)
                    range_search(result, currentNode->left, searchRange);
                if(currentNode->right.get() != NULL)
                    range_search(result, currentNode->right, searchRange);
            }
        }
    }
};

} // namespace spatial_index

#endif /* KDTREE_H_ */
