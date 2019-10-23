//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_PARTITION_HPP
#define ADLBIRREG_PARTITION_HPP

#include <algorithm>
#include <GeometricUtils.hpp>
#include <Utils.hpp>
#include <Communicator.hpp>
#include <set>

namespace lb{

double compute_mu(double grid_size, double max_normalized_load);
int get_rank_from_vertices(const std::array<int, 4>& vertices_id, const std::map<int, Communicator>& neighborhoods);

class Partition {
public:

    using DataIndex = int;
    using ProcRank = int;
    using VertIndex = int;


    const Box3 d;
    double grid_size;
    Point_3 coord;

    std::set<int> neighbor_list;
    std::array<Point_3, 8>       vertices;
    std::array<int, 8>           vertices_id;
    std::array<Tetrahedron_3, 6> tetrahedra;
    std::map<int, Communicator>  vertex_neighborhood;
    std::array<Plane_3, 12> planes;

    Partition(Point_3 coord, const Box3& domain, double grid_size,
              Point_3 v1, Point_3 v2, Point_3 v3, Point_3 v4,
              Point_3 v5, Point_3 v6, Point_3 v7, Point_3 v8)
            : coord(coord), d(domain), grid_size(grid_size),
              vertices({std::move(v1), std::move(v2), std::move(v3), std::move(v4),
                        std::move(v5), std::move(v6), std::move(v7), std::move(v8)}){
        get_MPI_worldsize(num_part);
        for(int j = 0; j < 8; ++j)
            vertices_id[j] = get_vertex_id<int>(vertices[j], num_part);
        update();
    }

    template<class InputIterator>
    Partition(Point_3 coord, const Box3& domain, double grid_size,
              InputIterator beg)
            : coord(coord), d(domain), grid_size(grid_size),
              vertices(deserialize(beg)) {
        get_MPI_worldsize(num_part);
        for(int j = 0; j < 8; ++j)
            vertices_id[j] = get_vertex_id<int>(vertices[j], num_part);
        update();
    }

    std::array<double, 24> serialize(){
        std::array<double, 24> buffer;
        for(int i = 0; i < 8; ++i){
            buffer[i*3+0] = vertices[i].x();
            buffer[i*3+1] = vertices[i].y();
            buffer[i*3+2] = vertices[i].z();
        }
        return buffer;
    }

    template<class InputIterator>
    std::array<Point_3, 8> deserialize(InputIterator begin){
        std::array<Point_3, 8> points;
        for(int i = 0; i < 8; ++i) points[i] = Point_3(*(begin + i*3), *(begin + i*3+1), *(begin + i*3 + 2));
        return points;
    }

    void update(){
        construct_tetrahedra();
        construct_planes();
    }

    bool isGeometryValid() {
        return lb_isGeometryValid(vertices, planes, grid_size);
    }

    void construct_tetrahedra() {
        tetrahedra = get_tetrahedra(vertices);
    }
/*
 *     6---------7
 *    /|        /|
 *   / |       / |
 *  /  |      /  |
 * 4-- 2 ----5 --3
 * |  /      |  /
 * Y Z       | /
 * |/        |/
 * 0----X----1
 */
    void construct_planes() {
        planes =  get_planes(vertices);
    }

    double get_volume() {
        return (vertices[1].x()-vertices[0].x()) * (vertices[4].y()-vertices[0].y()) * (vertices[2].z()-vertices[0].z());
    }

    bool is_inside(const Point_3& p) {
        return in_tetrahedron(tetrahedra[0], p) || in_tetrahedron(tetrahedra[1], p) || in_tetrahedron(tetrahedra[2], p) || in_tetrahedron(tetrahedra[3], p) || in_tetrahedron(tetrahedra[4], p) || in_tetrahedron(tetrahedra[5], p);
    }

    template<class CartesianPointTransformer, class A>
    bool is_ghost(const CartesianPointTransformer& c, const A& p) {
        return is_ghost(c.transform(p));
    }

    bool is_ghost(const std::array<Plane_3, 12>& planes, const Point_3& p){
        std::array<double, 12> distances_to_plane;
        std::transform(planes.cbegin(), planes.cend(), distances_to_plane.begin(), [&p](Plane_3 plane){return CGAL::squared_distance(p, plane);});
        auto min_distance = std::sqrt(*std::min_element(distances_to_plane.cbegin(), distances_to_plane.cend()));
        return min_distance <= (std::sqrt(3) * grid_size);
    }

    std::vector<Plane_3> find_planes_closer_than(const std::array<Plane_3, 12>& planes, const Point_3& p, double cut_off) {
        std::array<double, 12> distances_to_plane;
        std::transform(planes.cbegin(), planes.cend(), distances_to_plane.begin(), [&p](Plane_3 plane){return std::sqrt(CGAL::squared_distance(p, plane));});
        std::vector<Plane_3> closest_planes;
        for(int i = 0; i < 12; ++i) {
            if(distances_to_plane[i] <= cut_off) closest_planes.push_back(planes[i]);
        }
        return closest_planes;
    }

    template<class CartesianPointTransformer, class A>
    void move_data(const std::vector<A>& elements) {
        CartesianPointTransformer c;
        std::array<std::vector<A>, 26> to_move;
        const auto size = elements.size();
        std::vector<std::pair<Plane_3 , int>> ghosts;
        for(int i = 0; i < size; ++i) {
            const auto& e = elements.at(i);
            auto closest_planes = find_planes_closer_than(planes, c.transform(elements[i]), sqrt_3*grid_size);

            if(closest_planes.size() > 0) {
                ghosts.push_back(std::make_pair(closest_planes, i));
            }

            if(this->is_inside(c.transform(e))) {

            } else {

            }
        }
    }

    template<class CartesianPointTransformer, class A>
    std::vector<std::pair<std::vector<Plane_3>, int> > get_ghosts_and_destination_planes(const std::vector<A>& elements) {
        CartesianPointTransformer c;
        std::vector<std::pair<std::vector<Plane_3>, int>> ghosts;
        const auto size = elements.size();
        for(int i = 0; i < size; ++i) {
            auto closest_planes = find_planes_closer_than(planes, c.transform(elements[i]), sqrt_3*grid_size);
            if(closest_planes.size() > 0)
                ghosts.push_back(std::make_pair(closest_planes, i));
        }
        return ghosts;
    }

    std::vector<std::pair<std::vector<int>, int>> compute_ghosts_destination_ranks(const std::vector<std::pair<std::vector<Plane_3>, int>>& ghosts_and_planes) {
        std::vector<std::pair<std::vector<int>, int>> ghost_and_ranks;
        for(auto const& ghost_and_planes : ghosts_and_planes) {
            std::vector<int> neighbors_of;
            for(auto& plane : ghost_and_planes.first) {
                auto vertices_on_plane = get_points_on_plane(vertices.begin(), vertices.end(), vertices_id.begin(), vertices_id.end(), plane);
                neighbors_of.push_back(get_rank_from_vertices(vertices_on_plane, vertex_neighborhood));
            }
            ghost_and_ranks.emplace_back(neighbors_of, ghost_and_planes.second);
        }
        return ghost_and_ranks;
    }

    std::array<std::pair<int, std::vector<Point_3>>, 8> spread_centers_of_load(const Point_3& cl, const std::map<int, Communicator>& neighborhood) {
        int N;
        std::array<std::pair<int, std::vector<Point_3>>, 8> all_cl;
        std::array<std::pair<int, std::vector<double>>,  8> buff_all_cl;
        std::array<MPI_Request, 8> requests;

        std::array<double, 3> my_cl = {cl.x(), cl.y(), cl.z()};
        int idx = 0;
        for(auto&  comm : neighborhood) {
            int i = comm.first;
            N     = comm.second.comm_size;
            all_cl[idx] = std::make_pair(i, std::vector<Point_3>());
            buff_all_cl[idx] = std::make_pair(i, std::vector<double>());
            auto& curr_buff_all_cl = (*search_in_linear_hashmap<8>(buff_all_cl, i)).second;
            auto& curr_all_cl      = (*search_in_linear_hashmap<8>(all_cl, i)).second;
            curr_buff_all_cl.resize(3*N);
            curr_all_cl.resize(N);
            comm.second.Allgather(my_cl.data(), 3, MPI_DOUBLE, curr_buff_all_cl.data(), 3, MPI_DOUBLE);
            for(int j = 0; j < N; ++j) {// reconstruct the points
                curr_all_cl[j] = Point_3(curr_buff_all_cl[j*3],curr_buff_all_cl[j*3+1],curr_buff_all_cl[j*3+2]);
            }
            idx++;
        }
        return all_cl;
    }

    template<class LoadComputer, class A>
    double get_load_imbalance(const std::vector<A>& elements, MPI_Comm neighborhood = MPI_COMM_WORLD){
        LoadComputer lc;
        double my_load = lc.compute_load(elements);
        int N; MPI_Comm_size(neighborhood, &N);
        auto all_loads = get_neighbors_load(my_load, neighborhood); //global load balancing with MPI_COMM_WORLD
        auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / N;
        return (*std::max_element(all_loads.cbegin(), all_loads.cend()) / avg_load) - 1;
    }

    template<class CartesianPointTransformer, class LoadComputer, class A>
    std::array<std::pair<ProcRank, std::vector<DataIndex>>, 26> move_vertices(std::vector<A>& elements, MPI_Comm neighborhood = MPI_COMM_WORLD) {

        LoadComputer lc;
        CartesianPointTransformer transformer;
        std::vector<Point_3> points;

        double my_load = lc.compute_load(elements);

        const std::vector<double>& weights = lc.get_weights(elements);

        std::transform(elements.cbegin(), elements.cend(), std::back_inserter(points), [&transformer](auto el){return transformer.transform(el);});

        int N; MPI_Comm_size(neighborhood, &N);

        auto all_loads = get_neighbors_load(my_load, neighborhood); //global load balancing with MPI_COMM_WORLD

        auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / N;

        std::vector<double> normalized_loads;

        std::transform(all_loads.cbegin(), all_loads.cend(), std::back_inserter(normalized_loads), [&avg_load](auto n){return n/avg_load;});

        auto max_normalized_load = *std::max_element(normalized_loads.cbegin(), normalized_loads.cend());

        double mu = compute_mu(grid_size, max_normalized_load);

        std::transform(all_loads.cbegin(), all_loads.cend(), std::back_inserter(normalized_loads), [&avg_load](auto v){return v / avg_load;});

        auto center_of_load = get_center_of_load(weights, points);

        std::array<std::pair<VertIndex, std::vector<Point_3>>, 8> all_cl = spread_centers_of_load(center_of_load, vertex_neighborhood);
        std::vector<int> vertex_local_id_to_move = {0, 1, 2, 3, 4, 5, 6, 7, 8};

        const auto original_mu = mu;
        const int MAX_TRIAL    = 1;
        unsigned int iteration = 0;

        // move the vertices, keep valid geometry among neighbors of a given vertex
        while(iteration < 8) {
            auto x = iteration&1; auto y = (iteration&2)>>1; auto z = (iteration&4)>>2;
            int com = (((unsigned int) coord.x()+x) & 1)  +
                      (((unsigned int) coord.y()+y) & 1)*2+
                      (((unsigned int) coord.z()+z) & 1)*4;

            for(int trial = 0; trial < MAX_TRIAL; ++trial) {
                auto new_vertices = vertices;
                auto v   = vertices[com];
                auto vid = vertices_id[com];
                auto cls = (*search_in_linear_hashmap<8, Point_3>(all_cl, vid)).second;
                auto f1  = get_vertex_force(v, cls, normalized_loads);
                auto f1_after = constraint_force(d, v, f1);
                new_vertices[com] = move_vertex(v, f1_after, mu);
                int valid = lb_isGeometryValid(new_vertices, get_planes(new_vertices), grid_size);
                const auto& communicator = vertex_neighborhood[vid];
                std::vector<int> allValid(communicator.comm_size);
                communicator.Allgather(&valid, 1, MPI_INT, allValid.data(), 1, MPI_INT);
                int are_all_valid = std::accumulate(allValid.begin(),allValid.end(),1,[](auto k,auto v){return k*v;});
                if(are_all_valid) {
                    vertices = new_vertices;
                    update();
                    iteration++;
                    mu = original_mu;
                    break;
                } else {
                    mu /= 2.0;
                    if(trial == MAX_TRIAL-1) iteration++;
                }
            }

        }
        //*******************************************************
        // at this point we are sure that:
        // 1) Vertices have moved
        // 2) Geometry is valid for all partitions

        // Now, to proceed to the next lb iteration we must:
        // 1) Tag our data (REAL, EMPTY, or GHOST)
        // 2) Move what does not belong to us anymore
        //*******************************************************

        // after moving all the vertices, move the data
        // create neighbors domain
        const auto my_domain = this->serialize();

        std::array<std::pair<ProcRank, std::vector<A>>, 26>  all_neighbors;
        std::array<std::pair<ProcRank, std::vector<DataIndex>>, 26> neighbor_ghosts;

        std::fill(all_neighbors.begin(), all_neighbors.end(), std::make_pair(-1, std::vector<A>()));
        std::transform(neighbor_list.begin(), neighbor_list.end(), all_neighbors.begin(),
                [](ProcRank i){ return std::make_pair(i, std::vector<A>()); });
        std::vector<double> recv_buff(26*8*3);
        std::vector<MPI_Request> srequests(26), rrequests(26);

        for(int i = 0; i < 26; ++i) {
            ProcRank neighbor_rank = all_neighbors[i].first;
            MPI_Isend(my_domain.data(), 24, MPI_DOUBLE, neighbor_rank, 9999, MPI_COMM_WORLD, &srequests[i]);
            MPI_Irecv((recv_buff.data() + i * 24), 24, MPI_DOUBLE, neighbor_rank, 9999, MPI_COMM_WORLD, &rrequests[i]);
        }
        MPI_Waitall(26, rrequests.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(26, srequests.data(), MPI_STATUSES_IGNORE);

        std::vector<Partition> neighborhoods;
        for(int i = 0; i < 26; ++i) {
            // building neighborhood of vertex
            neighborhoods.emplace_back(coord, d, grid_size, (recv_buff.begin() + i * 24));
        }

        //TAG MY DATA
        const auto nb_elements = elements.size();
        size_t data_id = 0;
        while(data_id < nb_elements) {
            auto point = transformer.transform(elements[data_id]);
            bool is_real_cell = is_inside(point);
            auto closest_planes = find_planes_closer_than(planes, point, sqrt_3*grid_size);
            if(closest_planes.size() > 0 && is_real_cell) { // IS A NEIGHBORING CELL, i will move it to them proc that share the closest planes
                //detect to whom I will periodically share my states
                for(int nid = 0; nid < 26; ++nid) {
                    if(all_neighbors[nid].first < 0) break;
                    const auto& neighborhood_domains = neighborhoods[nid];
                    bool found = false;
                    for(const Plane_3& close_plane : closest_planes){
                         found = std::any_of(neighborhood_domains.planes.cbegin(), neighborhood_domains.planes.cend(),
                                 [&close_plane](const Plane_3& p){return p == close_plane;});
                         if(found) break;
                    }
                    if(found) {
                        //add to ghost data for this proc

                    }
                }
            } else if(!is_real_cell){ // transfer data to neighbor
                for(int nid = 0; nid < 26; ++nid) {
                    if(all_neighbors[nid].first < 0) break;
                    const auto& neighborhood_domains = neighborhoods[nid];
                    const auto neighbor_rank  = all_neighbors[nid].first;
                    auto& neighbor_data = all_neighbors[nid].second;
                    if(neighborhood_domains.is_inside(point)) {
                        neighbor_data.push_back(elements[data_id]);
                        std::iter_swap(elements.begin() + data_id, elements.end() - 1);
                        elements.pop_back();
                    }
                }
            } else {
                // those are my data
            }
            data_id++;
        }

    }

    template<class NumericalType>
    const NumericalType get_vertex_id(const Point_3& v, int N) const {
        auto proc_per_row = (NumericalType) std::cbrt(N) + 1;
        const double step = (d.xmax - d.xmin) / (proc_per_row-1);
        return position_to_cell<NumericalType>(v, step, proc_per_row, proc_per_row);
    }

    void init_communicators(int world_size) {
        int N = world_size, x;
        get_MPI_rank(my_rank);
        std::array<int, 8>  sbuff;
        std::vector<int> buff(8*N);
        std::map<int, std::set<int>> neighbors;

        std::copy(vertices_id.begin(), vertices_id.end(), sbuff.begin());
        std::sort(sbuff.begin(), sbuff.end());

        MPI_Allgather(sbuff.data(), 8, MPI_INT, buff.data(), 8, MPI_INT, MPI_COMM_WORLD);
        for(int j = 0; j < 8; ++j) {
            int vid = vertices_id[j];
            for(int i = 0; i < N; ++i) {
                auto beg = buff.begin() +  i * 8;
                auto end = buff.begin() + (i+1) * 8;
                if(std::binary_search(beg, end, vid)) {
                    neighbors[vid].insert(i);
                     neighbor_list.insert(i);
                }
            }
            neighbors[vid].insert(my_rank);
            std::vector<int> n_list(neighbors[vid].begin(), neighbors[vid].end());
            Communicator c(n_list);
            vertex_neighborhood[vid] = c;
        }


    }

    friend std::ostream &operator<<(std::ostream &os, const Partition &partition) {
        os << " ID(" << partition.vertices_id[0] << "): " << partition.vertices[0]
           << " ID(" << partition.vertices_id[1] << "): " << partition.vertices[1]
           << " ID(" << partition.vertices_id[2] << "): " << partition.vertices[2]
           << " ID(" << partition.vertices_id[3] << "): " << partition.vertices[3]
           << " ID(" << partition.vertices_id[4] << "): " << partition.vertices[4]
           << " ID(" << partition.vertices_id[5] << "): " << partition.vertices[5]
           << " ID(" << partition.vertices_id[6] << "): " << partition.vertices[6]
           << " ID(" << partition.vertices_id[7] << "): " << partition.vertices[7];
        return os;
    }

};
}

#endif //ADLBIRREG_PARTITION_HPP
