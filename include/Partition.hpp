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
        {
            vertices_id[j] = get_vertex_id<int>(vertices[j], num_part);
            //std::cout <<vertices_id[j] << std::endl;
        }
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

    bool contains(const Point_3& p) const {
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

    inline std::vector<Plane_3> find_planes_closer_than(const std::array<Plane_3, 12>& planes, const Point_3& p, double cut_off) {
        std::array<double, 12> distances_to_plane;
        std::transform(planes.cbegin(), planes.cend(), distances_to_plane.begin(), [&p](Plane_3 plane){return CGAL::squared_distance(p, plane);});
        std::vector<Plane_3> closest_planes;
        for(int i = 0; i < 12; ++i) {
            if(distances_to_plane[i] <= cut_off*cut_off) closest_planes.push_back(planes[i]);
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

            if(this->contains(c.transform(e))) {

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

    LinearHashMap<int, std::vector<Point_3>, 8> spread_centers_of_load(const Point_3& cl, const std::map<int, Communicator>& neighborhood) {
        int N;
        LinearHashMap<int, std::vector<Point_3>, 8> all_cl;
        LinearHashMap<int, std::vector<double>,  8> buff_all_cl;
        std::array<MPI_Request, 8> requests;

        std::array<double, 3> my_cl{cl.x(), cl.y(), cl.z()};
        int idx = 0;
        for(auto&  comm : neighborhood) {
            int i = comm.first;
            N     = comm.second.comm_size;
            all_cl[idx] = std::make_pair(i, std::vector<Point_3>());
            buff_all_cl[idx] = std::make_pair(i, std::vector<double>());
            auto& curr_buff_all_cl = (*search_in_linear_hashmap<int, std::vector<double>,  8>(buff_all_cl, i)).second;
            auto& curr_all_cl      = (*search_in_linear_hashmap<int, std::vector<Point_3>, 8>(all_cl, i)).second;
            curr_buff_all_cl.resize(3*N);
            curr_all_cl.resize(N);
            comm.second.Allgather(my_cl.data(), 3, MPI_DOUBLE, curr_buff_all_cl.data(), 3, MPI_DOUBLE, i);
            for(int j = 0; j < N; ++j) { // reconstruct the points
                curr_all_cl[j] = Point_3(curr_buff_all_cl[j*3],curr_buff_all_cl[j*3+1],curr_buff_all_cl[j*3+2]);
            }
            idx++;
        }
        return all_cl;
    }
    using ElementCount   = int;
    using Load           = double;
    using Imbalance      = double;
    struct LoadStatistics : std::tuple<ElementCount, Load, Load, Imbalance, Imbalance> {
        ElementCount e;
        Load max_load, avg_load, my_load;
        Imbalance global, mine;
        LoadStatistics(ElementCount e, Load maxl, Load avgl, Load myl, Imbalance global, Imbalance mine):
            e(e), max_load(maxl), avg_load(avgl), my_load(myl), global(global), mine(mine), std::tuple<ElementCount, Load, Load, Imbalance, Imbalance>(e, maxl, avgl, global, mine){}
    };

    template<class LoadComputer, class A>
    LoadStatistics get_load_statistics(const std::vector<A>& elements, MPI_Comm neighborhood = MPI_COMM_WORLD){
        LoadComputer lc;
        int count = elements.size();
        int buf;
        double my_load = lc.compute_load(elements);
        int N; MPI_Comm_size(neighborhood, &N);
        auto all_loads = ::get_neighbors_load(my_load, neighborhood); //global load balancing with MPI_COMM_WORLD
        auto avg_load  = std::accumulate(all_loads.cbegin(), all_loads.cend(), 0.0) / N;
        auto max_load  =*std::max_element(all_loads.cbegin(), all_loads.cend());
        MPI_Allreduce(&count, &buf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        return {buf, max_load, avg_load, my_load, (max_load / avg_load) - 1.0, (my_load / avg_load) - 1.0};
    }

    template<class CartesianPointTransformer, class LoadComputer, class A>
    std::vector<std::pair<ProcRank, std::vector<DataIndex>>> move_vertices(
            std::vector<A>& elements, MPI_Datatype datatype, double avg_load, LinearHashMap<int, double, 8> vertex_mu, MPI_Comm neighborhood = MPI_COMM_WORLD) {

        get_MPI_rank(my_rank);
        LoadComputer lc;
        CartesianPointTransformer transformer;
        std::vector<Point_3> points;

        double my_load = lc.compute_load(elements);

        const std::vector<double>& weights = lc.get_weights(elements);

        std::transform(elements.cbegin(), elements.cend(), std::back_inserter(points), [&transformer](auto el){return transformer.transform(el);});

        int N; MPI_Comm_size(neighborhood, &N);

        auto all_loads = get_neighbors_load(my_load, neighbor_list, vertex_neighborhood); //global load balancing with MPI_COMM_WORLD

        auto center_of_load = get_center_of_load(weights, points);

        LinearHashMap<VertIndex, std::vector<Point_3>, 8> all_cl = spread_centers_of_load(center_of_load, vertex_neighborhood);

        std::vector<int> vertex_local_id_to_move = {0, 1, 2, 3, 4, 5, 6, 7, 8};

        const int MAX_TRIAL    = 3;
        unsigned int iteration = 0;
        // move the vertices, keep valid geometry among neighbors of a given vertex
        while(iteration < 8) {
            auto x = iteration&1; auto y = (iteration&2)>>1; auto z = (iteration&4)>>2;
            int com = (((unsigned int) coord.x()+x) & 1)  +
                      (((unsigned int) coord.y()+y) & 1)*2+
                      (((unsigned int) coord.z()+z) & 1)*4;

            auto v   = vertices[com];
            auto vid = vertices_id[com];
            const auto& communicator = vertex_neighborhood[vid];
            std::vector<double> normalized_loads, loads = (*search_in_linear_hashmap<int, std::vector<double>, 8>(all_loads, vid)).second;
            double mu = (*search_in_linear_hashmap<int, double, 8>(vertex_mu, vid)).second;
            std::transform(loads.cbegin(), loads.cend(), std::back_inserter(normalized_loads), [&avg_load](auto n){return n/avg_load;});
            auto cls = (*search_in_linear_hashmap<int, std::vector<Point_3>, 8>(all_cl, vid)).second;
            auto f1  = -get_vertex_force(v, cls, normalized_loads);
            auto f1_after = constraint_force(d, v, f1);
            if(communicator.comm_size > 1)
                for(int trial = 0; trial < MAX_TRIAL; ++trial) {
                    auto new_vertices = vertices;
                    new_vertices[com] = move_vertex(v, f1_after, mu);
                    int valid = lb_isGeometryValid(new_vertices, get_planes(new_vertices), grid_size);
                    std::vector<int> allValid(communicator.comm_size);
                    communicator.Allgather(&valid, 1, MPI_INT, allValid.data(), 1, MPI_INT, vid);
                    int are_all_valid = std::accumulate(allValid.begin(),allValid.end(),1,[](auto k,auto v){return k*v;});
                    if(are_all_valid) {
                        vertices = new_vertices;
                        update();
                        break;
                    } else {
                        mu /= 2.0;
                    }
                }
            iteration++;
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
        int number_of_neighbors = neighbor_list.size();

        std::vector<std::pair<ProcRank, std::vector<A>>>  migration_and_destination(number_of_neighbors);
        std::vector<std::pair<ProcRank, std::vector<DataIndex>>> neighbor_ghosts(number_of_neighbors);

        std::transform(neighbor_list.begin(), neighbor_list.end(), migration_and_destination.begin(),
                [](ProcRank i){ return std::make_pair(i, std::vector<A>()); });
        std::transform(neighbor_list.begin(), neighbor_list.end(), neighbor_ghosts.begin(),
                       [](ProcRank i){ return std::make_pair(i, std::vector<DataIndex>()); });

        std::vector<double> recv_buff(neighbor_list.size()*8*3);
        std::vector<MPI_Request> srequests(neighbor_list.size()), rrequests(neighbor_list.size());

        {
            int i = 0;
            for(ProcRank neighbor_rank : neighbor_list) {
                MPI_Isend(my_domain.data(), 24, MPI_DOUBLE, neighbor_rank, 9999, MPI_COMM_WORLD, &srequests[i]);
                MPI_Irecv((recv_buff.data() + i * 24), 24, MPI_DOUBLE, neighbor_rank, 9999, MPI_COMM_WORLD, &rrequests[i]);
                i++;
            }
        }

        MPI_Waitall(neighbor_list.size(), rrequests.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(neighbor_list.size(), srequests.data(), MPI_STATUSES_IGNORE);

        std::vector<Partition> neighborhoods;
        for(int i = 0; i < number_of_neighbors; ++i) {
            // building neighborhood of vertex
            neighborhoods.emplace_back(coord, d, grid_size, (recv_buff.begin() + i * 24));
        }

        //TAG MY DATA
        int DBG_RANK = 7;
        //const auto nb_elements = elements.size();
        size_t data_id = 0;
        size_t nb_migrate = 0;
        double sent_load = 0;
        while(data_id < elements.size()) {
            auto point = transformer.transform(elements[data_id]);
            bool is_real_cell = this->contains(point);

            auto closest_planes = find_planes_closer_than(planes, point, sqrt_3*grid_size);
            /*if(closest_planes.size() > 0 && is_real_cell) { // IS A NEIGHBORING CELL, i will move it to the proc that share the closest planes
                for(int nid = 0; nid < number_of_neighbors; ++nid) {
                    const auto& neighborhood_domains = neighborhoods[nid];
                    bool found = false;
                    for(const Plane_3& close_plane : closest_planes) {
                         found = std::any_of(neighborhood_domains.planes.cbegin(), neighborhood_domains.planes.cend(),
                                 [&close_plane](const Plane_3& p){return p == close_plane;});
                         if(found) break;
                    }
                    if(found) { // add to ghost data for this proc
                        neighbor_ghosts[nid].second.push_back(data_id);
                    }
                }
            } else*/ if(!is_real_cell){ // transfer data to neighbor
                //throw std::runtime_error("?");
                for(int nid = 0; nid < number_of_neighbors; ++nid) {
                    const auto  neighbor_rank = migration_and_destination[nid].first;
                    if(neighbor_rank < 0 || neighbor_rank == my_rank) continue;
                    const auto& neighborhood_domains = neighborhoods[nid];
                    if(neighborhood_domains.contains(point)) {
                        //std::cout << "Transferring " << point << " to " << neighbor_rank << std::endl;
                        migration_and_destination[nid].second.push_back(elements[data_id]);
                        std::iter_swap(elements.begin() + data_id, elements.end() - 1);
                        elements.pop_back();
                        nb_migrate++;
                        break;
                    }
                }
            } else {
                // those are my data
            }
            data_id++;
        }

        //std::cout << "Transfer " << sent_load << " to "<< DBG_RANK << std::endl;
        //std::cout << my_rank << " transferred " << nb_migrate << std::endl;

        {
            //migrate the data
            std::vector<MPI_Request> srequests(number_of_neighbors);
            for(int nid = 0; nid < number_of_neighbors; ++nid) {
                int dest_rank = migration_and_destination[nid].first;
                assert(migration_and_destination[nid].second.size() == 0 | dest_rank != my_rank);
                MPI_Isend(migration_and_destination[nid].second.data(),
                          migration_and_destination[nid].second.size(), datatype, dest_rank, 101010, MPI_COMM_WORLD, &srequests[nid]);
            }

            std::vector<A> buf;
            double recv_load = 0;
            for(int nid = 0; nid < number_of_neighbors; ++nid) {
                int dest_rank = migration_and_destination[nid].first;
                MPI_Status status; MPI_Probe(dest_rank, 101010, MPI_COMM_WORLD, &status);
                int count;
                MPI_Get_count(&status, datatype, &count);
                buf.resize(count);
                MPI_Recv(buf.data(), count, datatype, dest_rank, 101010, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                recv_load += lc.compute_load(buf);
                std::move(buf.begin(), buf.end(), std::back_inserter(elements));
            }
            //if(my_rank == DBG_RANK) std::cout << DBG_RANK<< " received " <<recv_load << std::endl;
            MPI_Waitall(number_of_neighbors, srequests.data(), MPI_STATUSES_IGNORE);
        }
        return neighbor_ghosts;
    }

    template<class NumericalType>
    const NumericalType get_vertex_id(const Point_3& v, int N) const {
        auto vertices_per_row = (NumericalType) std::cbrt(N)+1;
        const double step = (d.xmax - d.xmin) / (vertices_per_row);
//        std::cout
//        << (v.x() / step) << " + "
//        << vertices_per_row * v.z() / step << " + "
//        << vertices_per_row * vertices_per_row * v.y() / step << std::endl;
        //std::cout << v.y() << std::endl;
        return (v.x() / step) + vertices_per_row * v.z() / step  + vertices_per_row * vertices_per_row * v.y() / step;
        //return position_to_cell<NumericalType>(v, step, proc_per_row, proc_per_row);
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
