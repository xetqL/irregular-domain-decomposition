//
// Created by xetql on 12/28/17.
//

#ifndef NBMPI_GEOMETRIC_ELEMENT_HPP
#define NBMPI_GEOMETRIC_ELEMENT_HPP

#include <array>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include "Utils.hpp"
#include "GeometricUtils.hpp"
#include "Communicator.hpp"
#include "Types.hpp"

namespace elements {

    using ElementRealType = type::Real;
    using Point3D = lb::Point_3;
    //using Box3D = boost::geometry::model::box<Point3D>;

    template<int N>
    struct Element {
        template<bool UseDoublePrecision = std::is_same<ElementRealType, double>::value>
        static CommunicationDatatype register_datatype() {
            MPI_Datatype element_datatype,
                    vec_datatype,
                    oldtype_range[1],
                    oldtype_element[1];

            MPI_Aint offset[1];

            int blockcount_element[1], blockcount_range[1];

            // register particle element type
            int array_size = N;
            auto mpi_raw_datatype = UseDoublePrecision ? MPI_DOUBLE : MPI_FLOAT;

            MPI_Type_contiguous(array_size, mpi_raw_datatype, &vec_datatype);

            MPI_Type_commit(&vec_datatype);

            blockcount_element[0] = 2; //position, velocity

            oldtype_element[0] = vec_datatype;

            //MPI_Type_extent(MPI_INT, &intex);
            offset[0] = static_cast<MPI_Aint>(0);

            //MPI_Type_struct(1, blockcount_element, offset, oldtype_element, &element_datatype);
            MPI_Type_create_struct(1, blockcount_element, offset, oldtype_element, &element_datatype);
            MPI_Type_commit(&element_datatype);

            blockcount_range[0] = N;
            oldtype_range[0] = mpi_raw_datatype;

            return {element_datatype, vec_datatype};
        }

        std::array<ElementRealType, N> position, velocity;

        static const int number_of_dimensions = N;

        constexpr Element(std::array<ElementRealType, N> p, std::array<ElementRealType, N> v) :
            position(p),
            velocity(v) {}

        constexpr Element(const std::array<ElementRealType, N>& p, const std::array<ElementRealType, N>& v, std::array<ElementRealType,N> a) :
            position(p.cbegin(), p.cend()),
            velocity(v.cbegin(), v.cend()) {}

        constexpr Element() = default;

        /**
         * Total size of the structure
         * @return The number of element per dimension times the number of characteristics (3)
         */
        static constexpr int size() {
            return N * 2;
        }

        static constexpr int byte_size() {
            return N * 2 * sizeof(ElementRealType);
        }

        static Element<N> create(std::array<ElementRealType, N> &p, std::array<ElementRealType, N> &v, int gid, int lid){
            Element<N> e(p, v, gid, lid);
            return e;
        }

        static Element<N> createc(std::array<ElementRealType, N> p, std::array<ElementRealType, N> v, int gid, int lid){
            Element<N> e(p, v, gid, lid);
            return e;
        }

        template<class Distribution, class Generator>
        static Element<N> create_random( Distribution& dist, Generator &gen, int gid, int lid){
            std::array<ElementRealType, N> p, v;
            std::generate(p.begin(), p.end(), [&dist, &gen](){return dist(gen);});
            std::generate(v.begin(), v.end(), [&dist, &gen](){return dist(gen);});
            return Element::create(p, v, gid, lid);
        }

        template<class Distribution, class Generator>
        static Element<N> create_random( Distribution& distx, Distribution& disty, Distribution& distz, Generator &gen, int gid, int lid){
            std::array<ElementRealType, N> p, v;
            p[0] = distx(gen);
            p[1] = disty(gen);
            if(N > 2) p[2] = distz(gen);

            return Element::create(p, v, gid, lid);
        }

        template<class Distribution, class Generator, class RejectionPredicate>
        static Element<N> create_random( Distribution& dist, Generator &gen, int gid, int lid, RejectionPredicate pred){
            std::array<ElementRealType, N> p, v;
            //generate point in N dimension
            int trial = 0;
            do {
                if(trial >= 1000) throw std::runtime_error("Could not generate particles that satisfy the predicate. Try another distribution.");
                std::generate(p.begin(), p.end(), [&dist, &gen](){return dist(gen);});
                trial++;
            } while(!pred(p));
            //generate velocity in N dimension
            std::generate(v.begin(), v.end(), [&dist, &gen](){return dist(gen);});
            return Element::create(p, v, gid, lid);
        }

        template<class Distribution, class Generator, int Cnt>
        static std::array<Element<N>, Cnt> create_random_n( Distribution& dist, Generator &gen ){
            std::array<Element<N>, Cnt> elements;
            int id = 0;
            std::generate(elements.begin(), elements.end(), [&]() mutable {
                return Element<N>::create_random(dist, gen, id, id++);
            });
            return elements;
        }

        template<class Container, class Distribution, class Generator>
        static void create_random_n(Container &elements, Distribution& dist, Generator &gen ) {
            int id = 0;
            std::generate(elements.begin(), elements.end(), [&]()mutable {
                return Element<N>::create_random(dist, gen, id, id++);
            });
        }

        template<class Container, class Distribution, class Generator, class RejectionPredicate>
        static void create_random_n(Container &elements, Distribution& dist, Generator &gen, RejectionPredicate pred ) {
            //apply rejection sampling using the predicate given in parameter
            //construct a new function that does the work
            int id = 0;
            std::generate(elements.begin(), elements.end(), [&]() mutable {
                return Element<N>::create_random(dist, gen, id, id++, [&elements, pred](auto const el) {
                       return std::all_of(elements.begin(), elements.end(), [&](auto p){return pred(p.position, el);});
                });
            });
        }

        /**
         * Equality of two elements regarding the VALUE of their properties
         * @param rhs another element
         * @return true if the position and the velocity of the two elements are equals
         */
        bool operator==(const Element &rhs) const {
            return position == rhs.position;
        }

        bool operator!=(const Element &rhs) const {
            return !(rhs == *this);
        }

        friend std::ostream &operator<<(std::ostream &os, const Element &element) {
            std::string pos = std::to_string(element.position.at(0));
            for(int i = 1; i < N; i++){
                pos += " " + std::to_string(element.position.at(i));
            }

            std::string vel = std::to_string(element.velocity.at(0));
            for(int i = 1; i < N; i++){
                vel += " " + std::to_string(element.velocity.at(i));
            }

            os << pos << ";" << vel;
            return os;
        }


    };



    template<int N>
    void import_from_file_float(std::string filename, std::vector<Element<N>>& particles) {

        std::ifstream pfile;
        pfile.open(filename, std::ifstream::in);
        if(!pfile.good()) throw std::runtime_error("bad particle file");

        std::string line;
        while (std::getline(pfile, line)) {
            auto parameters = split(line, ';');
            auto str_pos = split(parameters[0], ' ');
            auto str_vel = split(parameters[1], ' ');
            auto str_gid = parameters[3];
            auto str_lid = parameters[4];
            Element<N> e;

            for(int i = 0; i < N; ++i)
                e.position[i] = std::stof(str_pos[i], 0);
            for(int i = 0; i < N; ++i)
                e.velocity[i] = std::stof(str_vel[i], 0);

            e.gid = std::stoi(str_gid);
            e.lid = std::stoi(str_lid);
            particles.push_back(e);
        }
    }

    template<int N>
    void import_from_file_double(std::string filename, std::vector<Element<N>>& particles) {

        std::ifstream pfile;
        pfile.open(filename, std::ifstream::in);
        if(!pfile.good()) throw std::runtime_error("bad particle file");

        std::string line;
        while (std::getline(pfile, line)) {
            auto parameters = split(line, ';');
            auto str_pos = split(parameters[0], ' ');
            auto str_vel = split(parameters[1], ' ');
            auto str_gid = parameters[3];
            auto str_lid = parameters[4];
            Element<N> e;

            for(int i = 0; i < N; ++i)
                e.position[i] = std::stod(str_pos[i], 0);
            for(int i = 0; i < N; ++i)
                e.velocity[i] = std::stod(str_vel[i], 0);

            e.gid = std::stoi(str_gid);
            e.lid = std::stoi(str_lid);
            particles.push_back(e);
        }
    }

    template<int N, class RealType, bool UseDoublePrecision = std::is_same<RealType, double>::value>
    void import_from_file(std::string filename, std::vector<Element<N>>& particles) {
        if(UseDoublePrecision) import_from_file_double<N>(filename, particles);
        else import_from_file_float<N>(filename, particles);
    }

    template<int N>
    void export_to_file(std::string filename, const std::vector<Element<N>> elements) {
        std::ofstream particles_data;
        if (file_exists(filename)) std::remove(filename.c_str());
        particles_data.open(filename, std::ofstream::out);
        for(auto const& e : elements) {
            particles_data << e << std::endl;
        }
        particles_data.close();
    }

    template<int N>
    const inline ElementRealType distance2(const std::array<ElementRealType, N>& e1, const std::array<ElementRealType, N>& e2)  {
        Point3D point1(e1.at(0),  e1.at(1), N > 2 ? e1.at(2) : 0.0);
        Point3D point2(e2.at(0),  e2.at(1), N > 2 ? e2.at(2) : 0.0);
        return CGAL::squared_distance(point1, point2);
    }

    template<int N>
    const inline ElementRealType distance2(const elements::Element<N> &e1, const elements::Element<N> &e2)  {
        return elements::distance2<N>(e1.position, e2.position);
    }

    template<int N, typename T>
    std::vector<Element<N>> build_elements(const int length, const T* positions, const T* velocities) {
        std::vector<Element<N>> elements;
        elements.reserve(length);
        for(int i=0; i < length; ++i){
            Element<N> e({positions[2*i], positions[2*i+1]}, {velocities[2*i],velocities[2*i+1]}, i, i);
            elements.push_back(e);
        }
        return elements;
    }

//    template<int N, typename T>
//    void build_elements(std::vector<Element<N>>& elements, const int length, const T* positions, const T* velocities) throw() {
//        if(elements.empty()) {
//            throw std::runtime_error("Can not transform data into an empty vector");
//        }
//        std::generate(elements.begin(), elements.end(), [i = 0, id=0, &positions, &velocities]() mutable {
//            Element<N> e({positions[i], positions[i+1]}, {velocities[i], velocities[i+1]}, id, id);
//            i = i+N;
//            id++;
//            return e;
//        });
//    }

    template<int N>
    void serialize_positions(const std::vector<Element<N>>& elements, ElementRealType* positions){
        size_t element_id = 0;
        for (auto const& el : elements){
            for(size_t dim = 0; dim < N; ++dim)
                positions[element_id * N + dim] = el.position.at(dim);
            element_id++;
        }
    }

    template<int N, typename T>
    void serialize(const std::vector<Element<N>>& elements, T* positions, T* velocities, T* acceleration){
        size_t element_id = 0;
        for (auto const& el : elements){
            for(size_t dim = 0; dim < N; ++dim){
                positions[element_id * N + dim] = (ElementRealType) el.position.at(dim);  positions[element_id * N + dim] = (ElementRealType) el.position.at(dim);
                velocities[element_id * N + dim] = (ElementRealType) el.velocity.at(dim); velocities[element_id * N + dim] = (ElementRealType) el.velocity.at(dim);
                acceleration[element_id * N + dim] = (ElementRealType) el.acceleration.at(dim);  acceleration[element_id * N + dim] = (ElementRealType) el.acceleration.at(dim);
            }
            element_id++;
        }
    }

//    template<int N>
//    bool is_inside(const Element<N> &element, const std::array<std::pair<ElementRealType, ElementRealType>, N> domain){
//        auto element_position = element.position;
//
//        for(size_t dim = 0; dim < N; ++dim){
//            if(element_position.at(dim) < domain.at(dim).first || domain.at(dim).second < element_position.at(dim)) return false;
//        }
//
//        return true;
//    }

    template <int N, typename RealType>
    void init_particles_random_v(std::vector<elements::Element<N>> &elements, RealType T0, int seed = 0) {
        auto n = elements.size();
        std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
        std::normal_distribution<ElementRealType> ndist(0.0, T0 * T0);
        for (int i = 0; i < n; ++i) {
            elements[i].velocity[0] = ndist(gen);
            elements[i].velocity[1] = ndist(gen);
            if (N == 3) elements[i].velocity[2] = ndist(gen);
        }
    }
}

template<int N>
struct MESH_DATA {
    std::vector<elements::Element<N>> els;
};

#endif //NBMPI_GEOMETRIC_ELEMENT_HPP
