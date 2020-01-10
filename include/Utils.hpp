//
// Created by xetql on 10/21/19.
//

#ifndef ADLBIRREG_UTILS_HPP
#define ADLBIRREG_UTILS_HPP

#include <algorithm>
#include <vector>
#include <array>
#include <mpi.h>
#include <fstream>
#include <Communicator.hpp>
#include <set>
#include <random>
#include "GeometricUtils.hpp"
#include "Types.hpp"

std::vector<std::string> split(const std::string& s, char delimiter);
bool file_exists(const std::string fileName);

template<class Numeric, class InputOp>
Numeric mean(InputOp begin, InputOp end, int N) {
    return std::accumulate(begin, end, (Numeric) 0.0) / N;
}

template<class Numeric, class Container>
Numeric mean(Container container, int N) {
    Numeric acc = 0;
    for(Numeric v : container) acc += v;
    return acc / N;
}

template<typename Realtype, typename ContainerA>
Realtype get_slope(const ContainerA& y) {
    //std::cout << y.size() << std::endl;
    std::vector<Realtype> x(y.size());
    std::iota(x.begin(), x.end(), 0);

    double n = x.size();

    Realtype avgX = std::accumulate(x.begin(), x.end(), 0.0) / n;
    Realtype avgY = std::accumulate(y.begin(), y.end(), 0.0) / n;

    Realtype numerator = 0.0;
    Realtype denominator = 0.0;

    for(int i=0; i < n; ++i) {
        numerator   += (x[i] - avgX) * (y[i] - avgY);
        denominator += (x[i] - avgX) * (x[i] - avgX);
    }

    return denominator == 0 ? (Realtype) 0 : numerator / denominator;
}

template<class RealType, class Iter>
inline RealType mean(Iter b, Iter e) {
    const long N = std::distance(b, e);
    if(N > 0)
        return std::accumulate(b, e, (RealType) 0.0) / N;
    else
        return (RealType) 0.0;
}

template<class Realtype, class Iter>
Realtype median(Iter begin, Iter  end) {
    std::vector<typename Iter::value_type> tmp(begin, end);
    std::sort(tmp.begin(), tmp.end());
    if(std::distance(begin, end) == 0)
        return (Realtype) 0.0;
    else
        return (Realtype) tmp[(int) (tmp.size() / 2)];
}

template<class K, class V, int N>
using LinearHashMap = std::array<std::pair<K, V>, N>;

template<class K, class V, int N>
typename LinearHashMap<K, V, N>::iterator
search_in_linear_hashmap(std::array<std::pair<K, V>, N> &linmap, K key) {
    //for(std::pair<int, std::vector<A>>& entry : linmap) {
    auto entry = linmap.begin();
    for (; entry != linmap.end(); entry++) {
        if ((*entry).first == key) return entry;
    }
    return entry;
}

template<class K, class V, int N>
typename LinearHashMap<K, V, N>::const_iterator
search_in_linear_hashmap(const std::array<std::pair<K, V>, N> &linmap, K key) {
    //for(std::pair<int, std::vector<A>>& entry : linmap) {
    auto entry = linmap.begin();
    for (; entry != linmap.end(); entry++) {
        if ((*entry).first == key) return entry;
    }
    return entry;
}
std::set<int> filter_active_neighbors(const std::array<type::VertexIndex, 8>& vertices_id,
                                      const LinearHashMap<type::VertexIndex, int, 8>& vertices_trial,
                                      const std::map<int, Communicator>& vertex_neighborhood);
/**
 * To compute the local average load use local communicator, otherwise use MPI_COMM_WORLD
 * @param my_load
 * @param average_load
 * @param neighborhood
 */
inline void compute_average_load(const type::Real my_load, type::Real *average_load, MPI_Comm neighborhood) {
    int N;
    MPI_Comm_size(neighborhood, &N);
    MPI_Allreduce(&my_load, average_load, 1, MPI_DOUBLE, MPI_SUM, neighborhood);
    *average_load /= N;
}

inline std::vector<double> get_neighbors_load(type::Real my_load, MPI_Comm neighborhood) {
    int N;
    MPI_Comm_size(MPI_COMM_WORLD, &N);
    std::vector<double> all_loads(N);
    MPI_Allgather(&my_load, 1, MPI_DOUBLE, all_loads.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    return all_loads;
}

std::pair<int, std::vector<double>> get_neighbors_load(type::Real my_load, int vid, const Communicator &v_comm);

LinearHashMap<int, std::vector<double>, 8>
get_neighbors_load(double my_load, const std::set<int> &neighbors, const std::map<int, Communicator> &v_neighborhood);

template<class A, int N=8>
LinearHashMap<int, std::vector<A>, N>
get_neighbors_info(A my_info, MPI_Datatype datatype, const std::set<int> &neighbors, const std::map<int, Communicator> &v_neighborhood) {
    const auto nb_neighbors = neighbors.size();
    LinearHashMap<int, std::vector<A>, N> neighborhoods_info;

    std::vector<MPI_Request> srequests(nb_neighbors);
    int cnt = 0;
    for (int neighbor_rank : neighbors){
        MPI_Isend(&my_info, 1, datatype, neighbor_rank, 40111, MPI_COMM_WORLD, &srequests[cnt]);
        cnt++;
    }

    std::vector<A> neighbors_load(nb_neighbors);
    cnt = 0;
    for (int neighbor_rank : neighbors){
        MPI_Recv(&neighbors_load[cnt], 1, datatype, neighbor_rank, 40111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cnt++;
    }

    MPI_Waitall(nb_neighbors, srequests.data(), MPI_STATUSES_IGNORE);

    cnt = 0;
    for (const auto &v_comm : v_neighborhood) {
        const auto &comm_ranks = v_comm.second.get_ranks();
        std::vector<A> infos;
        int id_neighbor = 0;
        for (int neighbor_rank : neighbors) {
            if (std::find(comm_ranks.cbegin(), comm_ranks.cend(), neighbor_rank) != comm_ranks.cend()) { //in here?
                infos.push_back(neighbors_load[id_neighbor]);
            }
            id_neighbor++;
        }
        neighborhoods_info[cnt] = std::make_pair(v_comm.first, infos);
        cnt++;
    }
    return neighborhoods_info;
}


namespace std {
template<class InputIt, class BinaryFunction, class T>
T reduce(InputIt begin, InputIt end, BinaryFunction f, T start) {
    T result = start;
    while (begin != end) {
        result = f(result, *begin);
        begin++;
    }
    return result;
}

template<class InputIt1, class InputIt2, class OutputIt>
void zip(InputIt1 begin1, InputIt1 end1, InputIt2 begin2, InputIt2 end2, OutputIt out) {
    while (begin1 != end1 && begin2 != end2) {
        *out = std::make_pair(*begin1, *begin2);
        ++out; ++begin1; ++begin2;
    }
}
}

namespace io {
/**
 * Format of output is:
 * PROCID; x y z
 * */
template<class CartesianPointTransformer, class A>
void scatterplot_3d_output(int rank, std::string filename, const std::vector<A> &elements) {
    CartesianPointTransformer cpf;
    std::ofstream f(filename, std::ofstream::app);
    for (const A &el : elements) {
        lb::Point_3 p = cpf.transform(el);
        f << rank << " " << p.x() << " " << p.y() << " " << p.z() << std::endl;
#ifdef DEBUG
        std::cout << rank << " " << p.x() << " " << p.y() << " " << p.z() << std::endl;
#endif
    }
    f.close();
}
}

template<typename T>
inline T dto(double v) {
    T ret = (T) v;

    if(std::isinf(ret)){
        if(ret == -INFINITY){
            ret = std::numeric_limits<T>::lowest();
        } else {
            ret = std::numeric_limits<T>::max();
        }
    }

    return ret;
}

namespace statistic {
template<class RealType>
std::tuple<RealType, RealType, RealType> sph2cart(RealType azimuth, RealType elevation, RealType r){
    RealType x = r * std::cos(elevation) * std::cos(azimuth);
    RealType y = r * std::cos(elevation) * std::sin(azimuth);
    RealType z = r * std::sin(elevation);
    return std::make_tuple(x,y,z);
}

template<int N, class RealType>
class UniformSphericalDistribution {
    const RealType sphere_radius, spherex, spherey, spherez;
public:
    UniformSphericalDistribution(RealType sphere_radius, RealType spherex, RealType spherey, RealType spherez):
            sphere_radius(sphere_radius), spherex(spherex), spherey(spherey), spherez(spherez) {}

    std::array<RealType, N> operator()(std::mt19937& gen) {
        /*
        r1 = (np.random.uniform(0, 1 , n)*(b**3-a**3)+a**3)**(1/3);
        phi1 = np.arccos(-1 + 2*np.random.uniform(0, 1, n));
        th1 = 2*pi*np.random.uniform(0, 1, n);
        x = r1*np.sin(phi1)*np.sin(th1) + X;
        y = r1*np.sin(phi1)*np.cos(th1) + Y;
        z = r1*np.cos(phi1) + Z;
        */
        RealType a = sphere_radius, b = 0.0;
        std::uniform_real_distribution<RealType> udist(0.0, 1.0);

        RealType r1 = std::pow((udist(gen) * (std::pow(b, 3) - std::pow(a, 3)) + std::pow(a, 3)), 1.0/3.0);
        RealType ph1 = std::acos(-1.0 + 2.0 * udist(gen));
        RealType th1 = 2.0 * M_PI * udist(gen);

        auto p = std::make_tuple<RealType, RealType, RealType>(
                r1*std::sin(ph1) * std::sin(th1),
                r1*std::sin(ph1) * std::cos(th1),
                r1*std::cos(ph1)
        );
        if(N > 2) return {(std::get<0>(p))+spherex, std::get<1>(p)+spherey, std::get<2>(p)+spherez};
        else return {std::get<0>(p)+spherex, std::get<1>(p)+spherey};
    }
};

template<int N, class RealType>
class NormalSphericalDistribution {
    const RealType sphere_size, spherex, spherey, spherez;
public:
    NormalSphericalDistribution(RealType sphere_size, RealType spherex, RealType spherey, RealType spherez):
            sphere_size(sphere_size), spherex(spherex), spherey(spherey), spherez(spherez) {}

    std::array<RealType, N> operator()(std::mt19937& gen) {
        std::array<RealType, N> res;
        std::normal_distribution<RealType> ndistx(spherex, sphere_size/2.0); // could do better
        std::normal_distribution<RealType> ndisty(spherey, sphere_size/2.0); // could do better
        if(N == 3) {
            RealType x,y,z;
            do {
                std::normal_distribution<RealType> ndistz(spherez, sphere_size/2.0); // could do better
                x = ndistx(gen);
                y = ndisty(gen);
                z = ndistz(gen);
                res[0] = x;
                res[1] = y;
                res[2] = z;
            } while( (spherex-x)*(spherex-x) + (spherey-y)*(spherey-y) + (spherez-z)*(spherez-z) <= (sphere_size*sphere_size/4.0) );
        } else {
            RealType x,y;
            do {
                x = ndistx(gen);
                y = ndisty(gen);
                res[0] = x;
                res[1] = y;
            } while((spherex-x)*(spherex-x) + (spherey-y)*(spherey-y) <= (sphere_size*sphere_size/4.0) );
        }
        return res;
    }
};


/**
 * From http://www.tangentex.com/RegLin.htm
 * @tparam ContainerA
 * @tparam ContainerB
 * @param x x data
 * @param y y data
 * @return (a,b) of ax+b
 */
template<typename Realtype, typename ContainerA, typename ContainerB>
std::pair<Realtype, Realtype> linear_regression(const ContainerA& x, const ContainerB& y) {
    int i; Realtype xsomme, ysomme, xysomme, xxsomme;

    Realtype ai, bi;

    xsomme = 0.0; ysomme = 0.0;
    xysomme = 0.0; xxsomme = 0.0;
    const int n = x.size();
    for (i=0;i<n;i++) {
        xsomme = xsomme + x[i]; ysomme = ysomme + y[i];
        xysomme = xysomme + x[i]*y[i];
        xxsomme = xxsomme + x[i]*x[i];
    }
    ai = (n*xysomme - xsomme*ysomme)/(n*xxsomme - xsomme*xsomme);
    bi = (ysomme - ai*xsomme)/n;

    return std::make_pair(ai, bi);
}

} // end of namespace statistic

#endif //ADLBIRREG_UTILS_HPP
