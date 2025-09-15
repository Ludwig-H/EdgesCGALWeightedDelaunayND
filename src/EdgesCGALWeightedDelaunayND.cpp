// EdgesCGALWeightedDelaunayND.cpp
// Important: include Eigen BEFORE CGAL NewKernel_d to ensure EIGEN macros are visible
#include <Eigen/Core>

#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/Triangulation.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_set>

namespace npy {

[[noreturn]] void die(const std::string& m){ std::fprintf(stderr,"[error] %s\n", m.c_str()); std::exit(1); }

struct Header {
  std::string descr;            // e.g. "<f8", ">f4", "=f8", "|f8"
  bool fortran=false;
  std::vector<int64_t> shape;   // tuple
};

static Header read_header(std::istream& in){
  Header h;
  char magic[6]; in.read(magic,6);
  if(!in || std::memcmp(magic,"\x93NUMPY",6)!=0) die("fichier .npy invalide (magic)");
  char ver[2]; in.read(ver,2);
  uint32_t hlen32=0; uint16_t hlen16=0; size_t hlen=0;
  if(ver[0]==1){ in.read(reinterpret_cast<char*>(&hlen16),2); hlen=hlen16; }
  else { in.read(reinterpret_cast<char*>(&hlen32),4); hlen=hlen32; }
  std::string hdr(hlen,'\0');
  in.read(hdr.data(), hlen);
  if(!in) die("header .npy tronqué");
  auto p_descr = hdr.find("descr");
  if(p_descr==std::string::npos) die("header invalide (descr)");
  auto colon = hdr.find(':', p_descr);
  if(colon==std::string::npos) die("header invalide (descr)");
  auto q1 = hdr.find('\'', colon);
  auto q2 = hdr.find('\'', q1+1);
  if(q1==std::string::npos || q2==std::string::npos) die("header invalide (descr)");
  h.descr = hdr.substr(q1+1, q2-q1-1);
  auto p_for = hdr.find("fortran_order"); if(p_for==std::string::npos) die("header invalide (fortran_order)");
  h.fortran = hdr.find("True", p_for)!=std::string::npos;
  auto p_shape = hdr.find("shape"); if(p_shape==std::string::npos) die("header invalide (shape)");
  auto lp = hdr.find('(', p_shape), rp = hdr.find(')', lp);
  if(lp==std::string::npos || rp==std::string::npos) die("tuple shape invalide");
  std::string tup = hdr.substr(lp+1, rp-lp-1);
  size_t i=0;
  while(i<tup.size()){
    while(i<tup.size() && (tup[i]==',' || tup[i]==' ')) ++i;
    size_t j=i; while(j<tup.size() && isdigit((unsigned char)tup[j])) ++j;
    if(j>i){ h.shape.push_back(std::stoll(tup.substr(i,j-i))); i=j; } else ++i;
  }
  return h;
}

static inline bool is_little_endian_host(){
  uint16_t x = 1; return *(reinterpret_cast<uint8_t*>(&x))==1;
}

template<class T> inline void bswap_inplace(void* p);
template<> inline void bswap_inplace<float>(void* p){
  uint8_t* b = reinterpret_cast<uint8_t*>(p);
  std::swap(b[0],b[3]); std::swap(b[1],b[2]);
}
template<> inline void bswap_inplace<double>(void* p){
  uint8_t* b = reinterpret_cast<uint8_t*>(p);
  std::swap(b[0],b[7]); std::swap(b[1],b[6]); std::swap(b[2],b[5]); std::swap(b[3],b[4]);
}

template<class OutT>
std::vector<OutT> load_array_real(const std::string& path, std::vector<int64_t>& shape_out){
  std::ifstream f(path, std::ios::binary); if(!f) die("impossible d'ouvrir "+path);
  Header h = read_header(f);
  if(h.fortran) die("fortran_order=True non supporté");
  if(h.descr.size() < 3) die("descr invalide: "+h.descr);
  const char endian = h.descr[0]; const char code = h.descr[1]; const char bytes = h.descr[2];
  const bool host_le = is_little_endian_host();
  const bool is_be = (endian == '>');
  const bool is_le = (endian == '<') || (endian == '='); // '=' = native endian, sur x86 little
  const bool endian_ok = is_le || (!host_le && endian=='=') || (is_be);
  if(!endian_ok) die("endianness non supportée");
  size_t count = 1; for(auto d: h.shape) count *= size_t(d);
  std::vector<OutT> out(count);
  if(code=='f' && bytes=='8'){
    std::vector<double> buf(count);
    f.read(reinterpret_cast<char*>(buf.data()), count*sizeof(double));
    if(!f) die("corpus tronqué "+path);
    if(is_be && host_le){
      for(size_t i=0;i<count;++i) bswap_inplace<double>(&buf[i]);
    }
    for(size_t i=0;i<count;++i) out[i] = OutT(buf[i]);
  } else if(code=='f' && bytes=='4'){
    std::vector<float> buf(count);
    f.read(reinterpret_cast<char*>(buf.data()), count*sizeof(float));
    if(!f) die("corpus tronqué "+path);
    if(is_be && host_le){
      for(size_t i=0;i<count;++i) bswap_inplace<float>(&buf[i]);
    }
    for(size_t i=0;i<count;++i) out[i] = OutT(buf[i]);
  } else {
    die("dtype non supporté: "+h.descr+" (attendu f4 / f8)");
  }
  shape_out = h.shape;
  return out;
}

void save_u64_2col(const std::string& path, const std::vector<std::pair<uint64_t,uint64_t>>& E){
  std::string dict = "{'descr': '<u8', 'fortran_order': False, 'shape': (" + std::to_string(E.size()) + ", 2), }";
  while( (10 + dict.size()) % 16 != 0 ) dict.push_back(' ');
  dict.back() = '\n';
  std::ofstream f(path, std::ios::binary); if(!f) die("impossible d'écrire "+path);
  f.write("\x93NUMPY",6); char ver[2]={1,0}; f.write(ver,2);
  uint16_t hlen = (uint16_t)dict.size(); f.write(reinterpret_cast<char*>(&hlen),2);
  f.write(dict.data(), dict.size());
  for(const auto& e: E){ uint64_t a=e.first,b=e.second; f.write((char*)&a,sizeof a); f.write((char*)&b,sizeof b); }
}

} // namespace npy

struct PairHash{
  size_t operator()(const std::pair<uint64_t,uint64_t>& p) const noexcept{
    uint64_t x = p.first ^ (p.second + 0x9E3779B97F4A7C15ULL + (p.first<<6) + (p.first>>2));
    return static_cast<size_t>(x);
  }
};

namespace WeightedND {

using K  = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;             // modèle RegularTriangulationTraits pour dD
using Adapter = CGAL::Regular_triangulation_traits_adapter<K>;
using Vb = CGAL::Triangulation_vertex<Adapter, std::uint64_t>;     // on stocke l'index d'origine
using Cb = CGAL::Triangulation_full_cell<Adapter>;
using TDS = CGAL::Triangulation_data_structure<CGAL::Dynamic_dimension_tag, Vb, Cb>;
using RT  = CGAL::Regular_triangulation<K, TDS>;

void run(const std::vector<double>& P, const std::vector<double>& W, size_t N, int dim,
         std::vector<std::pair<uint64_t,uint64_t>>& edges)
{
  if(dim<2) npy::die("dimension d>=2 exigée");
  if(W.size()!=N) npy::die("weights.npy: longueur N exigée");

  RT rt(dim);
  using BP = typename K::Point_d;
  using WP = typename K::Weighted_point_d;
  using FT = typename K::FT;

  for(size_t i=0;i<N;++i){
    const double* r = &P[i*size_t(dim)];
    BP p(dim, r, r+dim);
    WP wpt(p, FT(W[i]));
    auto vh = rt.insert(wpt);
    if(vh != RT::Vertex_handle()) vh->data() = (uint64_t)i;
  }

  int curd = rt.tds().current_dimension();
  if(curd < 1){ edges.clear(); return; }
  const int verts_per_cell = curd + 1;

  std::unordered_set<std::pair<uint64_t,uint64_t>, PairHash> E;
  E.reserve(N*3);

  for(auto it = rt.finite_full_cells_begin(); it != rt.finite_full_cells_end(); ++it){
    for(int a=0;a<verts_per_cell;++a){
      uint64_t ia = it->vertex(a)->data();
      for(int b=a+1;b<verts_per_cell;++b){
        uint64_t ib = it->vertex(b)->data();
        if(ia==ib) continue;
        if(ia<ib) E.emplace(ia,ib);
        else      E.emplace(ib,ia);
      }
    }
  }
  edges.assign(E.begin(), E.end());
  std::sort(edges.begin(), edges.end());
}

} // namespace WeightedND

int main(int argc, char** argv){
  if(argc!=4){
    std::fprintf(stderr,"Usage: %s points.npy weights.npy out_edges.npy\n", argv[0]);
    return 64;
  }
  std::vector<int64_t> shpP, shpW;
  auto P = npy::load_array_real<double>(argv[1], shpP);
  auto W = npy::load_array_real<double>(argv[2], shpW);

  if(!(shpP.size()==2)) npy::die("points.npy: forme (N,d) exigée");
  size_t N = (size_t)shpP[0]; int dim=(int)shpP[1];
  size_t NW = (shpW.size()==1? (size_t)shpW[0] : ((shpW.size()==2 && shpW[1]==1)? (size_t)shpW[0] : 0));
  if(N==0) npy::die("nuage vide");
  if(NW!=N) npy::die("weights.npy: longueur N exigée");

  std::vector<std::pair<uint64_t,uint64_t>> edges;
  WeightedND::run(P, W, N, dim, edges);
  npy::save_u64_2col(argv[3], edges);
  std::fprintf(stderr,"[info] N=%zu d=%d edges=%zu\n", N, dim, edges.size());
  return 0;
}
