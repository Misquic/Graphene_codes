#ifndef UTILITIES_H
#define UTILITIES_H

#ifdef DEBUG
#include <stacktrace>
#define assertVerbose(cond, mess) do \
{ \
  if (!(cond)) \
  { \
    std::cerr << mess << "\nStacktrace:\n" << std::stacktrace::current() << '\n'; \
    throw std::runtime_error(#cond " is false\n"); \
  } \
} \
while (false) \

#else
#define assertVerbose(cond, mess)
#endif

#ifdef DEBUG
#define dmsg(x) std::cerr << x << '\n'
#else
#define dmsg(x)
#endif

#define RED(r) "\033[38;2;" #r ";0;0m"
#define GREEN(g) "\033[38;2;0;" #g ";0m"
#define BLUE(b) "\033[38;2;0;0;" #b "m"
#define RESET "\033[0m"

#ifdef DEBUG
#define LOG(x) std::cerr << RED(200) "File: " << __FILE_NAME__ << ":"\
<<__LINE__<<" line:\n" GREEN(200) <<#x<<RESET"\n";\
x
#else
#define LOG(x)
#endif

#ifdef DEBUG
#define LOG_VAR_STR(x) RED(200) #x ": " GREEN(200) +  str(x) + "\n" RESET
#define LOG_VAR(x) std::cerr << RED(200) #x ": " GREEN(200)<< x << "\n" RESET;
#else
#define LOG_VAR_STR(x)
#define LOG_VAR(x)
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <random>
#include <cassert>
#include <type_traits>
#include <algorithm>
#include <functional>
#include <iomanip>

////////////////// always implemented: //////////////////

double newton(int n, int k);
size_t newton_size_t(int n, int k);
double newton_double(int n, int k);
double time();

template<typename T>
inline constexpr T clamp(const T value, const T min, const T max)
{
  if (value < min) return min;
  if (value > max) return max;
  return value;
}

template<typename T>
double progressBar(const T current, const T bound1, const T bound2)
{
  const T min = std::min(bound1, bound2);
  const T max = std::max(bound1, bound2);
  const T range = max - min;
  const T dist = current - min;
  const double percent = dist / range * 100;

  std::cout << "\rProgress = " << std::setw(5) << std::setprecision(4) <<
               percent << " %" << std::flush;
  return percent;
}

template<typename T>
std::vector<T> linspace(const T start, const T end, const T step)
{
  const size_t size = (end - start)/step + 1;

  std::vector<T> numbers(size);

  size_t index = 0;
  for (T& value: numbers)
  {
    value = start + index * step;

    index++;
  }

  return numbers;
};

template<typename T>
std::vector<T> linspaceN(const T start, const T end, const size_t num)
{
  std::vector<T> numbers(num);

  T step = (end - start)/(num - 1);

  size_t index = 0;
  for (T& value: numbers)
  {
    value = start + index * step;

    index++;
  }

  return numbers;
};

template<typename T>
constexpr inline T pow2(T x){
  return x*x;
}

constexpr inline double signum(double x)
{
  return x > 0. ? 1. : -1.;
};

constexpr inline int signum(int x)
{
  return x > 0 ? 1 : -1;
};

constexpr inline long signum(long x)
{
  return x > 0l ? 1l : -1l;
};

template<typename T>
constexpr inline bool isBetween(const T value, const T bound1, const T bound2)
{
  return (bound1 <= value && value <= bound2) ||
         (bound2 <= value && value <= bound1);
};

class Rnd{
protected:
  std::mt19937 mt_gen;
  std::uniform_real_distribution<double> rnd_dist;
  std::normal_distribution<double> norm_dist;
public:
  Rnd(); // mt19937 has default seed;
  Rnd(unsigned seed); //our own seed
  double operator()();
  double operator()(const double& min, const double& max);
  double norm(double mu = 0, double sigma = 1);
};

template<class T>
double avg(const std::vector<T>& vec){
  size_t size = vec.size();
  double sum = 0;
  for(size_t i = 0; i < size; i++){
    sum+=vec[i];
  }
  return sum/size;
};

template<class T>
double avg(const std::vector<std::vector<T>>& mat){
  size_t size_i = mat.size();
  size_t size_j = mat[0].size();
  double sum = 0;
  for(size_t i = 0; i < size_i; i++){
    for(size_t j = 0; j < size_j; j++){
      sum+=mat[i];
    }
  }
  return sum/(size_i*size_j);
};

////////////////// operators on std::vector //////////////////

template<class T>
std::vector<T> operator*(const std::vector<T>& vec, double x){
  const size_t size = vec.size();
  std::vector<T> result;
  result.reserve(size);
  for(size_t i = 0; i < size; i ++){
    result.push_back(vec[i] * x);
  }
  return result;
};

template<class T>
std::vector<T> operator+(const std::vector<T>& vec, double x){
  const size_t size = vec.size();
  std::vector<T> result;
  result.reserve(size);
  for(size_t i = 0; i < size; i ++){
    result.push_back(vec[i] + x);
  }
  return result;
};

template<class T>
std::vector<T> operator-(const std::vector<T>& vec, double x){
  const size_t size = vec.size();
  std::vector<T> result;
  result.reserve(size);
  for(size_t i = 0; i < size; i ++){
    result.push_back(vec[i] - x);
  }
  return result;
};

template<class T>
std::vector<T> operator*(double x, const std::vector<T>& vec){
  return vec*x;
};

template<class T>
std::vector<T> operator+(double x, const std::vector<T>& vec){
  return vec+x;
};

template<class T>
std::vector<T> operator/(const std::vector<T>& vec, double x){
  assert(x != 0);
  return vec*(1.0/x);
};

template<class T>
std::vector<T>& operator*=(std::vector<T>& vec, double x){
  const size_t size = vec.size();
  for(size_t i = 0; i < size; i ++){
    vec[i]*=x;
  }
  return vec;
};


template<class T>
std::vector<T>& operator+=(std::vector<T>& vec, double x){
  const size_t size = vec.size();
  for(size_t i = 0; i < size; i ++){
    vec[i]+=x;
  }
  return vec;
};


template<class T>
std::vector<T>& operator-=(std::vector<T>& vec, double x){
  const size_t size = vec.size();
  for(size_t i = 0; i < size; i ++){
    vec[i]-=x;
  }
  return vec;
};

template<class T>
std::vector<T>& operator/=(std::vector<T>& vec, double x){
  return vec*=(1./x);
};


////////////////// end of (vec, val) //////////////////
////////////////// (vec, vec) //////////////////

template<class T>
std::vector<T> operator-(const std::vector<T>& vec1, const std::vector<T>& vec2){
  assert(vec1.size() == vec2.size());
  const size_t size = vec1.size();
  std::vector<T> result;
  result.resize(size);
  for(size_t i = 0; i < size; i ++){
    result.push_back(vec1[i] - vec2[i]);
  }
  return result;
};

template<class T>
std::vector<T> operator+(const std::vector<T>& vec1, const std::vector<T>& vec2){
  assert(vec1.size() == vec2.size());
  const size_t size = vec1.size();
  std::vector<T> result;
  result.resize(size);
  for(size_t i = 0; i < size; i ++){
    result.emplace_back(vec1[i] + vec2[i]);
  }
  return result;
};

// template<class T>
// std::vector<T> operator*(const std::vector<T>& vec1, const std::vector<T>& vec2){
  // assert(vec1.size() == vec2.size());
  // const size_t size = vec1.size();
  // std::vector<T> result;
  // result.reserve(size);
  // for(size_t i = 0; i < size; i ++){
    // result.emplace_back(vec1[i] * vec2[i]);
  // }
  // return result;
// };

template<class T>
std::vector<T>& operator+=(std::vector<T>& vec1, const std::vector<T>& vec2){
  assert(vec1.size() == vec2.size());
  const size_t size = vec1.size();
  for(size_t i = 0; i < size; i ++){
    vec1[i] += vec2[i];
  }
  return vec1;
};

template<class T>
std::vector<T>& operator-=(std::vector<T>& vec1, const std::vector<T>& vec2){
  assert(vec1.size() == vec2.size());
  const size_t size = vec1.size();
  for(size_t i = 0; i < size; i ++){
    vec1[i] -= vec2[i];
  }
  return vec1;
};
template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec){
  if(vec.size() == 0){
    return out;
  }
  out << vec[0];
  for(size_t i = 1; i < vec.size(); i++){
    out << "," << vec[i];
  }
  return out;
}
//////////////////////////////////////////////////////////////////////
// Array2D
//////////////////////////////////////////////////////////////////////
template<class data_type>
class Array2D{
private:
  size_t m_Ni;
  size_t m_Nj;
  std::vector<data_type> m_data;
public:
  explicit Array2D(size_t Ni, size_t Nj) noexcept;
  explicit Array2D(size_t Ni, size_t Nj, data_type val) noexcept;
  explicit Array2D(const std::vector<std::vector<data_type>>& vec);

  Array2D(const Array2D& other);
  Array2D(Array2D&& other);

  Array2D& operator=(const Array2D& other);
  Array2D& operator=(Array2D&& other);

  inline data_type& operator()(size_t i, size_t j);
  inline data_type& operator[](size_t index);
  inline const data_type& operator()(size_t i, size_t j) const;
  inline const data_type& operator[](size_t index) const;

  //methods
  inline size_t size() const;
  inline size_t sizeI() const;
  inline size_t sizeJ() const;
  typename std::vector<data_type>::iterator begin();
  typename std::vector<data_type>::iterator end();
  typename std::vector<data_type>::const_iterator begin() const;
  typename std::vector<data_type>::const_iterator end() const;
  void fill(std::function<data_type()> func);
  void apply(std::function<data_type(data_type a)> func); // lambda that takes value of data and returns what to input there
  void apply_ref(std::function<data_type(data_type &a)> func); // lambda that takes value of data and returns what to input there
  data_type max() const;
  data_type min() const;
  void transpose_inplace();
  Array2D<data_type> transpose() const;
  std::vector<data_type> mean(int axis) const;
  data_type mean() const;

  //getters
  std::vector<data_type> getVector() const;
  std::vector<data_type>& getVectorRef();

  //arithmetic operators
  std::vector<data_type> operator*(const std::vector<data_type>& vec) const;
  Array2D<data_type> operator*(const data_type val) const;
  Array2D<data_type> operator+(const Array2D<data_type>& other) const;


  //logic operators
  bool operator==(const Array2D<data_type>& other) const;
  bool operator!=(const Array2D<data_type>& other) const;

  template <class T>
  friend std::ostream& operator<<(std::ostream& out, const Array2D<T>& array);
};
template<class data_type>
Array2D<data_type> operator*(const data_type val, const Array2D<data_type>& arr);

template<class data_type>
Array2D<data_type>::Array2D(size_t Ni, size_t Nj) noexcept: m_Ni{Ni}, m_Nj{Nj}{
  assert(Ni > 0 && Nj > 0);
  m_data.resize(Ni*Nj);
};
template<class data_type>
Array2D<data_type>::Array2D(size_t Ni, size_t Nj, data_type val) noexcept: m_Ni{Ni}, m_Nj{Nj}{
  assert(Ni > 0 && Nj > 0);
  m_data.resize(Ni*Nj, val);
};
template<class data_type>
Array2D<data_type>::Array2D(const std::vector<std::vector<data_type>>& vec): m_Ni{vec.size()}, m_Nj{vec[0].size()}{
  // #ifdef DEBUG
  for(size_t i = 0; i < m_Ni; i++){
    assert(m_Nj == vec[i].size());
    if(m_Nj != vec[i].size()) throw std::runtime_error("Cannot convert non square vec<vec<T>> to Array2D<T>!\n");
  }
  // #endif
  m_data.resize(m_Ni*m_Nj);
  for(size_t i = 0; i < m_Ni; i++){
    for(size_t j = 0; j < m_Nj; j++){
      m_data[i*m_Nj + j] = vec[i][j];
    }
  }
};

template<class data_type>
Array2D<data_type>::Array2D(const Array2D<data_type>& other): m_Ni{other.m_Ni}, m_Nj{other.m_Nj}{
  dmsg("Array2D copying constructor\n");
  m_data = other.m_data;
};
template<class data_type>
Array2D<data_type>::Array2D(Array2D<data_type>&& other): m_Ni{other.m_Ni}, m_Nj{other.m_Nj}{
  dmsg("Array2D moving constructor\n");
  m_data = std::move(other.m_data);
};

template<class data_type>
Array2D<data_type>& Array2D<data_type>::operator=(const Array2D<data_type>& other){
  dmsg("Array2D copying operator\n");
  assert(m_Ni == other.m_Ni && m_Nj == other.m_Nj);
  if(this!= &other){
    if(m_Ni != other.m_Ni || m_Nj != other.m_Nj){
      throw std::runtime_error("assigning different sizes Array2D is prohibited!\n");
    }
    m_data = other.m_data;
  }
  return (*this);
};

template<class data_type>
Array2D<data_type>& Array2D<data_type>::operator=(Array2D<data_type>&& other){
  dmsg("Array2D moving operator\n");
  assert(m_Ni == other.m_Ni && m_Nj == other.m_Nj);
  if(this!= &other){
    if(m_Ni != other.m_Ni || m_Nj != other.m_Nj){
      throw std::runtime_error("assigning different sizes Array2D is prohibited!\n");
    }
    m_data = std::move(other.m_data);
  }
  return (*this);
};


template<class data_type>
data_type& Array2D<data_type>::operator()(size_t i, size_t j){
  assertVerbose(i < m_Ni && j < m_Nj, i << ' ' << m_Ni << ' ' << j << ' ' << m_Nj);
  return m_data[i*m_Nj+j];
}
template<class data_type>
data_type& Array2D<data_type>::operator[](size_t index){
  assert(index < m_Ni*m_Nj);
  return m_data[index];
};

template<class data_type>
const data_type& Array2D<data_type>::operator()(size_t i, size_t j) const{
  assertVerbose(i < m_Ni && j < m_Nj, i << ' ' << m_Ni << ' ' << j << ' ' << m_Nj);
  return m_data[i*m_Nj+j];
}

template<class data_type>
const data_type& Array2D<data_type>::operator[](size_t index) const{
  assert(index < m_Ni*m_Nj);
  return m_data[index];
};
template<class data_type>
size_t Array2D<data_type>::size() const{
  return m_data.size();
};
template<class data_type>
size_t Array2D<data_type>::sizeI() const{
  return m_Ni;
};
template<class data_type>
size_t Array2D<data_type>::sizeJ() const{
  return m_Nj;
};
template<typename data_type>
typename std::vector<data_type>::iterator Array2D<data_type>::begin(){
  return m_data.begin();
};
template<typename data_type>
typename std::vector<data_type>::const_iterator Array2D<data_type>::begin() const{
  return m_data.begin();
};
template<typename data_type>
typename std::vector<data_type>::iterator Array2D<data_type>::end(){
  return m_data.end();
};
template<typename data_type>
typename std::vector<data_type>::const_iterator Array2D<data_type>::end() const{
  return m_data.end();
};

template<class data_type>
void Array2D<data_type>::fill(std::function<data_type()> func){
  for(size_t i = 0; i < m_Ni*m_Nj; i++){
    m_data[i] = func();
  }
};
template<class data_type>
void Array2D<data_type>::apply(std::function<data_type(data_type a)> func){
  auto it_begin = this->begin();
  auto it_end   = this->end();
  for(; it_begin != it_end; it_begin++){
    *it_begin = func(*it_begin);
  }
};
template<class data_type>
void Array2D<data_type>::apply_ref(std::function<data_type(data_type &a)> func){
  auto it_begin = this->begin();
  auto it_end   = this->end();
  for(; it_begin != it_end; it_begin++){
    *it_begin = func(*it_begin);
  }
};

template<class data_type>
data_type Array2D<data_type>::max() const{
  return std::max_element(m_data.begin(), m_data.end())[0];
};
template<class data_type>
data_type Array2D<data_type>::min() const{
  return std::min_element(m_data.begin(), m_data.end())[0];
};

size_t permute_transpose_inplace(size_t a, size_t new_first_dim, size_t new_second_dim);
size_t inverse_permute_transpose_inplace(size_t a, size_t new_first_dim, size_t new_second_dim);
bool is_minimal_in_cycle(size_t a, size_t new_first_dim, size_t new_second_dim);
template<class data_type>
void Array2D<data_type>::transpose_inplace(){
  if(m_Ni == m_Nj){
    for(size_t i = 0; i < m_Ni; i++){
      for(size_t j = i+1; j < m_Nj; j++){
        std::swap(m_data[i*m_Nj + j], m_data[j*m_Nj + i]);
      }
    }
  }else{
    // std::cout << (*this) << "\n";
    // std::cout << m_data << "\n";
    Array2D<data_type> transposed_good = this->transpose();
    Array2D<int> moved(m_Nj, m_Ni);
    size_t elemets_switched = 2; //start with 2 because fist and last element always stays in its place.
    for(size_t a = 1; a < m_Ni*m_Nj-1; a++){
      size_t next_comes_from = permute_transpose_inplace(a, m_Nj, m_Ni);
      if(next_comes_from == a) {
        elemets_switched++;
        continue;
      };
      if(!is_minimal_in_cycle(a,m_Nj,m_Ni)) continue;

      size_t current_switching_index = a;
      data_type value_that_goes = m_data[a]; // we have to remember that, because on its place we now take other value
      moved[a]++;
      size_t goes_to = inverse_permute_transpose_inplace(a, m_Nj, m_Ni);
      do{
        m_data[current_switching_index] = std::move(m_data[next_comes_from]);
        moved[next_comes_from]++; //debug
        if(moved[next_comes_from] >1){
          std::cout << "wait!";
          std::cout << "moved: \n" << moved << "\n";
          std::cout << "m_data" << m_data << "\n";
          std::cout << "(*this): " << (*this) << "\n";
        }
        current_switching_index = next_comes_from;
        next_comes_from = permute_transpose_inplace(next_comes_from, m_Nj, m_Ni);
        elemets_switched++;
      }
      while(a!=next_comes_from);
      m_data[goes_to] = std::move(value_that_goes);

      elemets_switched++;

      if(elemets_switched == m_Ni*m_Nj){
        break;
      }
    }

    size_t temp = m_Nj;
    m_Nj = m_Ni;
    m_Ni = temp;
    if((*this)!=transposed_good){
      std::cout << "moved:\n" << moved << "\n";
      std::cout << "elements_switched: " << elemets_switched << " N*M: " << m_Ni*m_Nj << "\n";
    }
    // std::swap(m_Ni, m_Nj);
    // std::cout << m_data << "\n";
    // std::cout << (*this) << "\n";
  }
};


template<class data_type>
Array2D<data_type> Array2D<data_type>::transpose() const{
  Array2D<data_type> result(m_Nj, m_Ni);
  for(size_t i = 0; i < m_Ni; i++){
    for(size_t j = 0; j < m_Nj; j++){
      result.m_data[j*m_Ni + i] = m_data[i*m_Nj + j];
    }
  }
  // std::cout << (*this) << "\n\n";
  // std::cout << result << "\n";
  return result;
};
template<class data_type>
std::vector<data_type> Array2D<data_type>::mean(int axis) const{
  if(axis != 0 && axis != 1) throw std::runtime_error("axis must be 0 or 1");

  size_t N_axis = axis==0?m_Ni:m_Nj;
  std::vector<data_type> result;
  result.reserve(N_axis);

  if(axis == 0){
    for(size_t i = 0; i < m_Ni; i++){
      data_type sum{};
      for(size_t j = 0; j < m_Nj; j++){
        sum+=operator()(i,j);
      }
      result.push_back(sum/m_Nj);
    }
  }else{
    for(size_t j = 0; j < m_Nj; j++){
      data_type sum{};
      for(size_t i = 0; i < m_Ni; i++){
        sum+=operator()(i,j);
      }
      result.push_back(sum/m_Ni);
    }
  }

  return result;
};

template<class data_type>
data_type Array2D<data_type>::mean() const{
  data_type sum{};
  for(int index = 0; m_Ni*m_Nj; index++){
    sum+=m_data[index];
  }
  return sum/(m_Ni*m_Nj);
};

//getters
template<class data_type>
std::vector<data_type> Array2D<data_type>::getVector() const
{
  return m_data;
};

template<class data_type>
std::vector<data_type>& Array2D<data_type>::getVectorRef()
{
  return m_data;
};

//arithmetic operators
template<class data_type> //matrix multiplication NxM * Mx1 vec
std::vector<data_type> Array2D<data_type>::operator*(const std::vector<data_type>& vec) const{
  std::vector<data_type> result;
  if(m_Nj!=vec.size()){
    throw std::runtime_error("matrix * vector multiplication invalid sizes!\n");
  }
  result.reserve(m_Ni);
  for(size_t i = 0; i < m_Ni; i++){
    data_type sum = 0;
    for(size_t j = 0; j < m_Nj; j++){
      sum+=operator()(i,j)*vec[j];
    }
    result.push_back(sum);
  }
  return result;
};

template<class data_type>
Array2D<data_type> Array2D<data_type>::operator*(const data_type val) const{
  Array2D<data_type> result((*this));

  for(size_t i = 0; i < m_Ni*m_Nj; i++){
    result[i]*=val;
  }
  return result;
};

template<class data_type>
Array2D<data_type> operator*(const data_type val, const Array2D<data_type>& arr){
  return arr*val;
};
template<class data_type>
Array2D<data_type> Array2D<data_type>::operator+(const Array2D<data_type>& other) const{
  if(m_Ni != other.m_Ni || m_Nj!=other.m_Nj){
    throw std::runtime_error("matrix + matrix invalid sizes!\n");
  }

  Array2D<data_type> result((*this)); // copy
  for(size_t i = 0; i < m_Ni*m_Nj; i++)
  {
    result[i]+=other[i];
  }
  return result;
};


// logic operators
template<class data_type>
bool Array2D<data_type>::operator==(const Array2D<data_type>& other)const {
  if(m_Ni!=other.m_Ni || m_Nj!=other.m_Nj) return false;
  for(size_t i = 0; i < m_Ni*m_Nj; i++){
    if(m_data[i] != other.m_data[i]) return false;
  }
  return true;
};
template<class data_type>
bool Array2D<data_type>::operator!=(const Array2D<data_type>& other)const{
  return !operator==(other);
};



template<class data_type>
std::ostream& operator<<(std::ostream& out, const Array2D<data_type>& array){
  for(size_t i = 0; i < array.m_Ni; i++){
    out << array(i,0);
    for(size_t j = 1; j < array.m_Nj; j++){
      out << ", " << array(i,j);
    }
    out << "\n";
  }
  return out;
};
////////////////// printing and strings //////////////////
template<typename T>
void set_colorbar(std::stringstream& ss, T val, T min, T max){
  constexpr int range = 255;
  double t = static_cast<double>(val-min)/(max-min);
  double tt = t*t;
  double t1 = 1.-t;
  double t1t1 = t1*t1;
  int r = t*range;
  int g = int(15  * t1t1 * tt     * range)%range;
  int b = t1*(int(8.5 * t1t1 * t1 * t * range)%range);
  ss << "\033[48;2;";
  ss << r << ';' << g << ';' << b << 'm';
};
template<typename T,
         typename = std::enable_if_t<std::is_integral_v<T>>>
void set_color_bcg(std::stringstream& ss, T color){
  ss << "\033[48;2;" << color << ';' << color << ';' << color << 'm';
};
template<typename T,
         typename = std::enable_if_t<std::is_integral_v<T>>>
void set_color_bcg(T color){
  std::cout << "\033[48;2;" << color << ';' << color << ';' << color << 'm';
};
template<typename T,
         typename = std::enable_if_t<std::is_integral_v<T>>>
void set_color_bcg(std::stringstream& ss, T r, T g, T b){
  ss << "\033[48;2;" << r << ';' << g << ';' << b << 'm';
};
template<typename T,
         typename = std::enable_if_t<std::is_integral_v<T>>>
void set_color_bcg(T r, T g, T b){
  std::cout << "\033[48;2;" << r << ';' << g << ';' << b << 'm';
};

template<typename T,
         typename = std::enable_if_t<std::is_integral_v<T>>>
void set_color_letter(std::stringstream& ss, T color){
  ss << "\033[38;2;" << color << ';' << color << ';' << color << 'm';
};
template<typename T,
         typename = std::enable_if_t<std::is_integral_v<T>>>
void set_color_letter(T color){
  std::cout << "\033[38;2;" << color << ';' << color << ';' << color << 'm';
};
template<typename T,
         typename = std::enable_if_t<std::is_integral_v<T>>>
void set_color_letter(std::stringstream& ss, T r, T g, T b){
  ss << "\033[38;2;" << r << ';' << g << ';' << b << 'm';
};
template<typename T,
         typename = std::enable_if_t<std::is_integral_v<T>>>
void set_color_letter(T r, T g, T b){
  std::cout << "\033[38;2;" << r << ';' << g << ';' << b << 'm';
};



void reset_color();
void reset_color(std::stringstream& ss);
template<typename T>
void print_hist(const std::vector<T>& hist){
  std::stringstream ss;
  T max = std::max_element(hist.begin(), hist.end())[0];
  T min = std::min_element(hist.begin(), hist.end())[0];

  ss << "\n";
  if(max < 100){
    for(size_t i = 0; i < hist.size(); i++){
      set_colorbar(ss, hist[i], min, max);
      for(int j = 0; j < hist[i]; j++){
        ss.put(' ');
      }
      ss << "\033[0m";
      ss << hist[i];
      ss.put('\n');
    }
  }else{
    for(size_t i = 0; i < hist.size(); i++){
      set_colorbar(ss, hist[i], min, max);
      int j_max = static_cast<double>(hist[i]-min)/(max-min)*100;
      for(int j = 0; j < j_max; j++){
        ss.put(' ');
      }
      ss << "\033[0m";
      ss << hist[i];
      ss.put('\n');
    }
  }
  // ss << "\033[0m";
  ss << std::endl;
  std::cout << ss.str();
};

template<typename T>
void print_hist(const std::vector<std::vector<T>>& hist){
  std::stringstream ss;
  T max = hist[0][0];
  T min = hist[0][0];

  for(size_t i = 0; i < hist.size(); i++){
    for(size_t j = 0; j < hist[i].size(); j++){
      T val = hist[i][j];
      if(val > max){
        max = val;
      }else if(val < min){
        min = val;
      }
    }
  }

  ss << "colorbar:\n" << min;
  for(int i = 0; i < 100; i++){
    set_colorbar(ss, i, 0, 100);
    ss.put(' ');
  }
  ss << "\033[0m";
  ss << max << "\n\n";

  // T abs_max = std::abs(max);
  for(size_t i = 0; i < hist.size(); i++){
    for(size_t j = 0; j < hist[i].size(); j++){
      T val = hist[i][j];

      set_colorbar(ss, val, min, max);
      ss << "  ";
    }
    ss << "\033[0m";
    ss.put('\n');
  }

  // ss << "\033[0m";
  ss << std::endl;
  std::cout << ss.str();
}

template<typename T>
void print_hist(const Array2D<T>& hist, bool cutting = false){
  std::stringstream ss;
  T max = std::max_element(hist.begin(), hist.end())[0];
  T min = std::min_element(hist.begin(), hist.end())[0];

  if(cutting){
    const double higher_bound = 0.99;
    const double lower_bound = 0.01;

    T range = max-min; // count how many below certain parts
    size_t count_below = 0;
    size_t count_above = 0;
    auto start = hist.begin();
    auto stop = hist.end();
    T lower_max = min;
    T higher_min = max;

    for(;start!=stop;++start){
      T val = start[0];
      if(val < min + higher_bound*range){
        count_below++;
        if(val > lower_max){
          lower_max = val;
        }
      }
      if(val > min + lower_bound*range){
        count_above++;
        if(val < higher_min){
          higher_min = val;
        }
      }
    }

    if(static_cast<double>(count_below)/hist.size() > higher_bound){
      max = lower_max;
      dmsg("seting lower_max\n");
    }
    if(static_cast<double>(count_above)/hist.size() > 1-lower_bound){
      min = higher_min;
      dmsg("seting higher_min\n");
    }
  }

  ss << "\ncolorbar:\n" << min;
  for(int i = 0; i < 100; i++){
    set_colorbar(ss, i, 0, 100);
    ss.put(' ');
  }
  ss << "\033[0m";
  ss << max << "\n\n";

  size_t Ni = hist.sizeI();
  size_t Nj = hist.sizeJ();

  // T abs_max = std::abs(max);
  for(size_t i = 0; i < Ni; i++){
    for(size_t j = 0; j < Nj; j++){
      T val = hist(i,j);
      if(cutting){
        val = std::max(val, min);
        val = std::min(val, max);
      }
      set_colorbar(ss, val, min, max);
      // ss << "  ";
      ss.put(' ');
      ss.put(' ');
    }
    ss << "\033[0m";
    ss.put('\n');
    // ss << "\033[B\r"; // to use it it needs making space
  }

  // ss << "\033[0m";
  ss << std::endl;
  std::cout << ss.str();
};




template<class T, class V>
std::ostream& operator<<(std::ostream& out, const std::pair<T,V>& pair){
  out << "(" << pair.first << "," << pair.second << ")";
  return out;
};



template <class T>
std::ostream& operator<<(std::ostream& out, const std::vector<std::vector<T>>& mat){
  if(mat.size() == 0){
    return out;
  }
  // if(mat[0].size() == 0){
  //     return out;
  // }
  for(size_t i = 0; i < mat.size(); i++){
    for(size_t j = 0; j < mat[i].size(); j++){
      out << mat[i][j];
      if(j != mat[i].size()-1){
        out << ",";
      }
    }
    out << "\n";
  }
  return out;
}

template<class T>
std::string str(const T& val){
  return std::to_string(val);
}

template<class T>
bool save(const T& data, const std::string& path = "save.csv", const char mode = 'n'){
  std::ofstream file;
  if(mode == 'a'){
    file.open(path, std::ios::app);
  }else{
    file.open(path);
  }

  if(!file.is_open()){
    std::cerr << "Failed to open file: \"" + path + "\"\n";
    return false;
  }

  try
  {
    file << data;
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return false;
  }
  return true;
}


bool save(double data, const std::string& path, const char mode);
bool save(size_t data, const std::string& path, const char mode);
bool save(int data, const std::string& path, const char mode);

template <class T>
bool save_scalar(T data, const std::string& path, const char mode){
  std::ofstream file;
  if(mode == 'a'){
    file.open(path, std::ios::app);
  }else{
    file.open(path);
  }

  if(!file.is_open()){
    std::cerr << "Failed to open file: \"" + path + "\"\n";
    return false;
  }

  file << data << ",";

  return true;
}


template<typename Func>
void testForBisection(
  double xLeft, double xRight, Func& func)
{
  constexpr unsigned int MAX_BISECTION_TEST_STEPS = 300;

  double yLeft = func(xLeft);
  double yRight = func(xRight);
  if (yLeft == 0 && yRight == 0) return;

  const double delta = (xRight - xLeft) / MAX_BISECTION_TEST_STEPS;

  if (yLeft < yRight)
  {
    assertVerbose(
      yLeft < 0 && yRight > 0,
      "function doesn't cross 0, yLeft=" << yLeft << " yRight=" << yRight << '\n');

    double yPrev = yLeft;
    // go from left to right, every next value must be bigger then previous one
    for (double x = xLeft + delta; x < xRight + delta; x += delta)
    {
      double yCurr = func(x);
      assertVerbose(
        yCurr > yPrev,
        "function is not monotonic: stops increasing on " << x <<
        " xLeft=" << xLeft << " xRight=" << xRight << '\n');

      yPrev = yCurr;
    }
  }
  else if (yLeft > yRight)
  {
    assertVerbose(
      yLeft > 0 && yRight < 0,
      "function doesn't cross 0, yLeft=" << yLeft << " yRight=" << yRight << '\n');

    double yPrev = yLeft;
    // go from left to right, every next value must be less then previous one
    for (double x = xLeft + delta; x < xRight + delta; x += delta)
    {
      double yCurr = func(x);
      assertVerbose(
        yCurr < yPrev,
        "function is not monotonic: stops decreasing on " + str(x) +
        " xLeft=" + str(xLeft) + " xRight=" + str(xRight) << '\n');

      yPrev = yCurr;
    }
  }
  else
  {
    assertVerbose(yLeft != yRight, "Both edge values must be different signs");
  }
}

// finds x for which func(x) is lowerr than tol, !!! make sure func crosses zero in range x_left xright at least and only one time !!!
template<typename Func>
double bisection(double xLeft, double xRight, Func&& func)
{
  #ifdef DEBUG
  testForBisection(xLeft, xRight, func);
  #endif
  constexpr int N_max = 32;
  constexpr double tol = 1e-7;

  double yLeft = func(xLeft);
  double yRight = func(xRight);
  // dmsg("xl: " << xLeft << " xr: " << xRight << " yl: " << yLeft << " yr: " << yRight);

  if (yLeft == 0 && yRight == 0) return (xLeft + xRight) * 0.5;

  if (yLeft < 0 && yRight > 0)
  {
    for(int i = 0; i < N_max; i++){
      double x_new = (xLeft + xRight)*0.5;
      double y_new = func(x_new);

      // dmsg( "y_new: " << y_new << "\n");
      if( std::abs(y_new) < tol )
      {
        return x_new;
      }

      if(y_new < 0)
      {
        xLeft = x_new;
      }
      else
      {
        xRight = x_new;
      }
    }
  }
  else if(yRight < 0 && yLeft > 0)
  {
    for(int i = 0; i < N_max; i++){
      double x_new = (xLeft + xRight)*0.5;
      double y_new = func(x_new);

      // dmsg( "y_new: " << y_new << "\n");
      if( std::abs(y_new) < tol )
      {
        return x_new;
      }

      if(y_new > 0)
      {
        xLeft = x_new;
      }
      else
      {
        xRight = x_new;
      }
    }
  }
  else{
    throw std::invalid_argument("Error in bisection! func(xLeft) == func(xRight)"
      " = 0 or both have same sign: " + str(yLeft) + " | " + str(yRight) + "\n");
  }


  dmsg("returning case 3");
  return (xLeft + xRight) * 0.5;
};

#endif//UTILITIES_H
